"""
18_test_body_lowfreq.py — test the frequency argument for the Calama deep-event body-wave
failure. Claim: the 10–33 s body band fails because source duration (~15 s) and depth-phase
smear (~8–10 s from 33 km of slip depth extent) violate the point-source approximation, while
at T ≳ 45 s a point source is valid again (but depth resolution is lost, pP unresolved).

Test at 25–60 s body band (120-s window holding the whole P+pP+sP group, ±10 s shifts),
origin at the best-constrained 107 km:
  (a) VR of our SW-only DC (14/75/-83), the USGS Mww (359/71/-97), and the flipped
      body-band solution (4/75/+83). Prediction: VR positive for ours/USGS, flip no longer wins.
  (b) body-only DC grid search — prediction: argmin mechanism is normal (small Kagan vs USGS),
      not the flip.

Run:  mpirun -n 20 python -u scripts/18_test_body_lowfreq.py
"""
import os, sys, json
import numpy as np

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.grid import DoubleCoupleGridRegular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util.cap import parse_station_codes, Trapezoid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
vr = import_module("14_mtuq_vr").vr

MAG, DEPTH_KM = 6.9, 107.0
ORIGIN = dict(time="2026-05-25T21:52:20", latitude=-22.3667, longitude=-68.6013)
BWF = (1.0 / 60.0, 1.0 / 25.0)      # 25-60 s
BWW, BWTS = 120.0, 10.0             # window holds P..sP group (sP-P ~ 33 s) + 15-s STF
M0 = 10 ** (1.5 * MAG + 9.1)


def mech(strike, dip, rake):
    return dict(M0=M0, v=0.0, w=0.0, kappa=float(strike), sigma=float(rake),
                h=float(np.cos(np.radians(dip))))


def main():
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    data_dir = os.path.join(C.DATA_DIR, "mtuq", "2026calama")
    wpath = os.path.join(data_dir, "weights.dat")
    origin = Origin({**ORIGIN, "depth_in_m": DEPTH_KM * 1000})

    process_bw = ProcessData(filter_type="Bandpass", freq_min=BWF[0], freq_max=BWF[1],
        pick_type="taup", taup_model="ak135", window_type="body_wave",
        window_length=BWW, capuaf_file=wpath)
    misfit_bw = Misfit(norm="L2", time_shift_min=-BWTS, time_shift_max=BWTS,
        time_shift_groups=["ZR"], normalize=False)

    sids = parse_station_codes(wpath)
    data = read(os.path.join(data_dir, "*.[zrt]"), format="sac", event_id="2026calama",
                station_id_list=sids, tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    if rank == 0:
        greens = download_greens(stations, origin, "ak135")
    comm.Barrier()
    if rank != 0:
        greens = download_greens(stations, origin, "ak135")
    greens.convolve(Trapezoid(magnitude=MAG))
    dbw, gbw = data.map(process_bw), greens.map(process_bw)

    if rank == 0:
        print(f"Body band {1/BWF[1]:.0f}-{1/BWF[0]:.0f} s | window {BWW:.0f} s | "
              f"shifts +/-{BWTS:.0f} s | depth {DEPTH_KM:.0f} km", flush=True)
        for name, m in [("ours SW-only DC (14/75/-83)", mech(14, 75, -83)),
                        ("USGS Mww (359/71/-97)", mech(359, 71, -97)),
                        ("flipped body sol (4/75/+83)", mech(4, 75, 83))]:
            print(f"  VR {name}: {100*vr(dbw, gbw, misfit_bw, m, ['Z','R']):.1f} %", flush=True)

    grid = DoubleCoupleGridRegular(npts_per_axis=25, magnitudes=[6.75, 6.9, 7.05])
    results = grid_search(dbw, gbw, misfit_bw, origin, grid)
    if rank != 0:
        return
    idx = results.source_idxmin()
    d = grid.get_dict(idx)
    print("Body-only DC argmin: Mw %.2f | strike %.0f dip %.0f rake %.0f" % (
        grid.get(idx).magnitude(), d["kappa"], np.degrees(np.arccos(d["h"])),
        d["sigma"]), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
