"""
15_run_mtuq_dc.py — DOUBLE-COUPLE-CONSTRAINED MTUQ grid search (same final setups as the
full-MT runs in 12_run_mtuq_sota.py). The DC constraint counters the shallow-source
indeterminacy (Kanamori & Given 1981): with vanishing Mrt/Mrp excitation at long periods,
a full-MT search can add "invisible" dip-slip components at zero misfit cost; a DC grid
cannot — orientation components are coupled.

Grid: DoubleCoupleGridRegular (npts_per_axis=40 -> 64k orientations x 5 magnitudes).
Figures: plot_misfit_dc (strike/dip/rake misfit surfaces), beachball, waveform fits
(surface-only figure for the 2026 VLP run), solution json + VR.

Run (mtuq env, MPI):
  mpirun -n 24 python -u scripts/15_run_mtuq_dc.py --event 2009puertocabello \
      --lat 10.709 --lon -67.927 --depth 14 --time "2009-09-12T20:06:25.47" --mag 6.4
  mpirun -n 24 python -u scripts/15_run_mtuq_dc.py --event tele_us6000t7zp \
      --lat 10.55 --lon -68.50 --depth 15 --time "2026-06-24T22:04:33" --mag 7.5 \
      --swband 200,450,900
"""
import argparse, os, sys
import numpy as np

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.graphics import (plot_data_greens1, plot_data_greens2, plot_beachball,
                           plot_misfit_dc)
from mtuq.grid import DoubleCoupleGridRegular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
bands = import_module("12_run_mtuq_sota").bands


def run(event_id, lat, lon, depth_km, time, mag, npts=40, model="ak135", swband=None, swts=None, nobody=False, bwband=None, onlybody=False, sh=False):
    data_dir = os.path.join(C.DATA_DIR, "mtuq", event_id)
    path_data = os.path.join(data_dir, "*.[zrt]")
    path_weights = os.path.join(data_dir, "weights.dat")
    out = os.path.join(C.RESULTS_DIR, f"mtuq_dc_{event_id}")
    os.makedirs(out, exist_ok=True)

    origin = Origin({"time": str(time), "latitude": float(lat), "longitude": float(lon),
                     "depth_in_m": float(depth_km) * 1000.0})
    b = bands(mag)
    if swband:
        t1, t2, win = [float(x) for x in swband.split(",")]
        b["swf"] = (1.0 / t2, 1.0 / t1); b["sww"] = win; b["swts"] = t1 / 4.
    if swts:
        b["swts"] = swts
    if bwband:   # "Tmin,Tmax,winlen[,ts]" — deep events: 25-60 s keeps the point source valid
        p = [float(x) for x in bwband.split(",")]
        b["bwf"] = (1.0 / p[1], 1.0 / p[0]); b["bww"] = p[2]
        b["bwts"] = p[3] if len(p) > 3 else p[0] / 4.
        b["use_bw"] = True   # explicit band = explicit intent (bands() disables bw for Mw>=7)
    if nobody:
        b["use_bw"] = False
    process_bw = ProcessData(filter_type="Bandpass", freq_min=b["bwf"][0], freq_max=b["bwf"][1],
        pick_type="taup", taup_model=model, window_type="body_wave",
        window_length=b["bww"], capuaf_file=path_weights)
    process_sw = ProcessData(filter_type="Bandpass", freq_min=b["swf"][0], freq_max=b["swf"][1],
        pick_type="taup", taup_model=model, window_type="surface_wave",
        window_length=b["sww"], capuaf_file=path_weights)
    misfit_bw = Misfit(norm="L2", time_shift_min=-b["bwts"], time_shift_max=b["bwts"],
        time_shift_groups=(["ZR","T"] if sh else ["ZR"]), normalize=True)
    misfit_sw = Misfit(norm="L2", time_shift_min=-b["swts"], time_shift_max=b["swts"],
        time_shift_groups=["ZR", "T"], normalize=True)

    station_id_list = parse_station_codes(path_weights)
    mags = [round(mag + d, 2) for d in (-0.3, -0.15, 0.0, 0.15, 0.3)]
    grid = DoubleCoupleGridRegular(npts_per_axis=npts, magnitudes=mags)
    wavelet = Trapezoid(magnitude=mag)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        print(f"{event_id} DC: Mw{mag} grid {len(grid)} | use_bw={b['use_bw']}", flush=True)

    data = read(path_data, format="sac", event_id=event_id,
                station_id_list=station_id_list, tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)

    if rank == 0:
        greens = download_greens(stations, origin, model)
    comm.Barrier()
    if rank != 0:
        greens = download_greens(stations, origin, model)
    greens.convolve(wavelet)
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)

    results_sw = None if onlybody else grid_search(data_sw, greens_sw, misfit_sw, origin, grid)
    results_bw = grid_search(data_bw, greens_bw, misfit_bw, origin, grid) if (b["use_bw"] or onlybody) else None
    if rank != 0:
        return
    if onlybody:
        results = results_bw
    else:
        results = results_sw if results_bw is None else results_sw + results_bw

    idx = results.source_idxmin()
    best = grid.get(idx)
    lune = grid.get_dict(idx)
    merged = merge_dicts(best.as_dict(), lune, {"M0": best.moment()},
                         {"Mw": best.magnitude()}, origin)
    print("Best DC: Mw %.2f | strike %.0f dip %.0f rake %.0f" % (
          merged["Mw"], merged["kappa"], np.degrees(np.arccos(merged["h"])),
          merged["sigma"]), flush=True)

    if b["use_bw"]:
        plot_data_greens2(os.path.join(out, f"{event_id}_waveforms.png"),
            data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw,
            misfit_bw, misfit_sw, stations, origin, best, lune)
    else:
        plot_data_greens1(os.path.join(out, f"{event_id}_waveforms_sw.png"),
            data_sw, greens_sw, process_sw, misfit_sw, stations, origin, best, lune)
    plot_beachball(os.path.join(out, f"{event_id}_beachball.png"), best, stations, origin)
    plot_misfit_dc(os.path.join(out, f"{event_id}_misfit_dc.png"), results)
    save_json(os.path.join(out, f"{event_id}_solution.json"), merged)
    results.save(os.path.join(out, f"{event_id}_misfit.nc"))
    print(f"Done -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--depth", type=float, required=True)
    ap.add_argument("--time", required=True)
    ap.add_argument("--mag", type=float, required=True)
    ap.add_argument("--npts", type=int, default=40)
    ap.add_argument("--swband", default=None)
    ap.add_argument("--onlybody", action="store_true", help="body waves only (shallow-event dip-slip is surface-wave-indeterminate)")
    ap.add_argument("--swts", type=float, default=None)
    ap.add_argument("--nobody", action="store_true", help="surface waves only")
    ap.add_argument("--bwband", default=None, help="body band override: Tmin,Tmax,winlen[,ts]")
    a = ap.parse_args()
    run(a.event, a.lat, a.lon, a.depth, a.time, a.mag, a.npts, swband=a.swband, swts=a.swts, nobody=a.nobody, bwband=a.bwband, onlybody=a.onlybody, sh=a.sh)


if __name__ == "__main__":
    main()
