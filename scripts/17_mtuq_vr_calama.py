"""
17_mtuq_vr_calama.py — variance reduction for the FINAL Calama MTUQ solution
(SW-only DC, results/mtuq_dc_calama_swonly) plus the USGS Mww reference mechanism
evaluated on the same data / same 85-km GFs.

Surface-wave VR is the inversion-relevant number (the final search was SW-only:
the 40-s body window at 109/85 km depth contains P+pP+sP and is not modeled by a
point-source direct-P Green's function — documented deep-event limitation).
Body-wave VR is printed as a diagnostic only.

Run:  conda run -n mtuq python -u scripts/17_mtuq_vr_calama.py
"""
import os, sys, json
import numpy as np

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.misfit import Misfit
from mtuq.misfit.waveform import calculate_norm_data
from mtuq.process_data import ProcessData
from mtuq.util.cap import parse_station_codes, Trapezoid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
bands = import_module("12_run_mtuq_sota").bands
vr = import_module("14_mtuq_vr").vr

MAG, DEPTH_KM = 6.9, 85.0
ORIGIN = dict(time="2026-05-25T21:52:20", latitude=-22.3667, longitude=-68.6013)
M0 = 10 ** (1.5 * MAG + 9.1)
SOURCES = {
    "MTUQ SW-only DC (14/75/-83)": None,  # loaded from solution json
    "USGS Mww (359/71/-97)": dict(M0=M0, v=0.0, w=0.0, kappa=359.0, sigma=-97.0,
                                  h=np.cos(np.radians(71.0))),
}


def main():
    data_dir = os.path.join(C.DATA_DIR, "mtuq", "2026calama")
    wpath = os.path.join(data_dir, "weights.dat")
    sol = json.load(open(os.path.join(C.RESULTS_DIR, "mtuq_dc_calama_swonly",
                                      "calama_swonly_solution.json")))
    SOURCES["MTUQ SW-only DC (14/75/-83)"] = sol

    origin = Origin({**ORIGIN, "depth_in_m": DEPTH_KM * 1000})
    b = bands(MAG)
    pbw = ProcessData(filter_type="Bandpass", freq_min=b["bwf"][0], freq_max=b["bwf"][1],
        pick_type="taup", taup_model="ak135", window_type="body_wave",
        window_length=b["bww"], capuaf_file=wpath)
    psw = ProcessData(filter_type="Bandpass", freq_min=b["swf"][0], freq_max=b["swf"][1],
        pick_type="taup", taup_model="ak135", window_type="surface_wave",
        window_length=b["sww"], capuaf_file=wpath)
    mbw = Misfit(norm="L2", time_shift_min=-b["bwts"], time_shift_max=b["bwts"],
                 time_shift_groups=["ZR"], normalize=False)
    msw = Misfit(norm="L2", time_shift_min=-b["swts"], time_shift_max=b["swts"],
                 time_shift_groups=["ZR", "T"], normalize=False)

    sids = parse_station_codes(wpath)
    data = read(os.path.join(data_dir, "*.[zrt]"), format="sac", event_id="2026calama",
                station_id_list=sids, tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    greens = download_greens(stations, origin, "ak135")
    greens.convolve(Trapezoid(magnitude=MAG))

    dsw, gsw = data.map(psw), greens.map(psw)
    dbw, gbw = data.map(pbw), greens.map(pbw)
    for name, src in SOURCES.items():
        print(f"\n=== {name} @ {DEPTH_KM:.0f} km ===")
        print(f"  surface waves ({1/b['swf'][1]:.0f}-{1/b['swf'][0]:.0f} s): "
              f"VR = {100*vr(dsw, gsw, msw, src, ['Z','R','T']):.1f} %")
        print(f"  body waves ({1/b['bwf'][1]:.0f}-{1/b['bwf'][0]:.0f} s) [diagnostic]: "
              f"VR = {100*vr(dbw, gbw, mbw, src, ['Z','R']):.1f} %")


if __name__ == "__main__":
    main()
