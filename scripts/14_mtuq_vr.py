"""
14_mtuq_vr.py — variance reduction (VR) of the best MTUQ moment tensors.

VR = 1 - ||d - s||^2 / ||d||^2, computed with MTUQ's own machinery:
un-normalized L2 misfit at the best source (same time-shift settings as the inversion)
over calculate_norm_data. Reported per wave type for both events.

Run:  conda run -n mtuq python -u scripts/14_mtuq_vr.py
"""
import os, sys, json
import numpy as np

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.grid import FullMomentTensorGridSemiregular
from mtuq.misfit import Misfit
from mtuq.misfit.waveform import calculate_norm_data
from mtuq.process_data import ProcessData
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.util.math import to_mij
from mtuq.event import MomentTensor

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
run_mod = import_module("12_run_mtuq_sota")
bands = run_mod.bands

EVENTS = [
    dict(tag="2009puertocabello", lat=10.709, lon=-67.927, depth=14.0,
         time="2009-09-12T20:06:25.47", mag=6.4),
    dict(tag="tele_us6000t7zp", lat=10.435, lon=-68.472, depth=10.0,
         time="2026-06-24T22:05:11", mag=7.5),
]


def best_source(res_dir, tag):
    return json.load(open(os.path.join(res_dir, f"{tag}_solution.json")))


def vr(data, greens, misfit, sol, components):
    """VR = 1 - res/||d||^2 via MTUQ's own misfit engine (un-normalized L2, same
    time-shift settings/groups as the inversion) evaluated at the best source."""
    from mtuq.grid import UnstructuredGrid
    rho = sol["M0"] * np.sqrt(2.0)
    g1 = UnstructuredGrid(dims=("rho", "v", "w", "kappa", "sigma", "h"),
                          coords=[np.array([rho]), np.array([sol["v"]]),
                                  np.array([sol["w"]]), np.array([sol["kappa"]]),
                                  np.array([sol["sigma"]]), np.array([sol["h"]])],
                          callback=to_mij)
    res = float(np.asarray(misfit(data, greens, g1)).min())
    norm = calculate_norm_data(data, "L2", components)
    return 1.0 - res / norm


def main():
    for ev in EVENTS:
        tag = ev["tag"]
        data_dir = os.path.join(C.DATA_DIR, "mtuq", tag)
        res_dir = os.path.join(C.RESULTS_DIR, f"mtuq_sota_{tag}")
        wpath = os.path.join(data_dir, "weights.dat")
        origin = Origin({"time": ev["time"], "latitude": ev["lat"],
                         "longitude": ev["lon"], "depth_in_m": ev["depth"] * 1000})
        b = bands(ev["mag"])
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
        data = read(os.path.join(data_dir, "*.[zrt]"), format="sac", event_id=tag,
                    station_id_list=sids, tags=["units:m", "type:displacement"])
        data.sort_by_distance()
        stations = data.get_stations()
        greens = download_greens(stations, origin, "ak135")
        greens.convolve(Trapezoid(magnitude=ev["mag"]))

        src = best_source(res_dir, tag)
        print(f"\n=== {tag} (Mw grid centre {ev['mag']}) ===")
        dsw, gsw = data.map(psw), greens.map(psw)
        print(f"  surface waves ({1/b['swf'][1]:.0f}-{1/b['swf'][0]:.0f} s): "
              f"VR = {100*vr(dsw, gsw, msw, src, ['Z','R','T']):.1f} %")
        if b["use_bw"]:
            dbw, gbw = data.map(pbw), greens.map(pbw)
            print(f"  body waves ({1/b['bwf'][1]:.0f}-{1/b['bwf'][0]:.0f} s):    "
                  f"VR = {100*vr(dbw, gbw, mbw, src, ['Z','R']):.1f} %")


if __name__ == "__main__":
    main()
