"""
16_plot_mtuq_fixed.py — regenerate the MTUQ waveform-fit / beachball figures for a
KNOWN solution (from a saved solution json) without re-running the grid search.

Needed because 15_run_mtuq_dc.py overwrites <event>_waveforms.png on every variant
run, so the figure on disk belongs to whichever variant ran last, not necessarily
the adopted solution.

Run (mtuq env, serial):
  python scripts/16_plot_mtuq_fixed.py --event 2026kumamoto \
      --solution results/mtuq_dc_2026kumamoto/solution_sw60.json --tag sw60
"""
import argparse, json, os, sys
import numpy as np

from mtuq import read, download_greens, MomentTensor
from mtuq.event import Origin
from mtuq.graphics import plot_data_greens2, plot_beachball
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util.cap import parse_station_codes, Trapezoid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
bands = import_module("12_run_mtuq_sota").bands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--solution", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--model", default="ak135")
    ap.add_argument("--sh", action="store_true")
    a = ap.parse_args()

    sol = json.load(open(a.solution))
    origin = Origin({"time": sol["time"], "latitude": sol["latitude"],
                     "longitude": sol["longitude"], "depth_in_m": sol["depth_in_m"]})
    best = MomentTensor(np.array([sol[k] for k in
                                  ("Mrr", "Mtt", "Mpp", "Mrt", "Mrp", "Mtp")]))
    lune = {k: sol[k] for k in ("rho", "v", "w", "kappa", "sigma", "h")}
    mag = round(sol["Mw"], 2)

    data_dir = os.path.join(C.DATA_DIR, "mtuq", a.event)
    path_weights = os.path.join(data_dir, "weights.dat")
    out = os.path.join(C.RESULTS_DIR, f"mtuq_dc_{a.event}")
    b = bands(mag)
    process_bw = ProcessData(filter_type="Bandpass", freq_min=b["bwf"][0], freq_max=b["bwf"][1],
        pick_type="taup", taup_model=a.model, window_type="body_wave",
        window_length=b["bww"], capuaf_file=path_weights)
    process_sw = ProcessData(filter_type="Bandpass", freq_min=b["swf"][0], freq_max=b["swf"][1],
        pick_type="taup", taup_model=a.model, window_type="surface_wave",
        window_length=b["sww"], capuaf_file=path_weights)
    misfit_bw = Misfit(norm="L2", time_shift_min=-b["bwts"], time_shift_max=b["bwts"],
        time_shift_groups=(["ZR", "T"] if a.sh else ["ZR"]), normalize=True)
    misfit_sw = Misfit(norm="L2", time_shift_min=-b["swts"], time_shift_max=b["swts"],
        time_shift_groups=["ZR", "T"], normalize=True)

    data = read(os.path.join(data_dir, "*.[zrt]"), format="sac", event_id=a.event,
                station_id_list=parse_station_codes(path_weights),
                tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)
    greens = download_greens(stations, origin, a.model)
    greens.convolve(Trapezoid(magnitude=mag))
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)

    plot_data_greens2(os.path.join(out, f"{a.event}_waveforms_{a.tag}.png"),
        data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw,
        misfit_bw, misfit_sw, stations, origin, best, lune)
    plot_beachball(os.path.join(out, f"{a.event}_beachball_{a.tag}.png"),
                   best, stations, origin)
    print(f"FIXED_PLOT_DONE {a.tag}", flush=True)


if __name__ == "__main__":
    main()
