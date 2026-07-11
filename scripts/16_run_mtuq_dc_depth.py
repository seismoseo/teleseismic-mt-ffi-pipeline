"""
16_run_mtuq_dc_depth.py — DC-constrained MTUQ grid search WITH DEPTH as a search dimension
(MTUQ origins-list API, per GridSearch.DoubleCouple+Magnitude+Depth example).

Same data/bands as the final DC runs (15_run_mtuq_dc.py); the origin depth is now searched
over a user-given list. Produces plot_misfit_depth (misfit vs depth with mechanism/moment
trade-offs) + the usual beachball/solution json for the best depth.

Run (mtuq env, MPI):
  mpirun -n 24 python -u scripts/16_run_mtuq_dc_depth.py --event 2009puertocabello \
      --lat 10.709 --lon -67.927 --time "2009-09-12T20:06:25.47" --mag 6.4 \
      --depths 8,11,14,17,20,24
  mpirun -n 24 python -u scripts/16_run_mtuq_dc_depth.py --event tele_us6000t7zp \
      --lat 10.55 --lon -68.50 --time "2026-06-24T22:04:33" --mag 7.5 \
      --depths 10,15,20,25 --swband 200,450,900
"""
import argparse, os, sys
import numpy as np

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.graphics import plot_beachball, plot_misfit_depth
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


def run(event_id, lat, lon, time, mag, depths_km, npts=25, model="ak135", swband=None,
        nobody=False, bwband=None, onlybody=False):
    data_dir = os.path.join(C.DATA_DIR, "mtuq", event_id)
    path_data = os.path.join(data_dir, "*.[zrt]")
    path_weights = os.path.join(data_dir, "weights.dat")
    out = os.path.join(C.RESULTS_DIR, f"mtuq_dcz_{event_id}")
    os.makedirs(out, exist_ok=True)

    base = Origin({"time": str(time), "latitude": float(lat), "longitude": float(lon),
                   "depth_in_m": depths_km[0] * 1000.0})
    origins = []
    for d in depths_km:
        o = base.copy(); setattr(o, "depth_in_m", d * 1000.0)
        origins.append(o)

    b = bands(mag)
    if swband:
        t1, t2, win = [float(x) for x in swband.split(",")]
        b["swf"] = (1.0 / t2, 1.0 / t1); b["sww"] = win; b["swts"] = t1 / 4.
    if nobody:
        b["use_bw"] = False
    if bwband:   # "Tmin,Tmax,winlen[,ts]" — deep events: 25-60 s keeps the point source valid
        p = [float(x) for x in bwband.split(",")]
        b["bwf"] = (1.0 / p[1], 1.0 / p[0]); b["bww"] = p[2]
        b["bwts"] = p[3] if len(p) > 3 else p[0] / 4.
    process_bw = ProcessData(filter_type="Bandpass", freq_min=b["bwf"][0], freq_max=b["bwf"][1],
        pick_type="taup", taup_model=model, window_type="body_wave",
        window_length=b["bww"], capuaf_file=path_weights)
    process_sw = ProcessData(filter_type="Bandpass", freq_min=b["swf"][0], freq_max=b["swf"][1],
        pick_type="taup", taup_model=model, window_type="surface_wave",
        window_length=b["sww"], capuaf_file=path_weights)
    misfit_bw = Misfit(norm="L2", time_shift_min=-b["bwts"], time_shift_max=b["bwts"],
        time_shift_groups=["ZR"], normalize=True)
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
        print(f"{event_id} DC+depth: depths {depths_km} km | grid {len(grid)} DCs x "
              f"{len(origins)} origins | use_bw={b['use_bw']}", flush=True)

    data = read(path_data, format="sac", event_id=event_id,
                station_id_list=station_id_list, tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)

    if rank == 0:
        print("Downloading GFs for all depths on rank 0...", flush=True)
        greens = download_greens(stations, origins, model)
    comm.Barrier()
    if rank != 0:
        greens = download_greens(stations, origins, model)
    greens.convolve(wavelet)
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)

    results_sw = None if onlybody else grid_search(data_sw, greens_sw, misfit_sw, origins, grid)
    results_bw = grid_search(data_bw, greens_bw, misfit_bw, origins, grid) \
        if (b["use_bw"] or onlybody) else None
    if rank != 0:
        return
    parts = [r for r in (results_sw, results_bw) if r is not None]
    results = parts[0] if len(parts) == 1 else parts[0] + parts[1]

    origin_idx = results.origin_idxmin()
    best_origin = origins[origin_idx]
    idx = results.source_idxmin()
    best = grid.get(idx)
    lune = grid.get_dict(idx)
    merged = merge_dicts(best.as_dict(), lune, {"M0": best.moment()},
                         {"Mw": best.magnitude()},
                         {"best_depth_km": best_origin.depth_in_m / 1000.0})
    print("Best DC+depth: depth %.0f km | Mw %.2f | strike %.0f dip %.0f rake %.0f" % (
          merged["best_depth_km"], merged["Mw"], merged["kappa"],
          np.degrees(np.arccos(merged["h"])), merged["sigma"]), flush=True)
    # misfit vs depth table — MTUQ layout is (sources, origins): origins is the LAST axis
    # (MTUQDataArray._get_shape). reshape(len(origins), -1) scrambles it.
    vals = results.values.reshape(-1, len(origins))
    for _i, o in enumerate(origins):
        print(f"  depth {o.depth_in_m/1000:4.0f} km: min misfit {vals[:, _i].min():.5f}", flush=True)

    plot_misfit_depth(os.path.join(out, f"{event_id}_misfit_depth.png"), results, origins,
                      title=event_id)
    plot_beachball(os.path.join(out, f"{event_id}_beachball.png"), best, stations, best_origin)
    save_json(os.path.join(out, f"{event_id}_solution.json"), merged)
    print(f"Done -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--time", required=True)
    ap.add_argument("--mag", type=float, required=True)
    ap.add_argument("--depths", required=True, help="comma-separated km")
    ap.add_argument("--npts", type=int, default=25)
    ap.add_argument("--swband", default=None)
    ap.add_argument("--nobody", action="store_true", help="surface waves only")
    ap.add_argument("--bwband", default=None, help="body band override: Tmin,Tmax,winlen[,ts]")
    ap.add_argument("--onlybody", action="store_true", help="body waves only")
    a = ap.parse_args()
    run(a.event, a.lat, a.lon, a.time, a.mag,
        [float(x) for x in a.depths.split(",")], a.npts, swband=a.swband,
        nobody=a.nobody, bwband=a.bwband, onlybody=a.onlybody)


if __name__ == "__main__":
    main()
