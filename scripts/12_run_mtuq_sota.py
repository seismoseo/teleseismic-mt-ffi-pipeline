"""
12_run_mtuq_sota.py — SOTA MTUQ full moment-tensor inversion with lune uncertainty.

Full moment-tensor grid search on the eigenvalue lune (Tape & Tape) with syngine ak135
Green's functions and magnitude-appropriate passbands, following MTUQ's own
GridSearch.FullMomentTensor.py + DetailedAnalysis.py guidance. Produces the SOTA figures:
  - lune MISFIT surface       (plot_misfit_lune)
  - lune LIKELIHOOD surface   (plot_likelihood_lune, with estimate_sigma variance)
  - observed-vs-synthetic waveform fits (plot_data_greens2)
  - beachball with station coverage (plot_beachball)

Works for events in the ComCat CSV (--event ID) or given explicitly
(--event TAG --lat --lon --depth --time --mag), so it serves both the 2026 mainshock
(regional data prepped by 10_prep) and the 2009 M6.4 (data prepped separately).

Run (MPI) in the mtuq env:
  conda run -n mtuq mpirun -n 24 python -u scripts/12_run_mtuq_sota.py --event us6000t7zp --mag 7.5
"""
import argparse, os, sys
import numpy as np, pandas as pd

from mtuq import read, download_greens
from mtuq.event import Origin
from mtuq.graphics import (plot_data_greens2, plot_beachball, plot_misfit_lune,
                           plot_likelihood_lune)
from mtuq.grid import FullMomentTensorGridSemiregular, DeviatoricGridSemiregular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.misfit.waveform import estimate_sigma
from mtuq.process_data import ProcessData
from mtuq.util import merge_dicts, save_json
from mtuq.util.cap import parse_station_codes, Trapezoid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C


def bands(mw):
    """(bw_fmin, bw_fmax, bw_win, sw_fmin, sw_fmax, sw_win, bw_ts, sw_ts, use_bw)
    Periods scale with magnitude; large events -> longer periods, surface-wave dominated.
    TELESEISMIC surface waves need LONG periods (>=50-60 s) or 1-D ak135 path-dispersion
    errors dominate; M7+ needs mantle-wave periods (100-250 s) to average over finiteness."""
    if mw < 5.0:
        return dict(bwf=(0.10, 0.333), bww=15., swf=(0.025, 0.0625), sww=150.,
                    bwts=2., swts=10., use_bw=True)
    if mw < 6.0:
        return dict(bwf=(0.05, 0.15), bww=25., swf=(0.02, 0.05), sww=150.,
                    bwts=3., swts=12., use_bw=True)
    if mw < 7.0:                                            # 2009 M6.4 teleseismic
        return dict(bwf=(0.03, 0.10), bww=40., swf=(0.01, 0.02), sww=300.,   # SW 50-100 s
                    bwts=4., swts=20., use_bw=True)
    return dict(bwf=(0.02, 0.05), bww=70., swf=(0.004, 0.01), sww=600.,   # M7+: 100-250 s mantle waves
                bwts=8., swts=40., use_bw=False)            # surface-wave only


def run(event_id, lat, lon, depth_km, time, mag, npts=10, model="ak135", swband=None, swts=None, deviatoric=False, bwband=None, onlybody=False):
    data_dir = os.path.join(C.DATA_DIR, "mtuq", event_id)
    path_data = os.path.join(data_dir, "*.[zrt]")
    path_weights = os.path.join(data_dir, "weights.dat")
    out = os.path.join(C.RESULTS_DIR, f"mtuq_{'dev' if deviatoric else 'sota'}_{event_id}")
    os.makedirs(out, exist_ok=True)

    origin = Origin({"time": str(time), "latitude": float(lat), "longitude": float(lon),
                     "depth_in_m": float(depth_km) * 1000.0})
    b = bands(mag)
    if swband:                       # override surface-wave band, e.g. "200,500,900"
        t1, t2, win = [float(x) for x in swband.split(",")]
        b["swf"] = (1.0/t2, 1.0/t1); b["sww"] = win; b["swts"] = t1/4.
    if swts:
        b["swts"] = swts
    if bwband:                       # enable + set body band, e.g. "35,80,150,12" (deep events)
        p = [float(x) for x in bwband.split(",")]
        b["bwf"] = (1.0/p[1], 1.0/p[0]); b["bww"] = p[2]
        b["bwts"] = p[3] if len(p) > 3 else p[0]/4.
        b["use_bw"] = True
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
    mags = [round(mag + d, 1) for d in (-0.3, -0.15, 0.0, 0.15, 0.3)]
    grid = (DeviatoricGridSemiregular(npts_per_axis=npts, magnitudes=mags) if deviatoric
            else FullMomentTensorGridSemiregular(npts_per_axis=npts, magnitudes=mags))
    wavelet = Trapezoid(magnitude=mag)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        print(f"{event_id}: Mw{mag} ({lat:.3f},{lon:.3f}) {depth_km:.0f} km | "
              f"grid {len(grid)} MTs | use_bw={b['use_bw']}", flush=True)

    data = read(path_data, format="sac", event_id=event_id,
                station_id_list=station_id_list, tags=["units:m", "type:displacement"])
    data.sort_by_distance()
    stations = data.get_stations()
    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)

    # rank 0 downloads Green's functions FIRST (populates the on-disk cache); the other
    # ranks wait, then read from cache — avoids the concurrent-download race that returns
    # a corrupt (non-zip) syngine response.
    if rank == 0:
        print(f"Downloading Green's functions (syngine {model}) for {len(stations)} "
              "stations on rank 0...", flush=True)
        greens = download_greens(stations, origin, model)
    comm.Barrier()
    if rank != 0:
        greens = download_greens(stations, origin, model)
    greens.convolve(wavelet)
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)

    if rank == 0:
        print("Grid search (surface waves)...", flush=True)
    results_sw = None if onlybody else grid_search(data_sw, greens_sw, misfit_sw, origin, grid)
    results_bw = None
    if b["use_bw"] or onlybody:
        if rank == 0:
            print("Grid search (body waves)...", flush=True)
        results_bw = grid_search(data_bw, greens_bw, misfit_bw, origin, grid)
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
    print("Best full MT: Mw %.2f | strike %.0f dip %.0f rake %.0f | "
          "gamma %.1f delta %.1f (CLVD/iso off DC)" % (
          merged["Mw"], merged["kappa"], np.degrees(np.arccos(merged["h"])),
          merged["sigma"], merged.get("gamma", 0), merged.get("delta", 0)), flush=True)

    print("Data variance estimation + plotting SOTA figures...", flush=True)
    sigma_sw = estimate_sigma(data_sw, greens_sw, best, misfit_sw.norm, ["Z", "R", "T"],
                              misfit_sw.time_shift_min, misfit_sw.time_shift_max)

    plot_data_greens2(os.path.join(out, f"{event_id}_waveforms.png"),
        data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw,
        misfit_bw, misfit_sw, stations, origin, best, lune)
    plot_beachball(os.path.join(out, f"{event_id}_beachball.png"), best, stations, origin)
    plot_misfit_lune(os.path.join(out, f"{event_id}_misfit_lune.png"), results,
                     title=f"{event_id}  Mw {merged['Mw']:.2f}  — misfit on the lune")
    plot_likelihood_lune(os.path.join(out, f"{event_id}_likelihood_lune.png"), results_sw,
                         var=sigma_sw**2, title=f"{event_id}  — likelihood (surface waves)")
    save_json(os.path.join(out, f"{event_id}_solution.json"), merged)
    results.save(os.path.join(out, f"{event_id}_misfit.nc"))
    print(f"Done -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True, help="event tag = data/mtuq/<tag>/ dir")
    ap.add_argument("--lat", type=float); ap.add_argument("--lon", type=float)
    ap.add_argument("--depth", type=float, help="km"); ap.add_argument("--time")
    ap.add_argument("--mag", type=float, required=True)
    ap.add_argument("--npts", type=int, default=10)
    ap.add_argument("--swband", default=None, help="Tmin,Tmax,winlen override for surface waves")
    ap.add_argument("--swts", type=float, default=None, help="override surface-wave time-shift range (s)")
    ap.add_argument("--deviatoric", action="store_true", help="constrain to zero trace (CLVD allowed, no isotropic)")
    ap.add_argument("--bwband", default=None, help="Tmin,Tmax,winlen[,ts] to ENABLE body waves (deep events)")
    ap.add_argument("--onlybody", action="store_true", help="body waves only (shallow-event dip-slip surface-wave-indeterminate)")
    a = ap.parse_args()
    if a.lat is None:                                  # pull from ComCat CSV
        ev = pd.read_csv(C.COMCAT_CSV).set_index("event_id").loc[a.event]
        a.lat, a.lon, a.depth, a.time = (float(ev["latitude"]), float(ev["longitude"]),
                                         float(ev["depth"]), str(ev["time"]))
    run(a.event, a.lat, a.lon, a.depth, a.time, a.mag, a.npts, swband=a.swband, swts=a.swts, deviatoric=a.deviatoric, bwband=a.bwband, onlybody=a.onlybody)


if __name__ == "__main__":
    main()
