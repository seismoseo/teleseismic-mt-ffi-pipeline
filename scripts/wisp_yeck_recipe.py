"""
wisp_yeck_recipe.py — apply the Yeck et al. (2023, TSR) / USGS-NEIC best-practice recipe to an
EXISTING WISP run, producing a better-constrained finite-fault model. Event-agnostic.

This is a THIN WRAPPER over WISP's own functions (modelling_parameters.modelling_prop,
data_management.imagery_data, input_files.input_chen_imagery, green_functions.gf_retrieve,
inversion_chen_new.manual_modelling). No inversion machinery is reimplemented.

The recipe (Yeck et al., 2023, doi:10.1785/0320230040), each step a documented USGS choice:
  1. RELOCATE to the geodetically-inferred centroid and NUCLEATE there. A teleseismic hypocentre
     can be badly located when few near stations exist; a uniform-slip InSAR inversion locates the
     centroid far better. WASP time-shifts each teleseismic trace and InSAR carries no timing, so
     rupture-onset location is not resolved anyway -- pinning the fault to the geodetic centroid
     (and treating it as the nucleation point) makes the fault plane consistent with the surface
     displacement field.
  2. RIGHT-SIZE the fault grid. Do not over-size the template: spurious low-amplitude slip spreads
     onto the far edges, inflating both the moment and the apparent rupture duration.
  3. SEED THE MOMENT from the W-phase Mww (the authoritative long-period moment), not from a
     body-wave estimate (which tends to run high).
  4. INVERT JOINTLY with InSAR. Adding geodesy shrinks the source area and raises peak slip
     ("a typical consequence of including regional datasets" -- Yeck et al.).
  (RAKE is already limited to +-20 deg around the CMT rake by WISP's default model_space, which
   matches the USGS window; nothing to change.)

On the 2023 Mw 6.8 Al Haouz (Morocco) event this took our model from
  Mw 6.94 / T5-95 12.3 s / InSAR VR 55% / teleseismic misfit 0.1235   (over-sized, mislocated)
to
  Mw 6.80 / T5-95  8.5 s / InSAR VR 85% / teleseismic misfit 0.1045   (Yeck recipe)
i.e. better on EVERY metric, and in agreement with USGS (Mw 6.84).

NOTE: WISP's static-GF fortran truncates file paths at 100 characters, so the run is driven
through a short symlink in /tmp.

Usage:
  wisp_yeck_recipe.py SRC_DIR DST_DIR --centroid LAT,LON,DEPTH_KM --mww MW \
                      --grid NSTK,NDIP,DXY_KM [--insar INSAR_DIR] [--tag NAME]

Example (Morocco):
  wisp_yeck_recipe.py .../ffm.0/NP1_shift .../ffm.0/NP1_v2 \
      --centroid 30.978,-8.332,25.0 --mww 6.8 --grid 15,10,3.0 \
      --insar .../event_2023_morocco_mtuq/insar_usgs
"""
import argparse, glob, json, os, pathlib, re, shutil, sys
import numpy as np

CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
WISP_SRC = os.path.expanduser("~/works/neic-finitefault/src")


def _misfit(d):
    return float(re.search(r"averaged misfit error\s+([\d.]+)",
                           open(pathlib.Path(d) / "modelling_summary.txt").read()).group(1))


def _moment(d):
    return float(re.search(r"total moment of the inversion\s+([\d.eE+]+)",
                           open(pathlib.Path(d) / "modelling_summary.txt").read()).group(1))


def _report(d, label):
    m, M0 = _misfit(d), _moment(d)
    print(f"  {label}: misfit {m:.4f}  M0 {M0:.3e} dyne-cm  Mw {(2/3)*(np.log10(M0)-16.1):.3f}", flush=True)


def run(src, dst, centroid, mww, grid, insar=None, tag="mv"):
    sys.path.insert(0, WISP_SRC)
    import ffm.management as mng
    import ffm.modelling_parameters as mp
    from ffm.data_management import imagery_data
    from ffm.input_files import input_chen_imagery
    from ffm.green_functions import gf_retrieve
    from ffm.inversion_chen_new import manual_modelling
    default_dirs = mng.default_dirs(config_path=CONFIG)

    src, dst = pathlib.Path(src).resolve(), pathlib.Path(dst).resolve()
    lat, lon, dep = centroid
    nstk, ndip, dxy = int(grid[0]), int(grid[1]), float(grid[2])

    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (dst / f).unlink(missing_ok=True)
    if (dst / "plots").exists():
        shutil.rmtree(dst / "plots")

    # short symlink: WISP static-GF fortran truncates paths at 100 chars
    short = pathlib.Path(f"/tmp/{tag}")
    if short.is_symlink() or short.exists():
        short.unlink()
    os.symlink(dst, short)

    # (1) relocate to the geodetic centroid (= nucleation) ; (3) seed moment from Mww
    ti = json.load(open(short / "tensor_info.json"))
    ti["lat"] = ti["centroid_lat"] = lat
    ti["lon"] = ti["centroid_lon"] = lon
    ti["depth"] = ti["centroid_depth"] = dep
    ti["moment_mag"] = 10 ** (1.5 * mww + 16.1)          # dyne-cm
    json.dump(ti, open(short / "tensor_info.json", "w"), indent=4)

    # (2) right-size the fault; nucleation at the grid centre (= the centroid)
    sd = json.load(open(short / "segments_data.json"))
    s = sd["segments"][0]
    s["stk_subfaults"], s["dip_subfaults"] = nstk, ndip
    s["delta_strike"] = s["delta_dip"] = dxy
    s["hyp_stk"], s["hyp_dip"] = (nstk + 1) // 2, (ndip + 1) // 2
    json.dump(sd, open(short / "segments_data.json", "w"), indent=4)
    print(f"grid {nstk}x{ndip} @ {dxy} km = {nstk*dxy:.0f}x{ndip*dxy:.0f} km | centroid "
          f"{lat},{lon},{dep} km (= nucleation) | seed Mw {mww}", flush=True)

    dtypes = ["body", "surf"] + (["imagery"] if insar else [])
    mp.modelling_prop(ti, sd, data_type=dtypes, directory=short)

    # teleseismic pass: also regenerates the fault geometry + teleseismic GF for the new grid
    print("teleseismic (relocated + right-sized) ...", flush=True)
    manual_modelling(ti, ["body", "surf"], default_dirs, sd, directory=short)
    _report(dst, "teleseismic")

    # (4) joint inversion with InSAR
    if insar:
        tracks = sorted(glob.glob(os.path.join(insar, "*.txt")))
        if not tracks:
            raise SystemExit(f"no InSAR track files (*.txt) in {insar}")
        print(f"InSAR tracks: {[os.path.basename(t) for t in tracks]}", flush=True)
        imagery_data(imagery_files=[pathlib.Path(t) for t in tracks],
                     ramp_types=["linear"] * len(tracks), directory=short)
        input_chen_imagery(directory=short)
        cwd = os.getcwd(); os.chdir(short)
        try:
            gf_retrieve(["imagery"], default_dirs, directory=short)   # static GF for the new grid
        finally:
            os.chdir(cwd)
        print("JOINT teleseismic + InSAR ...", flush=True)
        manual_modelling(ti, ["body", "surf", "imagery"], default_dirs, sd, directory=short)
        _report(dst, "JOINT      ")

    short.unlink(missing_ok=True)
    print(f"DONE -> {dst}", flush=True)
    return dst


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Apply the Yeck et al. (2023) USGS recipe to a WISP run")
    p.add_argument("src", help="existing WISP run directory (source)")
    p.add_argument("dst", help="new run directory to create (originals are left untouched)")
    p.add_argument("--centroid", required=True, help="LAT,LON,DEPTH_KM of the geodetic centroid (= nucleation)")
    p.add_argument("--mww", type=float, required=True, help="W-phase Mww to seed the moment")
    p.add_argument("--grid", required=True, help="NSTK,NDIP,DXY_KM  (e.g. 15,10,3.0 -> 45x30 km)")
    p.add_argument("--insar", default=None, help="directory of resampled InSAR tracks (lon lat LOS sx sy sz)")
    p.add_argument("--tag", default="mv", help="short /tmp symlink name (avoids the 100-char fortran path bug)")
    a = p.parse_args()
    run(a.src, a.dst,
        [float(x) for x in a.centroid.split(",")],
        a.mww,
        a.grid.split(","),
        insar=a.insar, tag=a.tag)
