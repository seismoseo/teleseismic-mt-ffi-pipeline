"""
rerun_morocco_insar.py — JOINT teleseismic + InSAR inversion, matching USGS's data set.

USGS's finite fault (Barnhart, run with this same WISP code) inverted teleseismic body+surface
waves AND three resampled Sentinel-1 interferograms (1 ascending, 2 descending). Our teleseismic-
only model has a larger fault, slower Vr, longer duration and ~1.6x moment. This adds the SAME
USGS InSAR tracks as a static constraint on the steep plane (NP1_shift) and re-inverts, to test
whether the near-field geodetic data pull our slip toward the compact, lower-moment USGS model
(and thereby shorten the apparent duration).

InSAR is STATIC: it constrains the final slip (extent, amplitude, location, total moment), not
timing directly. A per-track linear ramp is co-estimated (orbital/atmospheric). Data files are the
USGS resampled interferograms (lon lat LOS[m] sx sy sz) -- WISP's native format.

Run:  /home/msseo/miniforge3/envs/ff-env/bin/python scripts/rerun_morocco_insar.py
"""
import glob, json, os, pathlib, re, shutil, sys

BASE = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026"))
FFM = BASE / "event_2023_morocco_mtuq/20230908221101/ffm.0"
INSAR = BASE / "event_2023_morocco_mtuq/insar_usgs"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")


def misfit(d):
    return float(re.search(r"averaged misfit error\s+([\d.]+)",
                           open(d / "modelling_summary.txt").read()).group(1))


def total_moment(d):
    return float(re.search(r"total moment of the inversion\s+([\d.eE+]+)",
                           open(d / "modelling_summary.txt").read()).group(1))


def main():
    sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
    import ffm.management as mng
    from ffm.data_management import imagery_data
    from ffm.input_files import input_chen_imagery
    from ffm.green_functions import gf_retrieve
    from ffm.inversion_chen_new import manual_modelling
    default_dirs = mng.default_dirs(config_path=CONFIG)

    SRC = FFM / "NP1_shift"
    ARCHIVE = FFM / "NP1_insar"          # real run dir, inside the project
    if ARCHIVE.exists():
        shutil.rmtree(ARCHIVE)
    shutil.copytree(SRC, ARCHIVE)
    # WISP static-GF fortran truncates paths at 100 chars -> drive it through a SHORT symlink in /tmp
    DST = pathlib.Path("/tmp/mi")
    if DST.is_symlink() or DST.exists():
        DST.unlink()
    os.symlink(ARCHIVE, DST)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (DST / f).unlink(missing_ok=True)
    if (DST / "plots").exists():
        shutil.rmtree(DST / "plots")

    tracks = sorted(glob.glob(str(INSAR / "s1-*.txt")))
    print(f"InSAR tracks ({len(tracks)}):")
    for t in tracks:
        print("  ", os.path.basename(t))
    ramps = ["linear"] * len(tracks)   # co-estimate a planar ramp per track

    # 1) register imagery + write WISP input files (imagery_data.json -> imagery_data.txt)
    imagery_data(imagery_files=[pathlib.Path(t) for t in tracks], ramp_types=ramps, directory=DST)
    input_chen_imagery(directory=DST)
    print("imagery_data.txt points:", open(DST / "imagery_data.txt").readline().strip())

    # 2) build STATIC Green functions for the imagery points on this fault geometry
    print("building static InSAR Green functions ...")
    gf_retrieve(["imagery"], default_dirs, directory=DST)

    # 3) joint inversion: teleseismic body + surface + InSAR (data already shifted)
    tensor_info = json.load(open(DST / "tensor_info.json"))
    seg_data = json.load(open(DST / "segments_data.json"))
    print("running JOINT body+surf+imagery inversion ...")
    manual_modelling(tensor_info, ["body", "surf", "imagery"], default_dirs, seg_data, directory=DST)

    m = misfit(DST); M0 = total_moment(DST)
    import numpy as np
    Mw = (2 / 3) * (np.log10(M0) - 16.1)   # M0 here is dyne-cm
    print(f"\nJOINT (tele+InSAR): misfit {m:.4f}  M0 {M0:.3e} dyne-cm  Mw {Mw:.3f}  -> {DST}")
    print("compare  tele-only NP1_shift: misfit 0.1235  Mw 6.94  (USGS joint: Mw 6.84)")


if __name__ == "__main__":
    main()
