"""
rerun_morocco_shortdur.py — diagnostic: does the teleseismic data REQUIRE the long (~16 s) source,
or is the late low-amplitude tail (~20-30% of moment after 10 s) an unresolved artifact?

Take the shifted steep-plane run (NP1_shift), shorten the allowed source duration
(max_source_dur 16 -> 10 s, i.e. USGS-like), re-invert on the SAME (already shifted) data, and
compare misfit + total moment + STF. If the misfit barely changes, the tail is not required by the
data -> reconciles our ~20 s / Mw 6.94 with USGS ~10-15 s / Mww 6.8.

Run:  /home/msseo/miniforge3/envs/ff-env/bin/python scripts/rerun_morocco_shortdur.py [MAXDUR]
"""
import json, os, pathlib, re, shutil, sys

BASE = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026"))
FFM = BASE / "event_2023_morocco_mtuq/20230908221101/ffm.0"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
MAXDUR = int(sys.argv[1]) if len(sys.argv) > 1 else 10


def misfit(d):
    return float(re.search(r"averaged misfit error\s+([\d.]+)",
                           open(d / "modelling_summary.txt").read()).group(1))


def total_moment(d):
    return float(re.search(r"total moment of the inversion\s+([\d.eE+]+)",
                           open(d / "modelling_summary.txt").read()).group(1))


def main():
    sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
    import ffm.management as mng
    from ffm.inversion_chen_new import manual_modelling
    default_dirs = mng.default_dirs(config_path=CONFIG)

    SRC = FFM / "NP1_shift"
    DST = FFM / f"NP1_shortdur{MAXDUR}"
    if DST.exists():
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (DST / f).unlink(missing_ok=True)
    if (DST / "plots").exists():
        shutil.rmtree(DST / "plots")

    # shorten the allowed source duration; everything else (shifted data, weights, geometry) identical
    ap = json.load(open(DST / "annealing_prop.json"))
    ap["max_source_dur"] = MAXDUR
    json.dump(ap, open(DST / "annealing_prop.json", "w"), indent=4)

    tensor_info = json.load(open(DST / "tensor_info.json"))
    seg_data = json.load(open(DST / "segments_data.json"))
    manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)

    m = misfit(DST); M0 = total_moment(DST)
    import numpy as np
    Mw = (2 / 3) * (np.log10(M0 * 1e7) - 16.1)
    print(f"maxdur {MAXDUR}s: misfit {m:.4f}  M0 {M0:.3e} Nm  Mw {Mw:.3f}  -> {DST}")
    print(f"(compare NP1_shift: misfit 0.1235  M0 3.263e+19  Mw 6.94)")


if __name__ == "__main__":
    main()
