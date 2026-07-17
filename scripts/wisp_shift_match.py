"""
wisp_shift_match.py — apply WISP's cross-correlation time-shift step (shift_match) to an existing
finite-fault run and re-invert. Event-agnostic; a thin wrapper over WISP's own functions.

`ffm model run auto_model` does NOT apply shift_match (it is a separate manual step in the NEIC
operational workflow). Without it, residual timing mismatch from the 1-D velocity model and centroid
mislocation is left in the fit. shift_match2 (auto) cross-correlates each observed trace against the
converged synthetic and shifts it (+-5 s body, +-12 s surface), then we re-invert. On Al Haouz this
cut the misfit ~12-13% and aligned the SH peaks without changing Mw.

For each nodal-plane folder under <run>/ffm.0/: reproduce baseline -> shift_match2 body+surf ->
re-invert. Writes <PLANE>_shift alongside the original (original left untouched).

Usage: wisp_shift_match.py <ffm.0 dir> [NP1 NP2 ...]
"""
import json, os, pathlib, re, shutil, sys
import numpy as np

CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
WISP_SRC = os.path.expanduser("~/works/neic-finitefault/src")


def misfit(d):
    return float(re.search(r"averaged misfit error\s+([\d.eE+-]+)",
                           open(pathlib.Path(d) / "modelling_summary.txt").read()).group(1))


def run_plane(ffm0, plane, manual_modelling, shift_match2, save_waveforms, default_dirs):
    src = ffm0 / plane
    dst = ffm0 / f"{plane}_shift"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (dst / f).unlink(missing_ok=True)
    if (dst / "plots").exists():
        shutil.rmtree(dst / "plots")
    ti = json.load(open(dst / "tensor_info.json"))
    sd = json.load(open(dst / "segments_data.json"))

    manual_modelling(ti, ["body", "surf"], default_dirs, sd, directory=dst)
    m0 = misfit(dst)
    cwd = os.getcwd(); os.chdir(dst)
    try:
        save_waveforms("body", shift_match2("body", directory=dst))
        save_waveforms("surf", shift_match2("surf", directory=dst))
    finally:
        os.chdir(cwd)
    manual_modelling(ti, ["body", "surf"], default_dirs, sd, directory=dst)
    m1 = misfit(dst)
    print(f"{plane}: pre-shift {m0:.4f}  post-shift {m1:.4f}  (delta {m0 - m1:+.4f})  -> {dst}", flush=True)


def main():
    sys.path.insert(0, WISP_SRC)
    import ffm.management as mng
    from ffm.inversion_chen_new import manual_modelling
    from ffm.shift_match import shift_match2, save_waveforms
    default_dirs = mng.default_dirs(config_path=CONFIG)
    ffm0 = pathlib.Path(sys.argv[1]).resolve()
    planes = sys.argv[2:] or [p.name for p in sorted(ffm0.glob("NP*")) if p.is_dir()
                              and (p / "tensor_info.json").exists() and not p.name.endswith("_shift")]
    print(f"planes: {planes}", flush=True)
    for p in planes:
        run_plane(ffm0, p, manual_modelling, shift_match2, save_waveforms, default_dirs)


if __name__ == "__main__":
    main()
