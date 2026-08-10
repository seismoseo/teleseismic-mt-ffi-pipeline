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
    # Deep faults exceed the ~126.5 km surface-wave GF cap; auto_model silently drops
    # surf there but manual_modelling hard-errors, so shift only what the plane can use.
    # surf_waves.json exists whenever data prep ran; synthetics_surf.txt exists only
    # if the plane actually MODELLED surface waves (deep faults silently drop them).
    types = ["body"] + (["surf"] if (src / "synthetics_surf.txt").exists() else [])

    def modelling():
        try:
            manual_modelling(ti, types, default_dirs, sd, directory=dst)
        except Exception as e:
            if "surface waves" in str(e) and "surf" in types:
                types.remove("surf")
                print(f"{plane}: fault below surface-wave GF cap -> body-only shift", flush=True)
                manual_modelling(ti, types, default_dirs, sd, directory=dst)
            else:
                raise

    modelling()
    m0 = misfit(dst)
    cwd = os.getcwd(); os.chdir(dst)
    try:
        for t in types:
            save_waveforms(t, shift_match2(t, directory=dst))
    finally:
        os.chdir(cwd)
    modelling()
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
