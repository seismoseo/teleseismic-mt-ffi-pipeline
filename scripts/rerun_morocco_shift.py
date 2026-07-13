"""
rerun_morocco_shift.py — add the WISP cross-correlation time-shift step that the Morocco run
originally SKIPPED (every other event — 2009 Puerto Cabello, Calama, Sanriku — used it).

Faithful mirror of rerun_2009_np1_sota.py: copy the converged plane dir, reproduce the current
inversion (baseline misfit = built-in faithfulness check vs the known 0.1402/0.1365), apply
shift_match2 auto cross-correlation to align observed<->synthetic for body (P+SH) and surf
(Rayleigh+Love), then re-invert with the aligned data. Nothing else changed: same station
weights, same segments, same annealing_prop (max_source_dur=16 preserved by NOT regenerating it).

Run:  /home/msseo/miniforge3/envs/ff-env/bin/python scripts/rerun_morocco_shift.py [NP1 NP2]
"""
import json, os, pathlib, re, shutil, sys

BASE = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026"))
FFM = BASE / "event_2023_morocco_mtuq/20230908221101/ffm.0"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")


def misfit(d):
    return float(re.search(r"averaged misfit error\s+([\d.]+)",
                           open(d / "modelling_summary.txt").read()).group(1))


def run_plane(plane, manual_modelling, shift_match2, save_waveforms, default_dirs):
    SRC = FFM / plane
    DST = FFM / f"{plane}_shift"
    if DST.exists():
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (DST / f).unlink(missing_ok=True)
    if (DST / "plots").exists():
        shutil.rmtree(DST / "plots")

    tensor_info = json.load(open(DST / "tensor_info.json"))
    seg_data = json.load(open(DST / "segments_data.json"))

    # 1) reproduce the current (un-shifted) inversion -> baseline misfit
    manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)
    m0 = misfit(DST)

    # 2) cross-correlation alignment (auto): align observed to the converged synthetic
    cwd = os.getcwd(); os.chdir(DST)
    try:
        save_waveforms("body", shift_match2("body", directory=DST))
        save_waveforms("surf", shift_match2("surf", directory=DST))
    finally:
        os.chdir(cwd)

    # 3) re-invert with the aligned data
    manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)
    m1 = misfit(DST)
    print(f"{plane}: pre-shift {m0:.4f}  post-shift {m1:.4f}  (delta {m0 - m1:+.4f})  -> {DST}")
    return m0, m1


def main():
    sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
    import ffm.management as mng
    from ffm.inversion_chen_new import manual_modelling
    from ffm.shift_match import shift_match2, save_waveforms
    default_dirs = mng.default_dirs(config_path=CONFIG)
    planes = sys.argv[1:] or ["NP1", "NP2"]
    for p in planes:
        run_plane(p, manual_modelling, shift_match2, save_waveforms, default_dirs)


if __name__ == "__main__":
    main()
