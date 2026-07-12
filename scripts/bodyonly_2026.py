"""Body-only control for the 2026 Calama finite-fault models. Re-invert the final 2026
runs (ffm.1/NP1_ref3 shallow, NP2_ref3 steep) using ONLY teleseismic body waves — i.e.
discarding exactly the surface waves that WISP could not use for the 2024 event. If the
slip model, peak, depth extent and plane ordering are unchanged, the surface-wave
asymmetry between the two events does not affect any finite-fault conclusion.

Usage: python bodyonly_2026.py NP1|NP2
"""
import json, os, pathlib, re, shutil, sys
sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import ffm.management as mng
from ffm.inversion_chen_new import manual_modelling
import ffm.plot_graphic_NEIC as plot

EV = pathlib.Path(os.path.expanduser(
    "~/works/17.Venezuela_2026/event_2026_calama/20260525215220/ffm.1"))
NP = sys.argv[1]
SRC, DST = EV / f"{NP}_ref3", EV / f"{NP}_bodyonly"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")

if DST.exists():
    shutil.rmtree(DST)
shutil.copytree(SRC, DST)
for f in ("Solution.txt", "modelling_summary.txt"):
    (DST / f).unlink(missing_ok=True)
if (DST / "plots").exists():
    shutil.rmtree(DST / "plots")

tensor_info = json.load(open(DST / "tensor_info.json"))
seg_data = json.load(open(DST / "segments_data.json"))
default_dirs = mng.default_dirs(config_path=CONFIG)
# only body waves — surface waves discarded exactly as for the 2024 event
manual_modelling(tensor_info, ["body"], default_dirs, seg_data, directory=DST)
m = re.search(r"averaged misfit error\s+([\d.]+)",
              open(DST / "modelling_summary.txt").read()).group(1)
print(f"{NP}_bodyonly (2026) body-only misfit {float(m):.4f}")

os.chdir(DST)
plot.plot_misfit(["body"], directory=pathlib.Path("."))
print("BODYONLYDONE")
