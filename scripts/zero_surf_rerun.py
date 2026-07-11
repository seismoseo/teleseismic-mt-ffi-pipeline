"""Zero the unmodelled surface-wave channels (EDA Rayleigh+Love, PTCN Love) in the final
Calama model (ffm.1/NP1_ref2) and re-invert -> NP1_ref3. WISP plots zero-weight channels grey."""
import json, os, pathlib, re, shutil, sys

sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import ffm.management as mng
from ffm.inversion_chen_new import manual_modelling

EV = pathlib.Path(os.path.expanduser(
    "~/works/17.Venezuela_2026/event_2026_calama/20260525215220/ffm.1"))
SRC, DST = EV / "NP1_ref2", EV / "NP1_ref3"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
ZERO = {("EDA", "BHZ"), ("EDA", "BHT"), ("PTCN", "BHT")}

if DST.exists():
    shutil.rmtree(DST)
shutil.copytree(SRC, DST)
for f in ("Solution.txt", "modelling_summary.txt"):
    (DST / f).unlink(missing_ok=True)
if (DST / "plots").exists():
    shutil.rmtree(DST / "plots")

d = json.load(open(DST / "surf_waves.json"))
n = 0
for w in d:
    if (w["name"], w["component"]) in ZERO:
        w["trace_weight"] = 0.0
        nw = len(str(w.get("wavelet_weight", "0 1 1 1 1 1 1 1")).split())
        w["wavelet_weight"] = " ".join(["0"] * nw) + "\n"
        n += 1
json.dump(d, open(DST / "surf_waves.json", "w"), indent=4)
print(f"zeroed {n} surf channels; active:",
      sum(w["trace_weight"] > 0 for w in d))

tensor_info = json.load(open(DST / "tensor_info.json"))
seg_data = json.load(open(DST / "segments_data.json"))
default_dirs = mng.default_dirs(config_path=CONFIG)
manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)
m = re.search(r"averaged misfit error\s+([\d.]+)",
              open(DST / "modelling_summary.txt").read()).group(1)
print(f"NP1_ref3 misfit {float(m):.4f}")
print("DONE")
