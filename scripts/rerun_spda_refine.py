"""Robustness re-inversion for the 2024 San Pedro de Atacama Mw 7.4 WISP models.

Both auto planes already used identical data (62 body channels, no drops; surface waves
impossible on either plane — hypocentre 127 km > 125-km surf GF cap), so the fair-comparison
requirement is met by construction. This re-anneals each plane from scratch (NP -> NP_ref)
to verify the misfit ordering is stable against annealing stochasticity, then regenerates
the body-wave fit plots.

Usage: python rerun_spda_refine.py NP1|NP2
"""
import json, os, pathlib, re, shutil, sys

sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import ffm.management as mng
from ffm.inversion_chen_new import manual_modelling
import ffm.plot_graphic_NEIC as plot

EV = pathlib.Path(os.path.expanduser(
    "~/works/17.Venezuela_2026/event_2024_spda/20240719015048/ffm.0"))
NP = sys.argv[1]
SRC, DST = EV / NP, EV / f"{NP}_ref"
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
manual_modelling(tensor_info, ["body"], default_dirs, seg_data, directory=DST)
m = re.search(r"averaged misfit error\s+([\d.]+)",
              open(DST / "modelling_summary.txt").read()).group(1)
print(f"{NP}_ref misfit {float(m):.4f}")

os.chdir(DST)
plot.plot_misfit(["body"], directory=pathlib.Path("."))
print("DONE")
