"""Reduced-smoothing control for the 2024 steep-plane WISP model. Re-invert NP2_ext with
the slip- and time-smoothness regularization relaxed (and the finest wavelet scale enabled),
to test whether the model is 'too smooth' by choice or by the resolution limit of teleseismic
1-D Green's functions. If relaxing smoothing does NOT improve variance reduction but DOES add
slip roughness, the smoothness is the honest data limit, not an over-regularization artifact.
"""
import json, os, pathlib, re, shutil, sys
sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import ffm.management as mng
from ffm.inversion_chen_new import manual_modelling
import ffm.plot_graphic_NEIC as plot

EV = pathlib.Path(os.path.expanduser(
    "~/works/17.Venezuela_2026/event_2024_spda/20240719015048/ffm.0"))
SRC, DST = EV / "NP2_ext", EV / "NP2_rough"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")

if DST.exists():
    shutil.rmtree(DST)
shutil.copytree(SRC, DST)
for f in ("Solution.txt", "modelling_summary.txt"):
    (DST / f).unlink(missing_ok=True)
if (DST / "plots").exists():
    shutil.rmtree(DST / "plots")


# relax the smoothing regularization: slip_weight + time_weight are the
# slip- and time-smoothness penalization coefficients (0.15 each by default)
ap = DST / "annealing_prop.json"
a = json.load(open(ap))
for k in ("slip_weight", "time_weight"):
    a[k] = round(a[k] * 0.2, 4)
    print("  reduced", k, "-> ", a[k])
json.dump(a, open(ap, "w"), indent=4)

tensor_info = json.load(open(DST / "tensor_info.json"))
seg_data = json.load(open(DST / "segments_data.json"))
default_dirs = mng.default_dirs(config_path=CONFIG)
manual_modelling(tensor_info, ["body"], default_dirs, seg_data, directory=DST)
m = re.search(r"averaged misfit error\s+([\d.]+)",
              open(DST / "modelling_summary.txt").read()).group(1)
print(f"NP2_rough (less smooth, slip+time wt x0.2) misfit {float(m):.4f}  (NP2_ext was 0.1127)")

os.chdir(DST)
plot.plot_misfit(["body"], directory=pathlib.Path("."))
print("ROUGHDONE")
