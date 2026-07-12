"""Extend the WISP fault domain and re-invert (body-only), to remove slip clipping at the
fault edges. Copies SRC -> DST, enlarges the subfault grid keeping the hypocentre at the
same physical position (indices shifted by the number of added rows/columns on the
up-dip/along-strike-negative side), recomputes GFs + inversion via manual_modelling.

Usage: python extend_fault_rerun.py <run_dir_parent> <NP> <add_stk_neg> <add_stk_pos> <add_dip_up> <add_dip_down>
e.g.:  python extend_fault_rerun.py event_2024_spda/20240719015048/ffm.0 NP2 2 2 3 3
creates NP2_ext from NP2_ref.
"""
import json, os, pathlib, re, shutil, sys

sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import ffm.management as mng
from ffm.inversion_chen_new import manual_modelling
import ffm.plot_graphic_NEIC as plot

parent = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026")) / sys.argv[1]
NP = sys.argv[2]
a_sn, a_sp, a_du, a_dd = (int(x) for x in sys.argv[3:7])
SRC, DST = parent / f"{NP}_ref", parent / f"{NP}_ext"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")

if DST.exists():
    shutil.rmtree(DST)
shutil.copytree(SRC, DST)
for f in ("Solution.txt", "modelling_summary.txt"):
    (DST / f).unlink(missing_ok=True)
if (DST / "plots").exists():
    shutil.rmtree(DST / "plots")

sd = json.load(open(DST / "segments_data.json"))
s = sd["segments"][0]
s["stk_subfaults"] += a_sn + a_sp
s["dip_subfaults"] += a_du + a_dd
s["hyp_stk"] += a_sn
s["hyp_dip"] += a_du
json.dump(sd, open(DST / "segments_data.json", "w"), indent=4)
print(f"{NP}: grid {s['stk_subfaults']}x{s['dip_subfaults']} "
      f"({s['stk_subfaults']*s['delta_strike']:.0f} x {s['dip_subfaults']*s['delta_dip']:.0f} km), "
      f"hyp index ({s['hyp_stk']},{s['hyp_dip']})")

tensor_info = json.load(open(DST / "tensor_info.json"))
default_dirs = mng.default_dirs(config_path=CONFIG)
manual_modelling(tensor_info, ["body"], default_dirs, sd, directory=DST)
m = re.search(r"averaged misfit error\s+([\d.]+)",
              open(DST / "modelling_summary.txt").read()).group(1)
print(f"{NP}_ext misfit {float(m):.4f}")

os.chdir(DST)
plot.plot_misfit(["body"], directory=pathlib.Path("."))
print("EXTDONE")
