"""
ablate_morocco_recipe.py — ABLATION: the Yeck recipe changed four things at once. Which one
actually shortened the source duration (T5-95: 12.3 s -> 8.5 s) and dropped the moment?

From the NP1_shift baseline (54x44 km fault, catalogue location, Mw 6.95 seed) we change exactly
ONE thing at a time and re-invert teleseismic-only, so nothing is confounded:

  A  resize    : right-size the grid to 45x30 km (15x10 @ 3 km).  Location + seed unchanged.
  B  relocate  : move to the InSAR centroid 30.978/-8.332/25 km.  Grid + seed unchanged.
  C  moment    : seed the moment from Mww 6.8 instead of 6.95.    Grid + location unchanged.

(The 4th ingredient, joint InSAR, is already known NOT to change the duration: NP1_insar kept
T5-95 = 12.3 s and Mw 6.94.)

Run:  /home/msseo/miniforge3/envs/ff-env/bin/python scripts/ablate_morocco_recipe.py
"""
import json, os, pathlib, re, shutil, sys
import numpy as np

BASE = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026"))
FFM = BASE / "event_2023_morocco_mtuq/20230908221101/ffm.0"
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
CEN = (30.978, -8.332, 25.0)


def stats(d):
    m = float(re.search(r"averaged misfit error\s+([\d.]+)", open(d / "modelling_summary.txt").read()).group(1))
    t, mr = [], []
    for ln in open(d / "STF.txt"):
        p = ln.split()
        if len(p) == 2 and p[0][0].isdigit():
            t.append(float(p[0])); mr.append(float(p[1]))
    t = np.array(t); mr = np.array(mr); dt = t[1] - t[0]
    cum = np.cumsum(mr) * dt; M0 = cum[-1]
    Mw = (2 / 3) * (np.log10(M0 * 1e7) - 16.1)
    t5 = t[np.searchsorted(cum, 0.05 * M0)]; t95 = t[np.searchsorted(cum, 0.95 * M0)]
    aft10 = (1 - cum[np.searchsorted(t, 10)] / M0) * 100
    return m, Mw, t95 - t5, aft10


def variant(name, mutate):
    sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
    import ffm.management as mng
    import ffm.modelling_parameters as mp
    from ffm.inversion_chen_new import manual_modelling
    default_dirs = mng.default_dirs(config_path=CONFIG)

    dst = FFM / f"ABL_{name}"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(FFM / "NP1_shift", dst)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (dst / f).unlink(missing_ok=True)
    if (dst / "plots").exists():
        shutil.rmtree(dst / "plots")

    ti = json.load(open(dst / "tensor_info.json"))
    sd = json.load(open(dst / "segments_data.json"))
    mutate(ti, sd["segments"][0])
    json.dump(ti, open(dst / "tensor_info.json", "w"), indent=4)
    json.dump(sd, open(dst / "segments_data.json", "w"), indent=4)

    mp.modelling_prop(ti, sd, data_type=["body", "surf"], directory=dst)
    manual_modelling(ti, ["body", "surf"], default_dirs, sd, directory=dst)
    m, Mw, dur, aft = stats(dst)
    print(f"RESULT {name:10s} misfit {m:.4f}  Mw {Mw:.2f}  T5-95 {dur:5.1f}s  moment>10s {aft:4.0f}%", flush=True)


def m_resize(ti, s):
    s["stk_subfaults"], s["dip_subfaults"] = 15, 10
    s["delta_strike"] = s["delta_dip"] = 3.0
    s["hyp_stk"], s["hyp_dip"] = 8, 5

def m_relocate(ti, s):
    ti["lat"] = ti["centroid_lat"] = CEN[0]
    ti["lon"] = ti["centroid_lon"] = CEN[1]
    ti["depth"] = ti["centroid_depth"] = CEN[2]

def m_moment(ti, s):
    ti["moment_mag"] = 10 ** (1.5 * 6.8 + 16.1)


if __name__ == "__main__":
    print("BASELINE   NP1_shift  misfit 0.1235  Mw 6.94  T5-95  12.3s  moment>10s  32%", flush=True)
    print("(known)    +InSAR     misfit 0.1236  Mw 6.94  T5-95  12.3s  moment>10s  32%  <- InSAR alone: no change", flush=True)
    for name, fn in (("resize", m_resize), ("relocate", m_relocate), ("moment", m_moment)):
        variant(name, fn)
    print("FULL RECIPE           misfit 0.1045  Mw 6.80  T5-95   8.5s  moment>10s  11%", flush=True)
