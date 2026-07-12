"""Resolution diagnostics for the 2024 teleseismic finite-fault inversion:
(a) first-pulse misfit vs P-radiation coefficient (near-nodal stations fit worst);
(b) per-station variance reduction below vs above 0.2 Hz (Jia's teleseismic cap)."""
import os, sys, json
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, rfftfreq
from obspy.taup import TauPyModel

# Helvetica per house style
for f in ("Helvetica", "Arial", "DejaVu Sans"):
    try:
        plt.rcParams["font.family"] = f; break
    except Exception:
        pass
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})

EV = os.path.expanduser("~/works/17.Venezuela_2026/event_2024_spda/20240719015048/ffm.0/NP2_ext")
DT = 0.2

def parse_obs(fn):
    out = {}; lines = open(fn).read().splitlines(); i = 0
    while i < len(lines):
        if lines[i].startswith("name:"):
            nm = lines[i].split(":", 1)[1].strip(); ch = lines[i+1].split(":", 1)[1].strip()
            L = int(lines[i+3].split(":", 1)[1]); j = i+6; vals = []
            while j < len(lines) and not lines[j].startswith("name:"):
                if lines[j].strip():
                    try: vals.append(float(lines[j]))
                    except Exception: pass
                j += 1
            out[(nm, ch)] = np.array(vals[:L]); i = j
        else:
            i += 1
    return out

def parse_syn(fn):
    out = {}; lines = open(fn).read().splitlines(); i = 0
    while i < len(lines):
        p = lines[i].split()
        if len(p) >= 4 and p[0].isdigit():
            n = int(p[0]); nm = p[2]; ch = p[3]
            vals = [float(lines[k].split()[0]) for k in range(i+1, i+1+n)]
            out[(nm, ch)] = np.array(vals); i = i+1+n
        else:
            i += 1
    return out

obs = parse_obs(f"{EV}/waveforms_body.txt"); syn = parse_syn(f"{EV}/synthetics_body.txt")
d = json.load(open(f"{EV}/tele_waves.json"))
strike, dip, rake = np.radians([341.06, 69.21, -95.0])
model = TauPyModel(model="ak135")

def radP(az_deg, toa_deg):
    az = np.radians(az_deg); ih = np.radians(toa_deg); sr = az - strike
    return (np.cos(rake)*np.sin(dip)*np.sin(ih)**2*np.sin(2*sr)
            - np.cos(rake)*np.cos(dip)*np.sin(2*ih)*np.cos(sr)
            + np.sin(rake)*np.sin(2*dip)*(np.cos(ih)**2 - np.sin(ih)**2*np.sin(sr)**2)
            + np.sin(rake)*np.cos(2*dip)*np.sin(2*ih)*np.sin(sr))

rad, fpm, lo, hi, names = [], [], [], [], []
for w in d:
    if w["component"] != "BHZ":
        continue
    k = (w["name"], "BHZ")
    if k not in obs or k not in syn:
        continue
    o = obs[k]; s = syn[k]; n = min(len(o), len(s)); o, s = o[:n], s[:n]
    if n < 100:
        continue
    amax = np.max(np.abs(o)); on = int(np.argmax(np.abs(o) > 0.15*amax)); fp = slice(on, on+40)
    fpm.append(np.sum((o[fp]-s[fp])**2)/(np.sum(o[fp]**2)+1e-30))
    try:
        arr = model.get_travel_times(source_depth_in_km=127, distance_in_degree=w["distance"], phase_list=["P"])
        toa = arr[0].takeoff_angle
    except Exception:
        toa = 25.0
    rad.append(abs(radP(w["azimuth"], toa)))
    O = rfft(o); S = rfft(s); f = rfftfreq(n, DT)
    for band, store in ((f < 0.2, lo), ((f >= 0.2) & (f < 1.0), hi)):
        store.append(1 - np.sum(np.abs(O[band]-S[band])**2)/(np.sum(np.abs(O[band])**2)+1e-30))
    names.append(w["name"])
rad = np.array(rad); fpm = np.array(fpm); lo = np.array(lo); hi = np.array(hi)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
ax1.scatter(rad, fpm, s=28, c="#1f4e79", zorder=3)
ax1.axvspan(0, 0.2, color="orange", alpha=0.15)
ax1.text(0.10, ax1.get_ylim()[1]*0.9 if False else max(fpm)*0.95, "near-nodal", color="#b35c00",
         ha="center", fontsize=9)
for nm, x, y in zip(names, rad, fpm):
    if y > 1.0 or x < 0.15:
        ax1.annotate(nm, (x, y), fontsize=7, xytext=(3, 2), textcoords="offset points")
r = np.corrcoef(rad, fpm)[0, 1]
ax1.set_xlabel("P-wave radiation coefficient |F$_P$|")
ax1.set_ylabel("First-pulse misfit (first ~8 s)")
ax1.set_title(f"(a) Unfit first pulses are near-nodal  (r = {r:+.2f})")

x = np.arange(len(names)); order = np.argsort(lo)[::-1]
ax2.bar(x-0.2, lo[order], width=0.4, color="#2e7d32", label="< 0.2 Hz (Jia band)")
ax2.bar(x+0.2, hi[order], width=0.4, color="#c62828", label="0.2–1.0 Hz")
ax2.axhline(0, color="k", lw=0.6)
ax2.set_xticks([])
ax2.set_xlabel("Stations (sorted by low-band VR)")
ax2.set_ylabel("Variance reduction")
ax2.set_title(f"(b) Low band fits (med {np.median(lo):+.2f}); high band does not (med {np.median(hi):+.2f})")
ax2.legend(loc="lower left", framealpha=1, facecolor="white")
fig.tight_layout()
out = os.path.expanduser("~/works/17.Venezuela_2026/event_2026_calama/report/fig_resolution_2024.png")
fig.savefig(out, dpi=190, bbox_inches="tight")
print("->", out, f"| corr {r:+.2f} | low med {np.median(lo):+.2f} | high med {np.median(hi):+.2f}")
