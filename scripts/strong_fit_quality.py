"""Did regional (strong-motion) waveform modelling actually succeed? Reports the per-station
and median variance reduction of the strong-motion fits in a joint WISP run — the decisive
test (the code running is not enough; the near-field waveforms must actually fit). Contrast
with the user's note that regional modelling *failed* for the shallow Venezuela doublet.

Usage: python strong_fit_quality.py <NP_solution_dir>
"""
import sys, os
import numpy as np

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

d = sys.argv[1]
obs = parse_obs(os.path.join(d, "waveforms_strong.txt"))
syn = parse_syn(os.path.join(d, "synthetics_strong.txt"))
vrs = []
for k in sorted(obs):
    if k not in syn:
        continue
    o = obs[k]; s = syn[k]; n = min(len(o), len(s)); o, s = o[:n], s[:n]
    if n < 50 or np.sum(o**2) == 0:
        continue
    vr = 1 - np.sum((o-s)**2)/np.sum(o**2)
    vrs.append((k[0], k[1], vr))
vrs.sort(key=lambda x: x[2])
print(f"strong-motion channels fit: {len(vrs)}")
print(f"{'sta':10s}{'ch':5s}{'VR':>7s}")
for nm, ch, vr in vrs:
    print(f"{nm:10s}{ch:5s}{vr:7.2f}")
v = np.array([x[2] for x in vrs])
print(f"\nmedian strong-motion VR: {np.median(v):+.2f} | fraction with VR>0.3: {np.mean(v>0.3):.0%}")
print("(VR>~0.4 median => regional modelling SUCCEEDED; VR~0 or negative => failed like Venezuela)")
