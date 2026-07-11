"""Slip figures + stats for the final all-channel refined Calama WISP models."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import fig_wisp_slipdist_5s as w5
import fig_fault_isochrone as fi

BASE = os.path.join(os.path.dirname(HERE), "event_2026_calama/20260525215220")
RUNS = {
    "ffm.0/NP2_ref2": "2026 Calama M6.9 — NP2 (359/71 steep, USGS mechanism), all channels",
    "ffm.0/NP1_ref2": "2026 Calama M6.9 — NP1 (201/20 shallow, USGS mechanism), all channels",
    "ffm.1/NP2_ref2": "2026 Calama M6.9 — steep plane (7/74, this study's mechanism), all channels",
    "ffm.1/NP1_ref2": "2026 Calama M6.9 — shallow plane (194/16, this study's mechanism), all channels",
    "ffm.1/NP1_ref3": "2026 Calama M6.9 — final: shallow plane (194/16, this study's mechanism)",
}
for rel, title in RUNS.items():
    d = os.path.join(BASE, rel)
    if not os.path.isdir(d):
        print("skip (missing):", rel); continue
    print("==", rel)
    w5.plot(d)
    fi.plot(d, dip_ylim="auto", title=title)
    fi.plot_risetime(d)
    s = fi.parse_solution(os.path.join(d, "Solution.txt"))
    slip, rake, depth, trup = s["slip"], s["rake"], s["depth"], s["trup"]
    w = slip / slip.sum()
    mrake = np.degrees(np.arctan2((w*np.sin(np.radians(rake))).sum(),
                                  (w*np.cos(np.radians(rake))).sum()))
    m = slip > 0.15*slip.max()
    zc = (slip*depth).sum()/slip.sum()
    print(f"  strike {s['strike']:.1f} dip {s['dip']:.1f} | slip-wt rake {mrake:.0f} | "
          f"peak {slip.max():.2f} m | depth(>15%) {depth[m].min():.0f}-{depth[m].max():.0f} km | "
          f"z_centroid {zc:.1f} km | t_end {trup[m].max():.0f} s")
