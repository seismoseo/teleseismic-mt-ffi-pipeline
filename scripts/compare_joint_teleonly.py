"""Compare the tele-only vs body+regional (joint) WISP slip models, for whichever event/plane
directories are passed. Reports peak slip, depth extent, slip-weighted centroid, roughness,
and along-strike/along-dip localization (the width of the slip patch) — the diagnostics that
reveal whether regional near-field data sharpens the rupture image.

Usage: python compare_joint_teleonly.py <teleonly_solution_dir> <joint_solution_dir> [label]
"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig_fault_isochrone as fi

def stats(d):
    s = fi.parse_solution(os.path.join(d, "Solution.txt"))
    slip, z, Dx, Dy = s["slip"], s["depth"], s["Dx"], s["Dy"]
    pk = slip.max()
    m = slip >= 0.15 * pk
    L = m.any(axis=0).sum() * Dx          # along-strike extent of >15% slip (km)
    W = m.any(axis=1).sum() * Dy          # along-dip extent (km)
    rough = (np.abs(np.diff(slip, 2, axis=0)).mean()
             + np.abs(np.diff(slip, 2, axis=1)).mean()) / pk
    return dict(peak=pk, zmin=z[m].min(), zmax=z[m].max(),
                zc=np.average(z, weights=slip), L=L, W=W, rough=rough)

tele, joint = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else ""
a, b = stats(tele), stats(joint)
print(f"\n=== {label} : tele-only vs body+regional ===")
print(f"{'metric':22s}{'tele-only':>12s}{'joint':>12s}")
for k, name in (("peak", "peak slip (m)"), ("zmin", "slip top (km)"), ("zmax", "slip bottom (km)"),
                ("zc", "slip centroid (km)"), ("L", "along-strike L (km)"),
                ("W", "along-dip W (km)"), ("rough", "roughness")):
    print(f"{name:22s}{a[k]:12.2f}{b[k]:12.2f}")
print(f"\narea(>15%): tele {a['L']*a['W']:.0f} km^2 -> joint {b['L']*b['W']:.0f} km^2 "
      f"({100*(b['L']*b['W']/(a['L']*a['W'])-1):+.0f}%)")
