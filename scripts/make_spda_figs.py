"""Slip figures + stats + stress drop for the 2024 San Pedro de Atacama Mw 7.4 WISP models."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
import fig_wisp_slipdist_5s as w5
import fig_fault_isochrone as fi

BASE = os.path.join(os.path.dirname(HERE), "event_2024_spda/20240719015048/ffm.0")
RUNS = {
    "NP1_ref": "2024 San Pedro de Atacama M7.4 — shallow plane (172/21, USGS FFM choice)",
    "NP2_ref": "2024 San Pedro de Atacama M7.4 — steep plane (341/69, Jia et al. choice)",
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
    Dx, Dy = s["Dx"], s["Dy"]
    w = slip / slip.sum()
    mrake = np.degrees(np.arctan2((w*np.sin(np.radians(rake))).sum(),
                                  (w*np.cos(np.radians(rake))).sum()))
    m = slip > 0.15*slip.max()
    zc = (slip*depth).sum()/slip.sum()
    print(f"  strike {s['strike']:.1f} dip {s['dip']:.1f} | slip-wt rake {mrake:.0f} | "
          f"peak {slip.max():.2f} m | depth(>15%) {depth[m].min():.0f}-{depth[m].max():.0f} km | "
          f"z_centroid {zc:.1f} km | t_end {trup[m].max():.0f} s")
    # stress drop (Eshelby circular crack on threshold-trimmed slip; mu from shear_model.txt)
    mu = float(open(os.path.join(d, "shear_model.txt")).read().split("\n")[2].split()[0]) / 10.
    cell = Dx * Dy * 1e6
    M0_tot = mu * slip.sum() * cell
    print(f"  mu {mu/1e9:.1f} GPa | M0(model) {M0_tot:.2e} N·m (Mw {(2/3)*(np.log10(M0_tot)-9.1):.2f})")
    print(f"  {'thr':>4s} {'A_eff':>7s} {'L':>5s} {'W':>5s} {'D_bar':>6s} {'M0_eff/M0':>9s} "
          f"{'dS_Eshelby':>10s} {'dS_muD/W':>9s}")
    for thr in (0.10, 0.15, 0.20, 0.30):
        mm = slip >= thr * slip.max()
        A = mm.sum() * cell
        L = (mm.any(axis=0).sum()) * Dx
        W = (mm.any(axis=1).sum()) * Dy
        Dbar = slip[mm].mean()
        M0e = mu * slip[mm].sum() * cell
        R = np.sqrt(A / np.pi)
        ds_esh = 7.0 / 16.0 * M0e / R**3 / 1e6
        ds_w = 2.5 * mu * Dbar / (W * 1e3) / 1e6
        print(f"  {thr:4.2f} {A/1e6:7.0f} {L:5.0f} {W:5.0f} {Dbar:6.2f} {M0e/M0_tot:9.2f} "
              f"{ds_esh:10.1f} {ds_w:9.1f}")
