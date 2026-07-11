"""Stress drop + source dimensions from the Calama WISP slip models.

Effective dimensions from a slip threshold (fraction of peak); stress drop from the
Eshelby circular crack, dSigma = (7/16) M0_eff / R^3 with R = sqrt(A_eff/pi), plus the
width-limited estimate c*mu*Dbar/W. Rigidity from WISP's shear_model.txt.
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig_fault_isochrone as fi

BASE = os.path.join(os.path.dirname(HERE), "event_2026_calama/20260525215220")
RUNS = {"shallow 194/16 (final, NP1_ref3)": "ffm.1/NP1_ref3",
        "steep 7/74 (NP2_ref2)": "ffm.1/NP2_ref2"}

for name, rel in RUNS.items():
    d = os.path.join(BASE, rel)
    s = fi.parse_solution(os.path.join(d, "Solution.txt"))
    slip, Dx, Dy = s["slip"], s["Dx"], s["Dy"]          # m, km, km
    mu = float(open(os.path.join(d, "shear_model.txt")).read().split("\n")[2].split()[0]) / 10.  # dyn/cm2 -> Pa
    cell = Dx * Dy * 1e6                                 # m^2
    M0_tot = mu * slip.sum() * cell
    print(f"\n== {name} | mu {mu/1e9:.1f} GPa | M0(model) {M0_tot:.2e} N·m "
          f"(Mw {(2/3)*(np.log10(M0_tot)-9.1):.2f})")
    print(f"{'thr':>4s} {'A_eff':>7s} {'L':>5s} {'W':>5s} {'D_bar':>6s} {'M0_eff/M0':>9s} "
          f"{'dS_Eshelby':>10s} {'dS_muD/W':>9s}")
    for thr in (0.10, 0.15, 0.20, 0.30):
        m = slip >= thr * slip.max()
        A = m.sum() * cell                               # m^2
        L = (m.any(axis=0).sum()) * Dx                   # km along strike
        W = (m.any(axis=1).sum()) * Dy                   # km along dip
        Dbar = slip[m].mean()
        M0e = mu * slip[m].sum() * cell
        R = np.sqrt(A / np.pi)
        ds_esh = 7.0 / 16.0 * M0e / R**3 / 1e6           # MPa
        ds_w = 2.5 * mu * Dbar / (W * 1e3) / 1e6         # MPa (c~2.5 buried dip-slip)
        print(f"{thr:4.2f} {A/1e6:7.0f} {L:5.0f} {W:5.0f} {Dbar:6.2f} {M0e/M0_tot:9.2f} "
              f"{ds_esh:10.1f} {ds_w:9.1f}")
