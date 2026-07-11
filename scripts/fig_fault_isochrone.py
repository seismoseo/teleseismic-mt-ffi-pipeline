"""
fig_fault_isochrone.py — On-fault slip + rupture-time isochrone map at TRUE 1:1 scale.

WISP's default SlipDist plot vertically exaggerates the fault (along-dip stretched vs
along-strike). This parses Solution.txt and plots the fault plane with along-strike and
along-dip in km at equal aspect (no vertical exaggeration), slip as colour, rupture-time
isochrones as contours, hypocentre at the origin, and a secondary depth axis.

Usage: python scripts/fig_fault_isochrone.py <NP_dir> [out.png]
"""

import os
import sys
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _helvetica():
    from matplotlib import font_manager
    for n in ("Helvetica", "Arial", "Nimbus Sans"):
        if any(n in f.name for f in font_manager.fontManager.ttflist):
            matplotlib.rcParams["font.family"] = n; break


def parse_solution(path):
    """Return dict with grids: along_strike(km), along_dip(km), slip(m), t_rup(s),
    depth(km), and geometry (Dx, Dy, nx, ny, hyp_stk, hyp_dip, dip)."""
    txt = open(path).read().splitlines()
    # geometry from header line 2
    m = re.search(r"nx\(Along-strike\)=\s*(\d+)\s*Dx\s*=\s*([\d.]+)km\s*ny\(downdip\)=\s*(\d+)\s*Dy\s*=\s*([\d.]+)km", txt[1])
    nx, Dx, ny, Dy = int(m.group(1)), float(m.group(2)), int(m.group(3)), float(m.group(4))
    mh = re.search(r"EQ in cell\s*\(\s*(\d+),\s*(\d+)\)", txt[2])
    hyp_stk, hyp_dip = int(mh.group(1)), int(mh.group(2))
    # data rows: after the line containing 'slip rake strike'
    i0 = next(i for i, l in enumerate(txt) if "slip" in l and "rake" in l) + 1
    rows = []
    for l in txt[i0:]:
        p = l.split()
        if len(p) >= 11:
            rows.append([float(x) for x in p[:11]])
    a = np.array(rows)  # lat lon depth slip rake strike dip t_rup t_ris t_fal mo
    # inner loop = along-strike (nx), outer = down-dip (ny) -> reshape (ny, nx)
    depth = a[:, 2].reshape(ny, nx)
    slip = (a[:, 3] / 100.0).reshape(ny, nx)      # cm -> m
    rake = a[:, 4].reshape(ny, nx)                # deg (Aki & Richards)
    trup = a[:, 7].reshape(ny, nx)
    tris = (a[:, 8] + a[:, 9]).reshape(ny, nx)    # rise = t_ris(rise) + t_fal(fall) = total slip duration
    dip = float(a[0, 6])
    strike = float(a[0, 5])
    # on-fault coordinates centred on the hypocentre. WISP's hyp_stk/hyp_dip
    # ("EQ in cell (4,4)") are 1-BASED, so convert to 0-based to place the origin
    # exactly on the nucleation cell (where t_rup = 0).
    js = (np.arange(nx) - (hyp_stk - 1)) * Dx     # along-strike (km), + = east
    ids = (np.arange(ny) - (hyp_dip - 1)) * Dy    # along-dip (km), + = down-dip
    X, Y = np.meshgrid(js, ids)
    return dict(X=X, Y=Y, slip=slip, rake=rake, trup=trup, tris=tris, depth=depth, Dx=Dx, Dy=Dy,
                nx=nx, ny=ny, hyp_stk=hyp_stk, hyp_dip=hyp_dip, dip=dip, strike=strike,
                ids=ids, js=js)


def _strike_label(strike):
    """Along-strike axis label with the +X compass direction derived from the strike."""
    dirs = ["north", "NNE", "NE", "ENE", "east", "ESE", "SE", "SSE",
            "south", "SSW", "SW", "WSW", "west", "WNW", "NW", "NNW"]
    d = dirs[int(round(strike / 22.5)) % 16]
    return f"Distance along strike (km)  — {d} →"


# fixed on-fault along-dip window so the depth axis is IDENTICAL across models with
# different dip_subfaults (the hypocentre row is at 20.3 km depth in all of them, and the
# depth-per-along-dip slope = sin(dip) is the same, so a fixed along-dip range -> fixed depth
# range). -22..32 km along-dip spans both the 10-row and 11-row widened geometries.
DIP_LO, DIP_HI = -22.0, 32.0


def plot(np_dir, out=None, title=None, dip_ylim=None):
    """dip_ylim: None -> fixed 2026 range (DIP_LO,DIP_HI) for cross-notebook consistency;
    'auto' -> from this model's along-dip span (for other events); or an explicit (lo,hi).
    title: None -> auto ('... — Mw X.XX') from the solution moment."""
    _helvetica()
    sol = parse_solution(os.path.join(np_dir, "Solution.txt"))
    out = out or os.path.join(np_dir, "plots", "SlipIsochrone_trueScale.png")
    X, Y, slip, trup, rake = sol["X"], sol["Y"], sol["slip"], sol["trup"], sol["rake"]
    if dip_ylim is None:
        dlo, dhi = DIP_LO, DIP_HI
    elif dip_ylim == "auto":
        dlo, dhi = sol["ids"].min() - sol["Dy"], sol["ids"].max() + sol["Dy"]
    else:
        dlo, dhi = dip_ylim

    # figure width:height set to the true km aspect (so set_aspect('equal') is faithful)
    span_x = sol["js"].ptp() + sol["Dx"]
    span_y = dhi - dlo
    fig_w = 13.0
    fig_h = max(2.2, fig_w * span_y / span_x + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    pcm = ax.pcolormesh(X, Y, slip, cmap="magma_r", shading="gouraud", vmin=0)
    # slip vectors: rake direction (Aki & Richards) in fault-plane coords, length ∝ slip.
    # +X = along strike (east); up-dip = +sin(rake) -> in down-dip-positive Y that is -sin(rake).
    U = np.cos(np.radians(rake)); V = -np.sin(np.radians(rake))
    sc = slip / max(slip.max(), 1e-9)
    sub = np.zeros_like(slip, bool); sub[:, ::2] = True          # every 2nd along-strike column
    mq = (slip > 0.15 * slip.max()) & sub
    ax.quiver(X[mq], Y[mq], (U * sc)[mq], (V * sc)[mq], color="0.2",
              angles="xy", scale_units="xy", scale=0.09, width=0.0022,
              headwidth=4, headlength=5, alpha=0.85, zorder=4)
    # rupture-front isochrones every 5 s, WHITE. A contour AT the max time can't be
    # drawn (nothing beyond it), so add a near-terminus level so the front reaches the
    # eastern slip edge.
    tend = float(trup[slip > 0.05].max())
    step = 5.0 if tend > 6 else 1.0                          # 5 s isochrones (1 s for very short ruptures)
    levels = list(np.arange(step, tend, step))
    if levels and tend - levels[-1] > 0.4 * step:
        levels.append(round(tend - 0.3 * step, 1))
    levels = sorted(set(round(x, 2) for x in levels if 0 < x < tend))
    if len(levels) >= 1:
        cs = ax.contour(X, Y, trup, levels=levels, colors="white", linewidths=0.9)
        ax.clabel(cs, fmt="%.0f s", fontsize=7, inline_spacing=2, colors="white")
    ax.plot(0, 0, marker="*", ms=18, mfc="yellow", mec="black", mew=1.2, zorder=5)

    # slip-weighted rupture velocity (distance from hypocentre vs rupture time)
    dist = np.sqrt(X**2 + Y**2)
    m = slip > 0.2
    Vr = float(np.sum(slip[m] * dist[m] * trup[m]) / np.sum(slip[m] * trup[m]**2))
    # caption placed BELOW the axes so it never overlaps the slip map
    fig.text(0.5, 0.005,
             f"White contours = rupture-front arrival, every {step:.0f} s (front spans 0–{tend:.0f} s);  "
             f"V$_r$ ≈ {Vr:.1f} km/s ({'subshear' if Vr < 3.5 else 'supershear'})",
             ha="center", fontsize=9,
             bbox=dict(fc="white", ec="0.5", alpha=0.92))

    ax.set_aspect("equal")                         # TRUE 1:1, no vertical exaggeration
    ax.set_xlabel(_strike_label(sol["strike"]))
    ax.set_ylabel("Distance along dip (km)")
    ax.set_ylim(dhi, dlo)                           # down-dip downward
    if title is None:
        i0 = next(i for i, l in enumerate(open(os.path.join(np_dir, "Solution.txt")))
                  if "slip" in l and "rake" in l)
        a = np.array([[float(x) for x in l.split()[:11]]
                      for l in open(os.path.join(np_dir, "Solution.txt")).read().splitlines()[i0 + 1:]
                      if len(l.split()) >= 11])
        mw = (2 / 3) * (np.log10(a[:, 10].sum()) - 16.1)
        title = (f"On-fault slip, slip vectors, and rupture-time isochrones "
                 f"(true 1:1 scale) — Mw {mw:.2f}")
    ax.set_title(title, fontsize=11)

    # secondary depth axis (depth increases down-dip)
    dmin, dmax = sol["depth"].min(), sol["depth"].max()
    secax = ax.secondary_yaxis(
        "right",
        functions=(lambda yd: dmin + (yd - sol["ids"].min()) / (sol["ids"].ptp()) * (dmax - dmin),
                   lambda d: sol["ids"].min() + (d - dmin) / (dmax - dmin) * sol["ids"].ptp()))
    secax.set_ylabel("Depth (km)")

    cb = fig.colorbar(pcm, ax=ax, orientation="vertical", pad=0.10, shrink=0.9)
    cb.set_label("Slip (m)")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"max slip {slip.max():.2f} m | along-strike {span_x:.0f} km x along-dip {span_y:.0f} km")
    print(f"-> {out}")
    return out


def plot_risetime(np_dir, out=None):
    """Map of per-subfault slip duration (rise + fall = source-time-function width),
    true 1:1 scale, slip>5cm contour overlaid. Rise time IS solved by the inversion
    (within the windowed slip-rate basis: min_rise 1 s, delta_rise 2 s, 8 windows)."""
    _helvetica()
    sol = parse_solution(os.path.join(np_dir, "Solution.txt"))
    out = out or os.path.join(np_dir, "plots", "RiseTime_trueScale.png")
    X, Y, slip, tris = sol["X"], sol["Y"], sol["slip"], sol["tris"]

    span_x = sol["js"].ptp() + sol["Dx"]
    span_y = sol["ids"].ptp() + sol["Dy"]
    fig_w = 13.0
    fig_h = max(2.2, fig_w * span_y / span_x + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    tr = np.ma.masked_where(slip < 0.05, tris)     # only show where there is slip
    pcm = ax.pcolormesh(X, Y, tr, cmap="viridis", shading="gouraud")
    cs = ax.contour(X, Y, slip, levels=[0.3, 0.8, 1.5], colors="white", linewidths=0.7)
    ax.clabel(cs, fmt="%.1f m", fontsize=6, colors="white")
    ax.plot(0, 0, marker="*", ms=18, mfc="yellow", mec="black", mew=1.2, zorder=5)

    med = float(np.median(tris[slip > 0.2]))
    fig.text(0.5, 0.005,
             f"Per-subfault slip duration (rise+fall); white = slip contours; "
             f"median ≈ {med:.1f} s where slip > 0.2 m",
             ha="center", fontsize=9, bbox=dict(fc="white", ec="0.5", alpha=0.92))

    ax.set_aspect("equal")
    ax.set_xlabel(_strike_label(sol["strike"]))
    ax.set_ylabel("Distance along dip (km)")
    ax.invert_yaxis()
    ax.set_title("Per-subfault rise time (data-driven, true 1:1 scale) — "
                 "Mw 7.57 combined model", fontsize=11)
    cb = fig.colorbar(pcm, ax=ax, orientation="vertical", pad=0.10, shrink=0.9)
    cb.set_label("Slip duration (s)")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"rise time: median {med:.1f}s | range {tris[slip>0.2].min():.1f}-{tris[slip>0.2].max():.1f}s")
    print(f"-> {out}")
    return out


if __name__ == "__main__":
    which = sys.argv[2] if len(sys.argv) > 2 else "slip"
    if which == "rise":
        plot_risetime(sys.argv[1])
    else:
        plot(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 and which != "rise" else None)
