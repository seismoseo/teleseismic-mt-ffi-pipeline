"""
usgs_compare.py — Compare a WISP finite-fault model (run auto-labelled) with the
published USGS NEIC model in /home/msseo/works/17.Venezuela_2026/USGS.

Parses:
  - USGS complete_inversion.fsp  (multi-segment slip model, SRCMOD/FSP format)
  - USGS moment_rate.mr          (moment-rate function, dt + 2-column table)
  - my Solution.txt              (WISP slip model: lat lon dep slip rake strike dip trup tris tfal mo)
  - my STF.txt                   (WISP moment-rate function, same 2-column format as .mr)

Figures: moment-rate overlay, map-view slip comparison, along-strike slip profile.
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _helvetica():
    from matplotlib import font_manager
    for n in ("Helvetica", "Arial", "Nimbus Sans"):
        if any(n in f.name for f in font_manager.fontManager.ttflist):
            matplotlib.rcParams["font.family"] = n
            break


# ---------- parsers ----------
def parse_mr(path):
    """Return (time[s], moment_rate[Nm/s]) from a USGS .mr or WISP STF.txt file."""
    t, mr = [], []
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            try:
                t.append(float(p[0])); mr.append(float(p[1]))
            except ValueError:
                continue
    return np.array(t), np.array(mr)


def parse_fsp(path):
    """Return dict of arrays for the USGS FSP slip model (all segments concatenated)."""
    lat, lon, xs, ys, z, slip, rake, trup = [], [], [], [], [], [], [], []
    meta = {}
    for line in open(path):
        if line.startswith("%"):
            m = re.search(r"Mw = ([\d.]+)\s+Mo = ([\d.eE+]+)", line)
            if m:
                meta["Mw"] = float(m.group(1)); meta["Mo"] = float(m.group(2))
            m = re.search(r"avTr =\s*([\d.]+) s\s+avVr =\s*([\d.]+)", line)
            if m:
                meta["avTr"] = float(m.group(1)); meta["avVr"] = float(m.group(2))
            continue
        p = line.split()
        if len(p) >= 10:
            try:
                vals = [float(x) for x in p[:10]]
            except ValueError:
                continue
            lat.append(vals[0]); lon.append(vals[1]); xs.append(vals[2]); ys.append(vals[3])
            z.append(vals[4]); slip.append(vals[5]); rake.append(vals[6]); trup.append(vals[7])
    return dict(lat=np.array(lat), lon=np.array(lon), x=np.array(xs), y=np.array(ys),
                z=np.array(z), slip=np.array(slip), rake=np.array(rake),
                trup=np.array(trup), **meta)


def parse_solution(path):
    """Return dict of arrays for the WISP Solution.txt slip model."""
    txt = open(path).read().splitlines()
    i0 = next(i for i, l in enumerate(txt) if "slip" in l and "rake" in l) + 1
    a = np.array([[float(x) for x in l.split()[:11]] for l in txt[i0:] if len(l.split()) >= 11])
    return dict(lat=a[:, 0], lon=a[:, 1], depth=a[:, 2], slip=a[:, 3] / 100.0,
                rake=a[:, 4], trup=a[:, 7])


EPI_LON, EPI_LAT = -68.5277, 10.436   # 2026 Venezuela doublet default
STRIKE_DEG = 90.0                     # along-strike profile axis (90 = east)


def _ew_km(lon, lat):
    """Signed along-strike distance (km) from the epicentre."""
    east = (lon - EPI_LON) * 111.32 * np.cos(np.radians(lat))
    north = (lat - EPI_LAT) * 111.32
    s = np.radians(STRIKE_DEG)
    return east * np.sin(s) + north * np.cos(s)


def _runtag(mine_run):
    """Derive the run label (e.g. 'ffm.6') from the NP directory path."""
    return os.path.basename(os.path.dirname(os.path.normpath(mine_run)))


# ---------- figures ----------
def plot_stf(mine_run, usgs_dir, out=None, label=None):
    _helvetica()
    tag = label or _runtag(mine_run)
    out = out or os.path.join(mine_run, "plots", "Compare_STF.png")
    tm, mm = parse_mr(os.path.join(mine_run, "STF.txt"))
    tu, mu = parse_mr(os.path.join(usgs_dir, "moment_rate.mr"))
    Mo_m = np.trapz(mm, tm); Mo_u = np.trapz(mu, tu)
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.fill_between(tu, mu / 1e19, color="0.75", label=f"USGS NEIC ($M_0$={Mo_u:.2e} Nm)")
    ax.plot(tm, mm / 1e19, color="#c0392b", lw=1.8,
            label=f"This study, {tag} ($M_0$={Mo_m:.2e} Nm)")
    ax.set_xlim(0, max(tu.max(), tm.max())); ax.set_ylim(bottom=0)
    ax.set_xlabel("Time after origin (s)")
    ax.set_ylabel(r"Moment rate ($\times10^{19}$ Nm/s)")
    ax.set_title("Moment-rate function — this study vs USGS", fontsize=11)
    ax.legend(framealpha=1, edgecolor="0.5")
    fig.tight_layout(); fig.savefig(out, dpi=190, bbox_inches="tight"); plt.close(fig)
    print(f"STF: mine dur~{tm[mm>0.02*mm.max()].max():.0f}s peak@{tm[mm.argmax()]:.0f}s | "
          f"USGS dur~{tu[mu>0.02*mu.max()].max():.0f}s peak@{tu[mu.argmax()]:.0f}s -> {out}")
    return out


def plot_slip_map(mine_run, usgs_dir, out=None, label=None):
    _helvetica()
    tag = label or _runtag(mine_run)
    out = out or os.path.join(mine_run, "plots", "Compare_slip_map.png")
    me = parse_solution(os.path.join(mine_run, "Solution.txt"))
    us = parse_fsp(os.path.join(usgs_dir, "complete_inversion.fsp"))
    vmax = max(me["slip"].max(), us["slip"].max())
    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, sharey=True)
    for ax, d, ttl in [(axes[0], me, f"This study ({tag})"),
                       (axes[1], us, "USGS NEIC")]:
        sc = ax.scatter(d["lon"], d["lat"], c=d["slip"], cmap="magma_r", vmin=0, vmax=vmax,
                        s=90, marker="s", edgecolor="none")
        ax.plot(EPI_LON, EPI_LAT, marker="*", ms=16, mfc="yellow", mec="black", mew=1, zorder=5)
        ax.set_title(ttl, fontsize=10); ax.set_ylabel("Latitude (°)")
        ax.set_aspect(1.0 / np.cos(np.radians(EPI_LAT)))
    axes[1].set_xlabel("Longitude (°)")
    cb = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.02); cb.set_label("Slip (m)")
    fig.suptitle("Map-view slip distribution — this study vs USGS", fontsize=12, y=0.98)
    fig.savefig(out, dpi=190, bbox_inches="tight"); plt.close(fig)
    print(f"slip map -> {out}")
    return out


def plot_slip_profile(mine_run, usgs_dir, out=None, label=None):
    """Dip-summed slip vs along-strike (E-W) distance from epicenter, both models."""
    _helvetica()
    tag = label or _runtag(mine_run)
    out = out or os.path.join(mine_run, "plots", "Compare_slip_profile.png")
    me = parse_solution(os.path.join(mine_run, "Solution.txt"))
    us = parse_fsp(os.path.join(usgs_dir, "complete_inversion.fsp"))
    lo = min(_ew_km(d["lon"], d["lat"]).min() for d in (me, us))
    hi = max(_ew_km(d["lon"], d["lat"]).max() for d in (me, us))
    step = max(2.0, round((hi - lo) / 16.0))
    bins = np.arange(np.floor(lo), np.ceil(hi) + step, step)
    ctr = 0.5 * (bins[:-1] + bins[1:])

    def profile(d):
        x = _ew_km(d["lon"], d["lat"])
        idx = np.digitize(x, bins) - 1
        prof = np.array([d["slip"][idx == i].mean() if np.any(idx == i) else 0.0
                         for i in range(len(ctr))])
        return prof

    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.fill_between(ctr, profile(us), step="mid", color="0.75", label="USGS NEIC")
    ax.plot(ctr, profile(me), drawstyle="steps-mid", color="#c0392b", lw=1.8, label=f"This study ({tag})")
    ax.axvline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel(f"Along-strike distance from epicentre, azimuth {STRIKE_DEG:.0f}° (km)")
    ax.set_ylabel("Mean slip over down-dip (m)")
    ax.set_title("Along-strike slip profile — this study vs USGS", fontsize=11)
    ax.legend(framealpha=1, edgecolor="0.5")
    fig.tight_layout(); fig.savefig(out, dpi=190, bbox_inches="tight"); plt.close(fig)
    print(f"slip profile -> {out}")
    return out


if __name__ == "__main__":
    import sys
    run, usgs = sys.argv[1], sys.argv[2]
    if len(sys.argv) > 3:   # optional: epi_lon epi_lat [strike_deg]
        EPI_LON, EPI_LAT = float(sys.argv[3]), float(sys.argv[4])
        if len(sys.argv) > 5:
            STRIKE_DEG = float(sys.argv[5])
    plot_stf(run, usgs); plot_slip_map(run, usgs); plot_slip_profile(run, usgs)
