"""
fig_wisp_slipdist_5s.py — reproduce WISP's NATIVE SlipDist plot exactly (same slipcpt colormap,
rake-vector quivers, white hypocentre star, dashed grey rupture-time contours, strike/dip labels,
depth axis), but with rupture-front isochrones every 5 s instead of WISP's hardcoded 10 s.

This VENDORS WISP's PlotSlipDistribution (src/ffm/plot_graphic_NEIC.py, lines ~567-719): the only
changes are the contour `levels` (range(5,800,5)) and the caption. WISP source is imported / copied,
never edited in place. Helpers (__several_axes, __unpack_plane_data, slipcpt) are imported from WISP.

Run:  /home/msseo/miniforge3/envs/ff-env/bin/python scripts/fig_wisp_slipdist_5s.py <run_dir>
"""
import os, sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

WISP = os.path.expanduser("~/works/neic-finitefault/src")
sys.path.insert(0, WISP)
import ffm.plot_graphic_NEIC as G
import ffm.plane_management as pl_mng
from ffm.fault_plane import point_sources_param
from ffm.get_outputs import read_solution_static_format

slipcpt = G.slipcpt
_several_axes = getattr(G, "__several_axes")            # module-level, imported not mangled
_unpack = getattr(pl_mng, "__unpack_plane_data")


def plot(run_dir, out_name="SlipDist_plane0_5s.png", interval=5):
    directory = pathlib.Path(os.path.expanduser(run_dir))
    seg_data = json.load(open(directory / "segments_data.json"))
    segments = seg_data["segments"]
    rise_time = seg_data["rise_time"]
    connections = seg_data.get("connections")
    tensor_info = json.load(open(directory / "tensor_info.json"))
    point_sources = point_sources_param(segments, tensor_info, rise_time, connections=connections)
    solution = read_solution_static_format(segments, data_dir=directory)

    slip = solution["slip"]; rake = solution["rake"]; rupt_time = solution["rupture_time"]
    max_slip = max(np.max(s.flatten()) for s in slip)
    (directory / "plots").mkdir(exist_ok=True)
    outpath = None
    for i_seg, (segment, slip_seg, rake_seg, rupttime_seg, ps_seg) in enumerate(
            zip(segments, slip, rake, rupt_time, point_sources)):
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip
        plt.rc("axes", titlesize=20); plt.rc("axes", labelsize=20)
        plt.rc("xtick", labelsize=16); plt.rc("ytick", labelsize=16); plt.rc("font", size=20)
        fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111)
        ax.set_ylabel("Distance Along Dip (km)", fontsize=16)
        ax.set_xlabel("Distance Along Strike (km)", fontsize=16)
        ax.xaxis.tick_top(); ax.xaxis.set_label_position("top")
        for sp in ("bottom", "top", "left", "right"):
            ax.spines[sp].set_linewidth(3)
        fig.subplots_adjust(right=0.85, top=0.85, bottom=0.3)
        stk, dip_n, dstk, ddip, hstk, hdip = _unpack(segment)
        x = np.arange(stk) * dstk - hstk * dstk
        y = np.arange(dip_n) * ddip - hdip * ddip
        X = np.linspace(min(x), max(x), 20 * len(x))
        Y = np.linspace(min(y), max(y), 20 * len(y))
        z = (u ** 2 + v ** 2) ** 0.5
        xcols, yrows = np.meshgrid(x, y); XCOLS, YROWS = np.meshgrid(X, Y)
        orig = np.transpose(np.array([xcols.flatten(), yrows.flatten()]))
        new = np.transpose(np.array([XCOLS.flatten(), YROWS.flatten()]))
        grid_rt = griddata(orig, rupttime_seg.flatten(), new, method="linear").reshape(XCOLS.shape)
        cont = ax.contour(XCOLS, YROWS, grid_rt, colors="0.75", linestyles="dashed",
                          levels=range(interval, 800, interval), linewidths=1.0)   # <-- 5 s
        plt.clabel(cont, fmt="%.0f", inline=True, fontsize=14, colors="k")
        grid_z = griddata(orig, z.flatten(), new, method="linear").reshape(XCOLS.shape)
        ax.quiver(x, y, u, v, scale=30.0, width=0.002, color="0.5", clip_on=False)
        ax.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        ax, im = _several_axes(grid_z, segment, ps_seg, ax, max_val=1.0, autosize=False)
        cbar_ax = fig.add_axes((0.125, 0.15, 0.5, 0.07))
        sm = plt.cm.ScalarMappable(cmap=slipcpt, norm=plt.Normalize(vmin=0.0, vmax=max_slip / 100.0))
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.outline.set_linewidth(3); cb.set_label("Slip (m)", fontsize=18)
        ax.text(-0.1, 1.22, "Strike = " + str(int(segment["strike"])), fontsize=15,
                fontweight="bold", transform=ax.transAxes, va="top", ha="left")
        ax.text(-0.1, 1.13, "Dip = " + str(int(segment["dip"])), fontsize=15,
                fontweight="bold", transform=ax.transAxes, va="top", ha="left")
        ax.text(0, -0.04, f"Rupture Front Contours Plotted Every {interval} s", fontsize=15,
                fontweight="bold", transform=ax.transAxes, va="top", ha="left")
        name = out_name if i_seg == 0 else out_name.replace("plane0", f"plane{i_seg}")
        outpath = directory / "plots" / name
        plt.savefig(outpath, dpi=300); plt.close()
        print("->", outpath)
    return str(outpath)


if __name__ == "__main__":
    plot(sys.argv[1])
