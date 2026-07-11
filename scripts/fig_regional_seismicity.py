"""Regional intermediate-depth seismicity context for the 2026 Calama event.
Three panels: (a) map of USGS M>=4.5 events (depth-coloured), historical M>=6.8 labelled;
(b) trench-perpendicular cross-section (21.5-24.5 S band) showing the slab and the
100-130 km nest; (c) magnitude-time history of M>=6 intermediate-depth events."""
import os
import pandas as pd
import pygmt

B = os.path.expanduser("~/works/17.Venezuela_2026")
ELAT, ELON, EDEP = -22.3667, -68.6013, 109.0
cat = pd.read_csv(os.path.join(B, "data/catalogs/nchile_intraslab_M45.csv"))
cat = cat[cat.depth >= 60]
big = pd.read_csv(os.path.join(B, "data/catalogs/nchile_intraslab_M6.csv"))
big = big[big.mag >= 6.8]
big["year"] = big.time.str[:4]

fig = pygmt.Figure()
pygmt.config(FONT="10p,Helvetica", FONT_TITLE="12p,Helvetica",
             MAP_FRAME_TYPE="plain", FORMAT_GEO_MAP="D")
region = [-71.5, -66.0, -26.0, -19.0]
pygmt.makecpt(cmap="roma", series=[60, 300], reverse=False)

# (a) map
fig.basemap(region=region, projection="M13c",
            frame=["af", "+tIntermediate-depth seismicity (USGS, M>=4.5, z>=60 km)"])
fig.grdimage("@earth_relief_02m", cmap="gray", shading="+a315+ne0.4", transparency=55)
fig.coast(shorelines="0.4p,black", borders="1/0.3p,gray40")
fig.plot(x=[-71.55, -71.45, -71.35, -71.30, -71.25], y=[-19.5, -21.0, -22.5, -24.0, -25.5],
         pen="1.4p,black", style="f1c/0.22c+r+t", fill="black")
pygmt.makecpt(cmap="roma", series=[60, 300])
fig.plot(x=cat.longitude, y=cat.latitude, size=0.03 * (1.6 ** cat.mag).clip(0, 60),
         style="cc", fill=cat.depth, cmap=True, pen="0.1p,gray30", transparency=25)
fig.plot(x=big.longitude, y=big.latitude, size=0.14 * (big.mag - 5.4),
         style="ac", fill=big.depth, cmap=True, pen="1p,black")
lab = big[big.mag >= 7.0]
fig.text(x=lab.longitude, y=lab.latitude + 0.22,
         text=[f"{y} M{m:.1f}" for y, m in zip(lab.year, lab.mag)],
         font="8p,Helvetica-Bold,black", fill="white@30")
fig.plot(x=[ELON], y=[ELAT], style="a0.7c", fill="yellow", pen="1.2p,black")
fig.text(x=ELON - 0.15, y=ELAT + 0.25, text="2026 M6.9", font="9p,Helvetica-Bold,black",
         fill="white@20", justify="MR")
fig.plot(x=[-68.93], y=[-22.46], style="s0.22c", fill="black")
fig.text(x=-68.93, y=-22.66, text="Calama", font="8p,Helvetica", fill="white@30")
# cross-section band
fig.plot(x=[-71.5, -66.0, -66.0, -71.5, -71.5], y=[-24.5, -24.5, -21.5, -21.5, -24.5],
         pen="1p,blue,--")
pygmt.makecpt(cmap="roma", series=[60, 300])
fig.colorbar(frame="af+lHypocentre depth (km)", position="JMR+o0.4c/0c+w8c")

# (b) cross-section: distance east of trench (-71.35) vs depth, band 21.5-24.5 S
fig.shift_origin(yshift="-9.0c")
band = cat[(cat.latitude <= -21.5) & (cat.latitude >= -24.5)]
xb = (band.longitude - (-71.35)) * 111.19 * 0.92   # km, cos(22.5 deg) ~ 0.92
bb = big[(big.latitude <= -21.5) & (big.latitude >= -24.5)]
xbb = (bb.longitude - (-71.35)) * 111.19 * 0.92
fig.basemap(region=[0, 520, 0, 320], projection="X13c/-6c",
            frame=["xaf+lDistance east of trench (km)", "yaf+lDepth (km)",
                   "+tCross-section 21.5-24.5 S (dashed box)"])
pygmt.makecpt(cmap="roma", series=[60, 300])
fig.plot(x=xb, y=band.depth, size=0.03 * (1.5 ** band.mag).clip(0, 40), style="cc",
         fill=band.depth, cmap=True, pen="0.1p,gray30", transparency=25)
fig.plot(x=xbb, y=bb.depth, size=0.14 * (bb.mag - 5.4), style="ac",
         fill=bb.depth, cmap=True, pen="1p,black")
bl = bb[bb.mag >= 7.0].reset_index()
off = [-16 if i % 2 == 0 else 22 for i in range(len(bl))]
fig.text(x=(bl.longitude - (-71.35)) * 111.19 * 0.92, y=bl.depth + pd.Series(off),
         text=[f"{y} M{m:.1f}" for y, m in zip(bl.year, bl.mag)],
         font="8p,Helvetica-Bold,black", fill="white@30")
fig.plot(x=[(ELON + 71.35) * 111.19 * 0.92], y=[EDEP], style="a0.7c", fill="yellow",
         pen="1.2p,black")

# (c) magnitude-time of M>=6, z>=60 km
fig.shift_origin(yshift="-6.6c")
allm6 = pd.read_csv(os.path.join(B, "data/catalogs/nchile_intraslab_M6.csv"))
yr = allm6.time.str[:4].astype(int) + 0.5
pygmt.makecpt(cmap="roma", series=[60, 300])
fig.basemap(region=[1900, 2031, 5.8, 8.6], projection="X13c/4c",
            frame=["xaf+lYear", "ya0.5f+lMagnitude", "+tM>=6 intermediate-depth history (region above)"])
fig.plot(x=yr, y=allm6.mag, style="c0.16c", fill=allm6.depth, cmap=True, pen="0.2p,black")
fig.plot(x=[2026.4], y=[6.9], style="a0.5c", fill="yellow", pen="1p,black")
for _, r in allm6[allm6.mag >= 7.3].iterrows():
    fig.text(x=int(r.time[:4]) + 0.5, y=r.mag + 0.17, text=f"{r.time[:4]} M{r.mag:.1f}",
             font="7p,Helvetica,black")

out = os.path.join(B, "event_2026_calama/report/fig_regional_seismicity.png")
fig.savefig(out, dpi=200)
print("->", out, "| map events:", len(cat), "| band:", len(band), "| M6+:", len(allm6))
