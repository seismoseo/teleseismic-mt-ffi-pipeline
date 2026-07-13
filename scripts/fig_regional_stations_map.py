"""Map of the regional strong-motion stations used in the WISP joint inversions, showing why
the 2026 M6.9 (10 stations) got more near-field data than the 2024 M7.4 (3 stations): the 2026
epicentre sits closer to the dense CX/IPOC forearc network, so more stations fall inside WISP's
2-degree near-field selection radius (dashed circles)."""
import os
import pandas as pd
import pygmt

B = os.path.expanduser("~/works/17.Venezuela_2026")
E24 = (-23.0791, -67.8404); E26 = (-22.3667, -68.6013)
st = pd.read_csv(os.path.join(B, "data/catalogs/regional_stations_map.csv"))

fig = pygmt.Figure()
pygmt.config(FONT="11p,Helvetica", FONT_TITLE="13p,Helvetica", MAP_FRAME_TYPE="plain",
             FORMAT_GEO_MAP="D")
region = [-71.2, -66.5, -25.0, -20.3]
fig.basemap(region=region, projection="M15c",
            frame=["af", "+tRegional strong-motion stations used in the joint inversions"])
fig.grdimage("@earth_relief_02m", cmap="gray", shading="+a315+ne0.4", transparency=45)
fig.coast(shorelines="0.4p,black", borders="1/0.3p,gray40")

# 2-degree WISP near-field selection radius around each epicentre
for (la, lo), col in ((E24, "blue"), (E26, "red")):
    fig.plot(x=[lo], y=[la], style="E-444.8", pen=f"1p,{col},--")   # 2 deg ~ 222.4 km radius -> diameter 444.8

# all acquired stations (small grey), then highlight used
fig.plot(x=st.longitude if "longitude" in st else st.lon, y=st.lat,
         style="c0.12c", fill="gray70", pen="0.2p,gray40")
u26 = st[st.used26 == 1]; u24 = st[st.used24 == 1]
fig.plot(x=u26.lon, y=u26.lat, style="t0.42c", fill="red", pen="0.6p,black")
fig.plot(x=u24.lon, y=u24.lat, style="i0.42c", fill="blue", pen="0.6p,black")
# stations used by BOTH get a small marker note is implicit (AF01/PB09/PB15 carry both symbols)
fig.text(x=st.lon, y=st.lat + 0.09, text=st.sta, font="6p,Helvetica,black", fill="white@40")

# epicentres
fig.plot(x=[E26[1]], y=[E26[0]], style="a0.8c", fill="red", pen="1.2p,black")
fig.plot(x=[E24[1]], y=[E24[0]], style="a0.8c", fill="blue", pen="1.2p,black")
fig.text(x=E26[1] - 0.15, y=E26[0] + 0.22, text="2026 M6.9", font="11p,Helvetica-Bold,red",
         justify="MR", fill="white@30")
fig.text(x=E24[1] + 0.15, y=E24[0] - 0.22, text="2024 M7.4", font="11p,Helvetica-Bold,blue",
         justify="ML", fill="white@30")

# legend
fig.plot(x=[-70.9], y=[-24.35], style="t0.42c", fill="red", pen="0.6p,black")
fig.text(x=-70.75, y=-24.35, text="used by 2026 (10 stations)", font="9p,Helvetica", justify="ML")
fig.plot(x=[-70.9], y=[-24.6], style="i0.42c", fill="blue", pen="0.6p,black")
fig.text(x=-70.75, y=-24.6, text="used by 2024 (3 stations)", font="9p,Helvetica", justify="ML")
fig.plot(x=[-70.9], y=[-24.85], style="c0.12c", fill="gray70", pen="0.2p,gray40")
fig.text(x=-70.75, y=-24.85, text="acquired, not selected (>2 deg)", font="9p,Helvetica", justify="ML")

fig.basemap(map_scale="jBR+w100k+o0.6c/0.6c+f+l")
out = os.path.join(B, "event_2026_calama/report/fig_regional_stations_map.png")
fig.savefig(out, dpi=200)
print("->", out, "| stations:", len(st), "| used26:", len(u26), "used24:", len(u24))
