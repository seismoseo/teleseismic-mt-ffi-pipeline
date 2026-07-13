"""Teleseismic station map for the 2023 Morocco Mw 6.8 study: the 20 GSN broadband stations
used in the MTUQ moment-tensor inversion, on an azimuthal-equidistant projection centred on
the epicentre so distance is radial and azimuthal coverage is faithful."""
import os
import numpy as np
import pandas as pd
import pygmt

B = os.path.expanduser("~/works/17.Venezuela_2026")
ELAT, ELON = 31.058, -8.3847
st = pd.read_csv(os.path.join(B, "data/catalogs/morocco_mtuq_stations.csv"))

fig = pygmt.Figure()
pygmt.config(FONT="10p,Helvetica", FONT_TITLE="13p,Helvetica", MAP_FRAME_TYPE="plain")
proj = f"E{ELON}/{ELAT}/150/15c"      # azimuthal equidistant, 150 deg radius
fig.coast(region="g", projection=proj, land="gray85", water="white",
          shorelines="0.2p,gray50", frame="g30")
# distance rings (30, 60, 90 deg)
for d in (30, 60, 90):
    fig.plot(x=[ELON], y=[ELAT], style=f"E-{2*d*111.19*2}", pen="0.5p,gray60,--", no_clip=True)
# great-circle paths epicentre -> station
for _, r in st.iterrows():
    fig.plot(x=[ELON, r.lon], y=[ELAT, r.lat], pen="0.4p,steelblue@50")
fig.plot(x=st.lon, y=st.lat, style="t0.34c", fill="steelblue", pen="0.5p,black")
fig.text(x=st.lon, y=st.lat, text=st.sta.str.split(".").str[-1], font="5p,Helvetica,black",
         justify="LM", offset="0.15c/0c")
fig.plot(x=[ELON], y=[ELAT], style="a0.6c", fill="red", pen="1p,black")
fig.text(x=ELON, y=ELAT, text="Morocco M6.8", font="10p,Helvetica-Bold,red",
         justify="LM", offset="0.3c/0.3c")
fig.text(position="TC", text="MTUQ teleseismic stations (20; azimuthal-equidistant, rings 30/60/90 deg)",
         font="12p,Helvetica", offset="0/0.4c", no_clip=True)
out = os.path.join(B, "event_2023_morocco/fig_morocco_stations.png")
fig.savefig(out, dpi=200)
print("->", out, "|", len(st), "stations")
