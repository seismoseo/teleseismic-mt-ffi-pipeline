"""
wisp_fetch_fallback.py — obspy-based replacement for `ffm get-data teleseismic` when the
CLI fails on IRIS http:// service discovery. Writes WISP's expected layout into CWD:
  <NET>_<STA>_<CHA>_<LOC>.sac   (raw counts)
  SAC_PZs_<NET>_<STA>_<CHA>_<LOC>  (SACPZ response incl. LATITUDE/LONGITUDE/AZIMUTH/DIP)
Then run wisp_fix_headers.py + `ffm model run . auto_model ...` as usual.

Usage: python wisp_fetch_fallback.py <time> <lat> <lon> [mindeg=30] [maxdeg=90]
"""
import sys, os
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import locations2degrees

TIME, LAT, LON = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
MINDEG = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0
MAXDEG = float(sys.argv[5]) if len(sys.argv) > 5 else 90.0
O = UTCDateTime(TIME)
cl = Client("IRIS")
print("station inventory...", flush=True)
inv = cl.get_stations(network="IU,II,G,GT,IC,CU,GE", channel="BH?", level="response",
                      starttime=O - 120, endtime=O + 3600,
                      latitude=LAT, longitude=LON, minradius=MINDEG, maxradius=MAXDEG)
n = 0
for net in inv:
    for sta in net:
        deg = locations2degrees(LAT, LON, sta.latitude, sta.longitude)
        try:
            st = cl.get_waveforms(net.code, sta.code, "*", "BH?", O - 120, O + 3600)
        except Exception:
            continue
        locs = {}
        for tr in st:
            locs.setdefault(tr.stats.location, set()).add(tr.stats.channel)
        loc = next((l for l, chs in locs.items() if len(chs) >= 3), None)
        if loc is None:
            continue
        st = st.select(location=loc).merge(method=1, fill_value=0)
        sinv = inv.select(network=net.code, station=sta.code, location=loc, time=O)
        for tr in st:
            tag = f"{net.code}_{sta.code}_{tr.stats.channel}_{loc}"
            try:
                tr.write(f"{tag}.sac", format="SAC")
                cinv = sinv.select(channel=tr.stats.channel)
                cinv.write(f"SAC_PZs_{tag}", format="SACPZ")
                n += 1
            except Exception as e:
                print("  skip", tag, type(e).__name__)
        print(f"  + {net.code}.{sta.code} ({deg:.0f} deg)", flush=True)
print(f"wrote {n} channel files")
