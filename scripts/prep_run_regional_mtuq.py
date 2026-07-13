"""Prepare regional broadband data (CX/C1/C/IU/GE, ~1.5-6 deg) for MTUQ and write CAP-format
input, for the deep Chile events. Regional MT is MTUQ's native regime (Zhu & Helmberger CAP).
Removes response to displacement, rotates to ZRT is left to MTUQ; here we write raw Z/N/E
displacement SAC + weights.dat, sorted by distance. syngine ak135 GFs (valid to 700 km depth)
serve the regional distances.

Usage: <mtuq python> prep_run_regional_mtuq.py <tag> <lat> <lon> <depth_km> <time> <maxdeg>
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.io.sac import SACTrace

tag, lat, lon, depth, tstr, maxdeg = (sys.argv[1], float(sys.argv[2]), float(sys.argv[3]),
                                      float(sys.argv[4]), sys.argv[5], float(sys.argv[6]))
t = UTCDateTime(tstr)
out = os.path.expanduser(f"~/works/17.Venezuela_2026/data/mtuq/{tag}")
os.makedirs(out, exist_ok=True)

rows = []
for cname in ("IRIS", "GEOFON"):
    try:
        cl = Client(cname, timeout=90)
    except Exception:
        continue
    try:
        inv = cl.get_stations(latitude=lat, longitude=lon, minradius=1.0, maxradius=maxdeg,
                              starttime=t-60, endtime=t+600, channel="BH?,HH?", level="response")
    except Exception as e:
        print(cname, "stations err", e); continue
    for netw in inv:
        for st in netw:
            key = f"{netw.code}.{st.code}"
            if any(r[0] == key for r in rows):
                continue
            d_m, az, baz = gps2dist_azimuth(lat, lon, st.latitude, st.longitude)
            dist = kilometers2degrees(d_m/1000.0)
            band = "HH" if any(c.code.startswith("HH") for c in st.channels) else "BH"
            try:
                stt = cl.get_waveforms(netw.code, st.code, "*", band+"?", t-60, t+600,
                                       attach_response=True)
                if len(stt) < 3:
                    continue
                stt.merge(method=1, fill_value=0)
                stt.remove_response(output="DISP", water_level=60,
                                    pre_filt=(0.004, 0.006, 0.5, 1.0))
                stt.resample(10.0)
                rows.append((key, dist, az, netw.code, st.code, st.latitude, st.longitude, stt))
            except Exception:
                continue
rows.sort(key=lambda r: r[1])
rows = rows[:24]
wl = []
for key, dist, az, net, sta, sla, slo, stt in rows:
    for tr in stt:
        c = tr.stats.channel[-1]
        comp = {"Z": "z", "1": "r", "2": "t", "N": "r", "E": "t"}.get(c)
        if comp is None:
            continue
        sac = SACTrace.from_obspy_trace(tr)
        sac.reftime = t
        sac.b = float(tr.stats.starttime - t)
        sac.stla, sac.stlo, sac.evla, sac.evlo, sac.evdp = sla, slo, lat, lon, depth
        sac.o = 0.0
        sac.kstnm, sac.knetwk = sta, net
        sac.khole = "00"
        sac.kcmpnm = tr.stats.channel
        b = tr.stats.channel[:2]
        sac.write(os.path.join(out, f"{tag}.{net}.{sta}.00.{b}.{comp}"))
    b = stt[0].stats.channel[:2]
    wl.append(f"   {tag}.{net}.{sta}.00.{b}   {dist*111.19:.1f}   1 1 1 1 1   0. 0. 0. 0. 0. 0.")
with open(os.path.join(out, "weights.dat"), "w") as f:
    f.write("# event.net.sta.loc.ch offset_km w1 w2 w3 w4 w5 P bw S sw rw lw\n" + "\n".join(wl) + "\n")
print(f"{tag}: wrote {len(rows)} regional stations to {out}")
print("dist(deg):", [round(r[1], 1) for r in rows])
print("az:", sorted(int(r[2]) for r in rows))
