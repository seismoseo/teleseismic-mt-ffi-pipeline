"""Convert broadband data already fetched for WISP (event dir: 3-component SAC) into MTUQ
CAP-format input, reusing a single waveform fetch. Response removal + rotation follow the
VALIDATED path of 10_prep_mtuq_data.py exactly (remove_response with the StationXML inventory,
output DISP; rotate ->ZNE then NE->RT) — the SAC pole-zero attach_paz route is avoided because
it introduced a polarity flip. Only the (small, fast) StationXML metadata is fetched fresh;
the bulky waveforms are reused from disk.

Usage: <mtuq python> wisp_data_to_mtuq.py <event_dir> <tag> <lat> <lon> <depth> <time> <mindeg> <maxdeg> <nmax>
"""
import glob, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from obspy import read, UTCDateTime, Stream, Inventory
from obspy.clients.fdsn import Client
from obspy.io.sac import SACTrace
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

evdir, tag, lat, lon, depth, tstr, mindeg, maxdeg, nmax = (
    sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]),
    sys.argv[6], float(sys.argv[7]), float(sys.argv[8]), int(sys.argv[9]))
origin = UTCDateTime(tstr)
out = os.path.expanduser(f"~/works/17.Venezuela_2026/data/mtuq/{tag}")
os.makedirs(out, exist_ok=True)
pre = (0.004, 0.006, 0.5, 1.0)

# group SAC by station
byst = {}
for f in glob.glob(os.path.join(evdir, "*_BH?*.sac")):
    st = read(f, format="sac")[0]
    byst.setdefault((st.stats.network, st.stats.station, st.stats.location), {})[st.stats.channel[-1]] = f

# fetch StationXML metadata for all stations (fast; not bulky waveforms)
nets = ",".join(sorted(set(k[0] for k in byst)))
inv = Inventory([], None)
for cn in ("IRIS", "GFZ", "ORFEUS"):
    try:
        inv += Client(cn, timeout=60).get_stations(
            network=nets, starttime=origin - 60, endtime=origin + 60,
            channel="BH?", level="response", latitude=lat, longitude=lon,
            minradius=mindeg, maxradius=maxdeg)
    except Exception as e:
        print(cn, "inv err", e)

rows = []
for (net, sta, loc), chans in byst.items():
    if not all(c in chans for c in "ZNE"):
        continue
    s0 = read(chans["Z"], format="sac")[0].stats.sac
    d_m, az, baz = gps2dist_azimuth(lat, lon, s0["stla"], s0["stlo"])
    dist = kilometers2degrees(d_m / 1000.0)
    if not (mindeg <= dist <= maxdeg):
        continue
    try:
        sinv = inv.select(network=net, station=sta, time=origin)
        if len(sinv.get_contents()["channels"]) == 0:
            continue
        stt = Stream()
        for c in "ZNE":
            stt += read(chans[c], format="sac")[0]
        stt.remove_response(inventory=sinv, output="DISP", pre_filt=pre, water_level=60)
        stt.rotate("->ZNE", inventory=sinv)
        stt.rotate("NE->RT", back_azimuth=baz)
        stt.resample(10.0)
        stt.trim(starttime=origin - 200, pad=True, fill_value=0.0)
        traces = {c: stt.select(component=c)[0] for c in "ZRT"}
        rows.append((net, sta, loc, dist, az, s0["stla"], s0["stlo"], traces))
    except Exception:
        continue

rows.sort(key=lambda r: r[4])
if len(rows) > nmax:
    idx = sorted(set(np.linspace(0, len(rows) - 1, nmax).round().astype(int)))
    rows = [rows[i] for i in idx]

wl = []
for net, sta, loc, dist, az, sla, slo, traces in rows:
    for c, comp in (("Z", "z"), ("R", "r"), ("T", "t")):
        tr = traces[c]
        sac = SACTrace.from_obspy_trace(tr)
        sac.reftime = origin
        sac.b = float(tr.stats.starttime - origin)
        sac.o = 0.0
        sac.khole = loc if loc else "00"
        sac.kstnm, sac.knetwk, sac.kcmpnm = sta, net, "BH" + c
        sac.stla, sac.stlo, sac.evla, sac.evlo, sac.evdp = sla, slo, lat, lon, depth
        sac.write(os.path.join(out, f"{tag}.{net}.{sta}.{loc or '00'}.BH.{comp}"))
    wl.append(f"   {tag}.{net}.{sta}.{loc or '00'}.BH   {dist*111.19:.1f}   1 1 1 1 1   0. 0. 0. 0. 0. 0.")
with open(os.path.join(out, "weights.dat"), "w") as f:
    f.write("# event.net.sta.loc.ch offset_km w1 w2 w3 w4 w5 P bw S sw rw lw\n" + "\n".join(wl) + "\n")
print(f"{tag}: wrote {len(rows)} stations (waveforms reused, response via inventory) to {out}")
az = sorted(int(r[4]) for r in rows)
gaps = np.diff(np.array(az + [az[0] + 360]))
print("dist(deg):", sorted(round(r[3], 1) for r in rows))
print(f"az: {az}\nmax azimuthal gap: {gaps.max()} deg")
