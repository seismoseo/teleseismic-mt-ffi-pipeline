"""
fetch_teleseismic.py — robust teleseismic broadband fetcher for the MTUQ/WISP pipeline.

Replaces wisp_fetch_fallback.py, which under-sampled badly: on the 2026 Mw 7.3 Mexico event the
GSN inventory had 81 stations at 30-90 deg but the old fetcher kept only 15 (USGS used 63). Its two
faults were (1) it REQUIRED all three components, silently dropping every BHZ-only (P) station, and
(2) per-station get_waveforms with a bare `except: continue`, so any timeout / not-yet-arrived data
dropped the station with no record. That left an 87 deg azimuthal gap and an over-compact model.

This version:
  * queries the GSN core + a curated FDSN backbone that is GLOBALLY distributed (excludes dense
    regional arrays like TA/N4/CI that would pile many stations onto one azimuth),
  * keeps a station if it has at least BHZ (vertical/P is usable on its own; horizontals are added
    when present, enabling SH),
  * fetches with a single bulk request (fast, reliable) and RETRIES,
  * de-duplicates to one location code per station,
  * LOGS every dropped station with the reason and reports azimuthal coverage / gaps.

Writes WISP's layout into CWD:
  <NET>_<STA>_<CHA>_<LOC>.sac            (raw counts)
  SAC_PZs_<NET>_<STA>_<CHA>_<LOC>        (SACPZ response, incl. coords/azimuth/dip)
Then run wisp_fix_headers.py + `ffm model run . auto_model ...` as usual.

Usage: fetch_teleseismic.py <time> <lat> <lon> [mindeg=25] [maxdeg=95] [extra_networks]
"""
import sys
import numpy as np
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

TIME, LAT, LON = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
MINDEG = float(sys.argv[4]) if len(sys.argv) > 4 else 25.0
MAXDEG = float(sys.argv[5]) if len(sys.argv) > 5 else 95.0
# GSN core + globally-distributed backbone (NOT dense arrays: no TA/N4/CI/BK-dense)
NETS = "IU,II,G,GE,GT,IC,CU,GR,US,CN,DK,MN,GO,PS,IM,AU,JP,PM,NA,CX,IV"
if len(sys.argv) > 6:
    NETS = NETS + "," + sys.argv[6]
O = UTCDateTime(TIME)
T1, T2 = O - 120, O + 3600           # 1 h window (long-period surface waves)
CHANS = "BHZ,BHN,BHE,BH1,BH2"

cl = Client("IRIS", timeout=120)
print(f"inventory query: nets={NETS} dist {MINDEG}-{MAXDEG} deg ...", flush=True)
inv = cl.get_stations(network=NETS, channel="BH?", level="response",
                      starttime=T1, endtime=T2, latitude=LAT, longitude=LON,
                      minradius=MINDEG, maxradius=MAXDEG)

# one entry per (net, sta): distance/azimuth + the response inventory
stations = {}
for net in inv:
    for sta in net:
        key = (net.code, sta.code)
        if key in stations:
            continue
        d_m, az, _ = gps2dist_azimuth(LAT, LON, sta.latitude, sta.longitude)
        stations[key] = dict(deg=kilometers2degrees(d_m / 1000.0), az=az)
print(f"inventory: {len(stations)} unique stations", flush=True)

# one bulk request for everything (retry once)
bulk = [(n, s, "*", "BH?", T1, T2) for (n, s) in stations]
st = Stream()
for attempt in (1, 2):
    try:
        st = cl.get_waveforms_bulk(bulk)
        break
    except Exception as e:
        print(f"  bulk attempt {attempt} failed: {type(e).__name__}", flush=True)
# fill gaps in any station that returned partial data individually (best-effort)
have = set((tr.stats.network, tr.stats.station) for tr in st)
missing = [k for k in stations if k not in have]
for (n, s) in missing:
    try:
        st += cl.get_waveforms(n, s, "*", "BH?", T1, T2)
    except Exception:
        pass

# group by station -> choose the location code that has BHZ; require BHZ (P) at minimum
by_sta = {}
for tr in st:
    by_sta.setdefault((tr.stats.network, tr.stats.station), []).append(tr)

n_written, kept, dropped = 0, [], []
for (net, sta), info in sorted(stations.items(), key=lambda kv: kv[1]["az"]):
    traces = by_sta.get((net, sta), [])
    if not traces:
        dropped.append((net, sta, "no data returned")); continue
    # per-location channel sets
    locs = {}
    for tr in traces:
        locs.setdefault(tr.stats.location, set()).add(tr.stats.channel)
    # prefer a location that has BHZ; among those, the one with the most components
    zlocs = [(l, chs) for l, chs in locs.items() if any(c.endswith("Z") for c in chs)]
    if not zlocs:
        dropped.append((net, sta, "no BHZ")); continue
    loc = max(zlocs, key=lambda lc: len(lc[1]))[0]
    sel = Stream([tr for tr in traces if tr.stats.location == loc]).merge(method=1, fill_value=0)
    sinv = inv.select(network=net, station=sta, location=loc, time=O)
    wrote_here = 0
    for tr in sel:
        tag = f"{net}_{sta}_{tr.stats.channel}_{loc}"
        try:
            tr.write(f"{tag}.sac", format="SAC")
            sinv.select(channel=tr.stats.channel).write(f"SAC_PZs_{tag}", format="SACPZ")
            n_written += 1; wrote_here += 1
        except Exception as e:
            print(f"  skip {tag}: {type(e).__name__}", flush=True)
    if wrote_here:
        kept.append((net, sta, info["az"], info["deg"], wrote_here))

# coverage report
az = sorted(int(k[2]) for k in kept)
gaps = np.diff(np.array(az + [az[0] + 360])) if az else np.array([360])
print(f"\nKEPT {len(kept)} stations ({n_written} channel files); dist "
      f"{min(k[3] for k in kept):.0f}-{max(k[3] for k in kept):.0f} deg; "
      f"max azimuthal gap {gaps.max()} deg", flush=True)
print(f"DROPPED {len(dropped)} stations:", flush=True)
for net, sta, why in dropped[:40]:
    print(f"   {net}.{sta}: {why}", flush=True)
print("stations kept (az):", [f"{k[0]}.{k[1]}({k[2]:.0f})" for k in kept], flush=True)
