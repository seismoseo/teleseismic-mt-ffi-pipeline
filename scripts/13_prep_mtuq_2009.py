"""
13_prep_mtuq_2009.py — download + CAP-format teleseismic 3-C broadband for the 2009 M6.4
Puerto Cabello event, for MTUQ moment-tensor inversion.

The 2009 event predates the 2026 regional dataset, and regional coverage in 2009 was sparse,
so we use TELESEISMIC GSN 3-C data (30-90 deg) — clean, no clipping — with syngine ak135 GFs.
Deconvolve to DISPLACEMENT, rotate to Z/R/T, common sampling rate, CAP naming + weights.dat.

Run in the mtuq env:
  python scripts/13_prep_mtuq_2009.py                       # default = 2009 M6.4
  python scripts/13_prep_mtuq_2009.py --event tele_us6000t7zp \
      --time 2026-06-24T22:05:11 --lat 10.435 --lon -68.472 --depth 10 --sr 1
Outputs -> data/mtuq/<event>/
"""
import argparse, os, sys
import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.io.sac import SACTrace
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees, locations2degrees

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C

_ap = argparse.ArgumentParser()
_ap.add_argument("--event", default="2009puertocabello")
_ap.add_argument("--time", default="2009-09-12T20:06:25.47")
_ap.add_argument("--lat", type=float, default=10.7090)
_ap.add_argument("--lon", type=float, default=-67.9270)
_ap.add_argument("--depth", type=float, default=14.0)
_ap.add_argument("--mindeg", type=float, default=25.0)
_ap.add_argument("--maxdeg", type=float, default=90.0)
_ap.add_argument("--sr", type=float, default=5.0)
_ap.add_argument("--nmax", type=int, default=24)
_a = _ap.parse_args()

EVENT_ID = _a.event
ORIGIN = UTCDateTime(_a.time)
EVLA, EVLO, EVDP = _a.lat, _a.lon, _a.depth
NETS = "IU,II,G,GT,IC,CU,GE"
MINDEG, MAXDEG = _a.mindeg, _a.maxdeg
COMMON_SR = _a.sr
NMAX = _a.nmax


def main():
    out_dir = os.path.join(C.DATA_DIR, "mtuq", EVENT_ID)
    os.makedirs(out_dir, exist_ok=True)
    cl = Client("IRIS")
    print("Fetching station inventory...", flush=True)
    inv = cl.get_stations(network=NETS, channel="BH?", level="response",
                          starttime=ORIGIN - 60, endtime=ORIGIN + 3600,
                          latitude=EVLA, longitude=EVLO, minradius=MINDEG, maxradius=MAXDEG)
    cand = []
    for net in inv:
        for sta in net:
            d = locations2degrees(EVLA, EVLO, sta.latitude, sta.longitude)
            cand.append((d, net.code, sta.code, sta.latitude, sta.longitude))
    cand.sort()
    print(f"{len(cand)} candidate stations {MINDEG}-{MAXDEG} deg", flush=True)

    rows = []
    for deg, net, sta, stla, stlo in cand:
        if len(rows) >= NMAX:
            break
        try:
            st = cl.get_waveforms(net, sta, "*", "BH?", ORIGIN - 60, ORIGIN + 3600,
                                  attach_response=False)
        except Exception:
            continue
        # pick a location code present on 3 components
        locs = {}
        for tr in st:
            locs.setdefault(tr.stats.location, set()).add(tr.stats.channel[-1])
        loc = next((l for l, comps in locs.items()
                    if {"Z"} <= comps and len(comps) >= 3), None)
        if loc is None:
            continue
        st = st.select(location=loc).merge(method=1, fill_value=0)
        try:
            sinv = inv.select(network=net, station=sta, location=loc, time=ORIGIN)
            st.detrend("demean"); st.detrend("linear"); st.taper(0.05)
            sr = st[0].stats.sampling_rate
            st.remove_response(inventory=sinv, output="DISP",
                               pre_filt=(0.001, 0.002, 0.45 * sr, 0.49 * sr), water_level=None)
            st.rotate("->ZNE", inventory=sinv)
            dist_m, az, baz = gps2dist_azimuth(EVLA, EVLO, stla, stlo)
            st.rotate("NE->RT", back_azimuth=baz)
            for tr in st:
                if abs(tr.stats.sampling_rate - COMMON_SR) > 1e-6:
                    tr.resample(COMMON_SR)
        except Exception as e:
            print(f"  skip {net}.{sta}: {type(e).__name__}: {str(e)[:50]}")
            continue

        dist_km = dist_m / 1000.0
        ok = 0
        for comp, suf in [("Z", "z"), ("R", "r"), ("T", "t")]:
            sel = st.select(component=comp)
            if not len(sel):
                continue
            tr = sel[0]
            sac = SACTrace.from_obspy_trace(tr)
            sac.reftime = ORIGIN; sac.o = 0.0; sac.iztype = "io"
            sac.knetwk, sac.kstnm, sac.khole, sac.kcmpnm = net, sta, loc, "BH" + comp
            sac.stla, sac.stlo = stla, stlo
            sac.evla, sac.evlo, sac.evdp = EVLA, EVLO, EVDP
            sac.dist, sac.az, sac.baz = dist_km, az, baz
            sac.write(os.path.join(out_dir, f"{EVENT_ID}.{net}.{sta}.{loc}.BH.{suf}"))
            ok += 1
        if ok == 3:
            rows.append((f"{EVENT_ID}.{net}.{sta}.{loc}.BH", dist_km, az))
            print(f"  + {net}.{sta} {deg:.0f} deg az{az:.0f}", flush=True)

    with open(os.path.join(out_dir, "weights.dat"), "w") as fh:
        fh.write("# event.net.sta.loc.ch  offset_km  w1 w2 w3 w4 w5  "
                 "P bw S sw rw lw\n")
        for code, dist_km, az in rows:
            fh.write(f"   {code}   {dist_km:.0f}   1 1 1 1 1   0. 0. 0. 0. 0. 0.\n")
    print(f"\nPrepared {len(rows)} stations for {EVENT_ID}")
    print(f"  azimuths: {sorted(int(r[2]) for r in rows)}")
    print(f"  -> {out_dir}")


if __name__ == "__main__":
    main()
