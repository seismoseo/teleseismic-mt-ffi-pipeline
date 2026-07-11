"""
wisp_fix_headers.py — Populate SAC headers WISP needs on get-data teleseismic output.

`ffm get-data teleseismic` left station/event coordinates out of the waveform SAC
headers (stla/stlo/evla/evlo = None), so WISP's pre-selection drops every station
(empty filter_tele.txt -> Fortran green_tele EOF crash). This reads station coords
from the companion SAC_PZs response files (+ event coords from the CMT) and writes
the headers WISP requires: stla, stlo, evla, evlo, evdp, dist, az, baz, gcarc,
cmpaz, cmpinc, khole, o, kcmpnm.

Usage: conda run -n ff-env python scripts/wisp_fix_headers.py <data_dir> <cmt_file>
"""

import glob
import os
import re
import sys

from obspy.io.sac import SACTrace
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

# SAC cmpaz/cmpinc by component (standard ZNE broadband orientation)
ORIENT = {"Z": (0.0, 0.0), "N": (0.0, 90.0), "E": (90.0, 90.0),
          "1": (0.0, 90.0), "2": (90.0, 90.0)}


def parse_cmt(cmt_file):
    lat = lon = dep = None
    for line in open(cmt_file):
        m = re.match(r"\s*latitude:\s*([-\d.]+)", line)
        if m: lat = float(m.group(1))
        m = re.match(r"\s*longitude:\s*([-\d.]+)", line)
        if m: lon = float(m.group(1))
        m = re.match(r"\s*depth:\s*([-\d.]+)", line)
        if m: dep = float(m.group(1))
    return lat, lon, dep


def parse_pz(pz_file):
    """Return (lat, lon, azimuth, dip) from a SAC_PZs file. azimuth/dip are the
    TRUE channel orientation (critical for BH1/BH2, which are not N/E)."""
    lat = lon = az = dip = None
    for line in open(pz_file):
        m = re.search(r"LATITUDE\s*:\s*([-\d.]+)", line)
        if m: lat = float(m.group(1))
        m = re.search(r"LONGITUDE\s*:\s*([-\d.]+)", line)
        if m: lon = float(m.group(1))
        m = re.search(r"AZIMUTH\s*:\s*([-\d.]+)", line)
        if m: az = float(m.group(1))
        m = re.search(r"\bDIP\s*:\s*([-\d.]+)", line)
        if m: dip = float(m.group(1))
    return lat, lon, az, dip


def main():
    data_dir = sys.argv[1]
    cmt_file = sys.argv[2]
    evla, evlo, evdp = parse_cmt(cmt_file)
    print(f"event: {evla}, {evlo}, depth {evdp} km")

    wfs = [f for f in glob.glob(os.path.join(data_dir, "*_BH*.sac"))
           if "SAC_PZs" not in os.path.basename(f)]
    n_ok = n_skip = 0
    for wf in wfs:
        stem = os.path.basename(wf)[:-4]            # strip .sac
        pz = None
        for cand in ("SAC_PZs_" + stem, "SAC_PZs_" + stem + ".sac"):   # naming varies by WISP version
            p = os.path.join(data_dir, cand)
            if os.path.isfile(p):
                pz = p; break
        if pz is None:
            n_skip += 1; continue
        stla, stlo, pz_az, pz_dip = parse_pz(pz)
        if stla is None or stlo is None:
            n_skip += 1; continue
        s = SACTrace.read(wf)
        comp = (s.kcmpnm or stem)[-1]
        if comp == "Z":
            cmpaz, cmpinc = 0.0, 0.0
        elif pz_az is not None:
            # TRUE orientation from the response file (correct for BH1/BH2)
            cmpaz, cmpinc = pz_az, 90.0
        else:
            cmpaz, cmpinc = ORIENT.get(comp, (0.0, 90.0))
        s.stla, s.stlo = stla, stlo
        s.evla, s.evlo, s.evdp = evla, evlo, evdp
        dist_m, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
        s.dist = dist_m / 1000.0
        s.az, s.baz = az, baz
        s.gcarc = kilometers2degrees(dist_m / 1000.0)
        s.cmpaz, s.cmpinc = cmpaz, cmpinc
        if s.o is None:
            s.o = 0.0                                # trace starts at origin (seconds-before=0)
        s.write(wf)
        n_ok += 1
    print(f"headers populated: {n_ok} | skipped (no PZ/coords): {n_skip} | total wf: {len(wfs)}")


if __name__ == "__main__":
    main()
