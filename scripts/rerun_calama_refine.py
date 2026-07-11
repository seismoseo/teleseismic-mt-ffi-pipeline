"""
rerun_sanriku_refine.py — refine the Sanriku WISP FFI with (a) Venezuela-wavetrain
contamination QC (the M7.2+M7.5 doublet 25 min earlier sweeps teleseismic windows),
(b) generic-response rejection (US-network A0~3.5e15 pattern), (c) shift_match, re-invert.

Usage: python scripts/rerun_sanriku_refine.py NP2   (or NP1)
"""
import glob, json, os, pathlib, re, shutil, sys

NP = sys.argv[1] if len(sys.argv) > 1 else "NP2"
KEEP_RESP = "keepresp" in sys.argv   # keep the A0~3.5e15 stations (fits look fine; user call)
FFM = next((a for a in sys.argv[2:] if a.startswith("ffm.")), "ffm.0")
BASE = pathlib.Path(os.path.expanduser("~/works/17.Venezuela_2026/event_2026_calama"))
EV = BASE / "20260525215220" / FFM
SRC, DST = EV / NP, EV / (f"{NP}_ref2" if KEEP_RESP else f"{NP}_ref")
CONFIG = os.path.expanduser("~/works/neic-finitefault/config.ini")
# Venezuela doublet timing relative to Japan origin (22:30:12.99)
VT0, VT1 = -1e9, -1e9+1
VLAT, VLON = 0.0, 0.0


def main():
    sys.path.insert(0, os.path.expanduser("~/works/neic-finitefault/src"))
    import ffm.management as mng
    from ffm.inversion_chen_new import manual_modelling
    from ffm.shift_match import shift_match2, save_waveforms
    from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
    from obspy.taup import TauPyModel
    import numpy as np
    m = TauPyModel("ak135")

    if DST.exists():
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    for f in ("Solution.txt", "modelling_summary.txt"):
        (DST / f).unlink(missing_ok=True)
    if (DST / "plots").exists():
        shutil.rmtree(DST / "plots")

    # WISP surf GF bank ends at 126.5 km: trim down-dip rows until the fault bottom fits,
    # else WISP hard-errors (manual) or silently drops ALL surface waves (auto).
    import math
    sd = json.load(open(DST / "segments_data.json"))
    seg = sd["segments"][0]
    hyp_depth = json.load(open(DST / "tensor_info.json"))["depth"]
    while (hyp_depth + (seg["dip_subfaults"] - seg["hyp_dip"] + 0.5) * seg["delta_dip"]
           * math.sin(math.radians(seg["dip"]))) > 126.5 and seg["dip_subfaults"] > 2:
        seg["dip_subfaults"] -= 1
    json.dump(sd, open(DST / "segments_data.json", "w"), indent=4)
    print(f"dip_subfaults = {seg['dip_subfaults']}")

    def zero(w):
        w["trace_weight"] = 0.0
        n = len(str(w.get("wavelet_weight", "0 1 1 1 1 1 1 1")).split())
        w["wavelet_weight"] = " ".join(["0"] * n) + "\n"

    def bad_response(name):
        g = glob.glob(str(BASE / f"SAC_PZs_*_{name}_BH*"))
        if not g:
            return False
        t = open(g[0]).read()
        a0 = re.search(r"A0\s*:\s*([\deE.+-]+)", t)
        return a0 and 3.0e15 < float(a0.group(1)) < 4.0e15

    def contaminated(w, kind):
        """kind: 'P','SH','SURF' — does the Venezuela train overlap this channel's window?"""
        stla, stlo = w["location"]
        dV = gps2dist_azimuth(VLAT, VLON, stla, stlo)[0] / 1000.0
        v = (VT0 + dV / 4.7, VT1 + dV / 2.9)
        degJ = kilometers2degrees(w["distance"] * 111.19) if w["distance"] < 200 else w["distance"] / 111.19
        # WISP stores distance in km? entries use degrees for tele. use distance as deg if <180
        degJ = w["distance"] if w["distance"] < 180 else w["distance"] / 111.19
        if kind == "P":
            tp = m.get_travel_times(46.0, degJ, ["P"])[0].time
            win = (tp - 15, tp + 65)
        elif kind == "SH":
            ts = m.get_travel_times(46.0, degJ, ["S"])[0].time
            win = (ts - 15, ts + 85)
        else:
            dkm = degJ * 111.19
            win = (dkm / 4.6, dkm / 3.2)
        return (win[0] < v[1]) and (v[0] < win[1])

    for fjson, kinds in [("tele_waves.json", {"BHZ": "P", "BHT": "SH"}),
                         ("surf_waves.json", {"BHZ": "SURF", "BHT": "SURF"})]:
        d = json.load(open(DST / fjson))
        nz = nr = 0
        for w in d:
            if w["trace_weight"] == 0:
                continue
            if not KEEP_RESP and bad_response(w["name"]):
                zero(w); nr += 1; continue
            if contaminated(w, kinds[w["component"]]):
                zero(w); nz += 1
        d = [w for w in d if w["component"] == "BHZ"] + [w for w in d if w["component"] != "BHZ"]
        json.dump(d, open(DST / fjson, "w"), indent=4)
        act = sum(w["trace_weight"] > 0 for w in d)
        print(f"{fjson}: zeroed {nz} contaminated + {nr} bad-response; {act} active")

    tensor_info = json.load(open(DST / "tensor_info.json"))
    seg_data = json.load(open(DST / "segments_data.json"))
    default_dirs = mng.default_dirs(config_path=CONFIG)
    manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)

    def misfit():
        return float(re.search(r"averaged misfit error\s+([\d.]+)",
                               open(DST / "modelling_summary.txt").read()).group(1))
    print(f"pre-shift misfit  {misfit():.4f}")
    cwd = os.getcwd(); os.chdir(DST)
    try:
        save_waveforms("body", shift_match2("body", directory=DST))
        save_waveforms("surf", shift_match2("surf", directory=DST))
    finally:
        os.chdir(cwd)
    manual_modelling(tensor_info, ["body", "surf"], default_dirs, seg_data, directory=DST)
    print(f"post-shift misfit {misfit():.4f}")
    print(f"DONE -> {DST}")


if __name__ == "__main__":
    main()
