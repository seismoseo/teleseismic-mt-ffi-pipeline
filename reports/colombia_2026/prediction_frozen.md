# us6000tjl2 — M7.4 San Jose del Palmar, Colombia (frozen 2026-08-11, BEFORE USGS FFM posted)

Origin  : 2026-08-10T12:34:28.125Z  4.8436 N  76.2422 W  depth 110.3 km  Mww 7.4
W-phase : M0 1.66e20 N-m   centroid depth 120.5 km   DC 86.6%
  NP1  strike 123.99  dip 58.63  rake 158.52
  NP2  strike 225.57  dip 71.78  rake  33.24   (oblique reverse/strike-slip pair)
  MT (x1e20 N-m): Mrr 0.4887  Mtt -1.5429  Mpp 1.0541  Mrt -0.1918  Mrp -0.9009  Mtp -0.2243
Products at freeze: MT, shakemap, PAGER, DYFI, ground failure. NO finite-fault yet.

## Bucaramanga nest? NO.
The Bucaramanga nest is at ~6.82 N 73.15 W, 140-160 km depth (NE Colombia) — this event is
~415 km to the SW at 4.84 N 76.24 W, 110-120 km deep: the CAUCA intermediate-depth cluster
of the subducting NAZCA slab (cf. Chang et al. 2019), a different slab segment from the
Bucaramanga slab/nest. Largest instrumental event in the Cauca segment vicinity.

## Deep-event protocol (from Calama/SPdA lessons)
- MTUQ body band UP: --bwband 25,60,120,10 (10-33 s mixes P/pP/sP and fails confidently)
- Depth from the body-only misfit-vs-depth curve (surface band has no depth resolution)
- M7.4 surface waves: mantle band --swband 200,450,900
- WISP surface-wave GF bank ends ~126 km; centroid 120.5 -> fault WILL straddle the cap:
  grep logs for "Maximum depth", compare planes body-only if SW dropped asymmetrically

## Our prospective prediction (frozen 2026-08-11, before any USGS finite fault)

Data: 120 stations / 316 channels, 30-90 deg, max gap 47 deg (Pacific).

KEY SOURCE-COMPLEXITY DISCOVERY (drives everything):
- The event begins with a WEAK EMERGENT NUCLEATION PHASE lasting ~35 s; the main moment
  release starts near the W-phase centroid TIME (12:35:02.6 = origin + 34.5 s) and PLACE
  (4.585 N 76.161 W, 29 km S of the epicentre).
- Apparent nucleation-to-asperity velocity ~0.8 km/s << any allowed Vr: origin-anchored
  finite-fault runs are kinematically impossible -> slip piles at the fault edge (86% of
  peak at the domain edge) and MTUQ windows miss the main energy (Mw pinned at grid edge).
  BOTH analyses re-anchored at the W-phase centroid (space AND time).

WISP finite fault (centroid-anchored, both W-phase planes, body-only above/below the 126.5-km
surface-wave GF cap -> fair two-plane comparison; shift_match applied):
- PREFERRED plane: 226/72/33 (NE-SW, steep, oblique reverse) misfit 0.1363
  vs conjugate 124/59/159: 0.1485 -> 8.2% margin, near-decisive; conjugate additionally
  violates the edge rule (40% of peak at its top edge; preferred plane 1-3%).
  NOTE: WISP directory labels are SWAPPED vs USGS (WISP NP1 = W-phase NP2), as for Kumamoto.
- Slip: single dominant asperity ~17 km SSW of the centroid (~45 km SSW of the epicentre), up-dip,
  85-125 km depth, peak 1.63 m; edge slip 1-3% (clean containment)
- M0 1.68e20 N-m (W-phase 1.66e20); main-phase T5-95 ~18.7 s (excludes the precursor)
- Effective area 2695 km2 -> radius 29 km -> stress drop ~2.9 MPa (physical)

MTUQ point source (honest failure, documented):
- At the centroid, joint 202/81/-65 (Kagan 95) and body 220/61/-52 (Kagan 84), Mw pinned
  at the lower grid edge; body-only depth curve minimum 70-90 km (shallower than centroid),
  mechanism unstable except ~90 km where it snaps to 225/71/4 (Kagan 28.6).
- Interpretation: a 40-km-offset late asperity + 35-s precursor + 40-km-deep slip zone
  violate the point-source approximation; the W-phase (its own centroid search) remains the
  mechanism authority. Pipeline rule: for complex sources, fix mechanism from W-phase and
  spend teleseismic data on the finite fault.

To compare when USGS posts their FFM: plane choice (226/72 vs 124/59), asperity position
(SSW + up-dip of centroid), depth extent 85-125 km, peak slip ~1.6 m, main-phase duration
~19 s, the precursor treatment (do they relocate the hypocentre?).

## Directivity (added 2026-08-11, same freeze window; corrects an axis-misread 'NE' in the
first push - geographic computation below is authoritative)
Two-stage, predominantly unilateral SOUTHWARD + UP-DIP directivity:
1) model-independent: epicentre -> W-phase centroid = 30 km toward azimuth 163 (SSE),
   110->120.5 km, ~34.5 s (apparent 0.9 km/s);
2) preferred plane: early->late moment centroid 27 km toward azimuth 220 (SW), depth
   110->101 km (up-dip), 75%/24% forward/backward moment = unilateral;
3) up-dip shallowing robust on BOTH planes; horizontal direction on the (rejected) conjugate
   is northward and inconsistent with (1), further supporting the 226/72/33 plane.

## Post-freeze comparison (USGS FFM posted 2026-08-14T19:14 UTC = 4 days AFTER our frozen 42b65a2)

VERIFIED: plane 225.57/71.78/33.23 (bit-for-bit our preferred plane = W-phase NP2); M0 1.72e20
vs ours 1.68e20 (2.5%); main asperity 31 km SSW (az 209) of the epicentre at 100 km depth
(ours ~44 km S, 85-125 km; theirs 91-108); unilateral SW directivity 46 km @ az 225 (ours
27 @ 220); main-phase duration 24 s (ours ~19); main-phase front speed 2.34 km/s by the same
metric (ours 2.52).
THE NUCLEATION QUESTION (frozen scorecard item): USGS ALSO could not use the cataloged
hypocentre - they relocated it (24 km NE of catalog, 125 km deep; neither catalog nor W-phase
centroid) and let a single front CRAWL silently (~30 s, avVr header 1.28 km/s) to the asperity:
their STF carries 0.2% of moment in the first 25 s - independent confirmation of our stage-1
moment bound. Two parameterizations of the same two-stage physics.
Real differences: peak slip 3.98 vs 1.63 m and dsigma 8.1 vs 2.9 MPa (their 6-km subfaults /
less smoothing = more compact asperity); initiation point location remains the genuinely
ill-determined element.
