# us6000kzv3 — M6.1 Mariana Islands region, 2023-08-14 (analysis started 2026-07-29)

Origin  : 2023-08-14T13:51:53.761Z  13.3470 N  147.5297 E  depth 8 km  Mww 6.1
W-phase : M0 1.74e18 N-m   centroid depth 13.5 km   DC 91.0%
  NP1  strike 253.28  dip 56.11  rake -28.72   (oblique: strike-slip + normal)
  NP2  strike 000.27  dip 66.49  rake -142.55  (oblique conjugate)
  MT (x1e18 N-m): Mrr -0.7986  Mtt 0.0078  Mpp 0.7908  Mrt 0.4953  Mrp -0.6972  Mtp 1.2952
Also USGS Mwb 6.0 (body-wave MT, 30% DC only — unstable, use Mww).
Tectonics: Mariana Islands region, east of Guam near the trench/outer rise; shallow oblique
normal + strike-slip faulting (intraplate/outer-rise style, NOT interface thrust).
USGS products: moment tensor, shakemap, PAGER, DYFI, ground failure. NO finite-fault product
(none expected at M6.1) — so no reference model exists.

## The research question
Can a RELIABLE finite-fault model be resolved for an M6.1 from teleseismic data alone?
Prior pipeline experience says this is at/beyond the resolution limit (Sanriku M6.9 FFI was
already marginal by the duration-vs-period criterion). Expected source scales for Mw 6.1:
duration ~3-6 s, rupture length ~8-12 km — comparable to ONE WISP subfault and far below the
10-30 s teleseismic body-wave passband. The deliverable is an honest quantified resolution
assessment, not a forced slip model.

## Precedent (user-supplied): Hartzell, Mendoza & Zeng (2013, GRL)
Teleseismic-P finite-fault of the 2011 Mw 5.8 Mineral, Virginia earthquake — proof that
M<6 FFI from teleseismic data IS possible when: P waves are used to ~1 Hz (t*~1 s passband),
resolution is carried by TIMING not amplitude, subfaults are ~1 km, and the fault is compact.
Strategy here follows that recipe: high-frequency body-wave band, compact fault with small
subfaults, resolution argued from P-wave timing; surface waves constrain moment/mechanism only.
