# us6000tgb9 — M6.8 Kumamoto (near Uto, Japan), USGS preliminary (frozen 2026-07-28, BEFORE USGS FFM posted)

Origin  : 2026-07-28T07:27:15.512Z  32.6817 N  130.7217 E  depth 10 km  Mww 6.8
W-phase : M0 1.73e19 N-m   centroid depth 11.5 km   DC 87.6%
  NP1  strike 212.45  dip 75.07  rake -163.76   (right-lateral, ~Futagawa-Hinagu NE-SW trend -> likely fault plane)
  NP2  strike 118.16  dip 74.32  rake -15.52    (left-lateral conjugate)
  MT (x1e19 N-m): Mrr -0.3273  Mtt 1.5717  Mpp -1.2445  Mrt -0.1197  Mrp 0.527  Mtp 0.7936
Tectonics: shallow crustal strike-slip, Kumamoto fault zone (2016 Mw7.0 sequence region).
Products at freeze time: moment-tensor (Mww), shakemap, PAGER, ground-failure. NO finite-fault yet.

## Our prospective prediction (frozen 2026-07-28, before any USGS finite fault)

Data: 241 stations / 641 channels, 32-89 deg, max azimuthal gap 30 deg (best coverage of any event in this project).

WISP finite fault (both W-phase planes, Mww-seeded, shift_match applied):
- PREFERRED plane: 212/75/-164 (right-lateral, Futagawa-Hinagu trend), misfit 0.1607 -> 0.1507 after shift
- Conjugate 118/74/-16: 0.1667 -> 0.1523 (plane preference weak ~1%, but agrees with tectonics)
- Mw 6.75 (both planes) vs W-phase 6.8
- Slip: SINGLE compact asperity, peak 1.88 m, ~2-16 km depth, centred at/slightly updip of the 10-km hypocentre; no shallow surface slip maximum
- Duration T5-95 = 13.6 s, moment-rate peak at ~3.5 s, Vr ~2.05 km/s
- Qualitatively: a smaller sibling of the 2016 Kumamoto Mw7.0 rupture

MTUQ point source (independent check):
- 25-station subset (gap 97 deg, 14/25 in one 45-deg band): surface-wave solutions aliased ~90 deg
  (four-lobe strike-slip radiation) -> Kagan 67-78 deg vs W-phase. Documented failure.
- 60-station azimuth-balanced set (gap 80 deg): joint DC 130/88/7, Mw 6.80, Kagan 27.3 deg vs
  W-phase -> snapped back to the correct lobe family, confirming the aliasing diagnosis.
  Body P+SH remains aliased (91 deg): teleseismic SH unusable at this geometry.
- Pipeline lessons frozen with the prediction: (1) for strike-slip, azimuthal BALANCE beats station
  count; (2) the shallow-event "body-only" rule is a dip-slip (Kanamori-Given) cure and does not
  transfer to vertical strike-slip.

To compare when USGS posts their FFM: plane choice, Mw, peak slip amplitude + depth extent,
single-vs-multiple asperities, duration/moment-rate shape, rupture velocity.
