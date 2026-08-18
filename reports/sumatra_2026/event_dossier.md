# us6000tkyk — M6.9 near Pematangsiantar, N Sumatra, INDONESIA (frozen 2026-08-15, BEFORE any USGS FFM)

Origin  : 2026-08-15T10:54:51.796Z  3.0856 N  99.0150 E  depth 172.5 km  Mww 6.9
W-phase : M0 2.62e19 N-m   centroid depth 180.5 km   DC 98.6%  (very high DC)
  NP1  strike 207.96  dip 58.15  rake   1.48   (left-lateral strike-slip on a moderately dipping plane)
  NP2  strike 117.18  dip 88.75  rake 148.14   (near-vertical, right-lateral oblique)
  MT (x1e19 N-m): Mrr 0.0716  Mtt -1.865  Mpp 1.7934  Mrt 1.2132  Mrp -0.6633  Mtp -1.2696
  W-phase centroid: 3.0802 N 99.0130 E, 10:55:02.7 = origin + 10.9 s (11 km... see note)
Other USGS MTs: Mwc 6.9 (M0 2.98e19, DC 74%, centroid 2.39N -> 77 km S: LOW-DC, suspect),
                Mwb 6.8 (M0 1.88e19, DC 81%, centroid ~origin time)
Tectonics: INTRASLAB, subducted Indo-Australian (Sunda) slab beneath N Sumatra, 170-180 km deep,
below the Toba caldera region. NOT the megathrust (which is offshore, <60 km here).
Products at freeze: MT, shakemap, PAGER, DYFI, ground failure. NO finite-fault yet.

## Deep-event protocol (Calama/SPdA/Colombia lessons)
- MTUQ body band UP: --bwband 25,60,120,10 (10-33 s mixes P/pP/sP -> confident ~90 deg failure)
- Depth from the BODY-ONLY misfit-vs-depth curve; surface band has no depth resolution
- 180 km is FAR below the ~126.5 km WISP surface-wave GF cap -> surface waves WILL be dropped
  on both planes -> body-only two-plane comparison, fair by construction (SPdA rule)
- Check hypocentre-vs-centroid consistency BEFORE inverting (Colombia lesson): W-phase centroid
  time is +10.9 s, spatially ~1 km from the epicentre -> apparent migration is negligible,
  so a normal origin-anchored inversion should be valid here (unlike Colombia).
- Mw 6.9 at 180 km: expect a compact rupture (~10-20 km); plane discrimination likely weak
  (Calama M6.9 at 109 km gave a 1.6% near-tie)
