# Teleseismic MT + FFI pipeline

A reproducible, battle-tested workflow for determining the **point-source moment tensor**
(MTUQ) and a **kinematic finite-fault model** (WISP, the USGS operational code) of any
moderate-to-large earthquake from openly available teleseismic data — from data fetching to
publication-grade figures and a PDF report.

Built and validated on nine real cases (2023–2026), each of which broke the standard workflow
in a different, instructive way. Every failure mode found along the way is encoded in the scripts
or documented in the cookbook, so the next event starts from here instead of rediscovering
them.

## What's here

```
scripts/      the pipeline (numbered MTUQ runners, WISP helpers, figure generators)
notebooks/    23 = the step-by-step COOKBOOK (start here)
              24 = case study: 2026 Calama M6.9 intraslab, 109 km deep (executed, all figures)
              25 = case study: 2024 San Pedro de Atacama M7.4, 127 km deep — teleseismic
                   plane discrimination on a disputed event (executed)
              22 = case study: 2026 Sanriku-oki M6.9 under wavetrain interference (executed)
              26 = case study: 2023 Al Haouz (Morocco) M6.8 retrospective vs published models
              27 = BLIND TEST: 2026 Mexico M7.3 — prediction frozen before the USGS finite
                   fault, then compared head-to-head (executed)
              28 = BLIND TEST: 2026 Kumamoto M6.8 strike-slip — prediction frozen 3.2 h before
                   the USGS finite fault, then compared head-to-head (executed)
              29 = RESOLUTION TEST: 2023 Mariana M6.1 — smallest event attempted; Hartzell-style
                   teleseismic FFI where USGS publishes no finite fault (executed)
              30 = BLIND TEST: 2026 Colombia M7.4, 110 km deep — emergent 35-s nucleation
                   discovered; prediction frozen before the USGS finite fault (executed)
              31 = deep dive: the 35-s two-stage nucleation of the Colombia M7.4 — raw-data
                   evidence, the diagnostic inversion failure, the centroid re-anchoring fix,
                   and the foreshock-cascade-vs-slow-slip discrimination (executed)
reports/      self-contained PDF reports for both case studies (tectonic context, methods
              primer for non-specialists, full inversion narrative, reproducibility appendix)
```

## Quickstart

1. Install [MTUQ](https://github.com/mtuqorg/mtuq) (conda env `mtuq`, with `mpi4py`) and the
   USGS [neic-finitefault / WISP](https://github.com/usgs/neic-finitefault) (env `ff-env`).
   Neither is vendored here.
2. Copy `scripts/config_template.py` to your project root as `config.py`; edit paths.
3. Open `notebooks/23_mtuq_wisp_cookbook.ipynb` and follow it top to bottom for your event:
   USGS metadata check → MTUQ data prep → full-MT / DC / DC+depth searches → VR + Kagan
   comparison → WISP data fetch → auto model both planes → fair refinement → figures.

## The rules the hard way taught us

| situation | rule (and where it came from) |
|---|---|
| any event | tag data `type:displacement`; remove response with `water_level=None`; keep data ≥2 Hz for syngine (nb19–21) |
| any event | compare mechanisms with the Kagan angle via `pyrocko.moment_tensor.kagan_angle` — never plane-by-plane angles, never hand-rolled Kagan code |
| any event | evaluate the *reference* (GCMT/USGS) mechanism's VR on your own data — if the accepted answer fits terribly, your protocol is broken, not the Earth |
| deep events (≳60–70 km) | the default 10–33 s body band mixes P/pP/sP and fails *confidently* (~90° mechanism rotation). Move the body band up: `--bwband 25,60,120,10` (validated: body VR −80% → +72%). Depth must be read from the body-only misfit-vs-depth curve; the surface-wave band has no depth resolution |
| deep events + WISP | the surface-wave GF bank ends near 126 km; a fault extending below it makes WISP **silently** drop all surface waves for that plane, faking a plane discrimination. The refine script auto-trims the fault; always grep the pipeline log for "Maximum depth" |
| plane discrimination | refine BOTH nodal planes with identical wave types, channels, QC and alignment before comparing misfits. Whether teleseismic data can discriminate is a function of rupture size: a ~10-km rupture gave a 1.6% near-tie (2026 Calama, nb24); a ~50-km rupture gave a decisive 10% preference (2024 M7.4, nb25). For compact ruptures argue planes from regional precedent / near-field data |
| near-nodal stations | unmodelled clean first motions at a few stations can discriminate mechanisms less than 10° apart — check the radiation coefficient before blaming the data |
| interfering events | a large earthquake within ~2 h contaminates long teleseismic windows; compute per-station group-velocity overlap windows and zero swept channels (nb22 method, empirically validated) |
| depth scans | a minimum on the EDGE of the searched range is not a minimum; MTUQ misfit arrays are (sources, origins) — origins LAST — when reshaping |
| MPI | `grid_search` returns None on non-root ranks; guard before combining results |
| Mw ≥ 7 | magnitude-dependent defaults silently disable body waves; an explicit `--bwband` now re-enables them (a "joint" run that was secretly SW-only was caught because the waveform figure had no body panels — always look at the fits; a missing panel is a missing dataset) |
| any FFM | check slip at the fault edges (>~15–20% of peak in an edge row = domain too small); re-invert on an enlarged fault — extended *identically* on both planes if a plane comparison is at stake (nb25: edge slip 41% → ≤16%, conclusion unchanged) |

## Case studies

**2026 Calama M6.9 (us6000t04s), 109 km deep intraslab, N Chile** — `notebooks/24`,
`reports/atacama_pair_2024_2026/report.pdf` (Part I; the same PDF covers both Chile events
plus the comparative Part III with the supershear-mechanics discussion). Final MT 194/16/−83 (Kagan 8.4° from USGS Mww, better VR
than the reference on both wave types); slip model: thin subhorizontal lens, 0.70 m peak at
105–113 km, ~15 s, subshear; centroid ~105–115 km from four depth-phase estimators; shallow
nodal plane weakly preferred (fit + Tarapacá-2005 precedent). No USGS finite-fault model
exists for this event.

**2024 San Pedro de Atacama M7.4 (us7000n05d), 127 km deep intraslab, N Chile** —
`notebooks/25`. The published rupture models disagree on the fault plane (USGS finite fault:
shallow 172/21; Jia et al. 2025 Nat. Comms., with strong motion + GNSS + aftershocks: steep
341/71). Our teleseismic-only pipeline: MT 176/20/−70 first-try under both protocol variants
(conjugate 7.8° from Jia's plane); with the hypocentre below the 125-km surface-wave GF cap
on *both* planes, the body-only two-plane comparison is fair by construction — and the
**steep plane wins by 10%** (0.1099 vs 0.1224), its slip at 81–170 km independently
reproducing Jia's 120–180 km subevent depths. Together with the 2026 near-tie, the pair
brackets the rupture size at which teleseismic plane discrimination turns on.

**2026 Sanriku-oki M6.9 (us6000t7zq), Japan Trench interface** — `notebooks/22`,
`reports/sanriku_2026/report.pdf`. Same-magnitude shallow counterpart: MT robust (centroid
depth 46 km recovered independently); FFI marginal exactly as the duration/period criterion
predicts; quantifies how the Venezuela M7.2+M7.5 doublet's wavetrain swept 69% of the
surface-wave windows — a plausible operational explanation for the missing USGS finite-fault
product.

**2026 Kumamoto M6.8 (us6000tgb9), shallow strike-slip, Japan — blind test, frozen
2026-07-28 13:19 UTC, 3.2 h before USGS posted their finite fault** — `notebooks/28`,
`reports/kumamoto_2026/report.pdf` + `prediction_frozen.md`. Prediction: right-lateral
212/75/−164 (Futagawa–Hinagu trend, ~1% misfit preference over the conjugate), Mw 6.75,
single compact asperity (peak 1.9 m, 2–16 km depth, no shallow slip maximum), T5-95 ≈ 14 s,
Vr ≈ 2.1 km/s — a smaller sibling of the 2016 Kumamoto rupture. **Outcome: every headline
verified** — USGS chose the identical plane (212/75/−164), Mw 6.754 vs our 6.75, peak slip
2.01 vs 1.88 m at the same place and depth, and the moment-rate functions overlay almost
perfectly including the secondary pulses (caveat: WISP is the USGS code family, so this
validates data handling/seeding/QC, not an independent method; Vr, weakly constrained, differs
2.05 vs 2.46 km/s). Also documents a new strike-slip failure mode: an azimuthally clustered
25-station subset aliased the four-lobed radiation pattern (surface-wave solutions locked ~90°
rotated, Kagan 67–78°); a 60-station azimuth-balanced set snapped the joint DC back to Kagan
27° — for strike-slip, azimuthal *balance* beats station count, and the shallow-event
"body-only" rule (a dip-slip cure) does not transfer.

**2023 Mariana Islands M6.1 (us6000kzv3), oblique normal/strike-slip — resolution test at a
magnitude where USGS publishes no finite fault** — `notebooks/29`,
`reports/mariana_2023/report.pdf`. Following Hartzell, Mendoza & Zeng (2013, GRL; teleseismic-P
FFI of the M5.8 Mineral, Virginia earthquake): P waves to 1 Hz, 2.4×2.0 km subfaults, 205
stations (gap 28°). **Resolved** (plane-stable, physically checked): a single compact ~10-km
asperity at 2–13 km depth centred on the hypocentre, the STF including its two-pulse structure
(cross-plane r=0.98), M0 to 1% of the W-phase, Vr≈2.0 km/s, stress drop 1.0–3.7 MPa. **Not
resolved**: the fault plane (2.2% margin — consistent with the rupture-size criterion),
independent point-source mechanism (MTUQ Kagan 52–86°) and centroid depth (misfit flat to 0.2%
over 6–33 km) — at M≲6.5 take both from the W-phase and spend the teleseismic data on the
finite fault. Extends the validated range from M6.8 down to M6.1.

**2026 San José del Palmar (Colombia) M7.4, 110 km deep — blind test, frozen 2026-08-11,
USGS posted 2026-08-14: VERIFIED** — `notebooks/30`, `reports/colombia_2026/report.pdf` +
`prediction_frozen.md`. Not a Bucaramanga-nest event: it belongs to the Cauca
intermediate-depth cluster of the Nazca slab, 407 km SW of the nest. Discovery: a **~35-s
weak emergent nucleation phase** (W-phase centroid time = origin + 34.5 s, 29 km S) makes
origin-anchored inversions kinematically impossible (slip piles at the fault edge; MTUQ
windows miss the main energy) — both analyses re-anchored at the W-phase centroid in space
and time. Frozen prediction: plane **226/72/33** (8.2% margin + conjugate edge violation =
near-decisive at this rupture size), single asperity ~45 km SSW of the epicentre, up-dip,
85–125 km depth, peak 1.6 m, main-phase T5-95 ≈ 19 s, Δσ ≈ 2.9 MPa, M0 to 1% of W-phase.
MTUQ point source documented as breaking on this complex source (mechanism stable nowhere
except ~90 km) — for complex sources, mechanism+depth from the W-phase, teleseismic data
into the finite fault. Also hardened `wisp_shift_match.py` for deep faults (body-only
fallback below the surface-wave GF cap). **Outcome (USGS posted 4 days after the freeze):
every physical headline verified — identical plane (225.57/71.78/33.23), M0 to 2.5%, same
SSW asperity at ~100 km, same unilateral SW directivity (46 km @ az 225 vs our 27 @ 220),
same main-phase duration and speed (2.34 vs 2.52 km/s, same metric) — and the methodological
prediction too: USGS also abandoned the cataloged hypocentre (relocated 24 km NE, then let a
single front crawl silently ~30 s; their STF has 0.2% of moment in the first 25 s,
independently confirming our stage-1 bound). Honest differences: their 6-km subfaults give a
compacter asperity (peak 3.98 vs 1.63 m; Δσ 8.1 vs 2.9 MPa), and the initiation point —
catalog vs NE relocation vs W-phase centroid — remains the ill-determined element.**

## Data & code availability

Waveform data are open (IRIS/EarthScope, GEOFON et al.) and fetched by the scripts
(`wisp_fetch_fallback.py` replaces the deprecated `ffm get-data teleseismic` discovery).
Green's functions: syngine (ak135) for MTUQ; WISP's own banks for the FFI. No waveforms or
Green's functions are stored in this repository.

MTUQ: Modrak et al.; uses the cut-and-paste strategy of Zhu & Helmberger (1996).
WISP / neic-finitefault: USGS NEIC operational finite-fault code (Goldberg et al.).
Please cite the upstream packages when using this pipeline.

## License

MIT (this repository's scripts and documents). Upstream packages carry their own licenses.
