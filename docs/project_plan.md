# Project plan: PBL height from EC data and its impact on footprint estimates

Status: agreed design, v0.1 (2026-10-08). Living document: update as decisions change.

## Research questions

1. How well can the PBL height be estimated from eddy-covariance data (PBLHmod) at ICOS stations?
2. How do PBLHmod and ERA5 BLH compare with measured PBLH, and with each other?
3. In which conditions, and why, do the estimates disagree?
4. What is the impact of PBLH errors on FFP (Kljun et al. 2015) footprint estimates, per half-hour and as climatologies?
5. Practical outcome: a recommendation on which PBLH source to use, and when.

## Scope

- Sites: ICOS ecosystem stations only, for now.
- Language: R (a Python port is possible later; keep the logic language-agnostic and well documented).
- Execution: developed locally (Windows, VS Code), run on the Unitus cluster.
- Sources compared (G, 2026-10-23): the evaluation (phases 3–4) compares **PBLHmod, ERA5 and measured PBLH only**.
  The operational ICOS PBLH (`icos_oper`) is used only in the phase-1 regression (operational vs PBLHmod v2,
  spec 03); it is not part of the three-way comparison, its figures or its tables. The FI-Hyy spec-07 results that
  include it are kept as a reference test.
- ERA5 in figures (G, 2026-10-23): tower cell only (`era5_blh_center`). The bilinear values (`era5_blh_bilinear`)
  appear in the numerical results (tables) only; at FI-Hyy the two differ very little (convective pair vs measured:
  bias +111 vs +113 m, median log ratio +0.32 both), which is itself a result to report.

## Phases

| Phase | Content | Output |
|---|---|---|
| 0 | Site metadata; standard PBLH format spec; readers for ICOS L2 (FLUXES, FLUXNET, METEO) and ERA5 | `R/read_*.R`, `docs/specs/standard_format.md` |
| 1 | PBLHmod v2 (Kljun App. B + B&G slab), with options and unit tests; regression vs legacy (PBLH column in L2 FLUXES) | `R/pblh_*.R`, `tests/` |
| 2 | Standardizers for reference data (ERA5, measured, PBLHmod) into standard files | `data/standard/` |
| 3 | Evaluation: PBLHmod vs ERA5 at all sites; three-way comparison where measurements exist; triple collocation | stats, figures |
| 4 | Error attribution: conditions and causes of disagreement | GAM models, figures |
| 5 | Footprint error propagation (Monte Carlo; half-hourly and climatologies) | sensitivity results |
| 6 | Paper (Quarto) and data/code release (Zenodo) | `paper/` |

## Phase 1: PBLHmod v2

- Base: Kljun et al. (2015) App. B (notes: `docs/references/kljun2015_appB_notes.md`). Diagnostic formulas for
  neutral and stable conditions; Batchvarova & Gryning (1991) slab model for convective conditions, integrated as a
  daily cycle from sunrise (no carry-over between days), with internal sub-steps.
- Selectable options (the comparison of options is itself a result):
  - stable/neutral formula: Nieuwstadt (1981) Eq. B1 (default, as in Kljun App. B) vs Zilitinkevich Eq. B4
    (c_s = 0.4 / 0.63 / 0.3; legacy used 0.4) vs neutral Eq. B2;
  - γ (free-atmosphere lapse rate, B&G): **optional input**. Three levels, from simplest to richest:
    1. constant (**default 0.01 K/m**, Batchvarova & Gryning 1991 as cited in Kljun 2015; legacy 0.0025 K/m run as
       sensitivity case) — EC-only, always available;
    2. climatological lookup table (site × month), computed once from ERA5 pressure levels and shipped with the code —
       operationally still EC-only;
    3. time-varying γ from ERA5 pressure levels (hourly/daily) — requires ERA5 download.
       Data note (2026-10-15): the pressure levels received are 6-hourly with only 1000/925/850 hPa below ~1.5 km.
       Sufficient for option 2 (e.g. θ gradient 850–700 hPa, site × month); option 3 would need hourly data and the
       intermediate levels (975, 950, 900, 875, 825, 800, 775, 750 hPa). Established recipe for γ from profiles:
       extrapolate the free-troposphere θ gradient above the mixed-layer top from the morning profile (`Wouters2019`).
    The function uses the richest γ provided and falls back to the constant. Comparing 1 vs 2 vs 3 is an explicit
    result: how much accuracy does the EC-only method lose without external data? (Supports the "simplicity" goal.)
  - convective trigger: H > 0 vs 1/L class (must be consistent);
  - H forcing: `H_F_MDS` (gap-filled, default; values with `H_F_MDS_QC` above a threshold, default > 1, treated as
    missing because poor gap-filling can be very bad) vs `H` (FLUXES) vs `H_CORR` (energy-balance corrected,
    FLUXNET). The closure-corrected run is an exploratory sensitivity test, not based on the assumption that EC
    underestimates H.
- B&G needs continuous forcing from sunrise: forcing gaps ≤ 1 h are interpolated; longer gaps stop that day's
  integration (flagged); restart at the next sunrise.
- Model output is **raw**: no smoothing, no clipping (in particular no 50 m floor), no gap-filling of PBLH.
  Post-processing (smoothing, clipping, filling) may be added in a later stage, only if the evaluation shows it is
  needed, as a separate step that keeps the raw output.
- Spec 02: the model. Spec 03: legacy regression — run the legacy script on BE-Bra, compare with the operational
  `PBLH` column (which has an unexplained 50 m floor, probably applied after `calc_PBL_height()`) and with PBLHmod,
  then switch each legacy fix on one at a time to quantify its effect (see `legacy_review.md`).

## Standard PBLH format (all sources: PBLHmod, ERA5, measurements)

Data table per site × source:

- Required: `timestamp_start_utc`, `timestamp_end_utc`, `pblh_m`, `source_id`
- Recommended (NA allowed): `qc_flag`, `uncertainty_m`, `n_samples`

Metadata sidecar (YAML) per file. **Every field except `site_id`, `source_id` and `method_class` may be `unknown`**:

- `method_class`: `turbulence` (Doppler lidar TKE/σw, sodar), `aerosol` (ceilometer, backscatter lidar), `thermodynamic` (radiosonde, Ri_b, microwave radiometer), `model` (ERA5), `parameterisation` (PBLHmod), `unknown`
- `instrument`, `retrieval_algorithm`, `height_reference` (agl / above instrument / asl), `instrument_altitude_m`, `time_reference`, `original_resolution`, `original_qc`, `processing_history`, `contact`, `licence`, `citation`
- `metadata_completeness`: computed automatically (fraction of known fields)

## Comparing different method classes

Different classes measure different "tops". The approach:

1. **Regime where definitions converge**: in a well-developed daytime convective BL, turbulence-, aerosol- and Ri-based tops roughly coincide. This is the main ground for direct comparison of all sources.
2. **Regimes where they diverge** (night, morning and evening transitions): compare only against matching classes (e.g. PBLHmod vs turbulence-based). Treat other differences as expected by definition (e.g. aerosol-based ≥ turbulence-based at night because of the residual layer), and report them as such.
3. **Co-located instruments of different classes** (e.g. Cloudnet sites with ceilometer and Doppler lidar): quantify the inter-method gap itself, as a "definitional uncertainty" band against which PBLHmod and ERA5 errors are judged.

The class is therefore not a stratum to compare against other classes; it tells which comparison is valid in which regime.

**Distance rule for measured PBLH** (G, 2026-10-22; applied in code via `distance_km` in `config/measured_sites.csv`
and `distance_class()`): < 5 km all regimes; 5–25 km daytime convective only; 25–50 km daytime climatologies only;
> 50 km not used. BE-Bra has no usable measurement within 50 km (Uccle ≈ 58 km, Diepenbeek ≈ 74 km, Cabauw ≈ 78 km).

**Censored measurements** (spec 07, G 2026-10-23). A turbulence-based MLH from Doppler-lidar ε has three kinds of
values: exact (qc 0/2), censored below the lowest reliable gate z_min (qc 1, value in (0, z_min]) and signal-limited
(qc 3, the search ended at a gate without signal: a lower bound, value in [MLH, ∞)). Metrics (bias, log ratio, …)
are computed on exact values only; qc 1 and qc 3 are summarised separately (share of model values ≤ z_min, share
below the bound and median model/bound ratio). Diurnal cycles use censored medians (Turnbull estimator,
`survival::survfit`, `Surv(type = "interval2")`) next to the plain median of exact values. Because the deep daytime
layers are mostly signal-limited in clean air, the exact values over-represent shallow layers: metrics on exact
values alone are biased towards "models too deep".

## ERA5 representativeness

1. Site metadata includes: ERA5 cell centre and distance to site, site elevation vs ERA5 model orography, land–sea fraction of the cell, subgrid orography std (ERA5 `sdor`), and a terrain-complexity class.
2. Neighbouring cells are available (3×3 around the tower cell, confirmed 2026-10-15): spatial spread of the 3×3 cells as a per-timestep representativeness metric. Bilinear interpolation of the 4 nearest cells to the tower position as a sensitivity case.
3. Use these as covariates in phase 4, and flag complex sites (mountain, coastal) so headline statistics are given with and without them.

## Evaluation and error attribution (phases 3–4)

- Errors: absolute and log-ratio (log(mod/ref)), since PBLH spans orders of magnitude.
- Metrics: bias, MAE, RMSE, correlation; Bland–Altman per regime and season; timing metrics (morning growth rate, time of peak, daily maximum).
- Triple collocation (log space) where three sources are available; state the "same quantity" assumption.
- Drivers: stability class, time of day and transitions, season, cloud/precipitation, wind, H magnitude, land cover, terrain, latitude, measurement height, time since last data gap.
  Candidate ERA5-based drivers (proposed 2026-10-15, to discuss): ERA5 cloud base height (cloud presence and
  low-cloud capping of the CBL); the difference between ERA5 and EC sensible heat flux (does an ERA5 BLH error follow
  an ERA5 surface-forcing error?); 3×3 BLH spread.
- Methods: binned statistics, then GAMs (`mgcv`).

## Footprint error propagation (phase 5)

- Model: official Kljun FFP R code (`calc_footprint_FFP.R`, `calc_footprint_FFP_climatology.R`).
- Monte Carlo: perturb h with error distributions conditional on stability and site, from phases 3–4.
- Express results using dimensionless zm/h and zm/L.
- Metrics: change in fetch (50/70/80/90 %), in the footprint fraction inside the ICOS target area, and in the number of records crossing FFP validity limits (zm < h, zm/L, u*).
- Show explicitly that a bias in h does not average out in climatologies.

## Inputs

- ICOS L2 FLUXES: `H`, `USTAR`, `MO_LENGTH`, `ZL`, `T_SONIC`, `WS`, `WD`, `V_SIGMA`, `PBLH` (legacy operational output), `FETCH_*`, `FOOTPRINT_*`, QC flags.
- ICOS L2 FLUXNET: `TA_F`, `PA_F` (kPa), `H_F_MDS`, `H_CORR`, `NIGHT`, `SW_IN_POT`, `WS_F`, `P_F`.
- ICOS L2 METEO: `TA` and the `TA_<n>` profile (station-specific number of levels, detected from the file header:
  BE-Bra 10, FI-Hyy 9), `P`, `SW_IN`. Optional; the temperature profile may help at tall towers.
- Measured PBLH (spec 07): ACTRIS/Cloudnet (CC BY 4.0, API https://cloudnet.fmi.fi/api), pairs in
  `config/measured_sites.csv`: FI-Hyy–hyytiala (0.56 km, Doppler lidar HALO 46, `epsilon-lidar` 2022-01 → 2024-04),
  FI-Lom–kenttarova (1.81 km, ε 2022-08 → 2023-12 and from 2026-04), DE-RuS–juelich (5.26 km: daytime convective only,
  ε from 2023-04), SE-Nor–norunda (0.06 km, CL51 ceilometer only: aerosol, no MLH product; primary route PI data),
  FR-Gri–palaiseau (23.56 km: daytime convective only, ε from 2025-01). Inventory: `scripts/07a_cloudnet_inventory.R`.
  PI data requested by G for the atmosphere stations next to SE-Htm, SE-Nor, SE-Svb, FI-Hyy, CZ-BK1 (list and contacts
  kept locally in `docs/local/sites_PBLH_measure.txt`).
- ERA5 (0.25°, tower cell + 8 neighbours, 2020-01-01 → 2026-10-07; formats and conventions in `CLAUDE.md`):
  - hourly `boundary_layer_height` (instantaneous; bulk Richardson diagnostic, Ri_c = 0.25), `cloud_base_height`,
    `surface_sensible_heat_flux` (hourly accumulation, J m-2, positive downward);
  - 6-hourly pressure levels (13 levels, 50–1000 hPa): temperature, geopotential;
  - static: surface geopotential, land–sea mask, sdor.
- Site metadata: `config/sites.csv` (static; `zm_m` = latest value, `canopy_height_m` = median of measurements)
  + `config/sites_timevarying.csv` (long format; `zm_m` as step periods, `canopy_height_m` as point measurements).
  Both derived from `config/sites_H_EC_canopy_formBIF.csv` (ETC/BIF source, do not edit by hand) and ETC ancillary files.
  Rules for the code: zm is a step function (before the first record use the first value, and flag it);
  canopy height is interpolated linearly between measurements and held constant outside the measured range.
- **Exclusions** (must be applied by every reader/analysis step, never silently skipped):
  - `config/exclusion_periods.csv`: hard, date-based exclusions (`site_id, date_start, date_end, scope, reason`).
    Currently: SE-Nor 2022-10-01 → 2023-04-20 (clearcut, tower lowered 36 → 12 m).
  - Tall canopy relative to zm (e.g. IT-BCi maize seasons): **phases 1–4** keep the data but add a per-half-hour flag
    and the ratio hc/zm as a covariate in the error analysis; **phase 5** applies the FFP validity criteria.
  - FFP validity criteria (Kljun et al. 2015): implemented with **configurable thresholds**, at three strictness levels:
    `strict` (as published), `relaxed` (as used operationally at ETC; in past work some criteria, notably the
    roughness-length and h/zm ones, had to be relaxed or nearly ignored), `none`. Report how many records each level
    keeps per site, and whether the phase-5 conclusions change with strictness. The exact published thresholds and the
    ETC relaxed values must be written in the phase-5 spec (to be provided by G. Nicolini / verified in the paper).
- ICOS target areas: KML polygons (ETC), one per station, stored in `data/raw/target_areas/XX-Yyy_TA.kml` (38 stations; some are multi-polygon, e.g. IT-Cp2, IT-Ren)
  (gitignored, not published for now). Read with `sf`; for FFP overlap, transform to a local
  metric projection centred on the tower (e.g. azimuthal equidistant), matching FFP x/y coordinates in metres.
  Needed in phase 5 only.
- Missing values: `-9999`.

## Open points

- Site metadata to verify: canopy-height jumps at FR-Fon, FR-Hes, SE-Svb (use current values for now);
  SE-Nor IGBP after clearcut (not used in calculations).

- Time zone and timestamp convention of ICOS L2 files (UTC vs local standard time) vs ERA5 (UTC, instantaneous at the hour).
  Spec 01 check (BE-Bra, 2026-10-12): after converting with UTC+1, the mean diurnal METEO `SW_IN` peaks at
  11:15–11:45 UTC (solar noon ≈ 11:42 UTC), consistent with local standard time. To confirm at other sites.
- ~~Units of `T_SONIC` in L2 (°C or K).~~ Resolved (spec 01, BE-Bra): `T_SONIC` is in K (median 290 K). Readers
  now rename only, and abort if a column's median is implausible for its declared unit, so other sites are checked
  automatically.
- ~~ERA5 for SE-Myc must be re-extracted with the corrected coordinates.~~ Done (confirmed 2026-10-16).
- Operational L2 FLUXES contain implausible USTAR (max 49 m/s) and PBLH (max 91 km) at BE-Bra: report to ETC QC.
- METEO `PA` unit: kPa at BE-Bra (median 101.1), as assumed in spec 01.
- Polar-day convective rule (e.g. start at solar midnight) for GL-ZaF, SE-Sto, FI-Lom: spec 02 runs diagnostic only on
  polar days (`flag_polar_day`).
- ~~Operational PBLH has a 50 m floor and the exact value 1654.6 m in ~13% of rows at BE-Bra: investigate in spec 03.~~
  Explained (2026-10-16): `legacy/proc_PBL_height.R` floors at 50 m and fills NA with the day/night mean of the
  processing run (reproduced exactly at BE-Bra; see `STATUS.md`). Implication for phase 5 and for the paper: the
  operational PBLH used for ICOS footprints contains ≈ 18 % run-mean fills; the comparison must flag these rows.
- ~~ERA5 extraction details (nearest cell? downscaling? QC?), and availability of neighbouring cells.~~ Resolved
  (2026-10-15): nearest 0.25° grid point, no downscaling, 3×3 neighbours delivered (BE-Bra: cell centre 6.5 km from
  the tower). SE-Myc corrected coordinates confirmed (2026-10-16).
- Spec 04, first PBLHmod vs ERA5 comparison (BE-Bra, 2026-10-17), points to discuss:
  - Operational fill detection uses a per-calendar-year rule (spec 04 amendment 7); to be checked at other sites,
    where processing runs may not align with calendar years.
  - EC H (PBLHmod forcing) exceeds ERA5 H at midday in MAM/JJA by ≈ 90–100 W m-2 (mean diurnal peak ≈ 215 vs
    ≈ 125 W m-2), yet PBLHmod convective h is close to ERA5 (median log ratio −0.07). With γ = 0.01 K m-1 the larger
    H is offset by a γ about twice the ERA5 850–700 hPa θ gradient (monthly medians 0.0044–0.0051 K m-1). The
    850–700 hPa layer lies mostly above the BE-Bra PBL, so it is only context for γ (options 2–3).
  - ~~`mad_log_ratio`: unscaled or scaled MAD?~~ Resolved (2026-10-21): scaled MAD (× 1.4826, sd-consistent;
    spec 04 amendment 8).
  - Pair 5 (PBLHmod vs operational PBLH) is compared on half-hours, not on ERA5 hours.
  - ~~PBLHmod values > 5 km at BE-Bra (179 half-hours, 70 > 10 km) come from implausible L2 USTAR (median 6.9 m/s,
    max 49 m/s; 1,726 rows > 2 m/s), across B1/B3/B5. Decide on an input plausibility filter for the forcing (e.g.
    u*/WS limit) — input QC, not output clipping.~~ Resolved (spec 05, 2026-10-21; spec 02 amendment 12): u* is
    rejected if u*/WS > 0.5 (WS ≥ 1 m s-1) or u* > 3 m s-1; u* and L set to NA, short gaps interpolated. BE-Bra:
    2,149 of 100,186 u* values rejected (ratio 1,328, cap 46, both 775; 618 interpolated, 1,531 left NA);
    PBLHmod > 5 km 179 → 11, > 10 km 70 → 0. Still to do: **report the implausible L2 USTAR to ETC QC** (on 99 % of
    the rejected rows L2 `TAU` is missing, so USTAR cannot be cross-checked there).
  - u* filter, low wind (spec 05, BE-Bra): 577 ratio rejections at WS 1–2 m s-1 in convective conditions. Kept as
    rejected: σw/u* 0.88 on those rows vs 1.21 on kept rows suggests an overestimated u*, not free convection.
    Revisit with measured PBLH.
  - ~~The 11 PBLHmod values > 5 km left after the u* filter: 10 come from B5 on 2021-08-18 and 2022-03-11, where L2 H
    is QC-rejected but MO_LENGTH (computed from H_UNCLEANED ≈ 0, |L| ≈ 3e6 m) is kept, while the forcing uses the
    gap-filled H_F_MDS. L2 L vs L from the forcing: median ratio 0.97, but 1,237 convective rows differ by more than
    a factor 3. To be addressed in spec 06 (L computed from the forcing).~~ Resolved (spec 06, option C with
    amendment 1, 2026-10-21): L2 L checked against the forcing L where the heat flux is positive; 26,652 rows checked,
    2,476 replaced (sign 1,349, magnitude 1,127) on 956 days; values > 5 km 11 → 1; pair 1 RMSE 349 → 338 m.
  - Remaining PBLHmod value > 5 km after spec 06 (BE-Bra): 2024-01-05 09:30, B1 stable, u* 2.67 m s-1 (u*/WS 0.43,
    WS 6.2 m s-1), H missing so L not checked. Known limitation of the spec 05 u* filter (passes both tests).
  - u* of 2–3 m s-1 passes the filter at BE-Bra but is suspect (V2 in spec 06 amendment 1: a cap of 2 m s-1 removes
    all outliers): include in the report to ETC QC.
  - ~~285 night-time L replacements under spec 06 amendment 1 (Nov 2022–Jan 2023), on rows with a positive
    gap-filled H_F_MDS (QC ≤ 1) while H_UNCLEANED is negative.~~ Handled in the model by spec 06 amendment 2 (L check
    in daytime only; 0 night replacements). Still in the ETC QC report list (below).
- Source of zm, d, z0 and target area per site.
- Measured PBLH datasets: Cloudnet ε MLH at FI-Hyy done (spec 07); other Cloudnet sites next (cluster); further open
  data candidates: E-PROFILE, IGRA radiosondes; PI data requested (see "Inputs").
- Spec 07 at FI-Hyy, points to discuss: (a) "shallow too deep / deep too shallow": where the lidar MLH is below
  105 m, PBLHmod is ≤ 105 m on only 26 % of half-hours (ERA5 24 %); where it is only a lower bound, PBLHmod is below
  the bound on 40 % (JJA 62 %, JJA stable 85 %; ERA5 33 %, JJA 45 %); candidate finding; (b) half of the measured
  profiles are signal-limited (clean air, low SNR at SMEAR II), so JJA daytime censored medians are not estimable;
  (c) Cloudnet classification "drizzle/rain" near the ground on 15 % of profiles vs `rain_detected` 3.5 %: default
  precipitation flag = `rain_detected` (G), the other rules kept as sensitivity; (d) ±31 m ground reference (Cloudnet
  altitude 150 m vs ICOS 181 m; instruments ≈ 200 m from the main station, Sinclair et al. 2022).
- Timing of EC vs ERA5 sensible heat (evaluation; possible error driver for the evening transition). EC timestamps
  are correct (SW_IN_POT and EC H centred on solar noon within 4 min, Cowork 2026-10-23). With every value at its
  interval midpoint (EC start + 15 min, ERA5 accumulation timestamp − 30 min), ERA5 H is earlier than EC H:
  best-fit lag FI-Hyy 25–37 min, BE-Bra 16–25 min (`output/tables/07_h_timing_<site>.csv`, 2020–2025); at FI-Hyy in
  JJA the afternoon decline is ≈ 63 min earlier, the morning rise ≈ 15–18 min (asymmetric: physical, consistent with
  canopy heat storage at the tall forest, not a time convention). The L2 air-storage term `SH` reduces the lag
  (FI-Hyy 12–27 min, BE-Bra −1 to 14 min) but not the FI-Hyy afternoon asymmetry (51–64 min); biomass/canopy storage
  is not in `SH`. No change to the PBLHmod forcing for now; to examine with the evening-transition errors (PBLHmod
  decay timing vs ERA5 and measured MLH).
- Cloudnet ε files at Hyytiälä are volatile/experimental (2026-10); before paper results, re-check checksums against
  the API, re-run if changed, and create a Cloudnet collection with a persistent identifier for the files used
  (spec 07).
- ~~Literature: the current bib is footprint-focused. A PBL collection is needed.~~ Done (2026-10-15):
  `pblh_collection.bib` (Kljun App. B references) + `pblh_literature_2026-10.bib` and annotated
  `literature_review_pblh.md`. Keep extending as papers are added to `C:\WORK\bibliography_new_win11\`.

## Points for the ETC QC report

Kept in the local (untracked) file `docs/local/etc_qc_report.md`: one entry per L2 data issue (site, period,
numbers, how PBLHmod handles it, status).

## Positioning in the literature (from the 2026-10-15 search; details in `docs/references/literature_review_pblh.md`)

- The Kljun App. B PBLH is already used network-wide without evaluation: ICOS ETC FFP tools (with optional ERA5
  override), FLUXNET2015 footprint climatologies (`Fang2024`), AmeriFlux footprint representativeness (`Chu2021`),
  Campbell EasyFlux (`Zhou2026`; it calls γ the "dry adiabatic lapse rate (commonly 9.8 × 10⁻³ K m⁻¹)", which is
  numerically close to our 0.01 K m⁻¹ default, but in B&G γ is the free-atmosphere potential-temperature gradient).
  `Rebmann2018` states ICOS uses Kljun 2015.
- Closest prior work: PBLH vs surface fluxes at a few US sites, daytime only (`Beamesderfer2023`, `Chu2026`); ERA5 vs
  radiosonde stratified by EC stability at one ICOS site (`Sinclair2022`, FI-Hyy). No study found evaluates EC-derived
  PBLH against ERA5 and measurements across a network and all stabilities. `Helbig2021` names this gap.
- Expectations to test: ERA5 too deep in very stable conditions and evening transitions (`Sinclair2022`,
  `Arganaraz2026`); the sign of the daytime ERA5 bias varies between studies in Europe; surface-flux stable formulas
  scatter widely (`Vickers2004`, `Steeneveld2007`); slab models overestimate without subsidence/advection
  (`Pietersen2015`); FFP sensitivity to h has only been tested as ±20 % (`Kljun2015`).
- Reference choice per class: turbulence-based (Doppler lidar, `Manninen2018`) for EC-derived PBLH; radiosonde bulk Ri
  for ERA5; ceilometer MLH only late morning to mid-afternoon (`Schween2014`: aerosol-based height +300 to +600 m vs
  Doppler lidar at transitions).
