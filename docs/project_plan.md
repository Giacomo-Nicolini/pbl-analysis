# Project plan: PBL height from EC data and its impact on footprint estimates

Status: agreed design, v0.1 (2026-10-08). Living document: update as decisions change.

## Research questions

1. How well can the PBL height be estimated from eddy-covariance data (PBLHmod) at ICOS stations?
2. How do PBLHmod and ERA5 BLH compare with measured PBLH, and with each other?
3. In which conditions, and why, do the estimates disagree?
4. What is the impact of PBLH errors on FFP (Kljun et al. 2015) footprint estimates, per half-hour and as climatologies?
5. Practical outcome: a recommendation for ICOS/ETC on which PBLH source to use, and when.

## Scope

- Sites: ICOS ecosystem stations only, for now.
- Language: R (a Python port is possible later; keep the logic language-agnostic and well documented).
- Execution: developed locally (Windows, VS Code), run on the CMCC cluster.

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

## ERA5 representativeness

1. Site metadata includes: ERA5 cell centre and distance to site, site elevation vs ERA5 model orography, land–sea fraction of the cell, subgrid orography std (ERA5 `sdor`), and a terrain-complexity class.
2. If neighbouring cells are available (file naming `lon±n__lat±n` suggests so; to check with the cluster provider): spatial spread of the 3×3 cells as a per-timestep representativeness metric. Otherwise bilinear interpolation of the 4 nearest cells.
3. Use these as covariates in phase 4, and flag complex sites (mountain, coastal) so headline statistics are given with and without them.

## Evaluation and error attribution (phases 3–4)

- Errors: absolute and log-ratio (log(mod/ref)), since PBLH spans orders of magnitude.
- Metrics: bias, MAE, RMSE, correlation; Bland–Altman per regime and season; timing metrics (morning growth rate, time of peak, daily maximum).
- Triple collocation (log space) where three sources are available; state the "same quantity" assumption.
- Drivers: stability class, time of day and transitions, season, cloud/precipitation, wind, H magnitude, land cover, terrain, latitude, measurement height, time since last data gap.
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
- ICOS L2 METEO: `TA` and `TA_1..10` profile, `P`, `SW_IN`. Optional; the temperature profile may help at tall towers.
- ERA5 `boundary_layer_height`: hourly, yearly CSV per station (nearest cell).
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
- ERA5 for SE-Myc must be re-extracted with the corrected coordinates.
- Operational L2 FLUXES contain implausible USTAR (max 49 m/s) and PBLH (max 91 km) at BE-Bra: report to ETC QC.
- METEO `PA` unit: kPa at BE-Bra (median 101.1), as assumed in spec 01.
- ERA5 extraction details (nearest cell? downscaling? QC?), and availability of neighbouring cells.
- Source of zm, d, z0 and target area per site.
- Measured PBLH datasets: to be collected; open data candidates: ACTRIS/Cloudnet, E-PROFILE, IGRA radiosondes.
- Literature: the current bib is footprint-focused. A PBL collection is needed (B&G 1991, Nieuwstadt 1981, Zilitinkevich, Seidel et al. 2010/2012, ERA5 BLH evaluations, ceilometer/lidar MLH methods).
