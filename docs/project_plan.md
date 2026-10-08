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

- Base: Kljun et al. (2015) App. B. Diagnostic formulas for neutral and stable conditions; Batchvarova & Gryning (1991) slab model for convective conditions.
- Selectable options (the comparison of options is itself a result):
  - stable formula: Nieuwstadt (1981) (verify it is the one in Kljun App. B) vs Zilitinkevich (as in legacy);
  - γ (free-atmosphere lapse rate): constant (default 0.0025 K/m as legacy; to be justified) vs site/season-varying from ERA5 pressure levels;
  - convective trigger: H > 0 vs 1/L class (must be consistent);
  - H forcing: `H` (FLUXES) vs `H_F_MDS` (gap-filled) vs `H_CORR` (energy-balance corrected, FLUXNET). The closure-corrected run is an exploratory sensitivity test, not based on the assumption that EC underestimates H.
- B&G needs continuous forcing from sunrise: define the gap rule (restart at next sunrise, flag integrated values after a gap).
- Each legacy bug fix behind a toggle, to quantify the effect of each fix (see `legacy_review.md`).
- Regression baseline: the operational `PBLH` column in ETC L2 FLUXES files.

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
- Site metadata: lat, lon, altitude, zm, canopy height, d, z0, target area (source: ICOS BADM; to confirm).
- Missing values: `-9999`.

## Open points

- Time zone and timestamp convention of ICOS L2 files (UTC vs local standard time) vs ERA5 (UTC, instantaneous at the hour).
- Units of `T_SONIC` in L2 (°C or K).
- ERA5 extraction details (nearest cell? downscaling? QC?), and availability of neighbouring cells.
- Source of zm, d, z0 and target area per site.
- Measured PBLH datasets: to be collected; open data candidates: ACTRIS/Cloudnet, E-PROFILE, IGRA radiosondes.
- Literature: the current bib is footprint-focused. A PBL collection is needed (B&G 1991, Nieuwstadt 1981, Zilitinkevich, Seidel et al. 2010/2012, ERA5 BLH evaluations, ceilometer/lidar MLH methods).
