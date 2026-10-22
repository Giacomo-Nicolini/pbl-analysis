# Spec 07: measured PBLH from ACTRIS/Cloudnet Doppler lidar; first three-way comparison at FI-Hyy

Phases 2–3. Branch: `spec-07-measured-cloudnet` (from `main`).
Read first: `CLAUDE.md`, `docs/project_plan.md` ("Comparing different method classes", "Evaluation"),
`docs/specs/standard_format.md`, `docs/specs/04_era5_standard_compare.md`, `docs/references/literature_review_pblh.md`
(section 3). Agreed with G on 2026-10-22.

## Goal

1. Inventory of open ACTRIS/Cloudnet data at the Cloudnet sites next to ICOS ecosystem stations of this project.
2. A turbulence-based mixing-layer height (MLH) from the Cloudnet Doppler-lidar dissipation rate (ε), after
   Manninen et al. (2018), written in the standard PBLH format (`method_class = turbulence`).
3. First three-way comparison at FI-Hyy: PBLHmod, ERA5 and the operational ICOS PBLH against the measured MLH.

FI-Hyy first (co-located with the Cloudnet site Hyytiälä, SMEAR II). Other sites: inventory only in this spec.

## 1. Site pairing and inventory

New tracked file `config/measured_sites.csv` (one row per ICOS site × measurement source):
`icos_site, provider, provider_site, provider_lat, provider_lon, distance_km, method_class, instruments, licence,
url, notes`. `distance_km` computed with `haversine_km()` from `config/sites.csv`. Initial rows (Cloudnet,
licence CC BY 4.0):

| ICOS site | Cloudnet site id | notes |
|---|---|---|
| FI-Hyy | `hyytiala` | co-located (SMEAR II) |
| FI-Lom | `kenttarova` | Kenttärova is the ICOS station FI-Kef (not in this project); FI-Lom (Lompolojänkkä) ≈ 2 km away |
| DE-RuS | `juelich` | JOYCE, 5.26 km (Cloudnet API coordinates): daytime convective comparisons only |
| SE-Nor | `norunda` | co-located (ICOS atmosphere NOR + ACTRIS) |
| FR-Gri | `palaiseau` | SIRTA, ≈ 22 km: daytime convective comparisons only |

Distance rule (project plan, G 2026-10-22): < 5 km all regimes; ≤ 25 km daytime convective only; 25–50 km daytime
climatologies only; > 50 km not used. Record this rule in `docs/project_plan.md` ("Comparing different method
classes") and apply it in the comparison code via `distance_km`.

Script `scripts/07a_cloudnet_inventory.R`: for each Cloudnet site above and each product `doppler-lidar`,
`epsilon-lidar`, `doppler-lidar-wind`, `lidar`, `mwr`, `classification`, query the Cloudnet API
(https://cloudnet.fmi.fi/api; check the API documentation for parameters and paging; be polite: one request at a
time, short pause) and write `output/tables/07_cloudnet_inventory.csv`: site, product, instrument (pid/name),
month, number of daily files, total size. Period 2020-01 → latest. No data download in this step.
Report: per site, the months with `epsilon-lidar` (and `doppler-lidar`) coverage.

## 2. Download (FI-Hyy only)

- New config key `cloudnet_dir` (default `data/raw/cloudnet/{site}`; cluster: `"LOCAL"`).
- `scripts/07b_cloudnet_download.R`: download the daily `epsilon-lidar` files for `hyytiala` over the ICOS FI-Hyy
  data period (and the `doppler-lidar` files only if section 3 needs a variable that `epsilon-lidar` lacks).
  Skip files already present with the same checksum (the API gives one); log failures and retry once.
- **Before downloading, report the total size from the API metadata; if it exceeds 10 GB, stop and ask G.**
- Data stays under `data/` (gitignored). Keep the API file metadata (uuid, checksum, processing version, PID of the
  instrument) in a small CSV next to the files, for provenance.

## 3. Reader and MLH (`R/measured_cloudnet.R`)

**First inspect** one `epsilon-lidar` file and report its variables, dimensions, units, time resolution, height
definition (above ground or above sea level; instrument altitude), the lowest valid gate, and any quality or
SNR/precipitation fields. Then implement. If the file content does not allow the method below, stop and report.

NetCDF reading: add `ncdf4` (or `RNetCDF` if `ncdf4` fails to install) to renv. Compiled package: if Windows blocks
its DLL (as happened with `cli`), stop and tell G.

Method (Manninen et al. 2018, Sect. 3.2.1; `Manninen2018`; thresholds as options):
- per profile, start at the lowest reliable gate `opts$z_min` (default: the lowest valid gate in the file, ≈ 105 m
  for the Halo lidar; report the value used) and go upward while ε ≥ `opts$eps_surface` (default 1e-4 m² s⁻³);
  the MLH is the height of the top of the last gate meeting the condition (surface-connected mixing);
- if the lowest reliable gate already has ε < `eps_surface` (or is missing): MLH not measurable, **censored below
  `z_min`** (`pblh_m = NA`, flagged; keep `z_min` in `x_censor_below_m`);
- if ε is missing in a gate inside the column (low SNR), stop the search there and flag `x_gap_in_column`;
- profiles flagged as precipitation or with insufficient signal in the file's own QC → `NA`;
- heights converted to metres above ground (`height_reference = agl`) with the instrument altitude and the site
  ground altitude from the file metadata; document the conversion.
- Cloud-driven mixing (Manninen's top-down step) is out of scope; note it as a limitation.

Aggregation to the ICOS half-hours: median of the profile MLHs whose time falls in [start, end); `n_samples` =
profiles used; `NA` unless at least `opts$min_coverage` (default 0.5) of the expected profiles are valid. Censored
profiles count as valid: if more than half the valid profiles are censored, the half-hour is censored
(`pblh_m = NA`, `qc_flag = 1`); otherwise the median of the non-censored ones with `qc_flag = 2` if any censored
profile was present. `qc_flag_meaning`: 0 all profiles measured; 1 censored below z_min; 2 partly censored.

Standard file: `source_id = "cloudnet_eps_mlh"`, `method_class = turbulence`, `time_reference = interval_mean`,
`instrument` = Cloudnet instrument name and PID, `retrieval_algorithm` = "surface-connected ε ≥ 1e-4 m2 s-3 from
the lowest reliable gate (Manninen et al. 2018)", `licence = CC BY 4.0`, `citation` = Cloudnet/ACTRIS data
citation as given by the API plus Manninen et al. 2018, `processing_history` with thresholds and z_min.

Tests (synthetic profiles): MLH at the right gate; censored case; gap inside the column; precipitation/QC → NA;
aggregation rules (coverage, censored majority, partly censored); height conversion.

## 4. FI-Hyy runs

Inputs G provides locally: ICOS L2 FI-Hyy (FLUXES, FLUXNET, METEO) under `data/raw/icos_l2/`, ERA5 FI-Hyy (three
families, from the cluster delivery) under the configured ERA5 dirs. If anything is missing, stop and list it.

- Spec-01 readers and checks for FI-Hyy (`scripts/01_check_inputs.R FI-Hyy`): report unit checks, coverage,
  exclusions.
- `scripts/02_run_pblh_model.R FI-Hyy` (all current cases) and `scripts/04_era5_compare.R FI-Hyy` (standard files,
  ERA5 cell metadata, PBLHmod vs ERA5). Report the operational fill detection at FI-Hyy (amendment 7 rule) and
  anything that looks site-specific (tall forest, zm ≈ 27 m, latitude 61.8° N: check polar-day flags and the
  convective window in June).

## 5. Three-way comparison (`scripts/07_measured_compare.R`)

Reference: `cloudnet_eps_mlh` (turbulence class). Pairs:
1. `pblhmod_default` vs measured — on half-hours;
2. `era5_blh_center` vs measured — measured matched to the ERA5 hour with the spec-04 rule (mean of the two
   half-hours adjacent to t, `require_both = TRUE`);
3. `icos_oper` (qc 0, and all) vs measured — on half-hours;
4. `era5_blh_center` vs `pblhmod_default` at FI-Hyy (as spec 04, for the record).

Metrics and strata as in spec 04 (`pblh_metrics()`, regime from PBLHmod, season, hour). In addition:
- the share of half-hours censored below z_min, by hour and season, reported for each pair (a censored reference
  cannot be compared; report how often each model value is itself below z_min on those half-hours);
- per the method-class rules in `project_plan.md`: the primary statistics are daytime convective (PBLHmod regime
  `convective`); stable and transition hours are reported separately and never pooled with daytime;
- ERA5 (thermodynamic/Ri definition) vs turbulence-based MLH: report, but label as cross-class.

Figures (PNG): `07_diurnal_by_season_FI-Hyy` (measured, PBLHmod, ERA5, operational qc 0); `07_scatter_FI-Hyy`
(each source vs measured, by regime, log axes); `07_censoring_FI-Hyy` (share censored by hour × season);
`07_logratio_by_hour_FI-Hyy`.
Tables: `07_metrics_FI-Hyy.csv` (long format), `07_censoring_FI-Hyy.csv`.

## 6. Docs

- `CLAUDE.md` "Input data": Cloudnet facts (product, variables, units, height reference, licence, API), the
  `cloudnet_dir` key.
- `docs/specs/standard_format.md`: example of a measured source (`cloudnet_eps_mlh`) and the censoring convention.
- `docs/STATUS.md`: spec 07 row and FI-Hyy results; `docs/project_plan.md`: distance rule, measured-data sources
  (Cloudnet pairs), and the PI requests in progress (G: SE-Htm, SE-Nor, SE-Svb, FI-Hyy, CZ-BK1 atmosphere
  stations; list kept locally in `docs/local/sites_PBLH_measure.txt`).

## Acceptance criteria

- Tests pass; scripts 07a, 07b, 01, 02, 04 and 07 run for FI-Hyy; standard files for all sources written and
  validated; nothing under `data/` or `output/` committed.
- Inventory table for the five Cloudnet sites.
- MLH method documented with thresholds, z_min and height conversion; censoring reported.

## Out of scope

Download and processing for sites other than FI-Hyy; radiosondes (IGRA) and ceilometer (aerosol) MLH; Manninen's
full ABL classification (cloud-driven mixing, skewness, shear); triple collocation; any change to PBLHmod.

## Summary to report back

As in `CLAUDE.md`, plus: inventory highlights; the ε file structure and z_min used; MLH coverage and censoring at
FI-Hyy; pair metrics (daytime convective first); FI-Hyy model/ERA5 results compared with BE-Bra; anything
unexpected.
