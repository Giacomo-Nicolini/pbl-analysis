# CLAUDE.md: pbl-analysis

Standing instructions for Claude Code. Read this first, every session.

## Project in one paragraph

Estimate Planetary Boundary Layer height (PBLH) from eddy-covariance (EC) data at ICOS ecosystem stations
(PBLHmod: Kljun et al. 2015 App. B diagnostics + Batchvarova & Gryning 1991 slab model), evaluate it against
ERA5 and measured PBLH (lidar, ceilometer, radiosonde...), explain when/why they disagree, and quantify the
impact of PBLH errors on FFP footprint estimates (Kljun et al. 2015). The work will become a peer-reviewed paper,
so code must be reproducible, documented and traceable to references.

- Full design: `docs/project_plan.md` (source of truth for scientific choices)
- Legacy code issues: `docs/legacy_review.md`
- Task specs: `docs/specs/NN_*.md`. **Implement only what the current spec asks.**

## Golden rules

1. **Never change scientific content silently.** Formulas, constants, thresholds, defaults: if a spec is ambiguous
   or you think something is physically wrong, stop and ask, or implement as specified and flag it clearly in your summary.
2. **Units everywhere.** Every function argument and returned column documents its unit in roxygen (`@param H Sensible heat flux [W m-2]`).
   Internally use SI: K for temperature, Pa for pressure, m, s. Convert at the reader boundary only.
3. **References in code.** Every formula carries a comment with its source (author, year, equation number when known).
4. **No hardcoded paths.** All paths come from `config.yml` via the `config` package (see Environment).
5. `legacy/` is read-only reference. Never edit, source, or import from it in new code (except in the explicit legacy regression test).
6. Never commit data, credentials or outputs (see `.gitignore`). Never print secrets.

## Environment

- R 4.4.1, both locally (Windows) and on the cluster (Linux, SLURM). Packages managed by `renv`:
  after adding a dependency run `renv::snapshot()`. Do not call `install.packages()` or `library()` inside functions.
- Paths: `config.yml` with profiles `default` (local PC) and `cluster`. Select with env var `R_CONFIG_ACTIVE=cluster`.
  Access with `cfg <- config::get()`.
- Development and testing happen locally on ONE site (BE-Bra); production loops over all ICOS stations on the cluster
  as a SLURM array job (one task per station). Code must run identically in both places.

## Input data (facts; do not guess beyond these)

- ICOS L2, per station `XX-Yyy`, in `<icos_l2_dir>/XX-Yyy_output_L2I_2026/`:
  - `ICOSETC_XX-Yyy_FLUXES_INTERIM_L2.zip`: H, USTAR, MO_LENGTH, ZL, T_SONIC, WS, WD, V_SIGMA, PBLH (legacy operational output), FETCH_*, FOOTPRINT_*, QC flags
  - `ICOSETC_XX-Yyy_FLUXNET_FLUXMET_HH_<years>_v*.zip`: TA_F, PA_F, H_F_MDS, H_CORR, NIGHT, SW_IN_POT, WS_F, P_F
  - `ICOSETC_XX-Yyy_METEO_INTERIM_L2.zip`: TA, TA_1..TA_10 (profile), P, SW_IN
  - Files are zipped CSVs. Read directly from the zip (no permanent unzipping). Use only needed columns.
  - Missing value: `-9999`. Timestamps: `TIMESTAMP_START`, `TIMESTAMP_END` as `YYYYMMDDHHMM`, **local standard time (no DST)**.
    Convert to UTC with the site UTC offset from site metadata. Store both start and end.
  - Units: `T_SONIC`, `TA*` in °C (convert to K); `PA_F` in kPa (convert to Pa); H in W m-2; USTAR in m s-1; MO_LENGTH in m.
- ERA5 BLH: yearly CSV per station, hourly, UTC, instantaneous at the hour, e.g.
  `BE-Bra__ERA5__reanalysis-era5-single-levels__boundary_layer_height__2025__lon+0__lat+0.csv`.
  Columns: `value` (BLH, m above ground), `time`, `latitude`, `longitude` (site), `location`, `pos_my_lat`, `pos_my_lon`,
  `my_lat`, `my_lon` (grid cell), `datetime`, `year`.
- Site metadata: `config/sites.csv` (lat, lon, altitude, UTC offset, zm, canopy height, d, z0...).
  Missing values allowed; functions must fail with an informative error only when a *required* value is missing.
  Time-varying zm / canopy height: `config/sites_timevarying.csv` (derived from `config/sites_H_EC_canopy_formBIF.csv`).
- **Exclusions**: `config/exclusion_periods.csv` lists periods that must be removed for a site (e.g. SE-Nor clearcut).
  Every pipeline step must apply them and report how many records were excluded. Never skip this silently.

## Code organisation

```
R/            functions only (no top-level code). One topic per file: read_icos.R, read_era5.R, pblh_model.R, ...
scripts/      thin entry points, numbered (01_..., 02_...)
_targets.R    pipeline (targets package) once phase 0 is done
tests/testthat/  unit tests
config/       sites.csv, other small metadata
config.yml    path profiles
cluster/      SLURM job scripts
docs/         plan, specs, references
paper/        Quarto manuscript (later)
data/         raw/, standard/, derived/ (gitignored)
output/       figures/, tables/ (gitignored)
```

## Style

- Tidyverse style, snake_case, explicit `pkg::fun()` for non-base functions (no `library()` in R/).
- Data handling: `dplyr`/`tidyr`/`readr`/`lubridate`; `data.table` allowed where speed matters (large I/O).
- Time: `POSIXct` with `tz = "UTC"` everywhere internally.
- Vectorised code; avoid growing vectors in loops. Sequential integration (B&G) may loop, but preallocate.
- Functions are pure: inputs in, data frame out. Plotting and file writing live in separate functions.
- Every exported function has roxygen docs: purpose, params with units, return value with units, references.
- Plots: `ggplot2`, publication style, saved via a helper to `output/figures/` (png + pdf). Never `windows()`.

## Testing

- `testthat` in `tests/testthat/`. Run with `testthat::test_dir("tests/testthat")`.
- For physics: known-value tests (hand-computed cases), limit tests (e.g. neutral limit of the stable formula),
  unit tests for conversions, and plausibility checks (PBLH > 0, finite).
- Readers: test on a small fixture file in `tests/testthat/fixtures/` (a few rows, no real data volume).
- Run the tests before declaring a task done, and report results.

## Git

- Small, focused commits with clear messages (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`).
- One branch per spec (`spec-01-readers`), merge to `main` when tests pass.
- Before committing, check `git status` for data, zip, credentials: they must never be staged.

## When you finish a task

Summarise: what was implemented, files touched, tests run and results, assumptions made,
open questions/flags for the scientist. Update `docs/project_plan.md` "Open points" if relevant.
