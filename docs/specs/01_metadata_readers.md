# Spec 01: Site metadata, exclusions and data readers

Phase 0. Branch: `spec-01-readers`.
Read `CLAUDE.md` and `docs/project_plan.md` first. This spec builds the input layer only:
**no PBLH modelling, no plots, no standard-format writer** (those are later specs).

## Goal

Reliable, tested functions that return, for one ICOS station, tidy half-hourly EC/meteo data in UTC and SI units,
plus hourly ERA5 BLH, together with site heights and exclusion flags. Tested locally on BE-Bra.

## Dependencies

Add with `renv::install()` and record with `renv::snapshot()`:
`config`, `here`, `readr`, `dplyr`, `tidyr`, `purrr`, `stringr`, `lubridate`, `cli`, `testthat`.
Do not add other packages without stating why in the summary.

## Files to create

```
R/utils_config.R
R/site_metadata.R
R/exclusions.R
R/read_icos.R
R/read_era5.R
scripts/01_check_inputs.R
tests/testthat/test-site_metadata.R
tests/testthat/test-exclusions.R
tests/testthat/test-read_icos.R
tests/testthat/test-read_era5.R
tests/testthat/fixtures/   (small synthetic files, see Tests)
tests/testthat/helper-load.R   (sources R/*.R for tests)
```

## 1. Configuration (`R/utils_config.R`)

- `get_cfg()`: wraps `config::get(file = here::here("config.yml"))`.
- `resolve_path(path, site = NULL)`: replaces the `{site}` placeholder with the station code; relative paths are
  resolved against the project root (`here::here()`), absolute paths are returned unchanged.
- Fail with an informative `cli::cli_abort()` if a configured path is `"TODO"` or does not exist.

## 2. Site metadata (`R/site_metadata.R`)

### `read_sites(file)`
- Reads `config/sites.csv`. `NA` strings are missing values.
- Column types: `site_id` character; `lat`, `lon`, `altitude_m`, `utc_offset_h`, `zm_m`, `canopy_height_m`, `d_m`,
  `z0_m`, `slope_perc` double; `igbp`, `slope_class`, `notes` character.
- Validation: `site_id` unique; `lat`, `lon`, `utc_offset_h` non-missing for every row (abort and list offending sites).
  Other columns may be missing.

### `read_sites_timevarying(file)`
- Reads `config/sites_timevarying.csv`: `site_id, variable, value, date_start, date_end, notes`.
- `variable` ∈ {`zm_m`, `canopy_height_m`}; abort on other values. Dates are calendar dates (`YYYY-MM-DD`);
  `date_end` may be `NA` (open-ended).

### `get_site_heights(site_id, timestamp_utc, sites, sites_tv)`
Returns a tibble with one row per input timestamp:
`timestamp_utc, zm_m, zm_source, canopy_height_m, hc_source`.

Rules (from `project_plan.md`; do not change without asking):
- Dates in `sites_timevarying.csv` are **local standard time** dates, starting at 00:00 local. Convert to UTC
  with the site `utc_offset_h` before comparing with `timestamp_utc`.
- **zm** is a step function: value valid from `date_start` (inclusive) to `date_end` (exclusive; `NA` = open-ended).
  Before the first period: use the first value, `zm_source = "extrapolated_before_first"`.
  Inside a period: `zm_source = "timevarying"`. No time-varying record for the site: use `sites.csv`,
  `zm_source = "static"`.
- **canopy height** (point measurements): linear interpolation in time between measurements
  (`hc_source = "interpolated"`), constant before the first / after the last measurement
  (`hc_source = "held_constant"`). Single measurement: constant everywhere. No record: `sites.csv` value,
  `hc_source = "static"`. All missing: `NA`, `hc_source = "missing"`.
- Do **not** compute displacement height or roughness length here (later spec).

## 3. Exclusions (`R/exclusions.R`)

### `read_exclusions(file)`
Reads `config/exclusion_periods.csv`: `site_id, date_start, date_end, scope, reason`.

### `flag_exclusions(data, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")`
- Adds `excl_flag` (logical) and `excl_reason` (character, `NA` if not excluded).
- Interval: from `date_start` 00:00 local standard time (inclusive) to the **end** of `date_end`
  (i.e. `date_end + 1 day` 00:00 local, exclusive). Both converted to UTC with `utc_offset_h`.
- Only `scope == "all"` exists now; abort on unknown scopes.
- Reports the number of flagged records with `cli::cli_inform()` (also when zero).
- Rows are **flagged, not dropped**. Provide `drop_excluded(data)` that removes flagged rows and reports how many.

## 4. ICOS L2 readers (`R/read_icos.R`)

### `find_icos_files(site_id, cfg)`
- Directory: `resolve_path(cfg$icos_l2_dir)/<site_id>_output_L2I_2026/`.
- Patterns: `ICOSETC_<site>_FLUXES_INTERIM_L2.zip`, `ICOSETC_<site>_METEO_INTERIM_L2.zip`,
  `ICOSETC_<site>_FLUXNET_FLUXMET_HH_*.zip` (the year range in the name varies between sites).
- FLUXES is required (abort if missing). METEO and FLUXNET optional (warning if missing).
  More than one match for a pattern: abort and list the files.
- Returns a named list of paths (`fluxes`, `fluxnet`, `meteo`; `NULL` if missing).

### Reading (common rules for the three files)
- Read the CSV directly from the zip with `readr::read_csv()` (each zip holds one CSV), using `col_select` to read
  only the needed columns. Use `tidyselect::any_of()` semantics: columns may be absent at some sites;
  list missing columns in a `cli` warning.
- `na = c("-9999", "-9999.0", "-9999.00", "NA", "")`.
- Read `TIMESTAMP_START` and `TIMESTAMP_END` as character (`YYYYMMDDHHMM`, local standard time).
  Create `timestamp_start_utc` and `timestamp_end_utc` (POSIXct, `tz = "UTC"`) =
  local clock time − `utc_offset_h` hours. Keep the original `TIMESTAMP_START` as character for traceability.
  Example: BE-Bra (UTC+1), `202001010000` → `2019-12-31 23:00:00 UTC`.
- Checks (report with `cli`, do not silently fix): duplicated timestamps (abort), non-monotonic order (sort and warn),
  step ≠ 30 min (warn with count). Return the data on a **complete regular 30-min grid** from first to last
  timestamp (missing steps become rows of `NA`), and report how many rows were added.
- **Units: convert at this boundary only**, and rename converted columns with a unit suffix so their unit is explicit:
  - °C → K (+273.15): `TA_F` → `TA_F_K`, `TA` → `TA_K`, `TA_1..TA_10` → `TA_1_K..TA_10_K`
  - Already K (rename only): `T_SONIC` → `T_SONIC_K` (*amended 2026-10-12: `T_SONIC` is in K in L2 FLUXES, verified on BE-Bra*)
  - kPa → Pa (×1000): `PA_F` → `PA_F_PA`, `PA` → `PA_PA` (METEO PA: verify on BE-Bra that it is kPa —
    values ≈ 100 — and state the result in the summary; if not, stop and ask)
  - All other columns keep their ICOS name and unit (document units in roxygen).
- No other modification of values (no filtering, no gap-filling).

### Columns to read
- **FLUXES**: `TIMESTAMP_START, TIMESTAMP_END, H, H_UNCLEANED, H_SSITC_TEST, USTAR, MO_LENGTH, ZL, T_SONIC, WS, WD,
  U_SIGMA, V_SIGMA, W_SIGMA, TAU, TAU_SSITC_TEST, PBLH, FETCH_MAX, FETCH_50, FETCH_70, FETCH_80, FETCH_90,
  FOOTPRINT_TA_CONTR, FOOTPRINT_FLAG`
- **FLUXNET**: `TIMESTAMP_START, TIMESTAMP_END, TA_F, TA_F_QC, PA_F, PA_F_QC, H_F_MDS, H_F_MDS_QC, H_CORR, NIGHT,
  SW_IN_POT, WS_F, P_F`
- **METEO**: `TIMESTAMP_START, TIMESTAMP_END, TA, PA, P, SW_IN, WS, TA_1 … TA_10`

### Functions
- `read_icos_fluxes(path, utc_offset_h)`, `read_icos_fluxnet(path, utc_offset_h)`, `read_icos_meteo(path, utc_offset_h)`.
- To avoid name clashes when merging, prefix FLUXNET and METEO columns (except timestamps):
  `fn_` for FLUXNET, `met_` for METEO. FLUXES columns keep their names.
- `read_icos_site(site_id, cfg, sites, sites_tv, exclusions)`:
  full outer join of the three tables on `timestamp_start_utc` (periods differ: e.g. FLUXNET ends in 2025 while
  FLUXES goes to 2026), then add `site_id`, the heights from `get_site_heights()` and the exclusion flags from
  `flag_exclusions()`. Output sorted by time, one row per half-hour.

### Plausibility report (flag, never modify)
- **Unit sanity check** (*amended 2026-10-12*), on raw values before conversion: abort with an informative message if the
  median of a column declared °C is outside [-40, 50], declared K outside [220, 330], declared kPa outside [50, 110].

`icos_range_report(data, include_excluded = FALSE)` returns a tibble
`variable, n_total, n_excluded, n_valid, pct_valid, n_out_of_range, min, max`. By default rows with `excl_flag` TRUE
are skipped (*amended 2026-10-12*); the number of excluded rows is always reported. Uses
using these ranges (only for reporting):
`H` [-200, 800] W m-2; `USTAR` [0, 3] m s-1; `T_SONIC_K` [223, 323] K; `PA_F_PA` [50000, 110000] Pa;
`PBLH` [10, 5000] m; `|MO_LENGTH|` > 0.
These ranges are not physical filters and must not be used to remove data.

## 5. ERA5 BLH reader (`R/read_era5.R`)

### `find_era5_files(site_id, cfg)`
- Directory: `resolve_path(cfg$era5_dir, site = site_id)`.
- Accept only file names matching
  `^<site>__(ERA5|ERA5T)__reanalysis-era5-single-levels__boundary_layer_height__(\d{4})__lon([+-]\d+)__lat([+-]\d+)\.csv$`.
  Other files in the folder are ignored, with an informative message listing them.
- Parse `product` (ERA5 / ERA5T), `year`, `lon_offset`, `lat_offset` (integers).
- If ERA5 and ERA5T exist for the same year and offsets: keep ERA5, warn.
- Returns a tibble of files with the parsed fields.

### `read_era5_blh(files, offsets = "center", sites = NULL, coord_tol_deg = 0.01)`
- `offsets = "center"` keeps `lon_offset == 0 & lat_offset == 0`; `offsets = "all"` keeps all
  (the 3×3 neighbours will arrive later; the code must already support them).
- Output columns: `timestamp_utc` (from `time`, POSIXct UTC, instantaneous value at the hour), `blh_m` (from `value`,
  m above model ground), `product`, `lon_offset`, `lat_offset`, `grid_lat` (`my_lat`), `grid_lon` (`my_lon`),
  `site_lat` (`latitude`), `site_lon` (`longitude`).
- Check that `time` and `datetime` are identical; if not, warn and use `time`.
- If `sites` is given (*amended 2026-10-12*): warn if `latitude`/`longitude` in the ERA5 files differ from `sites.csv`
  by more than 0.01° (ERA5 for SE-Myc was extracted with old, wrong coordinates).
- Checks: duplicated timestamps per offset (abort), hourly continuity (report gaps), `blh_m` < 0 or missing (report).
- Add `grid_dist_km`: great-circle (haversine) distance between site and grid-cell centre (R = 6371 km).

## 6. Script `scripts/01_check_inputs.R`

Thin entry point. For each site in `cfg$test_sites` (locally: BE-Bra):
1. load config, sites, time-varying heights, exclusions;
2. `read_icos_site()` and `read_era5_blh()`;
3. print a concise report: period covered by each file, rows added to complete the grid, `icos_range_report()`,
   `zm_source` / `hc_source` counts, number of excluded records, ERA5 years and products found, `grid_dist_km`,
   and the METEO `PA` unit check;
4. write the range report to `output/tables/01_check_inputs_<site>.csv`.

## 7. Tests

Use small **synthetic** fixtures in `tests/testthat/fixtures/` (a few rows; never real data volume):
- a zipped FLUXES CSV, a zipped FLUXNET CSV, a zipped METEO CSV (with `-9999` values, one missing timestamp,
  one column absent from METEO);
- ERA5 CSVs: one ERA5 and one ERA5T file for the same year, one file for another year, and one non-BLH file
  (to be ignored);
- small `sites.csv`, `sites_timevarying.csv`, `exclusion_periods.csv`.

Required test cases:
- timestamp conversion: `202001010000` at UTC+1 → `2019-12-31 23:00:00 UTC`; at UTC+2 and UTC+0;
- `-9999` → `NA`; °C → K and kPa → Pa conversion and column renaming;
- grid completion adds exactly the missing row; duplicated timestamps abort;
- missing optional column produces a warning, not an error;
- zm step function, including `extrapolated_before_first` and `static`; canopy height interpolation, held constant
  outside range, single measurement, static fallback;
- exclusions: SE-Nor-like case, with start and end boundaries (first and last half-hour of the period flagged,
  the half-hours just outside not flagged);
- ERA5: file-name parsing, ERA5 preferred over ERA5T, non-BLH file ignored, `offsets` selection, haversine distance
  against a hand-computed value.

## Acceptance criteria

- `testthat::test_dir("tests/testthat")` passes with no failures.
- `Rscript scripts/01_check_inputs.R` runs on BE-Bra without errors and prints the report.
- No hardcoded paths; no `library()` calls in `R/`; roxygen with units on every function.
- `renv::status()` reports a consistent state.
- Nothing under `data/` or `output/` is committed.

## Out of scope (do not implement)

PBLH model, displacement height / roughness length, standard-format files, footprint code, plots, `targets`
pipeline, SLURM scripts.

## Summary to report back

As in `CLAUDE.md`, plus: the METEO `PA` unit check result, the BE-Bra coverage numbers (valid % for `H`, `USTAR`,
`MO_LENGTH`, `T_SONIC_K`, `fn_H_F_MDS`, `PBLH`), the ERA5 files found, and anything in the real files that
differs from this spec.
