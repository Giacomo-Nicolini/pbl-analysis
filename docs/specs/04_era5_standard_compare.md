# Spec 04: ERA5 readers, standard PBLH format, first PBLHmod vs ERA5 comparison (BE-Bra)

Phases 0/2/3 (first step). Branch: `spec-04-era5-compare` (from `main`).
Read first: `CLAUDE.md` (section "Input data", ERA5), `docs/project_plan.md` (Standard PBLH format; ERA5
representativeness; Evaluation), `docs/references/literature_review_pblh.md` (sections 2–3).

## Goal

1. Read the new ERA5 delivery (hourly single levels, 6-hourly pressure levels, static fields; 3×3 cells) with
   correct cell positions, units and time conventions.
2. Describe the ERA5 cell of each site (position, orography, land–sea fraction, sub-grid orography) and the
   3×3 spread of BLH per hour.
3. Implement the standard PBLH format (data table + YAML sidecar) and write three sources for BE-Bra:
   PBLHmod v2, ERA5 BLH, operational ICOS PBLH (with its floor and fill rows flagged).
4. A first, **exploratory** comparison of PBLHmod vs ERA5 (and vs operational PBLH) at BE-Bra: matched pairs,
   metrics by regime/season/hour, diurnal cycles, a few figures. No inference or modelling yet.

Decisions agreed with G (2026-10-16), do not change without asking:
- Time matching: ERA5 BLH (instantaneous at the hour t) vs the **mean of the two half-hours adjacent to t**
  (the one ending at t and the one starting at t).
- ERA5 cell: tower cell = main series; bilinear interpolation to the tower position = sensitivity; 3×3 spread =
  per-hour representativeness metric.
- Period/site: BE-Bra, all overlapping data (PBLHmod valid from Oct 2020; ICOS L2 ends 2026-08-31; ERA5 ends
  2026-10-07).
- ERA5 sensible heat flux vs EC H and ERA5 cloud base height: descriptive figures only, no conclusions.

## Files

```
R/read_era5.R            rewrite (new format); keep haversine_km()
R/era5_cells.R           cell metadata, bilinear interpolation, 3×3 spread
R/standard_format.R      write/read/validate standard PBLH files
R/standardize_sources.R  PBLHmod, ERA5, operational PBLH -> standard tables
R/compare_pblh.R         matching to ERA5 hours, metrics, diurnal summaries
R/plot_utils.R           save_figure() helper (png + pdf to output/figures/)
R/plot_compare.R         figure functions (return ggplot objects; no file writing)
scripts/04_era5_compare.R
docs/specs/standard_format.md   (the format, as implemented; short)
tests/testthat/test-read_era5.R           (rewrite)
tests/testthat/test-era5_cells.R
tests/testthat/test-standard_format.R
tests/testthat/test-standardize_sources.R
tests/testthat/test-compare_pblh.R
tests/testthat/fixtures/make_fixtures.R   (extend; new synthetic ERA5 fixtures for XX-Tst)
```
Dependencies: add `ggplot2` (and `yaml` if not already a direct dependency); `renv::snapshot()`. No other new
packages without asking.

## 1. Configuration

Replace `era5_dir` in `config.yml` (both profiles) with:

| key | default (local) | cluster |
|---|---|---|
| `era5_sl_dir` | `data/raw/era5/{site}` | `/bigdata/c.trotta/era2oneflux_claude_zarr_giacomo/ERA5_surface_single_level/{site}` |
| `era5_pl_dir` | `data/raw/era5_pressure_levels/{site}` | `/bigdata/c.trotta/era2oneflux_claude_zarr_giacomo/ERA5_pressure_levels/{site}` |
| `era5_static_dir` | `data/raw/era5_single_levels` | `/bigdata/c.trotta/era2oneflux_claude_cds_giacomo_grib/out1/{var}` |

Extend `resolve_path()` to substitute `{var}` (same mechanism as `{site}`; abort if a placeholder is left
unresolved). Locally all static variables are in one folder (no `{var}`), on the cluster one folder per variable:
both must work. Also set `standard_dir` (cluster) if still `TODO`: leave `TODO`, do not invent paths.

## 2. ERA5 readers (`R/read_era5.R`, rewrite)

The old reader (yearly `value`/`time` files, `read_era5_blh()`, `find_era5_files()` with `lon±n__lat±n`) is
superseded: remove it, its fixtures and its tests. Keep `haversine_km()` and the idea of
`check_era5_site_coords()`.

### 2.1 Header

`parse_era5_header(path)` → named list of all `# key: value` lines between the first line and
`# ==== FINE METADATA (dati sotto) ====` (split at the first `": "`; trim), plus `n_header_lines`.
Numeric keys converted (`latitude_*`, `longitude_*`, `grid_*`, `pixel_*`, `n_rows`, `distance_km`); `levels`
→ integer vector. Abort if the end marker is missing.

### 2.2 File discovery

`find_era5_files(site_id, cfg, family = c("single_level", "pressure_level", "static"), variable = NULL)`

File-name patterns (`<site>` escaped):
- single level: `^<site>__reanalysis_era5_single_level__(<var>)__(\d{8}T\d{4})__(\d{8}T\d{4})__north_([+-]\d+)__south_([+-]\d+)__east_([+-]\d+)__west_([+-]\d+)\.csv$`
- pressure level: same with `reanalysis_era5_pressure_levels-pl`
- static: `^<site>__single-levels__(z|lsm|sdor)__(\d{12})__(\d{12})__north_([+-]\d+)__south_([+-]\d+)__east_([+-]\d+)__west_([+-]\d+)\.csv$`

Return a tibble, one row per file: `site_id, family, variable, path, file, grid_lat, grid_lon, lat_requested,
lon_requested, cell_row, cell_col` (+ `time_first`, `time_last` from the name).

**Cell position comes from the header**, never from the name:
- the tower cell is the file whose header has the smallest distance between (`latitude_grid`, `longitude_grid`)
  and (`latitude_requested`, `longitude_requested`);
- `cell_row = round((grid_lat - tower_grid_lat) / 0.25)` (+1 = north), `cell_col = round((grid_lon -
  tower_grid_lon) / 0.25)` (+1 = east); `grid_resolution_deg` from the header (abort if not 0.25 or inconsistent).
- Cross-check with the name (guards against a third convention): single/pressure level files use direction flags
  (`north_+1` → row +1, `south_+1` → row −1, `east_+1` → col +1, `west_+1` → col −1); static files use signed
  offsets in `north_` (row) and `east_` (col) with `south_+0__west_+0`. If name and header disagree → abort with
  both values.
- Per family × variable expect the 9 cells of a 3×3 block (rows and cols in −1..1). Fewer → warning listing the
  missing cells; duplicates → abort. Other files in the folder → ignored, listed in a message.

### 2.3 Variable registry

`era5_variables()` → tibble (single source of truth for names, units, conventions; cite ECMWF parameter
database / ERA5 documentation in comments):

| variable (file) | column | unit in → out | conversion | time reference |
|---|---|---|---|---|
| `boundary_layer_height` | `blh_m` | m → m | none | instantaneous at `timestamp_utc` |
| `cloud_base_height` | `cbh_m` | m → m | none; `NA` = no cloud (ERA5 sets it missing when there is no cloud) | instantaneous |
| `surface_sensible_heat_flux` | `H_era5_Wm2` | J m-2 (hourly accumulation, positive downward) → W m-2 positive upward | `-x / 3600` | interval (`t − 1 h`, `t`] |
| `temperature` (pl) | `T_K` | K | none | instantaneous |
| `geopotential` (pl) | `z_m2s2`, `zg_m` | m2 s-2 → m | `zg_m = z / g`, g = 9.80665 (geopotential height; differs from geometric height by < 0.5 % below 10 km — note in the docs) | instantaneous |
| `z` (static) | `orog_m` | m2 s-2 → m | `/ g` | – |
| `lsm`, `sdor` (static) | `lsm`, `sdor_m` | 0–1, m | none | – |

### 2.4 Readers

- `read_era5_single_level(files, variable, cells = c("all", "center"))` → long tibble:
  `site_id, cell_row, cell_col, grid_lat, grid_lon, grid_dist_km, timestamp_utc, <column>`;
  for `surface_sensible_heat_flux` also `timestamp_start_utc = timestamp_utc − 1 h`, `timestamp_end_utc =
  timestamp_utc`.
- `read_era5_pressure_levels(files, cells)` → long: `site_id, cell_row, cell_col, timestamp_utc, p_hPa, T_K,
  z_m2s2, zg_m, theta_K` with `theta_K = T_K * (1000 / p_hPa)^(r_d / cp)` (Poisson's equation; r_d and cp from
  `pblh_constants()`, r_d/cp = 0.2856). Temperature and geopotential files joined on cell × time × level (abort on
  mismatch). Column names `temperature_<p>hPa` / `geopotential_<p>hPa` parsed to `p_hPa`; check against the
  header `levels`.
- `read_era5_static(site_id, cfg)` → one row per cell: `cell_row, cell_col, grid_lat, grid_lon, orog_m, lsm,
  sdor_m` (header `units` checked against the registry; abort on mismatch).

Checks (report, never modify values):
- rows read = header `n_rows`; first/last time = header `first_time`/`last_time` (abort if not);
- continuity (hourly / 6-hourly), duplicated times per cell (abort);
- unit plausibility (abort if implausible, like the ICOS readers): median `blh_m` in [20, 3000]; median
  `T_K` at 1000 hPa in [230, 320]; `zg_m` increasing with decreasing pressure;
- **sign check for `H_era5_Wm2`**: the mean diurnal cycle (by hour of local solar time, from longitude) must have
  its maximum within ±3 h of solar noon; abort otherwise (catches a missing sign flip);
- coordinates: `lat_requested`/`lon_requested` vs `config/sites.csv` (tolerance 0.001°; warn).

## 3. ERA5 cells (`R/era5_cells.R`)

- `era5_cell_metadata(site_id, cfg, sites)` → one row: `site_id, site_lat, site_lon, site_altitude_m,
  cell_lat, cell_lon, cell_dist_km, orog_m (tower cell), dalt_m = site_altitude_m − orog_m, lsm, sdor_m`, and over
  the 3×3 block: `orog_min_m, orog_max_m, lsm_min, lsm_max, sdor_max_m, n_cells`.
- `era5_bilinear(era5_long, site_lat, site_lon, value_col)`: bilinear interpolation at the tower position from
  the 4 cells surrounding it (chosen from the sign of `site − cell centre` in lat and lon). Weights from the
  regular 0.25° grid (cite: standard bilinear interpolation). Returns `timestamp_utc, <value>_bilin`.
- `era5_spread(era5_long, value_col)` per time: `mean9, sd9, min9, max9, range9, cv9 = sd9/mean9, n9`
  (n9 < 9 when cells are missing; statistics computed on available cells and `n9` reported).

## 4. Standard PBLH format (`R/standard_format.R`, `docs/specs/standard_format.md`)

As in `docs/project_plan.md`, with these implementation rules:

- One data file per site × source: `<standard_dir>/<site>/<site>__<source_id>.csv.gz`, plus sidecar
  `<site>__<source_id>.yml`.
- Required columns: `timestamp_start_utc`, `timestamp_end_utc` (ISO 8601 UTC, `YYYY-MM-DDTHH:MM:SSZ`),
  `pblh_m` [m above ground], `source_id`. Recommended (may be all `NA`): `qc_flag` (integer; meaning in the
  sidecar), `uncertainty_m`, `n_samples`. Source-specific extra columns allowed only with prefix `x_`
  (e.g. `x_regime`).
- **Instantaneous values** (ERA5): `timestamp_start_utc == timestamp_end_utc`, and sidecar
  `time_reference: instantaneous`; interval values: `time_reference: interval_mean`.
- Sidecar fields: `site_id`, `source_id`, `method_class` (required; one of the classes in the plan), and
  `instrument`, `retrieval_algorithm`, `height_reference`, `instrument_altitude_m`, `time_reference`,
  `original_resolution`, `original_qc`, `qc_flag_meaning` (map), `processing_history`, `contact`, `licence`,
  `citation`, `created_utc`, `code_version` (git short hash if available) — any may be `unknown`.
  `metadata_completeness` = fraction of the optional fields not `unknown`, computed on write.
- Functions: `validate_standard_pblh(df, meta)` (required columns and types; start ≤ end; no duplicated
  start times; `pblh_m` finite or `NA`; `method_class` valid; abort with all problems listed),
  `write_standard_pblh(df, meta, dir)` (validates, then writes both files), `read_standard_pblh(path)`
  (reads data + sidecar, validates, returns the tibble with the sidecar as attribute `meta`).

## 5. Standardizers (`R/standardize_sources.R`)

Pure functions returning a standard table + meta list (writing is done by the script):

- `standardize_pblhmod(out, label)`: from `run_pblh_model()` output (or the derived csv). `source_id =
  "pblhmod_<label>"`, `method_class = parameterisation`, `time_reference = interval_mean`. `qc_flag`: 0 no flag;
  1 forcing interpolated (`flag_*_interp`); 2 h held by the model (`flag_sign_inconsistent`); 3 below zm
  (`flag_below_zm`); highest applicable value wins; meaning written in the sidecar. Extra: `x_regime`,
  `x_method`, `x_h_source`. Excluded rows: dropped (count reported).
- `standardize_era5_blh(era5_center, era5_bilin, spread, variant = c("center", "bilinear"))`:
  `source_id = "era5_blh_center"` / `"era5_blh_bilinear"`, `method_class = model`, `retrieval_algorithm =
  "IFS Cy41r2 diagnostic: bulk Richardson number = 0.25 from the lowest model level, no parcel excess (ECMWF 2016,
  Part IV Sect. 3.10.1, Eq. 3.90)"`, `height_reference = "above model surface"`,
  `time_reference = instantaneous`. Extra: `x_blh_sd9_m`, `x_blh_range9_m`, `x_blh_cv9`, `x_n9`. Rows inside
  exclusion periods (site exclusions, converted with the site UTC offset) are dropped and counted.
- `standardize_icos_oper(d)`: the operational `PBLH` column of FLUXES L2 (`read_icos_site()` output),
  `source_id = "icos_oper"`, `method_class = parameterisation`, `processing_history` describing
  `calc_PBL_height()` + `proc_PBL_height()` (see `docs/legacy_review.md`). `qc_flag`: 0 as computed;
  1 floored (`pblh_m == 50`); 2 run-mean fill. Fill detection, generic for all sites (the run boundaries are not
  known elsewhere): a value (rounded to 0.01 m) other than 50 that occurs in ≥ 50 rows **and** ≥ 0.5 % of the
  non-missing rows is a fill value; report the values found and their counts. At BE-Bra this must find exactly
  1654.58, 843.07, 1591.68, 682.86 m (test on the real data in the script, not in unit tests). Excluded rows dropped.

## 6. Comparison (`R/compare_pblh.R`)

- `match_halfhours_to_hours(std_hh, hours_utc, require_both = TRUE)`: for each ERA5 instant `t`, mean of the
  half-hour ending at `t` and the half-hour starting at `t`; `n_samples` = number of valid half-hours used; with
  `require_both = TRUE` the value is `NA` unless both are valid. Regime at the hour (`x_regime`): the common regime
  if equal, otherwise `"mixed"`. `qc_flag` at the hour: max of the two.
- `match_halfhours_to_accum(hh, hours_utc)`: for ERA5 accumulations over (`t − 1 h`, `t`]: mean of the two
  half-hours **ending** at `t − 30 min` and `t` (used only for EC H vs ERA5 H).
- `pblh_metrics(mod, ref)`: `n`, `bias = mean(mod − ref)`, `mae`, `rmse`, `r_pearson`, `r_spearman`, and in log
  space (only `mod > 0 & ref > 0`, number excluded reported): `median_log_ratio = median(log(mod/ref))`,
  `mad_log_ratio`, `frac_within_factor2`. Errors in absolute and log-ratio form, as in the plan.
- `pblh_metrics_by(df, mod, ref, by)`: metrics per stratum.
- `diurnal_summary(df, value_cols, by = c("season", "hour_utc"))`: median, p25, p75, n.
- Strata: all; PBLHmod regime at the hour (`stable`, `unstable_diag`, `convective`, `mixed`); season
  (DJF, MAM, JJA, SON); hour UTC; operational PBLH with `qc_flag == 0` only vs all rows.

Comparisons to compute (pairs `mod` vs `ref`):
1. `pblhmod_default` vs `era5_blh_center` (main);
2. `pblhmod_default` vs `era5_blh_bilinear` (cell sensitivity);
3. `pblhmod_gamma0025` vs `era5_blh_center` (γ sensitivity);
4. `icos_oper` (qc 0, and all) vs `era5_blh_center`;
5. `pblhmod_default` vs `icos_oper` (qc 0).

## 7. Exploratory extras (descriptive only)

- **EC H vs ERA5 H**: EC `H` as used by PBLHmod, recomputed from the model output as
  `H = wtheta * rho * cp` with `rho = p_Pa / (r_d * T_K)` (same constants as spec 02), matched with
  `match_halfhours_to_accum()`; ERA5 `H_era5_Wm2` of the tower cell. Metrics overall and
  by season; mean diurnal cycles.
- **Cloud base**: ERA5 `cbh_m` at the hour, classes `no_cloud` (NA), `cbh_le_2000` (≤ 2000 m), `cbh_gt_2000`.
  Distribution of `log(PBLHmod/ERA5)` per class and regime.
- **Free-atmosphere θ gradient (γ context)**: from the 6-hourly pressure levels of the tower cell,
  `gamma_850_700 = (theta_700 − theta_850) / (zg_700 − zg_850)` [K m-1]; monthly median and IQR at 12 UTC
  (and 00 UTC), as a table and a figure. Compare visually with the constant 0.01 and 0.0025 K m-1 used in
  PBLHmod. Not used as model input in this spec.

## 8. Figures (`R/plot_compare.R`, `R/plot_utils.R`)

`save_figure(p, name, width_cm, height_cm)` → `output/figures/<name>.png` (300 dpi) and `.pdf`. Publication style
(theme_bw or similar, readable fonts, units on axes, colour-blind-safe palette; one palette for the three sources
across all figures). Figure functions return ggplot objects; no `windows()`.

1. `04_diurnal_by_season_BE-Bra`: median and IQR ribbon per hour, 4 season panels; PBLHmod default, ERA5 centre,
   operational (qc 0) and operational (all, dashed).
2. `04_scatter_mod_era5_BE-Bra`: 2-D density (`geom_bin2d`), log10 axes, 1:1 and factor-2 lines, panels by regime.
3. `04_logratio_by_hour_BE-Bra`: boxplots of `log(PBLHmod/ERA5)` per hour UTC, by season.
4. `04_era5_spread_BE-Bra`: distribution of `x_blh_cv9` per hour; and centre vs bilinear scatter.
5. `04_H_ec_vs_era5_BE-Bra`: scatter and mean diurnal cycles.
6. `04_logratio_by_cloud_BE-Bra`: boxplots per cloud class × regime.
7. `04_gamma_850_700_BE-Bra`: monthly γ climatology with reference lines at 0.01 and 0.0025 K m-1.

## 9. Script `scripts/04_era5_compare.R`

For each site (`cfg$test_sites` or arguments):
1. read site metadata, exclusions, ICOS site data (spec 01), ERA5 (all three families);
2. cell metadata → `<derived_dir>/era5/era5_cell_metadata.csv` (one row per site; replace that site's row);
3. read PBLHmod outputs from `<derived_dir>/pblhmod/<site>_pblhmod_<label>.csv.gz` (`default`, `gamma0025`);
   if missing, abort with a message to run `scripts/02_run_pblh_model.R` first;
4. standardize and write `pblhmod_default`, `pblhmod_gamma0025`, `era5_blh_center`, `era5_blh_bilinear`,
   `icos_oper` to `<standard_dir>/<site>/`;
5. comparisons (section 6) → `output/tables/04_metrics_<site>.csv` (long format: pair, stratum, metric, value)
   and `04_diurnal_<site>.csv`; extras (section 7) → `04_H_metrics_<site>.csv`, `04_gamma_850_700_<site>.csv`;
6. figures (section 8);
7. print a summary: periods of each source and of the overlap; rows dropped by exclusions per source; number of
   matched hours per pair (and lost to `require_both`); the fill values detected in `icos_oper` with counts; the
   main metrics table (pair 1 by regime and season); the cell metadata row; run time.

## 10. Tests

Synthetic fixtures for `XX-Tst` built by `make_fixtures.R` (small, no real data): single-level files for the three
variables × 9 cells × 48 h; pressure-level T and z × 9 cells × 8 times × 3 levels; static z/lsm/sdor × 9 cells;
headers in the real format. Use the real naming conventions for each family (direction flags vs signed offsets).

- Header parsing (keys, numeric conversion, missing end marker aborts).
- Cell positions from header for both conventions; a fixture whose name contradicts its header aborts;
  missing cell warns; duplicate cell aborts.
- `n_rows` / first/last time mismatch aborts.
- sshf: `-360000 J m-2` → `H_era5_Wm2 = +100 W m-2`, interval (`t−1h`, `t`]; a fixture with flipped sign fails
  the diurnal sign check.
- θ known value: `T = 280 K`, `p = 850 hPa`, `r_d/cp = 287.05/1005 = 0.28562` →
  `θ = 280 × (1000/850)^0.28562 = 293.30 K` (show the calculation; tolerance 0.01 K); `zg = z / 9.80665`.
- Bilinear: field linear in lat and lon reproduced exactly at an arbitrary point; weights sum to 1.
- Spread: hand-computed sd/cv on 9 values; `n9 < 9` handled.
- Standard format: write/read round trip (values, timezone, sidecar); validation errors (missing required column,
  start > end, duplicated times, invalid `method_class`); `metadata_completeness` hand-computed;
  instantaneous rows with start == end accepted.
- Standardizers: PBLHmod qc mapping; operational fill detection on a synthetic series (floor rows, two constants,
  realistic noise) finds exactly the constants; excluded rows dropped and counted.
- Matching: hour mean of the two adjacent half-hours; `require_both`; `"mixed"` regime; accumulation matching uses
  the two half-hours ending at `t − 30 min` and `t`.
- Metrics: hand-computed case for every metric; log metrics skip non-positive values and report them.
- Old ERA5 reader tests and fixtures removed.

## Acceptance criteria

- `testthat::test_dir("tests/testthat")` passes (all tests).
- `Rscript scripts/04_era5_compare.R` runs on BE-Bra locally, writes the standard files, tables and figures, and
  prints the summary; the operational fill detection finds the four BE-Bra constants.
- Units documented in roxygen for every argument and returned column; conversions only at the reader boundary;
  every formula with its reference.
- `renv::snapshot()` done; nothing under `data/` or `output/` committed.

## Out of scope

γ from ERA5 as a model input (options 2–3); measured PBLH; other sites and cluster runs; statistical modelling
(GAMs), triple collocation, uncertainty; any change to PBLHmod or to its raw output; footprints.

## Summary to report back

As in `CLAUDE.md`, plus: the ERA5 cell metadata for BE-Bra; matched hours per pair; pair-1 metrics by regime and
season; the fill values detected; the monthly γ_850–700 range; anything physically odd (e.g. sign of the ERA5 vs
PBLHmod bias by regime compared with `literature_review_pblh.md` section 2).
