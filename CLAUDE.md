# CLAUDE.md: pbl-analysis

Standing instructions for Claude Code. Read this first, every session.

## Project in one paragraph

Estimate Planetary Boundary Layer height (PBLH) from eddy-covariance (EC) data at ICOS ecosystem stations
(PBLHmod: Kljun et al. 2015 App. B diagnostics + Batchvarova & Gryning 1991 slab model), evaluate it against
ERA5 and measured PBLH (lidar, ceilometer, radiosonde...), explain when/why they disagree, and quantify the
impact of PBLH errors on FFP footprint estimates (Kljun et al. 2015). The work will become a peer-reviewed paper,
so code must be reproducible, documented and traceable to references.

- Current status and next steps: `docs/STATUS.md`
- Full design: `docs/project_plan.md` (source of truth for scientific choices)
- Legacy code issues: `docs/legacy_review.md`
- L2 data issues for the ETC QC report: `docs/local/etc_qc_report.md` (local working file, not tracked)
- Task specs: `docs/specs/NN_*.md`. **Implement only what the current spec asks.**
- Literature: `docs/references/` (`pblh_collection.bib`, `pblh_literature_2026-10.bib`, annotated notes in
  `literature_review_pblh.md`). Cite by BibTeX key in code comments where a formula comes from a paper.

## Golden rules

1. **Never change scientific content silently.** Formulas, constants, thresholds, defaults: if a spec is ambiguous
   or you think something is physically wrong, stop and ask, or implement as specified and flag it clearly in your summary.
2. **Units everywhere.** Every function argument and returned column documents its unit in roxygen (`@param H Sensible heat flux [W m-2]`).
   Internally use SI: K for temperature, Pa for pressure, m, s. Convert at the reader boundary only.
3. **References in code.** Every formula carries a comment with its source (author, year, equation number when known).
4. **No hardcoded paths.** All paths come from `config.yml` (+ untracked `config.local.yml`) via `get_cfg()` (see Environment).
5. `legacy/` is read-only reference. Never edit, source, or import from it in new code. Exception (spec 03): only
   `scripts/03_legacy_regression.R` and `tests/testthat/test-legacy_port.R` may source `calc_PBL_height_operational.R`,
   `calc_PBL_height_developing.R` and `proc_PBL_height.R`, into dedicated environments. Kljun's FFP code in `legacy/`
   is local only (gitignored).
6. Never commit data, credentials or outputs (see `.gitignore`). Never print secrets.

## Environment

- R 4.4.1, both locally (Windows) and on the cluster (Linux, SLURM). Packages managed by `renv`:
  after adding a dependency run `renv::snapshot()`. Do not call `install.packages()` or `library()` inside functions.
  EBImage (Bioconductor) is optional, used only by `scripts/08v_ffp_verify.R` (comparison with the FFP code's
  smoothing; skipped if not installed), and ignored by renv (`ignored.packages` in `renv/settings.json`, G
  2026-10-28): run that check on the cluster.
- Paths: `config.yml` with profiles `default` (local PC) and `cluster`. Select with env var `R_CONFIG_ACTIVE=cluster`.
  Access with `cfg <- get_cfg()` (`R/utils_config.R`).
- **The repository is public: no absolute machine paths or user names in tracked files.** Machine-specific values
  are `"LOCAL"` in `config.yml` and are set in the untracked `config.local.yml` (gitignored; template
  `config.local.yml.example`), which `get_cfg()` merges over `config.yml` (its `default:` section, then the active
  profile). `resolve_path()` aborts on `"LOCAL"` and `"TODO"`.
- Development and testing happen locally on ONE site (BE-Bra); production loops over all ICOS stations on the cluster
  as a SLURM array job (one task per station). Code must run identically in both places.

## Input data (facts; do not guess beyond these)

- ICOS L2, per station `XX-Yyy`, in `<icos_l2_dir>/XX-Yyy_output_L2I_2026/`:
  - `ICOSETC_XX-Yyy_FLUXES_INTERIM_L2.zip`: H, USTAR, MO_LENGTH, ZL, T_SONIC, WS, WD, V_SIGMA, PBLH (legacy operational output), FETCH_*, FOOTPRINT_*, QC flags
  - **Reference period (G, 2026-10-22): the site-specific time range of the FLUXES file** (first to last row). Stations
    were labelled at different times, so the range differs by site; never hard-code years. FLUXNET and METEO files
    may cover a different range (e.g. FI-Hyy FLUXNET 2018–2026, BE-Bra 2020–2026): rows outside the FLUXES range are
    dropped. Inside it, the forcing takes `TA_F`/`PA_F`/`H_F_MDS` from FLUXNET where FLUXNET covers the row, and
    METEO `TA`/`PA` after the end of FLUXNET (at BE-Bra: from 2026-08-22 23:00 UTC). Comparisons with other sources
    use the overlap with this period (ERA5 starts 2020-01-01).
  - **BE-Bra FLUXNET (G, 2026-10-26)**: `ICOSETC_BE-Bra_FLUXNET_FLUXMET_HH_2020-2026_v1.3_r1.zip`, a new ONEFlux run
    matching the BIF file. It replaced `..._HH_2020-2025_v1.3_r1.zip` (kept in `data/archive/`, outside `data/raw/`,
    not read by the pipeline): 2026 now comes from FLUXNET, not METEO; 2020–2025 differs only on gap-filled values.
    Pre-replacement outputs in `data/derived/pre_bebra_fluxnet2026/` and `output/pre_bebra_fluxnet2026/`.
  - `ICOSETC_XX-Yyy_FLUXNET_BIF_<years>_v1.3_r1.csv` (ONEFlux BIF, replaces AUXNEE; next to the L2 zips): long format
    (`SITE_ID, GROUP_ID, VARIABLE_GROUP, VARIABLE, DATAVALUE`); u* thresholds of `NEE_VUT_REF` in `GRP_FLUX_REF`, one
    per `REF_YEAR`. Its `PRODUCT_NAME` years must match the HH file (`read_oneflux_bif()` aborts otherwise).
  - **L2 `ZL` is not (EC height − 2/3 hc) / `MO_LENGTH`** (spec 08 step 0; G, 2026-10-26): the implied z − d is
    17–18 m at BE-Bra 2020–2025 (ours 24.1 m) and steps between periods at both sites (probably a constant hc in the
    EddyPro metadata; `docs/local/etc_qc_report.md` item 11). Stability for footprints is always zm/L with
    zm = EC height − 2/3 hc, never L2 `ZL`.
  - `ICOSETC_XX-Yyy_FLUXNET_FLUXMET_HH_<years>_v*.zip`: TA_F, PA_F, H_F_MDS, H_CORR, NIGHT, SW_IN_POT, WS_F, P_F
  - `ICOSETC_XX-Yyy_METEO_INTERIM_L2.zip`: TA, TA_<n> profile (number of levels differs by station, G 2026-10-23:
    detect from the file header, never a fixed list), P, SW_IN
  - Files are zipped CSVs. Read directly from the zip (no permanent unzipping). Use only needed columns.
  - Missing value: `-9999`. Timestamps: `TIMESTAMP_START`, `TIMESTAMP_END` as `YYYYMMDDHHMM`, **local standard time (no DST)**.
    Convert to UTC with the site UTC offset from site metadata. Store both start and end.
  - Units: `TA*` in °C (convert to K); **`T_SONIC` already in K** (verified on BE-Bra, spec 01; rename only);
    `PA_F`, `PA` in kPa (convert to Pa); H in W m-2; USTAR in m s-1; MO_LENGTH in m.
    Readers check the median of each declared-unit column and abort if implausible.
  - `H2O` (FLUXES) is a water vapour mole fraction in mmol mol-1 (median 9.75 at BE-Bra; spec 03).
  - The QC of the FLUXNET files differs from that of the FLUXES file: FLUXNET additionally applies u* filtering, to
    NEE only (G, 2026-10-22); H and LE, and so H_F_MDS used by PBLHmod, are not affected.
  - ETC cleanset (G, 2026-10-22): EddyPro −9999 values in `air_density` / `air_heat_capacity` / `air_temperature`
    are converted to NA before `calc_PBL_height()` is called; the operational `PBLH` was computed from the same
    cleanset that became the published L2.
- **ERA5 (current delivery, C. Trotta, Oct 2026; verified on BE-Bra 2026-10-15).** 0.25° grid. For every variable,
  9 files per station: the cell containing the tower plus its 8 neighbours. All times UTC.
  - Hourly single levels, 2020-01-01 00:00 → 2026-10-07 23:00, one file per variable × cell (all years in one file):
    `XX-Yyy__reanalysis_era5_single_level__<var>__20200101T0000__20261007T2300__north_<n>__south_<s>__east_<e>__west_<w>.csv`
    with `<var>` = `boundary_layer_height`, `cloud_base_height`, `surface_sensible_heat_flux`.
    Local: `data/raw/era5/XX-Yyy/`; cluster path: set in `config.local.yml` (see `config.local.yml.example`).
  - Pressure levels, **6-hourly** (00/06/12/18 UTC), 2020-01-01 00:00 → 2026-10-07 18:00, 13 levels
    (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000 hPa), wide format, one column per level
    (`temperature_<p>hPa` [K], `geopotential_<p>hPa` [m2 s-2]):
    `XX-Yyy__reanalysis_era5_pressure_levels-pl__<temperature|geopotential>__20200101T0000__20261007T1800__north_..._west_<w>.csv`.
    Local: `data/raw/era5_pressure_levels/XX-Yyy/`; cluster path: set in `config.local.yml` (see
    `config.local.yml.example`).
  - Static fields, one row each (valid time 1979-01-20 12:00): `z` (surface geopotential [m2 s-2]; height = z / 9.80665),
    `lsm` (land–sea mask [0–1]), `sdor` (std of sub-grid orography [m]):
    `XX-Yyy__single-levels__<z|lsm|sdor>__197901201200__197901201200__north_<r>__south_+0__east_<c>__west_+0.csv`.
    Local: `data/raw/era5_single_levels/` (all stations in one folder); cluster: one folder per variable
    (`.../<z|lsm|sdor>/`), path set in `config.local.yml` (see `config.local.yml.example`).
  - **File format**: a `#`-commented metadata header (`# ==== METADATA ... ====` … `# ==== FINE METADATA (dati sotto) ====`),
    then a CSV with header. Time column: `time` (hourly and pressure-level files) or `valid_time_utc` (static files).
    Header keys used by readers: `point_id`, `variable`, `latitude_requested`, `longitude_requested` (tower),
    `latitude_grid`, `longitude_grid` (cell centre), `n_rows`, `first_time`, `last_time`, `levels`. Static files also
    carry `units` and `distance_km`; **hourly and pressure-level files carry no units** (units are stated here).
  - **Neighbour naming differs between file families. Never infer the cell position from the file name: read
    `latitude_grid` / `longitude_grid` from the header.** Hourly and pressure-level files use direction flags
    (`north_+1__south_+0__east_+0__west_+1` = NW cell); static files use signed offsets in `north_` (row) and `east_`
    (column), with `south_+0__west_+0` always (`north_-1__..__east_+1` = SE cell). In both families the tower cell is
    `north_+0__south_+0__east_+0__west_+0`.
  - Variables and conventions:
    - `boundary_layer_height` [m above model surface]: instantaneous at the timestamp. ERA5 diagnostic = lowest level
      where the bulk Richardson number (from the lowest model level, virtual dry static energy, 2 m wind = 0, no
      surface friction term, no parcel excess) reaches 0.25 (`ECMWF2016`, IFS Cy41r2 Part IV, Sect. 3.10.1, Eq. 3.90;
      after Seidel et al. 2012).
    - `surface_sensible_heat_flux`: **hourly accumulation [J m-2] over the hour ending at the timestamp, ECMWF sign
      (positive downward)**, not a mean. Convert at the reader: `H [W m-2, positive upward] = -sshf / 3600`, assigned to
      the interval (timestamp − 1 h, timestamp]. Verified on BE-Bra: the June–July midday maximum of upward H
      (≈ 130 W m-2) falls in the hours ending 11–12 UTC.
    - `cloud_base_height` [m above surface]: NaN when there is no cloud (≈ 13 % of hours at BE-Bra). Not an error.
  - Expected completeness: 59,328 hourly rows and 9,888 pressure-level rows per cell (no gaps at BE-Bra).
  - **Superseded format**: the earlier yearly files
    (`BE-Bra__ERA5__reanalysis-era5-single-levels__boundary_layer_height__2025__lon+0__lat+0.csv`) are no longer the
    source. Spec 04 replaced the reader (`R/read_era5.R`) and `config.yml` (`era5_sl_dir`, `era5_pl_dir`,
    `era5_static_dir`; `{var}` placeholder for the cluster static folders).
- **ACTRIS/Cloudnet (spec 07; checked on all FI-Hyy files 2026-10-22/23).** Licence CC BY 4.0. API
  `https://cloudnet.fmi.fi/api` (`/files?site=&product=&dateFrom=&dateTo=&showLegacy=true`, no paging documented,
  queried one year at a time; `/sites`); file metadata: `uuid`, `checksum` (SHA-256), `size`, `downloadUrl`,
  `instrument` (`name`, `pid`), `volatile` (no file PID/version yet), `legacy`. No citation endpoint. Pairs ICOS site ×
  Cloudnet site: `config/measured_sites.csv` (`distance_km` from the API coordinates; optional `z_min_agl_m`).
  Files under `cloudnet_dir` (`data/raw/cloudnet/{site}`; cluster `"LOCAL"`), with a provenance table
  `cloudnet_files_provenance_<product>.csv` written by `scripts/07b_cloudnet_download.R`.
  - `epsilon-lidar` (daily NetCDF4, one file per instrument): `epsilon` [m2 s-3] (height × time; missing where the
    signal is too low; no QC/SNR/precipitation field), `height` = gate centre [m **above mean sea level**], `altitude`
    [m a.s.l.] (Hyytiälä 150 m = API value; ICOS 181 m), `time` [hours since the UTC day], `rolling_window_period`
    180 s, `ray_accumulation_time` 6 s; global attributes `file_uuid`, `instrument_pid`, `references` (O'Connor et al.
    2010). Heights above ground = height − file `altitude`; gate top = centre + half the gate (30 m gates at
    Hyytiälä). Hyytiälä HALO 46: profiles every 14.23 s in vertical stares of ≈ 14–15 min every 30 min (≈ 63 profiles
    per half-hour, stare phase drifts between days); lowest valid gate 105 m agl.
  - `classification` (daily, 30 s): `target_classification` (0 clear … 2 drizzle/rain, 3 drizzle/rain + droplets …
    10), `cloud_base_height_amsl`/`_agl` [m], `rain_detected` (0/1), `altitude` [time]; height a.s.l., irregular gates
    (Hyytiälä lowest 102 m agl).
- Site metadata: `config/sites.csv` (lat, lon, altitude, UTC offset, zm, canopy height, d, z0...).
  Missing values allowed; functions must fail with an informative error only when a *required* value is missing.
  Time-varying zm / canopy height: `config/sites_timevarying.csv` (derived from `config/sites_H_EC_canopy_formBIF.csv`).
- **Exclusions**: `config/exclusion_periods.csv` lists periods that must be removed for a site (e.g. SE-Nor clearcut).
  Every pipeline step must apply them and report how many records were excluded. Never skip this silently.

## Footprints (phase 5, spec 08)

- **FFP code**: Kljun's `calc_footprint_FFP.R` (v1.3) and `calc_footprint_FFP_climatology.R` (v1.31), © N. Kljun,
  in `ffp_code_dir` (`config.yml`; default `legacy`, cluster `"LOCAL"`). Local only and gitignored: source it with
  `load_ffp_code()` into a dedicated environment; never copy, modify, commit or publish it; cite `Kljun2015`.
  Acceptance check: `git grep -n "calc_footprint_FFP <-" -- ':!docs' ':!CLAUDE.md'` returns nothing. All calls pass umean only
  (z0 = NaN; deviation from ETC, which passes both). Known code issues: `docs/local/ffp_code_notes.md`.
- **Heights**: d = 2/3 hc, z0 = 0.15 hc, zm = EC height − d (ETC convention). **Stability always from zm/L with this
  zm** and L2 `MO_LENGTH`, never from L2 `ZL` (see Input data).
- **Validity sets** (`footprint_validity()`; thresholds in `docs/project_plan.md`, phase 5): `kljun` (Kljun et al.
  2015, Eq. 27), `kljun_norsl` (Eq. 27 without 20 z0 < zm), `etc` (ETC operational), `none`; each also includes the
  FFP code's own checks and zm < h. `kljun` keeps no row at forest sites (BE-Bra, FI-Hyy): a result, not a bug.
- **Subsets** (`footprint_subsets()`): `valid`; `valid_qc` (FLUXES `NEE_DATA_FLAG == 0`); `valid_qc_used` (and FLUXNET
  `NEE_VUT_REF_QC == 0`: removes the u*-filtered rows, ONEFlux VUT thresholds from the BIF file, and the other FLUXNET
  NEE removals); `valid_qc_h` (valid and H measured and QC-good: `H_F_MDS_QC == 0`, or L2 `H` present where FLUXNET has
  no QC value; keeps the stable nights that the NEE u* filter removes).
- **Target-area fraction** = integral of the footprint over the ICOS target area relative to the full footprint (the
  FFP function integrates to 1 over the plane); `f_domain` = part inside the ±1000 m domain.

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
- Plots: `ggplot2`, publication style, saved via `save_figure()` to `output/figures/`: png by default; pdf only for paper figures (`formats = c("png", "pdf")`). Never `windows()`.

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
  Only exception: the tiny **synthetic** test fixtures in `tests/testthat/fixtures/` (zips included),
  rebuilt by `tests/testthat/fixtures/make_fixtures.R`. Never put real data there.

## When you finish a task

Summarise: what was implemented, files touched, tests run and results, assumptions made,
open questions/flags for the scientist. Update `docs/project_plan.md` "Open points" if relevant.
