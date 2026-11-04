# Spec 09: cluster deployment — all stations, reference-year footprints, batch infrastructure

Phases 1–5 (network scale). Branch: `spec-09-cluster` (from `main`).
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/project_plan.md` ("Scope"), `docs/specs/08_footprint_sensitivity.md`
(sections 3, 4b, 4c and the storage analysis), `docs/specs/07_measured_pblh_cloudnet.md`,
`docs/specs/11_gamma_era5_davos.md`. Agreed with G on 2026-10-28; γ strategy added 2026-10-30.

## Goal

Run the pipeline at the 40 class 1–2 ICOS ecosystem stations (`config/sites.csv`; 39 since amendment 1: IT-SR2
not run) on the Unitus cluster (SLURM): inputs check, PBLHmod, ERA5 comparison, footprint validity and coverage (full period), and the footprint sensitivity
on one reference year per station. Claude Code implements and tests on G's PC (light work only, no heavy runs:
RAM issue on the PC, see STATUS). **G runs everything on the cluster** with Cowork's runbook
(`docs/local/cluster_runbook.md`, untracked). The repo is public: no cluster paths, user names, host names or e-mail
addresses in tracked files (config.local.yml and an untracked `slurm/local.env` hold them).

## Cluster facts (G, 2026-10-28)

SLURM, partition set by G (local setting); R 4.4.1, renv works, packages can be installed; storage sufficient;
code in a home `dev/pbl-analysis` clone, data under a separate big-data project folder.
Partition facts (`sinfo`/`scontrol`, 2026-10-29): 5 nodes × 48 CPUs, ≈ 240 GB memory per node, **no wall-time
limit**, shared with other users. G usually requests 6 CPUs per job; `local.env` sets CPUs and memory per job
(default 6 CPUs, 4 GB per CPU). Chunks are kept short anyway (target ≤ 2 h, option `chunk_target_min`) so that a
failure costs little and the array spreads well across nodes; restartable as in spec 08.

## Part A — infrastructure (implemented on the PC, tested without a cluster)

### A1. Configuration

- New keys in `config.yml`:
  - `oneflux_bif_dir`: default `"{icos_l2_dir}/{site}_output_L2I_2026"` (BIF next to the L2 zips, as on the PC);
    cluster `"LOCAL"`. On the cluster the BIF files live in one ONEFlux run folder per site and processing run, so the
    value may contain a glob (`*`). Resolution: list the matching `ICOSETC_{site}_FLUXNET_BIF_*.csv`; keep those whose
    `PRODUCT_NAME` years match the HH FLUXNET zip in use (`read_oneflux_bif()` check); if several remain, take the
    newest processing timestamp (from the folder name or `PRODUCT_PROCESSING_DATE`) and report the choice. Abort if
    none matches.
  - `measured_dir`: default `"data/raw/measured"` (spec 10 inputs); cluster `"LOCAL"`.
  - cluster profile: `target_areas_dir`, `standard_dir`, `derived_dir`, `output_dir` from `"TODO"` to `"LOCAL"`.
- Update `config.local.yml.example` (placeholders only) and the reader of the BIF to use `oneflux_bif_dir`.
- Readers must find FLUXNET HH zips with any year span (`..._HH_<y1>-<y2>_v*.zip`); several matches → abort with a
  clear message.

### A2. Per-site driver and run status

- `scripts/run_site.R SITE --steps=01,11a,02,04,08a [--force] [--with-oper]`: runs the steps in order for one site, each
  in a fresh R process (`callr` or `Rscript`), logging to `{output_dir}/logs/{site}/{step}_{timestamp}.log`.
- Status file per site × step in `{derived_dir}/run_status/{site}__{step}.json`: start, end, exit code, git commit,
  R version, host name (written only to the untracked derived folder), inputs hash. A step whose status is "ok" with
  the same inputs hash is skipped unless `--force`.
- A failing step stops the site and is reported; other sites continue (one SLURM array task per site).

- **γ strategy (G, 2026-10-30).** Two PBLHmod variants at every site: **PBLHmod** (`pblhmod_default`, EC data
  only, constant γ: the code default) and **PBLHmod-clim** (`pblhmod_clim`, site × month γ climatology from ERA5,
  spec 11). Step `11a` computes the site's γ climatology before step 02; step 02 writes both variants (plus the
  existing sensitivity cases). The EC-only constant is **revised after the network γ climatology** (decision point
  in Part B), so the constant is an option read from `config.yml` (`pblh_gamma_default`, now 0.01), not hard-coded.
- **Reference cases for the supplement (G, 2026-10-30):** the two constants found in the literature and in practice
  are kept as named standard sources at every site, written by step 02/04 and compared like the main variants:
  `pblhmod_gamma001` (γ = 0.01 K m⁻¹, Batchvarova & Gryning 1991 as cited by Kljun et al. 2015) and
  `pblhmod_gamma0025` (γ = 0.0025 K m⁻¹, the value of the legacy ETC code; check whether it also has a literature
  source and cite it if so). Tables always; supplementary figures (same layout as the main ones, labelled "reference
  γ values"); not in the main figures.

### A3. SLURM templates (tracked, no machine data)

- `slurm/sites_array.sbatch`: array over the row index of a site list (`config/run_sites.txt`: 39 stations, all of
  `config/sites.csv` except `config/sites_not_run.csv`); calls `scripts/run_site.R`. `slurm/footprint_chunks.sbatch`: array over the rows of a chunk list (A5).
- Machine data (partition, mail address, R setup line, memory, wall time) from an untracked `slurm/local.env`,
  template `slurm/local.env.example`; `.gitignore` updated. Templates use `--cpus-per-task` from it (default 6).
- `scripts/slurm_submit.R --what=sites|chunks [--dry-run]`: writes the list, checks `local.env`, prints or runs the
  `sbatch` command. `--dry-run` renders the final script for review.

### A4. Network collection

`scripts/09_collect.R`: reads the per-site outputs and writes network tables `{output_dir}/tables/09_network_*.csv`:
input availability and reference period (01), γ climatology per site × month (11a) and the network summary
(median, p25, p75 per month and overall, sites with only 2 fit levels), PBLHmod and PBLHmod-clim summary per regime (02), u* filter / L check counts, PBLHmod vs
ERA5 metrics (04, tower cell and bilinear), ERA5 cell distance and orography, footprint coverage per case × set ×
subset (08a), reference-year table (spec 08 section 4c), run status. Flags sites with problems (missing inputs,
unit-check failures, > 5 km values, low coverage). Figures: network maps are out of scope; a dot plot of median ln
ratio PBLHmod/ERA5 by site and regime.

### A5. Footprint run on the cluster

Settings as spec 08 (ETC grid ±1000 m at 1 m, umean only, smoothing EBImage-compatible), with:
- **Storage `c5`** (G, 2026-10-28): summed grids stored as 5 m × 5 m block sums (mass-conserving, tower cell at a
  block centre), float32 + gzip, for all validity sets; per-half-hour metrics computed at 1 m and stored as tables.
  Monthly sums kept, annual climatologies built from them (sum, upsample, smooth, normalise).
- **Cases:** `pblhmod` (with the revised EC-only constant), `pblhmod_clim` and `era5` at every site; the reference
  cases `pblhmod_gamma001` and `pblhmod_gamma0025` and the operational ICOS PBLH (`oper`, i.e. `--with-oper`) at the
  extra sites below, for the supplement only (tables and supplementary figures; never in the main figures); `measured` where a measured source exists; `era5_bilinear` and the
  ±20 % anchors (`pblhmod_x0.8`, `pblhmod_x1.2`) at the sites listed in `config/footprint_extra_sites.txt` (default:
  BE-Bra, FI-Hyy, FR-Gri, CH-Dav, SE-Svb, IT-Cp2: one per climate/canopy type; G may edit). `oper` runs by default
  at the extra sites (supplement) and elsewhere only with `--with-oper` (internal outputs).
- **Sets/subsets:** `etc`, `kljun`, `kljun_norsl`, `none`; `valid`, `valid_qc`, `valid_qc_used`, `valid_qc_h`.
- **Reference-year rule** (proposal written by `09_collect.R`): as spec 08 section 4c — the calendar year with the most
  half-hours usable in both main cases (`pblhmod` and `era5` h present, FFP inputs present, valid under `etc`, and in
  `valid_qc_used`), within the FLUXES reference period and ERA5 (2020 →), with ≥ 9 months of data; within 2 %,
  prefer a year overlapping the site's measured source, otherwise the most recent.
- **Measured case:** the site's primary measured source from `config/measured_sites.csv` (not only Cloudnet). Aerosol
  sources (SIRTA, MeteoSwiss, later Lanzhot) only on their comparable set (daytime window, PBLHmod convective), as in
  specs 10–11; turbulence sources (Cloudnet) exact values plus the `bound` case, as in spec 08.
- **Chunks:** `scripts/08c_make_chunks.R` writes the list site × case × month for each site's reference year
  from the approved reference years only: `09_collect.R` writes the proposal (`output/tables/08_reference_years.csv`);
  after review, the approved list is committed as `config/reference_years.csv` (site, year, reason, approved); 08c
  aborts if a site has no approved year. Each chunk carries an estimated duration from the spec-08 timings;
  months are split into half-months when the estimate exceeds `chunk_target_min`. Each array task runs one chunk with
  the job's cores as workers; chunk outputs keyed by inputs hash (restartable, as spec 08).
- After all chunks of a site: `08_footprint_compare.R SITE` (reference year), then `09_collect.R`.
- Memory per worker from the spec-08 measurements (≈ 2.5 GB with c5 sums); `local.env` sets `--mem`.

### A6. Tests (no cluster needed)

BIF glob resolution (several candidates, product mismatch, newest processing run); FLUXNET HH year-span patterns;
run-status skip logic (same / changed hash, `--force`, failure stops the site); chunk list generation (reference
year approved / not approved, splitting); sbatch rendering in `--dry-run` with a fake `local.env` (no machine data in
the rendered output beyond what `local.env` provides); c5 block sums (mass conservation, tower-centred blocks,
annual = sum of months) — reuse spec-08 tests where they exist.

## Part B — runs on the cluster (G, with the runbook)

1. Setup (runbook): pull, `renv::restore()`, `config.local.yml` and `slurm/local.env` on the cluster, Kljun code and
   target-area KML copied to untracked folders, data folders created, test suite as a small job.
2. **Regression:** BE-Bra and FI-Hyy steps 01/02/04 on the cluster; compare with the PC outputs (copied over):
   PBLHmod max |diff| < 1e-9 m, identical NA/regime/method, identical pair-1 metrics. Same inputs required (check
   the md5 of the L2 zips and BIF on both machines; if they differ, report and stop).
3. All sites, steps 01, 11a, 02, 04, 08a (array). Review `09_network_*.csv` with Cowork: input problems, per-site
   specifics, the reference-year table → G approves the years.
3b. **Decision point (before any footprint run):** from the network γ table, G chooses the recommended EC-only
   constant (e.g. the network median, rounded) and the γ tables to publish (ERA5-derived: Copernicus attribution).
   Set `pblh_gamma_default`, commit `config/gamma_era5_clim.csv` (site × month), rerun step 02 and 04 at all sites
   (cheap; γ = 0.01 and 0.0025 stay as the reference sources `pblhmod_gamma001` / `pblhmod_gamma0025`), collect
   again.
4. Footprint chunks (array) on the reference years, then per-site compare and collection.
5. Measured (after spec 10 is merged, which generalises the comparison code): Cloudnet for FI-Lom, DE-RuS, FR-Gri
   (07b download on the login node, size first; 07c; comparison), FI-Hyy Cloudnet re-downloaded on the cluster
   (checksums) rather than copied; SIRTA (FR-Gri) and later PI data copied into `measured_dir`.

## Docs

`CLAUDE.md`: cluster workflow (driver, status files, SLURM templates, local.env, chunks), the γ step and the two
PBLHmod variants, the new config keys, the
"no machine data in tracked files" rule extended to `slurm/`. One colour and label scheme for all sources
(PBLHmod, PBLHmod-clim, ERA5, measured, reference γ cases, operational) defined once in `R/plot_utils.R` and used by
every figure, footprint figures included. **Figure style (G, 2026-10-30):** in `R/plot_utils.R`, (1) one named table
`pblh_sources_style` (source id → colour, label, plotting order, line type) and (2) one `theme_pblh()` function (fonts,
sizes, grid, legend position, panel spacing); every figure takes colours, labels and theme only from these two, so G
can restyle all figures later by editing them alone. Current colours are placeholders (colour-blind-safe, e.g. the
Okabe–Ito set: PBLHmod blue, PBLHmod-clim bluish green, ERA5 vermillion, measured black, reference γ values sky blue /
orange, operational reddish purple); G will replace them with the final theme. `docs/STATUS.md`: spec 09 row. `docs/project_plan.md`:
execution on the cluster, cases and extra-site list. The runbook stays in `docs/local/` (Cowork).

## Amendment 1 (2026-11-03): fixes from the first network run (job 11449, 2026-11-02)

1. **ERA5 longitudes.** Static files (and possibly others) give 0-360 longitudes at sites west of 0.125 E (FR-Bil,
   FR-Lus, UK-AMo, GL-ZaF): the tower cell was not found, and at FR-Lus the western neighbour (359.75) became col
   1439. `find_era5_files()` converts header longitudes to -180..180 (`era5_lon_180()`: values already in range are
   unchanged, bit-identical) and takes column offsets as wrapped differences, for every family; the ERA5 provenance
   text in the standard-file metadata uses the same convention.
2. **Polar sites (11a, 02).** Days without sunrise have no morning profile (`gamma_daily()` flag `polar_night` /
   `polar_day`; counts per month in the monthly table, `n_polar_night`, `n_polar_day`, and in the 11a log). A
   month without a valid daily γ is `NA` in the climatology (reported by 11a); PBLHmod-clim then uses the constant
   γ (`pblh_gamma_default`, `gamma_source_used = "constant_fallback"`, spec 11 B). With `gamma_scope = "all"` this
   constant enters N of Eq. B3 on those rows; B&G is not run on polar days (spec 02 Amendment 2). The SE-Sto failure
   was the missing-months warning (cli quantity from a vector of months).
3. **Exit status.** `scripts/run_site.R` exits with `run_site_exit_status()` (1 if any step failed) after the
   summary table; `slurm/local.env.example`: `SLURM_MAIL_TYPE=FAIL,ARRAY_TASKS` (per-task mails in an array).
4. **IT-SR2 not run** (no data in the L2I 2026 release, station work in progress; G 2026-11-03): removed from
   `config/run_sites.txt` (39 stations), kept in `config/sites.csv`; reason in `config/sites_not_run.csv` (config
   key `sites_not_run_file`). `slurm_submit.R` refuses a listed station; `09_collect.R` lists it in
   `09_network_flags.csv` as `not_run`, "not run (no L2 data)", never as a missing step.

## Amendment 2 (2026-11-03): after the first network collection (job 11546)

A. **Network figure** `09_network_logratio_*` (G 2026-11-03): PBLHmod as filled dots with a bar from p25 to p75 of
   the half-hourly ln(PBLHmod / ERA5 tower cell) behind each median (site × regime); PBLHmod-clim as open circles, no
   bar, drawn on top. Caption: stable values do not depend on γ (PBLHmod-clim = PBLHmod; max difference 0 at all 39
   sites in job 11546). `pblh_metrics()` gains `p25_log_ratio`, `p75_log_ratio` (type-7 quantiles, appended columns),
   so step 04 must rerun (it does anyway: `R/` changed). Functions `network_logratio_table()`,
   `plot_network_logratio()`; old metrics without the quartiles give a figure without bars and a warning.
B. **GL-ZaF** (0 `valid_qc_used` of 12,363 `valid_qc` half-hours). Cause (G, 2026-11-03, from the cluster file):
   `NEE_VUT_REF_QC` = −9999 on all 64,272 rows of the FLUXNET file (ONEFlux could not compute the VUT NEE: missing
   critical data), while `NEE_CUT_REF_QC` is filled (QC 0: 11,524; 1: 28,376; 2: 5,683; 3: 4,614; −9999: 14,075).
   **Rule for `valid_qc_used`** (G 2026-11-03; general, not site-specific): if `NEE_VUT_REF_QC` is entirely missing at
   a site (column absent, or −9999 on every FLUXNET row read), use `NEE_CUT_REF_QC == 0` instead (`nee_ref_qc()`; the
   reader now also keeps `NEE_CUT_REF_QC`). A site with any VUT QC value is unchanged (test: identical
   `valid_qc_used` on a VUT fixture whose CUT column differs). Step 08a writes `08_nee_qc_<site>.csv`
   (`nee_qc_status()`, one row for the whole period, `year` NA, and one per UTC year): `nee_qc_source` (`VUT` / `CUT` /
   `none`), VUT and CUT columns present, number of VUT QC values, and for the source the QC values, QC 0, overlap with
   FLUXES NEE QC good, `status` `ok` / `column_absent` / `all_missing` / `none_measured` / `no_overlap` /
   `no_fluxnet`. GL-ZaF FLUXNET covers 2023–2026 only, so its 2022 row keeps `no_fluxnet` (no `valid_qc_used` in
   2022). `09_collect.R` collects these as `09_network_nee_qc.csv`; `09_network_flags.csv` flags `nee_qc_cut` (CUT
   source) and `nee_qc` (whole-period status not `ok`). The u* check of 08a (`08_ustar_check`) still compares the VUT
   thresholds with `NEE_VUT_REF_QC`, so it has no rows at a CUT site. A reader test covers a GL-ZaF-like FLUXNET
   file (from 2023, offset 0) for the CUT fallback and each status.
   **Reference year** (G): GL-ZaF stays in the footprint run. `config/reference_year_exceptions.csv` (config key
   `reference_year_exceptions_file`) lists `site_id, rule, reason, decided`. Rule `coverage` applies only when the
   standard rule finds no year: the year with the most `valid_qc` half-hours (both main cases, set `etc`, any number
   of months; ties go to the most recent year), with the reason written in the proposal (`rule` column).
   `footprint_year_counts()` gains `n_valid_qc` and `n_months_valid_qc` (needs 08a).
C. **Measured data on the cluster.** New `run_site.R` steps `meas_std` (standard files of the primary measured
   source) and `meas_cmp` (per-site comparison, upstream 04 and `meas_std`), both `scripts/09m_measured.R`. The
   scripts depend on the provider (`measured_site_scripts()`): Cloudnet 07c, then 07 + 10; SIRTA 10a, then 10;
   MeteoSwiss 11b, with no per-site comparison (CH-Dav is in `11_gamma_evaluation.R`). The sub-scripts are part of the
   inputs hash (`extra_scripts`), and the raw measured files form a new input group `measured_raw`. The site list is
   `config/measured_run_sites.txt`. Order: Cloudnet download (07b: size, `--confirm`, re-check of the checksums),
   then `--steps=meas_std,meas_cmp,08a`, then `11_gamma_evaluation.R` and `09_collect.R`.

## Amendment 3 (2026-11-04)

1. **FI-Hyy reference year.** The Cloudnet MLH is counted: the cluster year counts give `n_measured` 1200 (2022), 657
   (2023), 231 (2024); only the exact values (qc 0/2) count, as in the footprint case `measured`; qc 1 (censored below
   z_min) has no value and qc 3 (signal-limited lower bound) is the case `bound`. The proposal showed 0 because the
   rule picked 2021 (3986 usable half-hours; 2022: 3510, −12 %, outside the 2 % tolerance within which a measured
   year was preferred). **Rule changed (G 2026-11-04):** among the years with ≥ 9 months, a year overlapping the
   site's measured source is preferred when its usable half-hours are ≥ 80 % of the best year's
   (`propose_reference_year(measured_min_share = 0.8)`); among those, the most measured half-hours (ties: most
   recent); otherwise the year with the most usable half-hours as before (2 % tolerance, most recent). FI-Hyy → 2022
   (88 % of 2021, 1200 measured). CH-Dav and FR-Gri may also move to another measured year with more measured
   half-hours (their year counts were not available on the PC): check the new proposal.
2. **Footprint cases.** `pblhmod_clim` was already a case at every site (spec 09 A5: 08a builds it from the standard
   file written by step 04, 08c takes every case of the input table), next to `pblhmod` (γ 0.01) and `era5`.
   Decision point 3b stays open: if the EC-only constant changes, the `pblhmod` chunks rerun (`pblhmod_clim` only in
   months that use the fallback constant). Estimate for the 39 proposed reference years of job 11546 (spec-08 timing
   2 s per FFP half-hour, 2 CPUs per chunk; half-hours = `n_valid` of set `none`, full period, scaled to one year):
   ≈ 1,840 chunks (12 per site × case; no month exceeds the 120 min split threshold, the largest ≈ 1,650 half-hours
   ≈ 28 min), ≈ 1.87 M FFP half-hours, ≈ 1,040 CPU-hours: `pblhmod`, `pblhmod_clim`, `era5` ≈ 470 chunks and
   ≈ 262 CPU-hours each; the six extra sites' cases (`era5_bilinear`, ±20 %, γ 0.01 / 0.0025, `oper`) 432 chunks,
   ≈ 261 CPU-hours; `measured`/`bound` at FI-Hyy (now 2022), FR-Gri, CH-Dav a few dozen chunks, a few CPU-hours.
   More than 1000 chunks: `slurm_submit.R --what=chunks` submits two arrays (MaxArraySize 1001).
3. **γ strategy (G 2026-11-04).** Default `pblhmod_default` = constant 0.01 K m-1 (EC only); PBLHmod-clim = ERA5 site
   × month climatology. Side results, **step 02 only, no footprints**: `gamma_med` = one constant for all stations,
   the median over the stations of each site's median daily ERA5 γ (`gamma_network_median()`, default fit layer; every
   station weighs the same; G 2026-11-04; the median of the daily values pooled over every site-day is kept as
   `gamma_med_pooled` for reference), and `gamma_site` = one constant per site, the median
   of the site's daily ERA5 γ (`gamma_site_constant()`). `09_collect.R` writes `gamma_med` to
   `<derived_dir>/gamma/network_gamma.csv` (only when every run site has its 11a table, and only if the value
   changed) and `09_network_gamma_site.csv`; step 02 reads the file (input group `gamma_network` of its inputs hash)
   and skips `gamma_med` while it is missing. Step 04 compares both with ERA5 (pairs `7_pblhmod_gamma_med_vs_era5_center`,
   `8_pblhmod_gamma_site_vs_era5_center`; no standard files); `09_collect.R` adds them to `09_network_era5_metrics.csv`,
   `09_network_regimes.csv` and the supplementary figure `09_network_logratio_supp_<n>sites` (PBLHmod with IQR bar,
   PBLHmod-clim, `gamma_med` filled triangles, `gamma_site` open triangles). Paper: a recommended constant and/or the
   γ tables for the ICOS stations as a side result.
   **Two-pass dependency and cluster order:** `gamma_med` needs the 11a results of every site. Order: (i) site array
   `--steps=01,11a`; (ii) `09_collect.R --gamma-only` (writes the network γ, seconds); (iii) site arrays with the
   default and measured steps (01 and 11a are skipped as current; 02 runs `gamma_med` and `gamma_site`); (iv)
   `11_gamma_evaluation.R`, `09_collect.R`. If the network γ changes later (e.g. a site added), the collection rewrites
   the file and the next site run reruns 02, 04 and 08a; footprint chunks are keyed by their content (rows, groups,
   code), so finished chunks stay valid.
   Check on BE-Bra (PC, run with the 4-site pooled γ 0.00463 before the switch to the median of the site medians,
   0.00469; BE-Bra 0.00473 K m-1): median ln(x / ERA5 tower cell), convective,
   default −0.07, `gamma_site` +0.18, `gamma_med` +0.19; stable identical in all cases.

## Reference years approved (G, 2026-11-04)

- `config/reference_years.csv`: the proposal of the cluster collection of 2026-11-04 (`08_reference_years.csv`, 39
  run sites), approved unchanged; columns `site_id`, `year`, `reason` (proposal text), `approved` (`G 2026-11-04`).
  All 39 years follow the standard rule (no `coverage` exception was needed, GL-ZaF included); the measured-overlap
  rule chose FI-Hyy 2022, FR-Gri 2023 and CH-Dav 2022. Years: 2020 (2 sites), 2022 (6), 2023 (8), 2024 (9),
  2025 (14).
- `read_reference_years()` (`R/footprint_chunks.R`) reads the list for `08c_make_chunks.R` and step `08`
  (`reference_year_of()`); it aborts if the file or a column (`approved` included) is missing, a site appears twice,
  a year or approval is empty, or a requested site is not listed.
- Footprint cases per site, from the code (08a `footprint_h_cases()`, 08c/08b case filter) and confirmed by the
  cluster 08a coverage table: `pblhmod`, `pblhmod_clim`, `era5` at all 39 sites; at the six extra sites also
  `era5_bilinear`, `pblhmod_x0.8`, `pblhmod_x1.2`, `pblhmod_gamma001`, `pblhmod_gamma0025` and `oper`; `measured` at
  FI-Hyy, FR-Gri, CH-Dav and `bound` at FI-Hyy (Cloudnet, turbulence). `oper` at the extra sites needs no
  `--with-oper` (08a, 08c, 08b and 08 include it there; `--with-oper` would add it at every site, internal only).
- Estimate (cluster year counts and coverage, 2 s per half-hour): about 1,880 chunks (no half-months: the longest
  month is far below `chunk_target_min`), 1.85 M FFP half-hours, about 1,030 CPU-hours at the nominal cost; up to
  about 2,000 CPU-hours if the PC cost per worker (4–6.6 s including the per-chunk climatologies) holds on the
  cluster. Two arrays (> 1000 tasks). Disk: c5 block sums ≈ 0.42 MB per group (measured on a BE-Bra test chunk),
  ≤ 47 GB in total (fewer where the `kljun` set is empty, as at forest sites), climatologies ≈ 0.7 GB, half-hour
  chunk tables ≈ 0.25 GB under `derived/`; step 08 tables and figures ≈ 0.5 GB under `output/`.

## Acceptance criteria

- Part A: tests pass on the PC; `--dry-run` renders valid sbatch scripts; `git grep` for `/bigdata`, `/home/`,
  `@unitus`, `g.nicolini`, `c.trotta` returns nothing outside `docs/local/` (untracked) and placeholders.
- Part B (reported by G): BE-Bra/FI-Hyy regression on the cluster passes; network tables for all sites; reference
  years approved; footprint outputs for every site's reference year; collection tables.

## Out of scope

Network maps; Monte Carlo (phase 5b); own MLH retrievals; Python port; changes to PBLHmod other than the γ constant.

## Summary to report back

Part A as in `CLAUDE.md`, plus: the dry-run sbatch scripts (rendered with the example env), estimated chunk count and
CPU-hours for all sites (from spec-08 timings), and anything that needs a cluster fact to be confirmed.
