# Spec 08: footprint sensitivity to the PBLH source (phase 5a) at BE-Bra and FI-Hyy

Phase 5 (first part). Branch: `spec-08-footprint-sensitivity` (from `main`).
Read first: `CLAUDE.md`, `docs/project_plan.md` ("Scope", "Footprint error propagation (phase 5)", FFP validity
criteria under "Inputs"/"Open points"), `docs/specs/standard_format.md`, `docs/specs/07_measured_pblh_cloudnet.md`
(qc convention, amendment 1), `docs/references/literature_review_pblh.md` (section 5). Agreed with G on 2026-10-23.

## Goal

How much does the choice of the PBLH change the footprint estimates, and does it matter for the data EC users keep?

1. Swap the PBLH in the Kljun et al. (2015) FFP model, all other inputs unchanged: PBLHmod vs ERA5 (vs measured MLH
   at FI-Hyy; optional operational ICOS PBLH, for internal use).
2. Effects per half-hour and on climatologies: fetch, footprint fraction inside the ICOS target area, and **coverage**
   (half-hours with a valid footprint).
3. Effects on all valid half-hours vs on the half-hours that survive QC and the u* filtering of the fluxes, i.e. the
   data actually used for fluxes and budgets (G, 2026-10-23: the largest h errors are expected at night, when most
   fluxes are removed by the u* filtering anyway).

This is a deterministic swap (phase 5a). The Monte Carlo with error distributions from all measured sites is phase 5b
(later spec).

## 1. Inputs

**FFP code.** Kljun's `calc_footprint_FFP.R` (v1.3) and `calc_footprint_FFP_climatology.R` (v1.31), © N. Kljun, local
and gitignored. New config key `ffp_code_dir` (default `legacy`, cluster `"LOCAL"`). Source them from there; never
copy, modify or commit them; cite `Kljun2015` in sidecars and docs. Read both files first and report: arguments,
internal validity checks (e.g. h > 10 m, zm < h, zm/L ≥ −15.5, u*), outputs (crosswind-integrated and 2D), `flag_err`
meaning, and `rslayer`.

**Site and turbulence inputs** (identical in every PBLH case; only h differs):
- heights from `config/sites_timevarying.csv`: EC height, canopy height hc; ETC convention d = 2/3 hc,
  z0 = 0.15 hc, zm = EC height − d (as `legacy/FFP_analysis_dev.R`, lines 336–353);
- L2 FLUXES: `WS` (umean), `MO_LENGTH` (L), `V_SIGMA` (σv), `USTAR` (u*), `WD`, `ZL`; as published (ETC practice).
  Rows rejected by the spec-05 u* plausibility rule are invalid in every case (reported, not replaced);
- target areas: `data/raw/target_areas/{site}_TA.kml` (read with `sf`, local metric projection centred on the tower,
  matching the FFP x/y in metres).

**PBLH cases** (half-hourly; same rows for all):

| case | source | half-hourly h |
|---|---|---|
| `pblhmod` | `pblhmod_default` standard file | as is |
| `era5` | `era5_blh_center` (tower cell) | linear interpolation in time between the two hourly instants to the half-hour midpoint (start + 15 min); NA unless both exist |
| `measured` (FI-Hyy only) | `cloudnet_eps_mlh` | exact values only (qc 0/2); qc 1 and qc 3 → NA here (reported; see section 5) |
| `oper` (optional, `--with-oper`) | `icos_oper` | as is; internal results only (G, 2026-10-23): outputs under `output/*/oper_internal/`, not in STATUS/project_plan figures or summaries |
| `pblhmod_x0.8`, `pblhmod_x1.2` | PBLHmod × 0.8 / × 1.2 | anchor: replicates the ±20 % test of Kljun et al. 2015 (Table 5) |

`era5_blh_bilinear` only in the tables of section 4 (`docs/project_plan.md` "Scope"), as `era5_bilinear`.

## 2. Validity, QC and u* subsets

**Validity sets** (configurable thresholds; one flag column per condition and case):

| condition | `kljun` (Kljun et al. 2015, Eq. 27) | `etc` (ETC operational, `legacy/FFP_analysis_dev.R` lines 1010–1115) |
|---|---|---|
| roughness sublayer | zm > 20 z0 | zm > m·z0, m from 5 down in steps of 0.1 until rejected ≤ 0.5 × retained (min 2) |
| entrainment layer | zm < 0.8 h | zm < 0.8 h |
| tower below PBL | – | EC height < h |
| stability | zm/L ≥ −15.5 | −15.5 < zm/L ≤ 1 |
| friction velocity | u* > 0.1 m s⁻¹ | u* ≥ 0.2 m s⁻¹ |
| inputs | missing → invalid | missing → invalid |

**u* conditions (wording clarified in amendment 2, G 2026-11-10):** the `etc` set always contains u* ≥ 0.2 m s⁻¹,
as legacy `FFP.invalid.cond.flag` (`legacy/FFP_analysis_dev.R` line 1094: `USTAR < 0.2` → invalid). This condition
does not depend on the legacy switch `DO.ustar.filter`, which is not used here (the optional u* threshold of the ETC
footprint script, see the subsets below). Independently of the set, the spec-05 u* plausibility rule
(`flag_ustar_implausible` of PBLHmod) marks rows as invalid inputs in every set, `none` included, and the FFP code's
own u* ≥ 0.1 m s⁻¹ check applies to every set.

**Stability (G, 2026-10-26):** zm is the height above the displacement height, zm = EC height − d with d = 2/3 hc, so
zm/L = (z − d)/L (Kljun's definition, confirmed by N. Kljun to G). One stability test for all sets, computed from this
zm and L2 `MO_LENGTH`. L2 `ZL` is **not** used: its implied z − d differs from the ICOS heights (BE-Bra 17–18 m vs
24.1 m in 2020–2025; steps between periods at both sites; `docs/local/etc_qc_report.md` item 11). The `etc` set
therefore deviates from the ETC pipeline, which used `ZL`.

Plus `none` (only the FFP code's own checks). Check the Kljun thresholds against the FFP code and the paper (`Kljun2015`, Eq. 27
and Sect. 2; PDF in G's bibliography folder, added 2026-10-26), and report any difference. The ETC z0 multiplier
`m` is data-dependent: compute it **once per site**, over the full period and on rows independent of h (all rows with
zm, z0), and use the same `m` in every case; report it. The FFP code's own checks apply on top of every set.

**QC and u* subsets** (published data only; no cleanset):
- `valid`: valid footprint under the chosen set;
- `valid_qc`: `valid` and the NEE (or FC) QC of L2 FLUXES says good (flag 1/2 removed, as ETC);
- `valid_qc_ustar`: `valid_qc` and the NEE was not removed by the u* filtering of the FLUXNET processing.
Wind-sector exclusions are already applied in L2 (G); the ETC optional u* threshold of the footprint script is not
used (G). **Step 0 must find how to identify the u*-filtered rows from the published files** (e.g. FLUXNET
`NEE_VUT_*_QC` vs FLUXES QC, or published u* thresholds). If it cannot be done unambiguously, stop and ask G.

## 3. Footprint computations

FFP settings as ETC operational (`legacy/FFP_calc_parallel_dev.R`): domain ±1000 m, dx = dy = 1 m, `rslayer = 1`,
`smooth_data = 1`, r levels 50/70/80/90 %; report any setting you take from there. Same settings for all cases.

- **A. Per half-hour, crosswind-integrated** (cheap): peak distance x_peak and distances x_50, x_70, x_80, x_90 along
  the wind.
- **B. Per half-hour, 2D** (expensive): footprint fraction inside the target area; 80 % contour area.
- **C. Climatologies** (expensive): monthly and annual, per case × validity set × subset (`valid`, `valid_qc`,
  `valid_qc_ustar`), day/night separately (day from `SW_IN_POT` > 0): fraction inside the target area, 80 % contour
  area and maximum distance, and the overlap between cases (Jaccard index of the 80 % contours, PBLHmod vs ERA5;
  FI-Hyy also vs measured on common rows).

Computational design (the run may go to the cluster):
- one script per step, `Rscript scripts/08x_... SITE [--case=] [--set=] [--with-oper]`, config profile from
  `R_CONFIG_ACTIVE`; chunks by site × case × year, written to `data/derived/footprint/{site}/` (gitignored); a chunk
  already written with the same inputs hash is skipped (restartable); a log per run;
- parallel over chunks with `parallel` (forking on Linux, PSOCK on Windows), number of workers from an option;
- **Runtime first (step 0):** time A, B and one monthly climatology on 1,000 half-hours of BE-Bra, extrapolate per
  site and in total, and report it with memory use. G then decides PC or cluster.

## 4. Analysis (`scripts/08_footprint_compare.R`)

Paired comparisons on common rows (rows valid in both cases); coverage separately.

- **Coverage tables** (`08_coverage_{site}.csv`): per case × validity set: rows with h; valid; valid_qc;
  valid_qc_ustar; rejected per condition (rows can fail several), by regime (PBLHmod), season and day/night. Includes
  the rows where h < zm (e.g. the 2,811 PBLHmod values at FI-Hyy) and h ≤ 10 m.
- **Half-hourly differences** (`08_halfhour_metrics_{site}.csv`): ln(x_r, case / x_r, pblhmod) and the elasticity
  d ln x_r / d ln h (per row and median by stratum); median, IQR, share > 10 % and > 25 %, by regime, zm/L class,
  zm/h class, season, hour (midpoint bins, `hour_bin_centred()`), per subset. Target-area fraction: difference in
  percentage points.
- **Which half-hours matter** (G's point): rank the rows by |ln x_80 ratio| (PBLHmod vs ERA5); for the top 10 % and
  25 %: share in each regime and day/night, and share removed by QC and by the u* filtering (i.e. not in
  `valid_qc_ustar`). Same for rows whose validity changes between cases.
- **Climatologies** (`08_climatology_{site}.csv`): per month and year, per subset: target-area fraction, 80 % area and
  distance, Jaccard overlap; and whether ICOS target-area criteria (`Rebmann2018`; check the exact criteria in the
  paper and write them in the docs) pass with PBLHmod and fail with ERA5 or vice versa.
- **Anchor**: the ±20 % cases vs Kljun et al. 2015 Table 5 (peak location and value in stable / neutral / convective).
- **FI-Hyy measured**: on rows with exact MLH, PBLHmod and ERA5 vs `measured`. State the selection (exact values
  over-represent shallow layers). Optional sensitivity: qc 3 lower bounds used as h, labelled "bound", to show the
  direction of the footprint change.

Figures (PNG): `08_coverage_{site}` (valid rows per case and subset, by month); `08_dx80_by_zmh_{site}` (ln x_80
ratio vs zm/h, by regime); `08_dx80_diurnal_{site}` (median ln x_80 ratio by hour and season, `valid` vs
`valid_qc_ustar`); `08_clim_contours_{site}` (annual 80 % contours, PBLHmod vs ERA5, over the target area, per
subset); `08_topaffected_{site}` (top-affected rows by regime and subset). ERA5 tower cell only in figures.

## 4b. PC test run (G, 2026-10-26)

The full run (both sites, whole period) is for the cluster, together with the other stations. On G's PC, run a
**test on one month** to get a first feeling for the size of the effect; no definitive conclusions are drawn from it.

- All scripts take `--start=YYYY-MM-DD --end=YYYY-MM-DD` (also needed for the cluster chunks).
- Test month: **September 2022** at both sites (day and night of similar length, so both convective and stable
  regimes; inside the FI-Hyy Cloudnet period, so the `measured` case is included). In step 0 check the coverage of
  each case in that month (rows with h, valid rows); if a case has < 70 % of the rows of the month, propose the
  nearest month with better coverage and ask G.
- Computed for the test month only: A, B and C (C = the monthly climatology, day/night; no annual climatology).
- Computed for the **full period** anyway, because they need no FFP run: the validity flags, the ETC z0 multiplier
  and the coverage tables (section 4, first item; the FFP code's own checks reproduced as conditions, not by
  running the model; report if any check can't be reproduced that way).
- Outputs of the test carry the period in the file name (e.g. `08_halfhour_metrics_FI-Hyy_2022-09.csv`) and are
  labelled "test, one month" in figures and in STATUS.
- The step-0 runtime estimate is still required (per month and extrapolated to the full period and to all
  stations), to plan the cluster run.

## 5. Tests

- ERA5 interpolation to half-hour midpoints (hand case; missing neighbour → NA).
- Validity flags per condition for `kljun`, `etc`, `none` (hand cases at the boundaries; u* 0.1 vs 0.2; zm/L 1).
- ETC z0 multiplier: data-dependent loop reproduces the legacy rule; computed once and reused across cases.
- Subset logic (`valid` ⊇ `valid_qc` ⊇ `valid_qc_ustar`).
- Elasticity and log-ratio metrics on synthetic values.
- FFP wrapper on a small synthetic set: skipped when `ffp_code_dir` is missing (CI has no Kljun code), with a clear
  message.

## 6. Docs

- `CLAUDE.md`: FFP code location and licence rule (`ffp_code_dir`), validity sets, u* subset derivation.
- `docs/project_plan.md` (phase 5): 5a/5b split; the validity table above (resolves the open point "write the
  published and ETC relaxed thresholds in the phase-5 spec"); the coverage question (G).
- `docs/STATUS.md`: spec 08 row and results (without the `oper` case); remove "Relaxed FFP validity thresholds" from
  "Waiting on others".

## Acceptance criteria

- Step 0 report (FFP code API and checks, u* subset derivation, runtime estimate) before any long run.
- Tests pass; on the PC: coverage over the full period and A/B/C for the test month (section 4b) at BE-Bra and
  FI-Hyy; scripts ready for the full run on the cluster; outputs in `output/` and `data/derived/` only.
- Kljun code untouched and untracked (`git status`, `git grep -n "calc_footprint_FFP <-"` returns nothing tracked).
- Coverage, half-hourly and climatology tables for both sites; anchor vs Kljun Table 5.

## Out of scope

Monte Carlo (5b); other sites; footprint-weighted flux attribution or land-cover analysis; changes to PBLHmod or to
the FFP code; the ETC footprint u* threshold option.

## Summary to report back

As in `CLAUDE.md`, plus: step-0 findings; coverage per case and subset; median and top-decile differences in x_80 and
target-area fraction, PBLHmod vs ERA5, by regime, for `valid` vs `valid_qc_ustar`; share of the most-affected rows
removed by QC/u*; climatology overlap and target-area criteria; anchor vs Kljun Table 5; FI-Hyy measured case;
anything unexpected.

## Amendment 1 (G, 2026-10-26/27): decisions during implementation

1. **Subsets.** `valid_qc_ustar` is implemented as `valid_qc_used`: NEE QC good in FLUXES and not gap-filled in
   FLUXNET (`NEE_VUT_REF_QC == 0`). It removes the u*-filtered rows (checked against the ONEFlux VUT thresholds of
   the BIF file: 100 % of the u*-filtered rows are gap-filled) and the other FLUXNET NEE removals (u* explains 74 % of
   the exclusions at FI-Hyy, 81 % at BE-Bra). New subset `valid_qc_h`: valid and H measured and QC-good
   (`H_F_MDS_QC == 0`, or L2 `H` present where FLUXNET has no QC value); the u* filter applies to NEE only, so this
   subset keeps the stable nights retained for H, LE and other gases.
2. **Validity sets.** `kljun` = Eq. 27 only (no u* condition); new set `kljun_norsl` = Eq. 27 without the roughness
   condition; `etc` as legacy; `none`. The FFP code's own checks apply to every set. Stability from zm/L with
   zm = EC height − 2/3 hc for every set (section 2).
3. **FFP calls.** umean only (z0 = NaN) in every call. A (crosswind-integrated distances) from the paper's
   equations, C (climatologies) from the summed B fields; both verified against Kljun's code
   (`scripts/08v_ffp_verify.R`).
4. **Monthly sums.** Each monthly chunk keeps the unsmoothed summed fields and counts per group
   (`<site>_<case>_<YYYY-MM>_sums.rds`, float32), so annual climatologies are built by summing months and smoothing
   at the end (`ffp_clim_from_sums()`; tested against a direct climatology).
5. **Target area.** The target-area fraction is relative to the full footprint (f_domain = part inside the
   ±1000 m domain). Reported as a distribution (p10, median, p90) by zm/L class and as the difference between cases
   [percentage points]; how often the 90 % region reaches the domain edge (analogue of the code's `flag_err = 2`).
   The ICOS criteria of `Rebmann2018` are kept as a site-design check: ≥ 90 % of the flux from the area of interest,
   also under stable conditions, for ≥ 50 % of the time (p. 475); cumulative footprint ≥ 80 % for an acceptable
   measurement (p. 487).
6. **Production run.** One reference year per station on the cluster (selection rule in the cluster spec). The
   EBImage check of the smoothing is done there.

## Amendment 2 (G, 2026-11-10): per-site override of the ETC z0 multiplier (DE-Tha), DE-Tha reference year 2025

Branch `spec08-m-override`.

**Why.** The roughness test of the `etc` set (zm ≥ m·z0, d = 2/3 hc, z0 = 0.15 hc, zm = EC height − d) depends only on
hc and the EC height, so it passes or fails whole height periods. The multiplier loop (section 2: m from 5 down by
0.1 until rejected ≤ 0.5 × retained, minimum 2, once per site over all rows with zm and z0) cannot keep a period that
holds about a third of the rows or less. At DE-Tha (EC height 42 m) hc rises from 27.74 m (2020-10-26) through
30.60 m (2024-04-15) to 32.59 m (2025-03-05; linear in between, constant after), so zm/z0 falls from 5.65 to 4.15.
The loop stops at **m = 4.6** (38,932 of 116,880 rows rejected, 0.499 × retained): every half-hour from 2024-06-11
on fails, and 2025 has 0 usable half-hours. (The amendment text said m = 5; the cluster table of 2026-11-04 and a
height-only recomputation on the PC both give 4.6, with the same 38,932 rows. The decision does not change: 4.6
still rejects all of 2025.)

**Decisions (G, 2026-11-10).**
- The heights in `config/sites_timevarying.csv` are correct and stay (the DE-Tha BIF is being checked by a
  colleague).
- **m = 4.1 at DE-Tha**, a documented per-site override: `config/etc_z0_multiplier.csv` (`etc_z0_multiplier_file` in
  `config.yml`), columns `site_id, z0_multiplier, decided_by, date, reason`; read by `read_etc_z0_multiplier()`
  (every field required, `site_id` in `config/sites.csv`, 2 ≤ value ≤ 5, one row per site). With 4.1 no DE-Tha row
  fails (minimum zm/z0 4.148); no row of 2020–2023 fails at 4.6 or at 4.1 (height-only check on the PC).
- The loop stays the default everywhere else, and is always computed: `etc_z0_multiplier_site()` returns `m_loop`,
  `m_used` and `m_source` (`loop` / `override`); 08a writes them to `08_etc_multiplier_<site>.csv` (also
  `n_rejected` with m_used, `n_rejected_loop`, `zm_over_z0_min`, `override_reason`) and prints them in its log;
  `09_collect.R` writes `09_network_etc_multiplier.csv` (tables written before the amendment: m_loop = m_used = m_etc,
  source loop). The proper fix, m per height period, belongs to the ETC pipeline (G's item 15 of
  `docs/local/etc_qc_report.md`) and is not done here.
- No change to the other validity conditions (u* wording clarified in section 2).
- **DE-Tha reference year 2025** (was 2023) in `config/reference_years.csv`, approved G 2026-11-10, so that the
  footprint year overlaps the ceilometer data requested from the site (TU Dresden SkyVUE raw files, later). The
  cluster checks the 9-month usable rule after the 08a re-run.

**Inputs hash.** `reference_years.csv` was never hashed as a file: only step 08 depends on it, through its argument
`--year=<the site's year>` (`pipeline_steps()`), which is part of its hash; so the DE-Tha change reruns step 08 at
DE-Tha only. Nothing changed there. 08a now hashes the site's own row of `etc_z0_multiplier.csv`
(`config_rows`, `site_config_row()`), added to the hash only when the site has a row: the hash of a site without a
row is computed exactly as before, and an edit of the DE-Tha row changes no other site's hash (tested). The chunk
hash of 08b / `08c --pending` (`footprint_task_hash()`) depends on the rows and groups of each chunk, not on these
files. As for every code change, the new code in `R/` changes every step's hash once (any step run again is
recomputed, with identical results where nothing else changed); finished footprint chunks of the other sites stay
valid (their task hashes do not change).

**Roughness condition of the `kljun` set.** zm > 20 z0 (Kljun et al. 2015, Eq. 27) fails at essentially every forest
site (zm/z0 well below 20 with d = 2/3 hc and z0 = 0.15 hc; medians BE-Bra 12.2, FI-Hyy 3.8, DE-Tha 4.9, range
4.1–5.7 at DE-Tha; six of the 39 run sites have a loop multiplier below 5): this is why
`kljun_norsl` (Eq. 27 without the roughness condition) exists. The ETC relaxation (m ≤ 5) is the same compromise.

**Collection.** `09_collect.R` reads only the per-site 08a tables (`08_coverage_<site>.csv`, full period,
`08_year_counts_`, `08_nee_qc_`, `08_etc_multiplier_`), never the step-08 outputs of a year
(`08_*_<site>_<year>.csv`, figures), so the DE-Tha 2023 outputs of the earlier step-08 run are not picked up. The
reference-year proposal `08_reference_years.csv` now shows the approved year next to it (`approved_year`), and
`09_network_etc_multiplier.csv` carries `reference_year` (approved).

**Methods note (paper).** The ETC roughness test removes whole height periods, because zm/z0 depends only on the
canopy and EC heights; the data-dependent multiplier is reported per site (loop value and value used), and the one
per-site override (DE-Tha, m = 4.1) is listed with its reason.
