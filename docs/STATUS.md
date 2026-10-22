# Project status (handover note)

Last update: 2026-10-21. Read this first in a new session, then `docs/project_plan.md` and `CLAUDE.md`.

## PBLH sources (names used in the code and in the standard files)

| Source | `source_id` | What it is | Class |
|---|---|---|---|
| PBLHmod v2 | `pblhmod_default`, `pblhmod_gamma0025` | Our EC-based model (Kljun et al. 2015 App. B): stable B1 Nieuwstadt, neutral B2/B3 auto, convective B5 Batchvarova–Gryning slab daily cycle. `run_pblh_model()`, `R/pblh_*.R`; label = options case (γ = 0.01 default, 0.0025 sensitivity) | parameterisation |
| ERA5 BLH | `era5_blh_center`, `era5_blh_bilinear` | ECMWF reanalysis diagnostic (bulk Ri = 0.25, IFS Cy41r2 §3.10.1); tower cell, or bilinear to the tower from 4 cells | model |
| Operational ICOS PBLH | `icos_oper` | `PBLH` column of L2 FLUXES, produced by the ETC pipeline with `legacy/calc_PBL_height_operational.R` + `legacy/proc_PBL_height.R`; `qc_flag` 1 = 50 m floor, 2 = run-mean fill | parameterisation |
| Legacy (spec 03, to come) | `legacy_<variant>` | Port of the legacy code with switches, for the regression | parameterisation |
| Measured PBLH (later) | – | Ceilometer / lidar / radiosonde (ACTRIS/Cloudnet, E-PROFILE, IGRA) | turbulence / aerosol / thermodynamic |

## Specs

| Spec | Content | Phase | State |
|---|---|---|---|
| 01 | Site metadata, time-varying heights, exclusions, ICOS L2 readers, unit checks | 0 | merged (`main`) |
| 02 | PBLHmod v2 | 1 | merged (`main`, 58d45ad) |
| 04 | ERA5 readers (3×3 cells, pressure levels, static), standard PBLH format, first PBLHmod vs ERA5 comparison at BE-Bra (amendments 1–8) | 0, 2, 3 (first step) | merged (`main`, aba3394) |
| 05 | u* plausibility filter in the PBLHmod forcing (spec 02 amendment 12); local config override (`config.local.yml`, no cluster paths in tracked files) | 1 | merged (`main`, 25ae55c) |
| 06 | Consistency check of L2 `MO_LENGTH` against the forcing (option C, G 2026-10-21): L2 L kept where consistent (same sign, within 3×), otherwise replaced by L from the forcing and flagged (spec 02 amendment 13); checked only where the heat flux is positive (spec 06 amendment 1) and in daytime (amendment 2) | 1 | merged (`main`, 2026-10-21; amendment 2 merged 2026-10-21) |
| 03 | Legacy regression: reproduce operational PBLH from L2 (operational code `calc_PBL_height_operational.R` + `proc_PBL_height.R`), port it, attribute differences to PBLHmod v2 (amendments 1–2) | 1 | merged (`main`, e224c4b) |

Phases (see `project_plan.md`): 0 done; 1 model done, legacy regression pending (spec 03); 2 standard format and
three sources done (spec 04), measured PBLH pending; 3 started at BE-Bra only; 4–6 not started.

- Tests: 787 passing (`main`, after spec 03). Environment: R 4.4.1 + renv (Posit Package
  Manager), same as the cluster. Machine-specific paths in the untracked `config.local.yml` (template
  `config.local.yml.example`).
- PBLHmod defaults: `stable_method = "nieuwstadt"` (B1), `neutral_method = "auto"` (B2 if 10 < N/|f| < 70, else
  B3), `combine = "bg_only"`, `gamma = 0.01` K/m, `h_source = "H_F_MDS"` with `h_qc_max = 1`; u* input filter
  `ustar_ws_max = 0.5` (where WS ≥ `ws_min_ratio_test = 1` m s-1), `ustar_max = 3` m s-1 (spec 02 amendment 12;
  u* and L set to NA, short gaps interpolated); L check `l_source = "L2_checked"` (`l_consistency_factor = 3`,
  daytime rows with wtheta0 > `l_check_wtheta_min = 0.005` K m s-1 only, `l_check_daytime_only = TRUE`;
  inconsistent L2 L replaced by the forcing L; spec 02 amendment 13, spec 06 amendments 1–2); raw output (no smoothing, clipping or filling). Forcing T and p from
  measured meteo (FLUXNET `TA_F`/`PA_F` 2020–2025, METEO `TA`/`PA` 2026).
- BE-Bra PBLHmod (default, spec 06 amendment 2): convective median 690 m (p5–p95 74–1575), stable 259 m, all rows
  412 m; γ = 0.0025 raises the convective median to 1095 m. Spec 06 amendment 1: 690 m (74–1575), 260 m, 413 m;
  1095 m. Spec 05: 698 m (72–1586), 261 m, 413 m; 1104 m. Spec 02: 704 m (72–1609), 271 m, 423 m; 1109 m.
- Spec 06 amendment 2 at BE-Bra (2026-10-21): L check in daytime only. Rows checked 26,198; replaced 2,191 (sign
  1,170, magnitude 1,021; convective 1,720, unstable_diag 471), night replacements 285 → 0 (all in Nov 2022–Jan
  2023, positive gap-filled H_F_MDS with L2 H QC-rejected); `flag_sign_inconsistent` 335; values > 5 km 1.
  Regression chain: `nofilter_l2` = spec 02, `l_L2` = spec 05, `l_check_allday` = amendment 1 (all < 1e-9 m).
  Pair 1 vs ERA5 practically unchanged vs amendment 1 (RMSE 337.6 → 337.5 m; unstable_diag RMSE 449 → 454 m).
- Spec 06 at BE-Bra (2026-10-21, amendment 1 = V1): 26,652 rows checked, 2,476 replaced (sign 1,349, magnitude
  1,127; 69 % convective, 31 % unstable_diag; 88 % daytime; on 956 days; winter clusters, e.g. 2023-01-19 → 01-28); L2 L used
  on 96.5 % of rows with L; `flag_sign_inconsistent` 1,228 → 335. PBLHmod > 5 km: 11 → 1 (2024-01-05 09:30, B1,
  u* 2.67 m s-1: spec 05 filter limitation). Regression chain: `nofilter_l2` = spec 02 and `l_L2` = spec 05
  (max |diff| ≤ 4e-12 m). Pair 1 vs ERA5, spec 05 → 06: RMSE 349 → 338 m, Pearson 0.725 → 0.737, within a
  factor 2 72.4 → 72.6 %, convective RMSE 385 → 353 m, stable unchanged (310 → 309 m).
- Spec 05 at BE-Bra (2026-10-21): the filter rejects 2,149 of 100,186 u* values (2.1 %: ratio 1,328, cap 46,
  both 775); 618 then interpolated, 1,531 left NA. PBLHmod > 5 km: 179 → 11; > 10 km: 70 → 0. `nofilter`
  reproduces the spec-02 default (max |diff| 3.6e-12 m). Pair 1 vs ERA5: RMSE 623 → 349 m, Pearson 0.49 → 0.73,
  bias −9 → −46 m, median log ratio −0.136 → −0.151, within a factor 2 71.8 → 72.4 %.
- Spec 04 at BE-Bra (2026-10-17): overlap 2020-10-25 → 2026-08-25; PBLHmod default vs ERA5 tower cell 50,269 matched
  hours, bias −9 m, median log ratio −0.14, 72 % within a factor 2, Spearman 0.78; stable −0.23, convective −0.07;
  JJA evenings (19–22 UTC) PBLHmod ≈ half of ERA5. Operational PBLH vs ERA5: median log ratio +0.57 (qc 0).
  Open points from spec 04 in `project_plan.md` ("Open points").

## Next steps

1. ~~Review and merge spec 04 (`spec-04-era5-compare`).~~ Done (2026-10-21, merged at aba3394; `mad_log_ratio`
   switched to the scaled MAD at review, amendment 8).
2. ~~Review and merge spec 05 (`spec-05-ustar-filter`).~~ Done (2026-10-21, merged at 25ae55c).
3. ~~Spec 06 — consistency of L2 `MO_LENGTH` with the forcing.~~ Done (2026-10-21). Why option C (G): at BE-Bra,
   10 of the 11 values > 5 km left after the u* filter came from B5 on days where L2 H is QC-rejected but
   `MO_LENGTH` (from H_UNCLEANED ≈ 0, |L| ≈ 3e6 m) is kept. A (L NA where L2 H is NA) would stop B5 early on 959 of
   1,770 convective days, because L2 H is QC-rejected on 47 % of rows with L; B (inconsistent L set to NA) on 268
   days; D (L always from the forcing) is kept as a sensitivity case (`l_forcing`). Checking stable rows too made
   results worse (spec 06 amendment 1).
4. ~~Review and merge spec 03 (`spec-03-legacy-regression`).~~ Done (2026-10-22, merged at e224c4b). BE-Bra results:
   - Part A: the operational code on inputs rebuilt from L2 (EddyPro engine formulas for AT, AP, rho, cp) reproduces
     the operational column: median |rel. diff| 0.009 %, p95 0.06 %, fill values within 0.006 %, 1,341 of 1,342 fill
     blocks; the local cleanset gives 0.0004 % (remaining differences = rounding of the published inputs). Runs:
     2020-10-26 → 2026-01-01 and 2026-01-01 → 2026-08-26 (local standard time). The development version is out of
     scope (G, 2026-10-22).
   - Part B: port = original (max |diff| 0 m on the full series, < 1e-8 m on the synthetic test).
   - Part C, on common rows (valid in every variant and in PBLHmod v2; spec 03 amendment 2 note 9). One at a
     time (n = 79,263; baseline median 707 m, PBLHmod v2 364 m): neutral formula −211 m (convective −740 m: the
     operational unstable bound −0.003 puts many daytime rows in the neutral class, 0.3 u*/|f| ≈ 1.3 km), smoothing
     −192 m, gap interpolation +176 m, stable formula −38 m, gamma −11 m (convective −101 m); fill, floor, mol_locf,
     forcing ≈ 0. Cumulative (n = 75,689; 724 → 366 m): smoothing −206 m, neutral formula −82 m (convective −970 m),
     stable formula −46 m; the last step equals pblhmod_default (residual 4.6e-13 m). On all valid rows the fill
     (−383 m) and forcing (±160–190 m) effects are row-set changes, not value changes. `at_meteo` (sonic temperature
     offset): negligible.

## Data

ERA5 from C. Trotta (Oct 2026), 0.25°, tower cell + 8 neighbours, all stations (details in `CLAUDE.md`):
hourly BLH, cloud base height, sensible heat flux (2020-01-01 → 2026-10-07); 6-hourly temperature and geopotential
on 13 pressure levels (enough for a γ climatology, not for time-varying γ); static orography, land–sea mask, sdor.
SE-Myc uses the corrected coordinates; all stations structured like BE-Bra. BE-Bra tower cell: 6.55 km from the
tower, orography 8.9 m, lsm 0.96, sdor 7.6 m.

Literature: `docs/references/literature_review_pblh.md`, `pblh_literature_2026-10.bib` (~50 verified papers),
`pblh_collection.bib`. Key papers and ECMWF IFS Cy41r2 Part IV downloaded; second-wave papers: G will notify.
Novelty check: no published network-scale evaluation of EC-derived PBLH against ERA5 and measurements.

L2 data issues found so far (implausible USTAR, MO_LENGTH kept where H is QC-rejected, positive night-time H_F_MDS
in winter 2022–23, sonic temperature offset) are collected for the ETC QC report in `docs/local/etc_qc_report.md`
(local working file, not tracked); all open.

## Waiting on others

- Second-wave papers (G).
- Relaxed FFP validity thresholds used at ETC (for the phase-5 spec).
- Verification of forest canopy-height jumps (FR-Fon, FR-Hes, SE-Svb).

## Open scientific points (details in `docs/project_plan.md`)

- B1 on windy stable nights tends to the B2 limit (p95 ≈ 1.5 km at BE-Bra); candidate: bound by B3. Decide after
  comparison with ERA5/measurements.
- Polar-day convective rule (GL-ZaF, SE-Sto, FI-Lom).
- γ options 2–3 (ERA5 climatology / time-varying) after phase 1.

## Local setup notes

- PC rebuild procedure: `docs/setup/pc_rebuild_checklist.md` (+ `backup_before_format.ps1`).
- Gitignored but needed: `data/raw/` (ICOS L2 zips for BE-Bra, ERA5 CSVs in `era5/`, `era5_pressure_levels/`,
  `era5_single_levels/`, `target_areas/*.kml`), `C:\WORK\secrets`, `docs/local/` (local working notes).
- Reference PDFs live outside the repo in `C:\WORK\bibliography_new_win11\` (not in git; the `.bib` `file` fields
  point there).
