# Project status (handover note)

Last update: 2026-10-26. Read this first in a new session, then `docs/project_plan.md` and `CLAUDE.md`.

## PBLH sources (names used in the code and in the standard files)

| Source | `source_id` | What it is | Class |
|---|---|---|---|
| PBLHmod v2 | `pblhmod_default`, `pblhmod_gamma0025` | Our EC-based model (Kljun et al. 2015 App. B): stable B1 Nieuwstadt, neutral B2/B3 auto, convective B5 Batchvarova–Gryning slab daily cycle. `run_pblh_model()`, `R/pblh_*.R`; label = options case (γ = 0.01 default, 0.0025 sensitivity) | parameterisation |
| ERA5 BLH | `era5_blh_center`, `era5_blh_bilinear` | ECMWF reanalysis diagnostic (bulk Ri = 0.25, IFS Cy41r2 §3.10.1); tower cell, or bilinear to the tower from 4 cells | model |
| Operational ICOS PBLH | `icos_oper` | `PBLH` column of L2 FLUXES, produced by the ETC pipeline with `legacy/calc_PBL_height_operational.R` + `legacy/proc_PBL_height.R`; `qc_flag` 1 = 50 m floor, 2 = run-mean fill. Used in the phase-1 regression only (spec 03), not in the evaluation (project plan "Scope", G 2026-10-23); script 04 compares it only with `--with-oper` | parameterisation |
| Legacy (spec 03) | `legacy_<variant>` | Port of the operational legacy code with switches, for the regression (in `data/derived`, not a standard file) | parameterisation |
| Measured MLH, Cloudnet (spec 07) | `cloudnet_eps_mlh` | Surface-connected mixing from the Doppler-lidar dissipation rate ε (ε ≥ 1e-4 m² s⁻³ upward from the lowest reliable gate; Manninen et al. 2018); half-hourly; `qc_flag` 0/2 exact, 1 censored below z_min (105 m at FI-Hyy), 3 signal-limited = lower bound; precipitation = Cloudnet `rain_detected` | turbulence |
| Other measured PBLH (later) | – | PI data (SE-Htm, SE-Nor, SE-Svb, FI-Hyy, CZ-BK1 requested), ceilometer / radiosonde (E-PROFILE, IGRA) | aerosol / thermodynamic / turbulence |

## Specs

| Spec | Content | Phase | State |
|---|---|---|---|
| 01 | Site metadata, time-varying heights, exclusions, ICOS L2 readers, unit checks | 0 | merged (`main`) |
| 02 | PBLHmod v2 | 1 | merged (`main`, 58d45ad) |
| 04 | ERA5 readers (3×3 cells, pressure levels, static), standard PBLH format, first PBLHmod vs ERA5 comparison at BE-Bra (amendments 1–8) | 0, 2, 3 (first step) | merged (`main`, aba3394) |
| 05 | u* plausibility filter in the PBLHmod forcing (spec 02 amendment 12); local config override (`config.local.yml`, no cluster paths in tracked files) | 1 | merged (`main`, 25ae55c) |
| 06 | Consistency check of L2 `MO_LENGTH` against the forcing (option C, G 2026-10-21): L2 L kept where consistent (same sign, within 3×), otherwise replaced by L from the forcing and flagged (spec 02 amendment 13); checked only where the heat flux is positive (spec 06 amendment 1) and in daytime (amendment 2) | 1 | merged (`main`, 2026-10-21; amendment 2 merged 2026-10-21) |
| 03 | Legacy regression: reproduce operational PBLH from L2 (operational code `calc_PBL_height_operational.R` + `proc_PBL_height.R`), port it, attribute differences to PBLHmod v2 (amendments 1–2) | 1 | merged (`main`, e224c4b) |
| 07 | Cloudnet inventory (FI-Hyy, FI-Lom, DE-RuS, SE-Nor, FR-Gri); turbulence-based MLH from Doppler-lidar ε (Manninen et al. 2018) in the standard format; FI-Hyy runs (PBLHmod, ERA5) and first comparison against measured MLH (amendment 1: no operational PBLH, ERA5 tower cell in figures); reference period = FLUXES range; METEO TA levels from the header | 2, 3 | merged (`main`, b22f7e2) |
| 08 | Footprint sensitivity to the PBLH source (FFP, Kljun et al. 2015) at BE-Bra and FI-Hyy; one-month PC test (FI-Hyy 2022-09); ONEFlux BIF reader and u* thresholds | 5a | in progress on `spec-08-footprint-sensitivity` |

Phases (see `project_plan.md`): 0 done; 1 done (model, input QC, legacy regression); 2 standard format, three
sources (spec 04) and the first measured source (spec 07, FI-Hyy) done; 3 BE-Bra (vs ERA5) and FI-Hyy (vs ERA5 and
measured MLH) done; 4–6 not started.

- Tests: 968 passing (`spec-07-measured-cloudnet`). Environment: R 4.4.1 + renv (Posit Package
  Manager), same as the cluster. Machine-specific paths in the untracked `config.local.yml` (template
  `config.local.yml.example`).
- PBLHmod defaults: `stable_method = "nieuwstadt"` (B1), `neutral_method = "auto"` (B2 if 10 < N/|f| < 70, else
  B3), `combine = "bg_only"`, `gamma = 0.01` K/m, `h_source = "H_F_MDS"` with `h_qc_max = 1`; u* input filter
  `ustar_ws_max = 0.5` (where WS ≥ `ws_min_ratio_test = 1` m s-1), `ustar_max = 3` m s-1 (spec 02 amendment 12;
  u* and L set to NA, short gaps interpolated); L check `l_source = "L2_checked"` (`l_consistency_factor = 3`,
  daytime rows with wtheta0 > `l_check_wtheta_min = 0.005` K m s-1 only, `l_check_daytime_only = TRUE`;
  inconsistent L2 L replaced by the forcing L; spec 02 amendment 13, spec 06 amendments 1–2); raw output (no smoothing, clipping or filling). Forcing T and p from
  measured meteo (BE-Bra: FLUXNET `TA_F`/`PA_F` 2020-01 → 2026-08-22, METEO `TA`/`PA` after; before 2026-10-26
  FLUXNET 2020–2025 and METEO for all of 2026).
- BE-Bra PBLHmod (default, FLUXNET 2020–2026, 2026-10-26): convective median 707 m (p5–p95 75–1581), stable 257 m,
  all rows 414 m; γ = 0.0025 raises the convective median to 1122 m. With FLUXNET 2020–2025 (spec 06 amendment 2):
  690 m (74–1575), 259 m, 412 m; 1095 m. Spec 06 amendment 1: 690 m (74–1575), 260 m, 413 m;
  1095 m. Spec 05: 698 m (72–1586), 261 m, 413 m; 1104 m. Spec 02: 704 m (72–1609), 271 m, 423 m; 1109 m.
- BE-Bra FLUXNET 2020–2026 replacement (G, 2026-10-26; old outputs in `*/pre_bebra_fluxnet2026/`):
  - Inputs, old vs new file on 2020–2025 (105,216 common rows): only gap-filled values differ. `TA_F` 14,318 rows (all
    QC 2, 99 % in 2020; max 0.033 K), `PA_F` 15,530 rows (all QC 2; max 28 Pa), `H_F_MDS` 69 rows (all gap-filled,
    2025-12-26 → 12-31; max 63 W m-2), `NEE_VUT_REF_QC` 754 rows (0.8 %); NA pattern unchanged; `H_CORR` differs on
    48 % (not used). 2026: FLUXNET `TA_F`/`PA_F` identical to METEO where both are valid; FLUXNET ends 2026-08-22
    22:30 UTC.
  - PBLHmod (13 cases, scripts 01, 02, 04, 03 rerun): before 2025-12-26 max |diff| ≤ 0.16 m (gap-filled TA/PA);
    2025-12-26 → 12-31 up to 534 m (default; 16 rows > 1 m, all on the days with changed `H_F_MDS`); 2026 changed
    because H now comes from `H_F_MDS` instead of L2 `H` (default: valid H rows 6,318 → 10,553, convective 1,860 →
    4,337, B&G stopped at gaps 8,959 → 2,148, convective median 594 → 783 m). `l_forcing` gains 3,387 valid 2026 rows.
  - Pair 1 vs ERA5 (49,400 h): RMSE 337.5 → 333.5 m, bias −49 → −42 m, median log ratio −0.153 → −0.146, within a
    factor 2 72.6 → 72.8 %, Spearman 0.797 → 0.801; convective n 15,242 → 16,442.
  - Spec 03: the reproduction of the operational PBLH is unchanged (all and per-year strata identical; regime strata
    change because regimes come from PBLHmod); variants with our forcing (`legacy_forcing`, `cum_06`–`cum_12`) change
    as expected. The `cleanset` variant is no longer produced (cleanset moved to `data/archive/`; old result kept in
    the backup).
  - ONEFlux u* check now runs at BE-Bra: VUT thresholds 0.18–0.28 m s-1 (CUT 0.2345); 100 % of u*-filtered rows
    gap-filled; u* explains 80.7 % of the QC-good NEE gap-filled in FLUXNET (day 71.7 %, night 84.2 %; FI-Hyy 74 %).
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
- Spec 07 at FI-Hyy (2026-10-23):
  - Inputs: FLUXES/FLUXNET/METEO 2018-01 → 2026-08 (data from 2018-11-06); unit checks pass; METEO TA_1–TA_9.
    PBLHmod default: 80 % valid, convective median 700 m (stable 315 m), u* rejected 1.8 %, L replaced 2,013, none
    > 5 km; 2,811 values below zm = 27 m (stable nights 2,282, early-morning B&G starts 430). ERA5 tower cell 11.1 km
    from the tower, 42 m lower. PBLHmod vs ERA5 (2020-01 → 2026-08, 51,584 h): median log ratio −0.10 (DJF +0.16,
    JJA −0.35), 73 % within a factor 2, Spearman 0.75 (BE-Bra −0.14 / 72 % / 0.78). Operational fills
    3316.87/2000.85 m (2018-11-06 → 2025-12-31) and 3365.22/1532.14 m (2026-01-01 → 2026-08-26), 28 % of rows.
  - Measured MLH (Cloudnet `epsilon-lidar`, 688 days 2022-01 → 2024-04, 390 MB; classification 682 days): z_min
    105 m agl; profiles: measured 67 %, censored ε < 1e-4 at 105 m 22 %, 105 m gate missing 6 %, no signal 1 %,
    precipitation (`rain_detected`) 3.5 %; half of the measured profiles end at a gate without signal (clean air).
    Half-hours (30,727): qc 0 3,171, qc 2 6,870, qc 3 (lower bound) 11,889, qc 1 (censored) 7,365, NA 1,432. JJA
    daytime: qc 3 on 54–62 % of half-hours, so the censored median is not estimable there.
  - Comparison on exact values (qc 0/2), convective regime: PBLHmod n 2,703, bias +202 m, median log ratio +0.44,
    56 % within a factor 2, Spearman 0.47; ERA5 tower cell n 880, +111 m, +0.32, 69 %, 0.53 (bilinear practically
    identical). Exact values over-represent shallow layers. Where the lidar MLH is below 105 m, PBLHmod is ≤ 105 m on
    26 % (ERA5 24 %); where it is only a lower bound, PBLHmod is below the bound on 40 % (JJA 62 %; ERA5 33 %, JJA
    45 %): "shallow too deep / deep too shallow", candidate finding. Outputs with the operational PBLH kept in
    `output/*/reference_spec07_with_oper/` (and `output/tables/reference_spec04_with_oper/`).
  - Sensible-heat timing (`scripts/07d_h_timing.R`, 2020–2025, values at interval midpoints): ERA5 H earlier than
    EC H, best-fit lag FI-Hyy 25–37 min (JJA 37), BE-Bra 16–25 min; FI-Hyy JJA afternoon crossings ≈ 63 min earlier,
    morning ≈ 15–18 min (asymmetric). Adding the L2 storage term `SH` (air below the EC system, ±12 W m-2 at
    FI-Hyy) reduces the lag (FI-Hyy 12–27 min, BE-Bra −1 to 14 min) but leaves the FI-Hyy afternoon difference at
    51–64 min. Diagnostic only (open point in `project_plan.md`).
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
5. Spec 07 — implemented on `spec-07-measured-cloudnet` (2026-10-23), to review and merge. Scripts: 07a inventory,
   07b download (size report first, `--confirm`; `--match=PRODUCT`), 07c MLH (profiles in `data/derived/cloudnet/`,
   standard file, sensitivity columns), 07 comparison, 07d sensible-heat timing; 04 with optional `--with-oper`.
   Diurnal figures place every value at its interval midpoint (hour bins centred on the full hour). Distance rule and censoring
   convention in `project_plan.md`; Cloudnet facts in `CLAUDE.md`; example in `docs/specs/standard_format.md`.
6. After spec 07: the other Cloudnet sites (FI-Lom, DE-RuS, FR-Gri; on the cluster), PI data into the standard
   format, ETC QC report. Before paper results: re-check the volatile Cloudnet files and create a Cloudnet collection
   (open point in `project_plan.md`).

## Data

ERA5 from C. Trotta (Oct 2026), 0.25°, tower cell + 8 neighbours, all stations (details in `CLAUDE.md`):
hourly BLH, cloud base height, sensible heat flux (2020-01-01 → 2026-10-07); 6-hourly temperature and geopotential
on 13 pressure levels (enough for a γ climatology, not for time-varying γ); static orography, land–sea mask, sdor.
SE-Myc uses the corrected coordinates; all stations structured like BE-Bra. BE-Bra tower cell: 6.55 km from the
tower, orography 8.9 m, lsm 0.96, sdor 7.6 m.

Literature: `docs/references/literature_review_pblh.md`, `pblh_literature_2026-10.bib` (~50 verified papers),
`pblh_collection.bib`. Key papers and ECMWF IFS Cy41r2 Part IV downloaded; second-wave papers checked against the PDFs (2026-10-22; Durre2008 from Crossref only, no PDF).
Novelty check: no published network-scale evaluation of EC-derived PBLH against ERA5 and measurements.

FI-Hyy (spec 07): ICOS L2 (FLUXES, FLUXNET 2018–2026, METEO) and ERA5 (three families) under `data/raw/`; Cloudnet
Hyytiälä `epsilon-lidar` and `classification` under `data/raw/cloudnet/FI-Hyy/` (with provenance tables).

L2 data issues found so far (BE-Bra: implausible USTAR, MO_LENGTH kept where H is QC-rejected, positive night-time
H_F_MDS in winter 2022–23, sonic temperature offset, latent pipeline bug; FI-Hyy: USTAR up to 189 m s-1 and WS
median 122 m s-1 on cap rows, T_SONIC 6–220 K in winter, L2 PBLH up to 309 km and 28 % run-mean fills) are
collected for the ETC QC report in `docs/local/etc_qc_report.md` (local working file, not tracked); all open.

## Waiting on others

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
- Reference PDFs live outside the repo in `C:\WORK\bibliography_new_win11\` (not in git; the `file` fields were
  removed from the tracked `.bib` files on 2026-10-26).
- BE-Bra: the replaced FLUXNET 2020–2025 file, the cleanset and the EddyPro outputs are in `data/archive/` (not read
  by the pipeline).
