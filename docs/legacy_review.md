# Legacy code review (`legacy/`)

Reviewed 2026-10-08. Reference only: new code goes in `R/` and `scripts/`.

## What's there

| File | Purpose | Reuse |
|---|---|---|
| `calc_PBL_height_operational.R` | **The version the ETC pipeline sources** (`FFP_analysis_dev.R` line 133) to produce the operational L2 `PBLH`. PBLH from EC data (L, u*, H, T) following Kljun et al. 2015 App. B: diagnostic formulas for neutral/stable, Batchvarova & Gryning (1991) slab model for convective conditions; stability classes from 1/L (Essa et al. 2006, with a modified bound). | Reference of spec 03 (reproduction, port, attribution). |
| `calc_PBL_height_developing.R` | G's later development version (formerly `calc_PBL_height.R`; the one reviewed below and used in spec 03 until 2026-10-22). Differences from `_operational`: table below. | **Not used in this project** (G, 2026-10-22: out of scope in every phase and spec; only the operational code counts). Kept as a read-only record. |
| `proc_PBL_height.R` | ETC operational post-processing of `calc_PBL_height()` output: negatives → NA, 50 m floor, NA → run day/night mean. Identical to the pipeline copy. | Reference for spec 03 (explains the operational `PBLH` column). |
| `FFP_analysis_dev.R` | ETC pipeline script: reads the cleanset, calls `calc_PBL_height()` (line 984) and `proc_PBL_height()` (line 1003), writes the L2 column as `PBLH = round(PBL, 2)` (lines 1376, 1442). | Reference for the operational call (spec 03 section 1). |
| `FFP_calc_parallel_dev.R` | ETC pipeline: parallel FFP computation. Line 155: `h = ifelse(is.na(PBL), 1000, PBL)`, i.e. a missing PBL is replaced by **1000 m** in the crosswind-integrated footprint. | Note for phase 5 (footprints). |
| `calc_footprint_FFP.R`, `calc_footprint_FFP_climatology.R` | Kljun's FFP R code (see below). **Kept locally only, not in git** (`.gitignore`). | Phase 5 (footprints). |
| `PBL_validation.R` | Validation at FI-Hyy (SMEAR II), 2016: Doppler lidar MLH (.mat), ceilometer cloud base, SMEAR fluxes and meteo. Diurnal cycles (year/season/month), Bland–Altman, error boxplots. | Reuse the analysis design. Rewrite the I/O. |
| `PBL_ERA_comparison_test.R` | BE-Bra 2024: PBLH from ETC L2 footprint output vs ERA5 `boundary_layer_height`. Unfinished, does not run. | Reuse the idea (ERA5 as a second reference). |
| `FFP_sens_to_PBL.R` | Sensitivity of the FFP (Kljun) 80% fetch to PBLH. Depends on external `calc_footprint_FFP*.R`. | Possible later module. |
| `ERA-PBL-stations.txt` | 85 stations (ICOS + NEON) with lat/lon for ERA5 extraction. | Keep as site metadata. |
| `STabilityClasses_on_MOL.*` | Alternative stability classes from L (a–h). | Reference. |

**Exception to the "no cluster paths in tracked files" rule (spec 05, agreed 2026-10-21):**
`PBL_ERA_comparison_test.R` contains two absolute cluster paths in comments (lines 7 and 13). `legacy/` is
read-only (CLAUDE.md rule 5), so the file is kept unchanged; it is an accepted hit of the cluster-path and
user-name checks of spec 05.

**Kljun's FFP code (local only).** Kljun, N., Calanca, P., Rotach, M. W., Schmid, H. P. (2015): A simple
two-dimensional parameterisation for Flux Footprint Prediction (FFP), Geosci. Model Dev., 8, 3695–3713,
doi:10.5194/gmd-8-3695-2015. R code `calc_footprint_FFP.R` v1.3 (08/12/2017) and
`calc_footprint_FFP_climatology.R` v1.31 (25/03/2018), © N. Kljun, from the FFP website (footprint.kljun.net).
Kept locally only, not redistributed (copyright, no licence stated). The local climatology file is dated 2022-04:
possible local modification by ETC, to be checked against the published v1.31 before phase 5.

**Legacy-only R package (spec 03).** `zoo` is in `renv.lock` only because the operational legacy code `require()`s
it; it is used only by `scripts/03_legacy_regression.R` and `tests/testthat/test-legacy_port.R`, never by `R/`.
`bigleaf` and `data.table` (needed only by the development version) were removed from renv on 2026-10-22; `crayon`
stays as a dependency of `readr`, `testthat` and others. `legacy/` is excluded from renv's dependency scan
(`.renvignore`).

**Latent bug in `calc_PBL_height_operational.R`** (pipeline code; listed in the local ETC QC report): the loop
attempts a height where L **or** u* is present (`!is.na(1/MOL) | !is.na(1/ustar)`); where L is missing and u* present
the class test `if (1/MOL[i] < -0.056)` evaluates `NA` and stops the function. In the pipeline this is masked
because L is filled by LOCF before the call, so only leading rows of a run could trigger it (they never do when the
run starts at its first data row). The spec 03 port returns `NA` there instead.

## `calc_PBL_height_operational.R` vs `calc_PBL_height_developing.R` (added 2026-10-22, spec 03 amendment 2)

Record only: the development version is not used in this project (G, 2026-10-22).

| Aspect | `_operational` (ETC pipeline) | `_developing` |
|---|---|---|
| Unstable class bound (1/L) | `-0.056 <= 1/L < -0.003`: the near-neutral class (`-0.016 <= 1/L < -0.004`) is **unreachable** (caught by "unstable"), and `[-0.004, -0.003)` is unstable too | `-0.056 <= 1/L < -0.016` |
| When a height is computed | if **L or u\*** is present (`!is.na(1/MOL) \| !is.na(1/ustar)`); with u\* missing the formulas give NA; with L missing the class test errors (never reached with L filled by LOCF) | if L **and** u\* are present |
| Nieuwstadt (1981) fallback for neutral / stable | **active** (used where the primary formula gives NA) | commented out |
| B&G previous height | the gap-interpolated h0 (`PBL.h.0.int[i-1]`), else the B&G value; **row 1 skipped** (`i != 1`, no legacy issue 5) | the h0 after interpolation (same when `GF.PBL = TRUE`); row 1 with H > 0 errors |
| Gap filling (`na.approx(maxgap = 2)`) and 5-point smoothing | **always on** (no options); `stats::filter` explicit | options `GF.PBL`, `smooth.PBL` (both TRUE in the ETC call); `filter()` can be masked by dplyr |
| Packages | `zoo` only | `zoo`, `data.table`, `bigleaf`, `crayon` |
| Other | `TIMESTAMP` may be POSIXlt; `PBL.dir` argument for plots / output; with neither `AP` nor altitude, `break()` (bug; not reached, AP passed) | `stop()` with a message |

Effect at BE-Bra (spec 03 part A, rebuilt L2 inputs): `_operational` reproduces the operational column (median
|relative difference| 0.009 %, p95 0.06 %, the four fill values within 0.006 %, 1,341 of 1,342 fill blocks), while
`_developing` leaves 75 % of the fill blocks unexplained (p95 49 %, fill values off by up to 3.3 %).

## `proc_PBL_height.R` (added 2026-10-16): ETC post-processing of `calc_PBL_height()` output

Steps, in order: (1) PBLH < 0 → NA; (2) PBLH < 50 m → 50 m; (3) each NA → mean of all "day" values
(hour of the timestamp 06–18 inclusive) or all "night" values of the same processing run. The fill is computed over the
vector already floored at 50 m; filling in a loop does not change the mean, so each run gets one day and one night
constant. Reproduced exactly at BE-Bra with the hour of `TIMESTAMP_END` (local standard time): 1654.58 / 843.07 m
(2020–2025 run), 1591.68 / 682.86 m (2026 run). Further issues:

11. **Run-mean gap filling.** The fill value depends on the processing period and on outliers (BE-Bra day mean
    1654.58 m vs median ≈ 1509 m of the same values; max 91 km), and ignores season and stability. ≈ 18 % of
    2020–2026 operational rows at BE-Bra are such fills, mostly in daytime (likely causes, to confirm in spec 03:
    B&G chains stopping when both `PBL.h.0[i-1]` and `PBL.h.BG[i-1]` are NA, and the 5-point smoothing returning NA
    near gaps).
12. **"Day" = clock hours 06–18**, not sun position or H sign; same at all latitudes and seasons.
13. **50 m floor** is applied before computing the fill means, and without regard to zm (the commented-out zm-based
    rule is not active).
14. **Upper bound missing.** No cap: implausible values (up to 91 km) pass and enter the fill means.
15. `H.EC` is a free (global) variable in the histogram; `windows()` devices; plots always produced.

## Issues found in `calc_PBL_height_developing.R` (formerly `calc_PBL_height.R`) and the validation script

1. **Temperature units (likely major).** The function expects `AT` in K, but `PBL_validation.R` passes `AT.24` in °C. In the B&G denominator the mechanical term `C·T·u*²` is then about 19× too small, so convective growth is overestimated. This may explain the TODO "PBL sovrastima molto".
2. **Order of gap filling and B&G.** `na.approx(maxgap = 2)` is applied to `PBL.h.0` *before* the B&G step. Short convective runs (≤ 2 steps) get linearly interpolated instead of modelled.
3. **Inconsistent convective trigger.** Stability classes use 1/L, but B&G is triggered by `H > 0`. The two criteria can disagree.
4. **Neutral formula `0.3 u*/|f|`.** At 61°N (f ≈ 1.28e-4 s⁻¹) and u* = 0.5 m/s it gives ≈ 1170 m. A likely source of the "higher in winter" TODO (windy, near-neutral).
5. **First time step.** At `i = 1`, `PBL.h.0[i-1]` is empty, so `if()` errors if `H[1] > 0`.
6. **`smooth.PBL = FALSE` breaks the function.** `PBL.h` is never assigned. The comment says 9-point smoothing while the code uses 5 points. `filter()` is masked if dplyr is loaded (use `stats::filter`).
7. **Pressure units.** Documented as Pa, passed as kPa (only harmless because `rho` is supplied).
8. **`gamma = 0.0025 K/m`.** The comment says ≈ 0.01. This needs a justified choice, ideally a parameter.
9. **Scripts don't run as-is.** `MLH.TS` is undefined (`MLH.df.TS` is the one created). dplyr and lubridate are used without being loaded. There are hardcoded `C:/WORK/DATA/...` paths, `install.packages()` inside a function, and `windows()` devices.
10. **Time conventions not handled.** GMT is assumed everywhere. The timestamp convention (start vs end of averaging period) and the 30-min vs 60-min ERA5 alignment are open TODOs.

## Open TODOs carried over from the legacy notes

- QC: exclude rain, low cloud and fog periods. MLH < 60 m is not resolvable by the lidar.
- Uncertainty of the method, of the lidar MLH, and of temporal aggregation.
- Bland–Altman per season. Error boxplots at fixed hours (Haefele 2014). 4-h means (Hellén et al. 2019).
- Is PBLH (from EC) comparable with lidar MLH at all? Address this in the design.
