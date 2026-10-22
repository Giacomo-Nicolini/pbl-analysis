# Spec 03: Legacy regression (operational ICOS PBLH vs legacy code vs PBLHmod v2)

Phase 1. Branch: `spec-03-legacy-regression` (from `main`). Independent of spec 04.
Read first: `CLAUDE.md`, `docs/legacy_review.md` (all issues, including 11–15 on `proc_PBL_height.R`),
`docs/specs/02_pblh_model.md`, `docs/STATUS.md` (spec 03 findings).

## Goal

1. **Reproduce** the operational `PBLH` column of the ICOS L2 FLUXES file at BE-Bra by running the legacy code
   (`legacy/calc_PBL_height_operational.R` + `legacy/proc_PBL_height.R`, Amendment 2) on the same inputs.
2. **Port** the legacy algorithm into new, testable code with switches, proven identical to the original.
3. **Attribute** the difference between operational PBLH and PBLHmod v2 to individual legacy choices, by switching
   them one at a time (and cumulatively) from the legacy baseline towards PBLHmod v2.

Already established (2026-10-16, see `STATUS.md`): `proc_PBL_height()` floors at 50 m and fills NA with the
day/night mean of the processing run (day = hour 06–18 of `TIMESTAMP_END`, local standard time); the fill values
1654.58 / 843.07 m (2020–2025 run) and 1591.68 / 682.86 m (2026 run) are reproduced exactly from the operational
column itself. This spec must reproduce them **from the inputs**.

## Exception to the `legacy/` rule

`CLAUDE.md` rule 5 allows sourcing `legacy/` only in the explicit legacy regression. In this spec, only
`scripts/03_legacy_regression.R` and `tests/testthat/test-legacy_port.R` may `source()` the legacy files
`calc_PBL_height_operational.R`, `calc_PBL_height_developing.R` and `proc_PBL_height.R` (Amendment 2), into dedicated
environments (`sys.source(..., envir = legacy_env)`), never into the global environment, and never from `R/`. The
legacy files are not edited.

Running the original requires packages it loads itself (`zoo`, `crayon`, `bigleaf`, `data.table`). Its code calls
`install.packages()` for missing ones: **the script must check that all four are installed and abort otherwise**
(never let the legacy code install packages). Add them to renv as legacy-only dependencies (note this in
`renv.lock` commit message and in `docs/legacy_review.md`).

Call the originals with all plotting/saving off: `calc_PBL_height(..., plot.atmc = FALSE, plot.PBL = FALSE,
save.PBL = FALSE)`; `proc_PBL_height(..., plot.PBL.h.freq = FALSE, PBL.dir = <tempdir>/)` (with `PBL.dir` set,
its final plot goes to a png in the temp dir instead of `windows()`; delete it afterwards).

## 1. Inputs of the operational run (confirmed by G, 2026-10-17)

The ETC pipeline call (G):

```r
calc_PBL_height(
  TIMESTAMP = set.cleaned.timestamp,  # LST, %Y-%m-%d %H:%M:%S
  LAT = st.lat, zm = H.EC, d = d.,    # d. = 2/3 canopy height
  MOL = data.table::nafill(set.cleaned$MO_LENGTH, "locf"),
  ustar = set.cleaned$USTAR, H = set.cleaned$H,
  AT = set.cleaned$AT,   # [K]  EddyPro full_output air_temperature
  AP = set.cleaned$AP,   # [Pa] EddyPro full_output air_pressure (no biomet: from site altitude)
  rho = set.cleaned$rho, # [kg m-3] EddyPro full_output air_density
  cp = set.cleaned$cp,   # [J kg-1 K-1] EddyPro full_output air_heat_capacity
  altit = st.alt)
```
then `proc_PBL_height()`. Implement as one function `legacy_inputs(d, site_meta)`.

**Data source: only the official, published ICOS L2 files** (FLUXES, plus METEO/FLUXNET where stated), read with
the spec-01 readers. Do not use the ETC cleanset or EddyPro full_output files, even though they hold the same data
(G, 2026-10-17).

| legacy argument | value | notes |
|---|---|---|
| `TIMESTAMP` | `TIMESTAMP_END`, local standard time | consistent with the fill reproduction (day/night from the hour of `TIMESTAMP_END`) |
| `LAT`, `altit` | `sites.csv` `lat`, `altitude_m` | – |
| `zm`, `d` | `zm_m` (time-varying), `2/3 × canopy_height_m` | in `calc_PBL_height()` they only enter `zdl`, computed but never used: verify and say so |
| `MOL` | FLUXES `MO_LENGTH` with **last-observation-carried-forward** filling (`data.table::nafill(type = "locf")` semantics; implement in base R/dplyr, no data.table in `R/`) | leading `NA` stay `NA` |
| `ustar`, `H` | FLUXES `USTAR`, `H` ("cleaned" = the published L2 values) | – |
| `AT`, `AP`, `rho`, `cp` | see 1.1 | – |
| options | `GF.PBL = TRUE`, `smooth.PBL = TRUE` | – |
| processing runs | 2020-01-01 → 2025-12-31 and 2026-01-01 → end of file, each processed separately (calc + proc) | control: one run over all years (one day and one night fill value; the data show two) |

Rows must form a continuous half-hourly series (the legacy aborts on jumps); use `NA` for missing values.

### 1.1 Ambient variables (`AT`, `AP`, `rho`, `cp`)

ETC pipeline (G, 2026-10-17): EddyPro → RFlux QC → "cleanset" → ICOS L2 FLUXES (plus storage and footprint).
No turbulence variable is recomputed for L2, so L2 `H`, `USTAR`, `MO_LENGTH`, `T_SONIC`, `H2O` are exactly the
values the operational run used. `AT`, `AP`, `rho`, `cp` (EddyPro ambient values) are **not published in L2**, so
they are rebuilt from L2 following the EddyPro formulas:

- Take the formulas from the EddyPro documentation (LI-COR, "Calculating ambient and cell parameters"; EddyPro
  manual) and, where the documentation is incomplete, from the open-source EddyPro engine
  (github.com/LI-COR/eddypro-engine); cite file and routine in comments:
  - `AP = 101300 × exp(−altitude / 8200)` [Pa] (no biomet: EddyPro uses site altitude, Campbell & Norman 1998);
    check in the engine whether the altitude includes the measurement height;
  - `AT`: EddyPro ambient temperature without biomet = sonic temperature corrected for humidity (Schotanus et al.
    1983); determine whether L2 `T_SONIC` is the raw or the corrected sonic temperature and apply the correction
    with L2 `H2O` (check its unit) only if needed;
  - `rho` (moist air density) and `cp` (moist air heat capacity) from `AT`, `AP` and `H2O` with EddyPro's formulas
    and constants.

  Unit tests for each reconstructed variable (hand-computed cases). If a formula cannot be established from the
  documentation or the engine source, stop and report (do not substitute other formulas silently).

## 2. Reproduction (`scripts/03_legacy_regression.R`, part A)

- Run the originals with the inputs of section 1 and compare with the operational `PBLH` on all rows where the
  operational value is not `NA`: share of rows with |legacy − operational| < 0.01 m; distribution of the relative
  difference (median, p95, max); by year and regime; and the four fill values.
- Acceptance: median |relative difference| < 0.5 %, p95 < 2 %, and the four fill values within 0.5 %. Exact
  reproduction is not required, because the rebuilt ambient variables may differ slightly from EddyPro's; report
  the share of rows within 0.01 m anyway (if the rebuild is exact, it should be close to 100 %).
- Also report the sensitivity of the reproduction to each ambient variable (replace one at a time by the simplest
  alternative: `AT = T_SONIC`, `AP` from altitude, `rho`/`cp` computed inside the legacy function). If the
  acceptance is not met, report the pattern of mismatches (by year, regime, hour, near gaps) and stop part C until
  discussed with G.

## 3. Port (`R/legacy_port.R`)

`legacy_pblh(t_local, lat, zm, d, L, ustar, H, T_K, p_Pa = NULL, altitude_m = NULL, sw = legacy_switches())`
reimplements `calc_PBL_height()` followed by `proc_PBL_height()` in vectorised/pure style (no plots, no package
loading). It returns `pblh_m` plus intermediate columns (`atmos_cond`, `h0_m`, `h_bg_m`, `h_smooth_m`,
`flag_floor`, `flag_fill`).

`legacy_switches()` returns a named list; defaults reproduce the legacy exactly:

| switch | legacy (default) | PBLHmod v2 value | legacy issue |
|---|---|---|---|
| `gamma` | 0.0025 K m-1 | 0.01 | 8 |
| `k` | 0.41 | 0.4 | – |
| `omega` | 2π/86400 | 7.2921159e-5 | – |
| `stability_classes` | Essa et al. 2006 1/L classes | regime from L sign + daily window | 3 |
| `stable_formula` | Zilitinkevich, c_s = 0.4 | Nieuwstadt B1 | – |
| `neutral_formula` | 0.3 u*/abs(f) | auto B2/B3 | 4 |
| `convective_trigger` | `H > 0` | daily window from sunrise (spec 02) | 3 |
| `bg_integration` | explicit Euler, one step per half-hour | RK4 with sub-steps | – |
| `gap_interp_before_bg` | TRUE (`na.approx(maxgap = 2)` on h0 before B&G) | FALSE | 2 |
| `smoothing` | 5-point centred moving mean | none | 6 |
| `floor_m` | 50 | none (`NA`) | 13 |
| `fill` | `"run_mean_day_night"` (06–18 h of the timestamp) | `"none"` | 11, 12 |
| `negative_to_na` | TRUE | – | – |
| `mol_locf` | TRUE (L filled by last observation carried forward before the call) | FALSE (spec 02: 1/L interpolated over ≤ 1 h, longer gaps NA) | – |
| `forcing` | ETC inputs (section 1: EddyPro/L2 `AT`, `AP`, `rho`, `cp`; L2 `H`) | PBLHmod v2 forcing (`T_K`, `p_Pa`, `rho`, `wtheta` from `prepare_pblh_forcing()`, H = `H_F_MDS`) | – |

Reproduce the legacy exactly, including its quirks, when all switches are at their defaults: e.g. B&G uses
`h0[i-1]` when not `NA`, else the B&G value of `i-1`; `i = 1` with `H > 0` (legacy issue 5: decide the port's
behaviour = return `NA` for that row and document it; check it never occurs in the BE-Bra runs, otherwise report);
smoothing returns `NA` where any of the 5 values is `NA`; the run-mean fill is computed after flooring.

Constants used only by the port go in `legacy_switches()` (documented as "legacy values"), not in
`pblh_constants()`.

## 4. Port equality (test + script part B)

- `tests/testthat/test-legacy_port.R`: on a synthetic 3-day series with gaps, NA runs, sign changes and stable /
  neutral / convective rows, `legacy_pblh()` with default switches equals the original legacy pipeline (sourced as
  described above) to 1e-8 m on every row (including `NA` positions). Also unit tests for each switch (one
  hand-checked case each, e.g. floor, fill value = mean of day rows of the run, smoothing window).
- Script part B: the same equality on the full BE-Bra series for the chosen input variant (report max difference;
  must be < 1e-6 m).

## 5. Attribution (script part C)

On BE-Bra, with the chosen inputs:

1. **One at a time**: from the legacy baseline, change each switch to its PBLHmod v2 value (others at legacy
   values). Variants named `legacy_<switch>`.
2. **Cumulative path**: switches changed in this order, each run keeping the previous changes:
   `fill → floor_m → smoothing → gap_interp_before_bg → mol_locf → forcing → k, omega → gamma → neutral_formula →
   stable_formula → convective_trigger + stability_classes → bg_integration`. The last step should approach PBLHmod v2; report the
   residual difference to `pblhmod_default` and explain it (remaining structural differences: daily window,
   sunrise h0, forcing preparation, H source, interpolation rules).
3. For every variant: quantiles (1, 5, 25, 50, 75, 95, 99 %) of `pblh_m` overall and by regime (regime from
   PBLHmod v2 on the same row), share of `NA`, share floored, share filled; median and MAD of
   `log(variant / pblhmod_default)` by regime; the same for `log(variant / operational)`.

Output tables: `output/tables/03_reproduction_<site>.csv`, `03_variants_<site>.csv` (long format: variant,
stratum, statistic, value). Output data: `<derived_dir>/legacy/<site>_legacy_<variant>.csv.gz`.

## 6. Figures

Using the `save_figure()` helper from spec 04 if merged; otherwise create `R/plot_utils.R` with the same
interface (coordinate with spec 04: one implementation only).

1. `03_reproduction_<site>`: operational vs legacy-reproduced (hexbin/2-D density) and difference time series.
2. `03_waterfall_<site>`: median `pblh_m` by regime along the cumulative path (one point per step, lines per
   regime), with PBLHmod v2 as a reference line.
3. `03_one_at_a_time_<site>`: median log-ratio to PBLHmod v2 per switch and regime (dot plot).
4. `03_diurnal_<site>`: median diurnal cycle (by season) of operational, legacy baseline, legacy without
   post-processing (fill, floor, smoothing off) and PBLHmod v2.

## Acceptance criteria

- Tests pass (old and new).
- Reproduction at BE-Bra with a documented input variant meets section 2 (median |relative difference| < 0.5 %,
  p95 < 2 %, the four fill values within 0.5 %), with the share of rows within 0.01 m reported; or a
  stop-and-report if not (Amendment 1; the earlier "≥ 99 % exact, fill values exact" criterion is withdrawn).
- Port equal to the original (< 1e-6 m) on the full BE-Bra series.
- Attribution tables and figures produced; summary printed.
- Legacy files unchanged (`git diff legacy/` empty); no legacy code sourced outside the two allowed files.

## Out of scope

Changing PBLHmod v2; fixing the ETC pipeline; other sites; ERA5 (spec 04).

## Summary to report back

As in `CLAUDE.md`, plus: the input variant that reproduces the operational column (to be confirmed by G); the
reproduction statistics; the ranking of legacy choices by their effect on the median per regime (one at a time and
along the cumulative path); the residual difference legacy-fully-switched vs PBLHmod v2 and its explanation; any
legacy behaviour that the port had to emulate deliberately (quirks).

## Amendment 1 (2026-10-21, agreed with G before implementation)

1. **Acceptance of the reproduction** is the section 2 criterion: median |relative difference| < 0.5 %, p95 < 2 %,
   the four fill values within 0.5 %; the share of rows within 0.01 m is reported. The "≥ 99 % exact, fill values
   exact" line in the acceptance criteria is corrected accordingly (the ambient variables are rebuilt from L2, so an
   exact match is not expected).
2. **Processing runs.** Each run starts at its first row with data (first run: `TIMESTAMP_END` 2020-10-26 00:30,
   the first FLUXES row with `MO_LENGTH`/`PBLH`; the FLUXES file starts on 2020-01-01 with empty rows, and the
   operational column is `NA` before that date, so those rows were not in the operational run). Split at
   `TIMESTAMP_START` 2026-01-01 00:00, local standard time (confirmed by G). The split must reproduce the four fill
   values; otherwise stop and report.
3. **Exception to "Out of scope" (changing PBLHmod v2):** `bg_step()` and `run_pblh_model()` get the option
   `bg_integration = c("rk4", "euler")` (default `"rk4"`). `"euler"` (one explicit step per half-hour, end-of-interval
   height, as the legacy code) is used only by the cumulative path of section 5, at the step where every other switch
   already has its PBLHmod v2 value. The default output is unchanged: the regression chain (`nofilter_l2` = spec 02,
   `l_L2` = spec 05, `l_check_allday` = spec 06 amendment 1) was rerun after the change and is identical
   (< 1e-9 m, same NA, regime, method).
4. **Input data facts found while planning (BE-Bra):** L2 `H2O` is in mmol mol-1 (median 9.75, range 4–17); L2
   `T_SONIC` is 5.2 K above METEO `TA` (median, IQR 3.6–6.8 K), more than the ≈ 0.9 K humidity effect (sonic
   temperature offset, inherited by the EddyPro ambient `AT`; item in the ETC QC report).
5. **G's answers (2026-10-22), recorded here because no Amendment 2 was needed** (the hypothesis test below did
   not change the inputs): (1) in the ETC cleanset, EddyPro −9999 values in `air_density`, `air_heat_capacity`,
   `air_temperature` are converted to NA before `calc_PBL_height()` is called; (2) the operational `PBLH` was
   computed from the same cleanset that became the published L2.
6. **Hypothesis test `amb_na` (2026-10-22), rejected.** Hypothesis: the operational `AT`/`rho`/`cp` were NA where
   EddyPro could not compute them, breaking the B&G chain and leaving the daytime NA blocks that
   `proc_PBL_height()` fills. Note: in the engine source (flux_params.f90) neither a missing H2O nor an RH outside
   0–130 % makes these three outputs −9999 (dry-air fallback; only RH becomes −9999), so the test is empirical.
   BE-Bra results: of the 15,136 rows filled in the operational column but computed in the reproduction, 1.3 % have
   L2 H2O missing and 1.2 % an engine RH outside 0–130 % (2.5 % either), vs 0.8 % / 0.7 % (1.4 %) on rows computed
   in both. Variant `amb_na` (AT, rho, cp NA there) vs `rebuilt`: fills 4,898 / 1,738 (day / night) vs 4,884 /
   1,732 (operational 18,298 / 3,166); operational fill blocks (1,342) reproduced ≥ 50 %: 336 vs 335, fully: 113 vs
   113; fill values 1601.81 / 846.11 / 1584.28 / 705.30 m (operational 1654.58 / 843.07 / 1591.68 / 682.86); p95
   |rel. diff| 49.4 %; within 0.01 m 9.5 %. Acceptance (section 2) still not met; part A stopped. Next: the ETC
   inputs for one month from G.
7. **Hypothesis test `tau_na` (2026-10-22), rejected.** Hypothesis: in the cleanset, u* was removed where the
   momentum flux failed QC (L2 `TAU` NA); in the legacy code a missing u* gives no h0 and no B5 step, so the
   convective chain stops. Overlap: L2 `TAU` is NA on 31.2 % of the 15,136 "operational fill, reproduction
   computes" rows vs 24.0 % of the 76,721 rows computed in both (29.6 % of all rows). Part A with `tau_na` (and
   `tau_amb_na`): fills 17,607 / 17,755 day / night (operational 18,298 / 3,166); operational fill blocks reproduced
   ≥ 50 %: 768 of 1,342 (57 %), fully: 383; share of operational fill rows filled 58 % (rebuilt 29 %); fill values
   1590.24 / 796.62 / 1649.76 / 734.82 m (operational 1654.58 / 843.07 / 1591.68 / 682.86); median |rel. diff|
   3.9 %, p95 1057 %; within 0.01 m 3.7 %. It reproduces more daytime blocks but creates about 14,600 night fills
   that the operational column does not have, and moves all four fill values further away. Acceptance not met; part
   A stopped. Next: the ETC cleanset inputs for one month from G.

## Amendment 2 (2026-10-22, agreed with G)

1. **Operational reference.** G reorganised `legacy/` (commit b020a93): `calc_PBL_height.R` →
   `calc_PBL_height_developing.R` (G's later development version, used in spec 03 until now), new
   `calc_PBL_height_operational.R` (the version the ETC pipeline sources), pipeline scripts `FFP_analysis_dev.R` and
   `FFP_calc_parallel_dev.R`, Kljun's FFP code (local only, gitignored); `proc_PBL_height.R` unchanged. The
   operational reference is `calc_PBL_height_operational.R` + `proc_PBL_height.R` (`FFP_analysis_dev.R` line 133;
   output written as `round(PBL, 2)`, no other change). Parts A–C target it. The differences from `_developing` are
   documented in `docs/legacy_review.md` (diff table): unstable class bound −0.003 instead of −0.016 (near-neutral
   class unreachable); a height is computed if L **or** u* is present; Nieuwstadt fallback active; B&G uses the
   gap-interpolated h0 for the previous step and skips row 1; gap filling and smoothing always on.
2. **Part A rerun (2026-10-22), accepted** (operational version, rebuilt L2 inputs, two runs): median |relative
   difference| 0.009 %, p95 0.06 %, rows within 0.01 m 13.3 %; fill values 1654.50 / 843.07 / 1591.59 / 682.86 m
   (operational 1654.58 / 843.07 / 1591.68 / 682.86; max deviation 0.006 %); fills 18,454 / 3,297 day / night
   (operational 18,298 / 3,166: the extra rows are after the end of the operational run 2, 2026-08-26 → 09-01,
   where the operational column is NA; they do not change the fill means); fill blocks reproduced 1,341 of 1,342.
   Ambient-variable variants change nothing material (`rho_cp_legacy`: within 0.01 m 16.0 %); one run over all
   years: p95 2.1 %, fill values off by 0.4–2 %, so the two runs are confirmed. The development version on the same
   inputs: p95 49 %, 75 % of the fill blocks not reproduced: the earlier part-A failure (Amendment 1, notes 5–7)
   came from the code version, not from the inputs; `amb_na` and `tau_na` are moot.
3. **Cleanset diagnostic** (local file `cfg$cleanset_file`, never committed; run 1 only, 2020-10-26 → 2026-01-01):
   the operational version on the cleanset gives median |relative difference| 0.0004 %, p95 0.006 %, 77.9 % within
   0.01 m, fill values 1654.51 / 843.07 m, 1,156 of 1,157 blocks. Cleanset vs L2 inputs: identical missing-value
   patterns for MO_LENGTH, USTAR, H, AT, rho, cp (AP: NA on 1,991 rows in the cleanset, constant in the rebuild);
   values differ only by rounding (median |diff|: MO_LENGTH 0.002 m, USTAR 0.00025 m s-1, H 0.0025 W m-2, AT
   0.0025 K, rho 1e-5 kg m-3, cp 0.003 J kg-1 K-1; AP 0.47 Pa). The remaining L2-vs-operational differences
   (13 % of rows within 0.01 m) are therefore rounding of the published inputs, not code.
4. **Port (part B)** emulates `_operational` by default, with switches for the five differences of item 1, so that
   part C can include "operational → developing" as a step.
5. **Development version out of scope** (G, 2026-10-22): `calc_PBL_height_developing.R` is not used in this project,
   in any phase or spec; only the operational code counts. Script 03 no longer sources it; there are no switches and
   no attribution step for the operational/developing differences. The file stays in `legacy/` as a read-only
   record (`docs/legacy_review.md`). Legacy-only packages: only `zoo` (needed by the operational version);
   `bigleaf` and `data.table` removed from renv; `crayon` kept (dependency of `readr`, `testthat`).
6. **Run 2 ends at the last operational row** (`TIMESTAMP_END` 2026-08-26 00:00, local standard time;
   `legacy_inputs(end_local)`). Part A rerun: fills 18,298 / 3,165 day / night (operational 18,298 / 3,166), fill
   values 1654.50 / 843.07 / 1591.59 / 682.86 m, median |rel. diff| 0.009 %, p95 0.06 %: accepted.
7. **Port (part B), implementation notes.** `legacy_pblh(inp, site_meta, sw, forcing, icos)` takes the
   `legacy_inputs()` table instead of the vectors of section 3 (the `forcing` switch needs the PBLHmod forcing and
   the daily window needs the site data). Defaults reproduce `calc_PBL_height_operational.R` + `proc_PBL_height.R`:
   synthetic 3-day series < 1e-8 m (test), full BE-Bra series max |diff| 0 m with the same NA pattern. Deviation
   from the original only where it would stop (missing L with u* present: `NA`). Switch notes: `k` and `omega`
   change together (`k_omega`); `convective_trigger` and `stability_classes` change together and, with the daily
   window, the core is `run_pblh_model()` (PBLHmod v2 with the switch's `gamma` and `bg_integration`), which
   requires the PBLHmod values of all other core switches, so they have no one-at-a-time variant;
   `bg_integration = "rk4"` inside the legacy loop (one-at-a-time) keeps the end-of-interval height.
8. **Attribution (part C), BE-Bra.** Cumulative path in the order of section 5; the last step equals
   `pblhmod_default` (max |diff| 4.6e-13 m, same NA pattern): the residual is zero because, with all switches at
   their PBLHmod values, the port delegates to `run_pblh_model()`; the structural differences (daily window, sunrise
   h0, holds and stops, regime from the L sign, interpolation rules) are attributed to the
   `convective_trigger_stability_classes` step. Results in `output/tables/03_variants_BE-Bra.csv` and the summary in
   `docs/STATUS.md`.
9. **Common-rows attribution (review 2026-10-22).** The medians of part C on all valid rows of each variant mix
   value changes with sample changes (e.g. the run-mean fill "−383 m" is rows dropping out). Added: medians (overall
   and by regime) on the rows valid in every variant of the analysis and in PBLHmod v2, with the number of valid rows
   per variant (both versions in `03_variants_<site>.csv`, strata `common_one_at_a_time|...`,
   `common_cumulative|...`; statistics `n_valid`, `n`). Figures: `03_waterfall_` and `03_one_at_a_time_` show the
   common rows; `03_waterfall_allrows_`, `03_one_at_a_time_allrows_` the all-rows version with n per step.
   BE-Bra, common rows (median PBLH, change in m; all / stable / unstable_diag / convective):
   - One at a time (n = 79,263): baseline 707 / 379 / 1317 / 1652 m; neutral formula −211 / −40 / −595 / −740;
     smoothing −192 / −71 / +179 / +50; gap interpolation +176 / +128 / +376 / +98; stable formula −38 / −50 / −9 /
     −3; gamma −11 / 0 / −11 / −101; forcing −8 / +1 / +20 / +5; k/omega −2; bg_integration −1; fill, floor,
     mol_locf, at_meteo 0. PBLHmod v2 on the same rows: 364 / 263 / 557 / 610 m.
   - Cumulative (n = 75,689): baseline 724 m → PBLHmod v2 366 m; steps smoothing −206, neutral formula −82
     (unstable_diag −1243, convective −970), stable formula −46, bg_integration −10, convective trigger + classes
     −8, gamma −8 (convective −99); fill, floor, gap interpolation, mol_locf, forcing, k/omega ≤ 2 m.
   On all valid rows the fill (−383 m) and forcing (±160–190 m) effects come from changed row sets; on common rows
   both are near zero. The neutral formula (operational class bound −0.003 → 0.3 u*/|f| on many daytime rows) and
   the 5-point smoothing remain the main legacy effects.
10. **`at_meteo`** (one at a time; AT = METEO TA instead of the EddyPro AT, 2,247 rows without METEO TA keep the
    EddyPro AT): negligible (22,341 rows change by > 0.01 m; relative change p95 0.0009 %, max 1.5 %; medians and
    fill values unchanged). The sonic temperature offset (local ETC QC report, item 4) does not affect the
    operational PBLH materially, because AT enters only the mechanical term of B&G.
11. **Figures** are written as png only (`save_figure(formats = "png")` default since 2026-10-22; pdf only for paper
    figures).
