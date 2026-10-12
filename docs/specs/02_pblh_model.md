# Spec 02: PBLH model (PBLHmod v2)

Phase 1. Branch: `spec-02-pblh-model` (from `main`).
Read `CLAUDE.md`, `docs/project_plan.md` (Phase 1) and `docs/references/kljun2015_appB_notes.md` first.

## Goal

A tested, documented R implementation of the PBL height model of Kljun et al. (2015, App. B): diagnostic formulas for
neutral/stable conditions and the Batchvarova & Gryning (1991) slab model for convective conditions, integrated as a
daily cycle. Input: the output of `read_icos_site()` (spec 01). Output: raw half-hourly PBLH with regime, method and
flags. Tested locally on BE-Bra.

**The output is raw**: no smoothing, no clipping (no minimum height, no 50 m floor), no gap-filling of PBLH.

## Files to create

```
R/pblh_constants.R      physical constants and model parameters (single source of truth)
R/solar.R               solar elevation and sunrise/sunset
R/pblh_forcing.R        preparation of the model forcing from read_icos_site() output
R/pblh_diagnostic.R     diagnostic formulas (B1, B2, B4)
R/pblh_bg.R             Batchvarova & Gryning integrator (B5)
R/pblh_model.R          pblh_options() and run_pblh_model() (orchestration)
scripts/02_run_pblh_model.R
tests/testthat/test-solar.R
tests/testthat/test-pblh_diagnostic.R
tests/testthat/test-pblh_bg.R
tests/testthat/test-pblh_forcing.R
tests/testthat/test-pblh_model.R
```

## 1. Constants (`R/pblh_constants.R`)

One function `pblh_constants()` returning a named list (with comment and reference for each value):

| Name | Value | Unit | Source |
|---|---|---|---|
| `k` | 0.4 | – | von Kármán constant, as in Kljun et al. 2015 App. B |
| `g` | 9.80665 | m s-2 | standard gravity |
| `omega` | 7.2921159e-5 | rad s-1 | Earth's rotation rate (sidereal day). Legacy used 2π/86400 (solar day): document the change |
| `r_d` | 287.05 | J kg-1 K-1 | gas constant of dry air |
| `cp` | 1005 | J kg-1 K-1 | specific heat of dry air at constant pressure |
| `bg_A`, `bg_B`, `bg_C` | 0.2, 2.5, 8 | – | Batchvarova & Gryning 1991, as in Kljun 2015 Eq. B5 |
| `nieuwstadt_a`, `nieuwstadt_b` | 3.8, 2.28 | – | Nieuwstadt 1981, Kljun 2015 Eq. B1 |
| `c_n` | 0.3 | – | Kljun 2015 Eq. B2 (Hanna & Chang 1993) |

## 2. Solar geometry (`R/solar.R`)

FLUXNET `SW_IN_POT` is not available for every period (e.g. at BE-Bra FLUXNET ends in 2025 while FLUXES runs to 2026),
so sunrise is computed from astronomy, always the same way.

- `solar_elevation(timestamp_utc, lat, lon)`: solar elevation angle [deg], NOAA solar-position algorithm
  (NOAA Global Monitoring Laboratory solar calculator equations; cite in the code). Vectorised.
- `sun_times(date_utc, lat, lon)`: sunrise and sunset [POSIXct UTC] per date, geometric (elevation = 0°,
  no refraction correction; document this choice).
- `is_day(timestamp_start_utc, timestamp_end_utc, lat, lon)`: TRUE if the solar elevation at the **midpoint** of the
  half-hour is > 0.
- Test against known values: e.g. sunrise/sunset for a fixed location and date compared with the NOAA calculator
  (tolerance ±3 min). Polar day/night must not crash (return `NA` sunrise/sunset and `is_day` all TRUE / all FALSE).

## 3. Forcing (`R/pblh_forcing.R`)

`prepare_pblh_forcing(data, site_meta, opts)` takes `read_icos_site()` output and returns, per half-hour:
`timestamp_start_utc, timestamp_end_utc, ustar, L, H, H_source, T_K, T_source, p_Pa, p_source, rho, wtheta,
f_cor, is_day, excl_flag`, and the flags `flag_ustar_interp, flag_L_interp, flag_H_interp, flag_H_qc_masked`.

Rules:
- `ustar` = `USTAR` [m s-1]; `L` = `MO_LENGTH` [m] (FLUXES; not gap-filled).
- `H` [W m-2] according to `opts$h_source`:
  - `"H_F_MDS"` (default): `fn_H_F_MDS`, with values where `fn_H_F_MDS_QC > opts$h_qc_max` (default 1) set to `NA`
    (`flag_H_qc_masked = TRUE`). QC meaning (FLUXNET convention): 0 measured, 1 good-quality gap-fill, 2 medium,
    3 poor — verify in the data and state the distribution in the summary.
    Where `fn_H_F_MDS` is unavailable (e.g. after the end of the FLUXNET file), fall back to `H` (FLUXES) if
    `opts$h_fallback = TRUE` (default). Record the source per row in `H_source`.
  - `"H"`: `H` from FLUXES only. `"H_CORR"`: `fn_H_CORR` only.
- `T_K`: `fn_TA_F_K` → fallback `met_TA_K` → fallback `T_SONIC_K`; record `T_source`.
- `p_Pa`: `fn_PA_F_PA` → fallback `met_PA_PA` → fallback barometric estimate from `altitude_m`
  (`p = 101325 * exp(-altitude_m / 8434)`, isothermal scale height ≈ 8.4 km; cite). Record `p_source`.
- `rho = p_Pa / (r_d * T_K)` [kg m-3] (ideal gas law, dry air).
- `wtheta = H / (rho * cp)` [K m s-1] (surface kinematic heat flux).
- `f_cor = 2 * omega * sin(lat)` [s-1]; use `abs(f_cor)` in the formulas (all sites are in the northern hemisphere,
  but do not rely on it).
- Short gaps: `ustar`, `L`, `wtheta`, `T_K`, `p_Pa` are linearly interpolated over gaps of at most
  `opts$max_gap_interp` consecutive half-hours (default 2 = 1 h); set the matching `flag_*_interp`.
  Longer gaps stay `NA`. Never interpolate across excluded rows.
- Excluded rows (`excl_flag = TRUE`): all forcing set to `NA`.

## 4. Diagnostic formulas (`R/pblh_diagnostic.R`)

Vectorised, pure functions; inputs in SI; return height [m]:
- `h_nieuwstadt(L, ustar, f)`: Eq. B1, `h = L/3.8 * (-1 + sqrt(1 + 2.28 * ustar / (|f| * L)))`, for `L > 0`.
  For `L <= 0` return `NA`. For very large `L` (|ustar / (f L)| < 1e-6) use the neutral limit (B2) to avoid loss of
  precision.
- `h_neutral(ustar, f, c_n = 0.3)`: Eq. B2, `c_n * ustar / |f|`.
- `h_zilitinkevich(L, ustar, f, c_s = 0.4)`: Eq. B4, `c_s * sqrt(ustar * L / |f|)`, for `L > 0`, else `NA`.
- `h_diagnostic(L, ustar, f, opts)`:
  - `L > 0`: `opts$stable_method` = `"nieuwstadt"` (default) → B1; `"zilitinkevich"` → B4 with `opts$c_s`;
    `"neutral"` → B2.
  - `L < 0` (unstable, outside the convective integration window, e.g. before the daily start or after the end):
    B2 (mechanical/neutral estimate).
  - `L` or `ustar` missing → `NA`.
  - Returns height and a `method` label (`"B1"`, `"B4"`, `"B2"`).

## 5. Batchvarova & Gryning integrator (`R/pblh_bg.R`)

- `bg_rate(h, wtheta, ustar, L, T_K, gamma, const)`: Eq. B5, dh/dt [m s-1]:
  `(wtheta/gamma) / ( h^2 / ((1+2A) h - 2 B k L) + C ustar^2 T / (gamma g ((1+A) h - B k L)) )`.
  If either denominator `(1+2A) h - 2 B k L` or `(1+A) h - B k L` is ≤ 0, return `NA` (the caller holds h constant
  and flags).
- `bg_step(h0, wtheta, ustar, L, T_K, gamma, dt = 1800, n_sub = 30)`: integrate over one half-hour with `n_sub`
  explicit sub-steps (default 60 s), forcing held constant over the half-hour. Returns the height at the end of the
  interval **and** the mean of the sub-step heights over the interval. Use 4th-order Runge–Kutta for each sub-step.
- `gamma` may be a scalar or a vector with one value per half-hour (interface for the later ERA5 options).

## 6. Daily cycle and orchestration (`R/pblh_model.R`)

### `pblh_options(...)`
Returns a named list with defaults (all overridable):
`stable_method = "nieuwstadt"`, `c_s = 0.4`, `gamma = 0.01` [K m-1], `h_source = "H_F_MDS"`, `h_qc_max = 1`,
`h_fallback = TRUE`, `max_gap_interp = 2` [half-hours], `n_sub = 30`, `n_neg_end = 2` [half-hours],
`combine = "max"`, `label = "default"`.
Validate values (abort on unknown options).

### `run_pblh_model(data, site_meta, opts = pblh_options())`
For each **day** (from one sunrise to the next, computed with `R/solar.R`):

1. **Night / before the daily start**: `pblh_m = h_diagnostic()`, `regime = "stable"` if `L > 0`,
   `"unstable_diag"` if `L < 0`.
2. **Daily start at sunrise**: the initial height `h0` is the diagnostic value of the last half-hour before sunrise;
   if missing, the nearest available diagnostic value in the 2 h before sunrise; if none, no integration that day
   (`flag_no_h0 = TRUE`, diagnostic values only).
3. **Convective window**: from sunrise, integrate with `bg_step()`:
   - while `wtheta <= 0` before the first positive value: hold `h` constant at `h0` (dh = 0);
   - when `wtheta > 0`: integrate B5;
   - short negative spells after the first positive value (fewer than `opts$n_neg_end` consecutive half-hours, e.g.
     clouds): hold `h` constant;
   - the window ends at the first run of `opts$n_neg_end` consecutive half-hours with `wtheta <= 0` after the first
     positive value, or at sunset, whichever comes first;
   - `wtheta > 0` but `L > 0` (sign inconsistency), or `bg_rate()` returns `NA`: hold `h` constant,
     `flag_sign_inconsistent = TRUE`;
   - a forcing gap (after interpolation) inside the window stops integration for the rest of that day
     (`flag_bg_stopped_gap = TRUE`; the remaining half-hours of the window get diagnostic values or `NA`).
   - Output in the window: `h_bg_m` = mean sub-step height over the half-hour; the state carried to the next
     half-hour is the end-of-interval height.
   - `pblh_m = max(h_bg_m, h_diag_m)` if `opts$combine = "max"` (default; mechanical mixing may exceed convective
     growth), or `h_bg_m` if `"bg_only"`. `regime = "convective"`; `method = "B5"` or `"B2"` depending on which
     term was taken.
4. **After the window** (evening, night): back to step 1.

Output: one row per input half-hour, columns
`site_id, timestamp_start_utc, timestamp_end_utc, pblh_m, regime, method, h_diag_m, h_bg_m, day_id (date of the
sunrise, local standard time), is_day, ustar, L, wtheta, T_K, p_Pa, H_source, T_source, p_source, zm_m,
flag_ustar_interp, flag_L_interp, flag_H_interp, flag_H_qc_masked, flag_no_h0, flag_bg_stopped_gap,
flag_sign_inconsistent, flag_below_zm (pblh_m < zm_m), excl_flag`.
Attribute `pblh_options` = the options used (for reproducibility). Excluded rows: `pblh_m = NA`, `regime = "excluded"`.

Do not modify `pblh_m` after this step (no smoothing, clipping or filling).

## 7. Script `scripts/02_run_pblh_model.R`

For each site (`cfg$test_sites` or command-line arguments):
1. read inputs with the spec-01 functions;
2. run `run_pblh_model()` with default options, and with `gamma = 0.0025` as a sensitivity case (`label = "gamma0025"`);
3. write each output to `<derived_dir>/pblhmod/<site>_pblhmod_<label>.csv.gz` (create folders as needed);
4. print a summary: valid % of `pblh_m`; counts per `regime` and `method`; quantiles (1, 5, 25, 50, 75, 95, 99 %) of
   `pblh_m` overall and per regime; the same quantiles for the operational `PBLH` column on the same rows (side by
   side, no statistics beyond that); flag counts; `H_source` / `T_source` / `p_source` counts; the `H_F_MDS_QC`
   distribution; run time.

No plots in this spec.

## 8. Tests

Known-value and limit tests (hand-computed in the test file, with the calculation shown in comments):
- `h_nieuwstadt`: one hand-computed case; tends to `h_neutral` for large L (relative difference < 1e-3 at L = 1e6 m);
  for large stability `h_nieuwstadt / sqrt(ustar L / |f|)` → `sqrt(2.28)/3.8 ≈ 0.397`.
- `h_zilitinkevich`, `h_neutral`: hand-computed cases; `NA` for invalid L.
- `bg_rate`: hand-computed case; `NA` when a denominator is ≤ 0.
- `bg_step`, analytic case: with `ustar = 0` and `L = 0`, B5 reduces to `dh/dt = (1+2A) wtheta / (gamma h)`, so
  `h(t)^2 = h0^2 + 2 (1+2A) wtheta t / gamma`. Check `bg_step` against this (relative error < 1e-4).
- Convergence: `n_sub = 30` vs `n_sub = 300` differ by < 0.5 % on a typical morning case (h0 = 100 m,
  wtheta = 0.1 K m s-1, ustar = 0.3, L = -20, T = 290 K, gamma = 0.01).
- Daily cycle on a synthetic 3-day series (sinusoidal H, constant ustar, L consistent with H): integration starts at
  sunrise, h is non-decreasing inside the window, window ends at the configured rule, and day 2 restarts from the
  diagnostic value (no carry-over).
- Gap rules: gap ≤ 2 half-hours interpolated and flagged; longer gap inside the window stops integration and flags.
- H QC masking and fallback; T and p fallbacks; excluded rows give `NA` and `regime = "excluded"`.
- Solar: sunrise/sunset against NOAA reference values; polar cases.
- Options validation: unknown option values abort.

## Acceptance criteria

- `testthat::test_dir("tests/testthat")` passes (old and new tests).
- `Rscript scripts/02_run_pblh_model.R` runs on BE-Bra and prints the summary.
- Every formula has its reference in a comment; constants only in `pblh_constants.R`.
- `renv::status()` consistent; nothing under `data/` or `output/` committed.

## Out of scope

Smoothing, clipping or filling of PBLH; γ from ERA5 (only the scalar/vector interface); legacy comparison (spec 03);
standard-format output (phase 2); plots; footprint.

## Summary to report back

As in `CLAUDE.md`, plus: the BE-Bra summary from the script (both γ values), the `H_F_MDS_QC` distribution, the
share of days with `flag_no_h0` / `flag_bg_stopped_gap`, the share of `flag_sign_inconsistent`, and anything that
looked physically odd.

## Amendments (2026-10-12, agreed before implementation)

1. **Sunrise reference / zenith.** `sun_times()` gets `zenith_deg = 90` (geometric, default; 90.833 reproduces the
   NOAA calculator with refraction). Unit tests use textbook checks (equation-of-time extremes, solstice declination,
   equinox day length). `scripts/02_run_pblh_model.R` additionally compares `is_day` with `fn_SW_IN_POT > 0` over the
   FLUXNET period and reports the agreement; mismatches should be at most one half-hour at each sunrise/sunset.
2. **Polar day / night.** Days without sunrise (polar day): diagnostic only, `flag_polar_day = TRUE`. Polar night:
   diagnostic only, `day_id` = local standard time date. A polar-day convective rule is an open point.
3. **L interpolation.** Short gaps in `L` are filled by interpolating `1/L` linearly (smooth through neutral), not `L`;
   still flagged `flag_L_interp`.
4. **H fallback scope.** With `h_source = "H_F_MDS"`, the fallback to `H` (FLUXES) applies only where `fn_H_F_MDS` is
   missing in the file; rows masked by `h_qc_max` stay `NA` (`flag_H_qc_masked`).
5. **Days without positive `wtheta`.** If `wtheta` is never > 0 between sunrise and sunset, there is no convective
   window that day: diagnostic only, `flag_no_positive_h = TRUE`. On normal days the rows before the first positive
   `wtheta` stay in the window (h held at `h0`, `regime = "convective"`).
6. **`h_bg_m`** (mean over the half-hour) is the trapezoidal mean of the `n_sub + 1` sub-step nodes.
7. **Monotonicity test** applies to `h_bg_m` (with `combine = "max"`, `pblh_m` can decrease when `h_diag_m` does).
8. **Window start** = first half-hour whose midpoint has solar elevation > 0 (consistent with `is_day`); `h0` from the
   half-hour before, else the nearest valid diagnostic value within the 2 h (4 half-hours) before.
9. **Barometric fallback** `p = 101325 exp(-z / 8434)`: isothermal atmosphere, scale height R_d T / g with T ≈ 288 K
   (e.g. Wallace & Hobbs 2006, Sect. 3.2).

10. *(Implementation decision, to confirm)* A day whose daylight `wtheta` has no positive value **and** contains gaps
    is a gap day (`flag_bg_stopped_gap`, diagnostic only), not a `flag_no_positive_h` day: the missing values could
    have been positive (e.g. days before the start of the H record).

11. **Defaults per Kljun et al. (2015) App. B** (decided after the first BE-Bra run, where B2 was taken in 85 % of
    convective half-hours with `combine = "max"`):
    - `combine = "bg_only"` is the default; `"max"` stays available.
    - New option `neutral_method = c("auto", "B2", "B3")`, default `"auto"`: `N = sqrt(g * gamma / T_K)` per
      half-hour (gamma scalar or vector); B2 (`c_n = 0.3`) where `10 < N/|f| < 70`, otherwise B3
      `h = c_nn * ustar / sqrt(|f| * N)`, `c_nn = 1.36` (Kljun Eq. B3, Zilitinkevich et al. 2012). Applied wherever B2
      was used (unstable rows outside the window, `combine = "max"`, `stable_method = "neutral"`). B1 for stable
      conditions (including its internal neutral limit) is unchanged. New functions `brunt_vaisala()`,
      `h_neutral_nn()`, `h_neutral_select()`; method label `"B3"`.
    - Script cases: `default` (bg_only, auto, gamma 0.01), `gamma0025`, `max_B2` (previous behaviour), `max_auto`.

Added output columns: `flag_polar_day`, `flag_no_positive_h`.
