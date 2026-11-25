# Spec 17: γ from the station's own data (EC-derived γ, "PBLHmod-γEC")

Extension after `analysis-v1` (reopens the stand-by of 2026-11-19). Branch: `spec-17-gamma-ec` (from `main`).
Agreed with G on 2026-11-25. Read first: `CLAUDE.md`, `docs/STATUS.md` (stand-by section),
`docs/specs/02_pblh_model.md` (B&G, amendments 11–14), `docs/specs/11_gamma_era5_davos.md` (+ amendment 1),
`docs/specs/09_cluster_deployment.md` (steps, hashes, two-pass order of `gamma_med`).

## Why

The EC-only default γ = 0.01 K m⁻¹ is too large: in the daytime window at the six aerosol sites the median |ln ratio|
to the measured MLH is 0.43; the ERA5-based γ (`clim`, `gamma_site`, `gamma_med`) brings it to about 0.20 (decision 3b,
2026-11-11). Those fixes need ERA5, while PBLHmod must be computable from EC data alone (G, 2026-10-30). This spec
estimates γ from the station's own sensible heat flux and tower air temperature, evaluates it, and prepares decision
D1: is the EC-derived γ the recommended EC-only γ (and possibly the default)?

Principle. In the zero-order slab model behind B&G (Batchvarova & Gryning 1991; Kljun et al. 2015 Eq. B5) the mixed
layer warms as dθ/dt = (1 + A) w'θ'_s / h (entrainment flux = −A × surface flux) while it grows into a layer of
gradient γ. The rise of θ at the tower, given the heat put in, therefore carries γ. A γ fitted with the same H that
drives the model also absorbs a systematic error of H (e.g. energy-balance non-closure): report this, do not correct
it (no H post-correction, G 2026-11-12).

## 1. Daily γ_EC (`R/gamma_ec.R`, `scripts/17a_gamma_ec.R`)

Inputs per site: the step-02 forcing (`prepare_pblh_forcing()`: `wtheta`, `T_K`, `T_source`, `p_Pa`, `ustar`, `L`,
`is_day`, flags) and the `default` PBLHmod output (the day's `h0`, convective window, flags). **No ERA5.**

- θ_obs = T_K (100000 / p_Pa)^(r_d / cp) [K] (dry θ, constants from `pblh_constants()`, as spec 11 A). T from the
  model's T source (`fn_TA_F_K` → `met_TA_K`); rows with `T_source = T_SONIC` are not used (count them).
- **Fit window** per day: the B&G window of the `default` run (first half-hour with `wtheta > 0` to the window end),
  cut at the half-hour of the day's maximum θ_obs if that comes earlier, and at local solar noon +
  `gamma_ec_end_offset_h` (default 3 h). At least `gamma_ec_min_n` valid half-hours (default 6); no forcing gap inside.
- **Method A (primary, model-consistent).** For a candidate γ, integrate h(t) with `bg_step()` from the day's `h0` with
  the day's forcing (exactly as `run_pblh_model()` in the window: same sub-steps, same rules for negative spells), and
  the mixed-layer temperature θ_mod(t) = θ0 + ∫ (1 + A) wtheta / h dt (same sub-steps; A = `bg_A`). Minimise
  Σ (θ_obs − θ_mod)² over θ0 (analytic: mean residual) and γ (Brent on log γ in [`gamma_min`, `gamma_max`] =
  [0.001, 0.02], the spec-11 bounds). Optimum at a bound → NA (flag). Report γ_A, θ0, RMSE, R², n.
- **Method B (closed form, sensitivity).** Encroachment with h0 = 0 gives h² = 2 (1 + 2A) Q / γ and
  θ − θ_s = γ (1 + A) / (1 + 2A) · h, with Q = ∫ wtheta dt from the window start [K m]. So θ_obs = a + s √Q, and
  γ_B = s² (1 + 2A) / (2 (1 + A)²) (factor 0.97 for A = 0.2: almost independent of A). Linear least squares; report
  s, R². Derivation in the code comment (Tennekes 1973; Batchvarova & Gryning 1991).
- **Day selection** (flags, nothing deleted; counts per flag in the summary): clear sky, mean `SW_IN / SW_IN_POT`
  over the window ≥ `gamma_ec_clear_min` (default 0.6; FLUXNET variables; where missing, flag `clear_unknown` and do
  not use the day); no `flag_no_h0`, `flag_bg_stopped_gap`, `flag_no_positive_h`; θ rise over the window ≥ 1 K;
  R²_A ≥ `gamma_ec_r2_min` (default 0.8); γ_A inside the bounds. A day passing all is `used`.
- **Climatology**: site × calendar month median, p25, p75 and n of γ_A over `used` days; months with fewer than
  `gamma_ec_min_days` (default 5) days are NA. Same layout as the ERA5 monthly table of script 11a, so the model reads
  either.
- Outputs: `{derived_dir}/gamma_ec/{site}_daily.csv` (date, window start/end, n, γ_A, γ_B, R²_A, R²_B, RMSE_A,
  θ rise [K], Q [K m], `h0`, flags, `used`) and `{site}_monthly.csv`. Parameters in `config.yml` (`gamma_ec_*`),
  documented.

## 2. PBLHmod option (amendment 15 to spec 02)

- `gamma_source` gains `"ec_clim"` (site × month γ_EC climatology) and `"ec_daily"` (the day's γ_A where `used`, else
  the climatology, flagged `ec_clim_fallback`). A month without a value uses the constant (`constant_fallback`,
  warning); no network fallback (as spec 11, decision 2). `gamma_scope` as for ERA5.
- Cases in step 02: `gamma_ec` (`ec_clim`), standard source `pblhmod_gamma_ec` (label "PBLHmod-γEC"), written at every
  site; `gamma_ec_daily` (tables only, no standard file).
- Order: `17a` needs the `default` output and step 02 needs the 17a table for the new cases: two passes, as
  `gamma_med` (step 02 skips the `gamma_ec*` cases while the table is missing; then `17a`; then 02 again). New
  `run_site.R` step `17a`. Propose the cheapest order and say which existing outputs it recomputes (every `R/` change
  reruns every requested step: report the cost before the cluster run).
- With `gamma_source = "constant"` the output is unchanged (< 1e-9 m) at BE-Bra, FI-Hyy, FR-Gri, CH-Dav.

## 3. Evaluation

**A. Measured sites** (BE-Bra ERA5 only, and `config/measured_run_sites.txt`): `scripts/11_gamma_evaluation.R` with
cases `default`, `clim`, `gamma_med`, `gamma_site`, `gamma_ec`, `gamma_ec_daily` (`--cases=`); same sets, metrics and
bootstrap as spec 11 amendment 1; outputs with prefix `17_` (no change to the existing `11_` tables). Clear-day growth
curves (`17_growth_clearday_{site}`) with PBLHmod, PBLHmod-clim, PBLHmod-γEC, ERA5, measured.

**B. γ benchmark: what the measurement asks for.** At each measured site, constant-γ cases on a grid 0.002, 0.003, …,
0.012 K m⁻¹ (step 02, tables only, comparison period only if that is simpler; no standard files). Median ln ratio in
the primary set against γ; γ_best = the γ where it crosses 0 (linear interpolation; NA if it does not cross), with a
day-block bootstrap interval if cheap. Figure `17_gamma_benchmark`: per site γ_best (with interval), the site median
γ_EC (p25–p75), the ERA5 site median (`gamma_site`), `gamma_med` and 0.01.

**C. Network** (39 sites): γ_EC monthly climatology vs the ERA5 climatology (11a): scatter and seasonal cycle per site
(`17_gamma_ec_vs_era5`), share of months with a value, used days per site; PBLHmod-γEC vs ERA5 tower cell (new pair in
step 04, `9_pblhmod_gamma_ec_vs_era5_center`, all strata as pairs 1–8).

**D. Diagnostic only, no correction:** per site, γ_EC / γ_ERA5 (median of the monthly ratios) against the site's
daytime energy-balance ratio (H + LE) / (NETRAD − G) on the reference year (L2 FLUXNET variables; report availability
per site; summer daytime sums). This is G's open question of decision 3b: does the residual spread follow the EBR?

**E. Footprint bracketing check (for D1):** on the step-08 rows of each site's reference year (`common_pe`, `etc`,
`valid_qc_used`), the share of half-hours where h of `gamma_ec` lies outside [min, max] of (ERA5, PBLHmod,
PBLHmod-clim), and the median distance outside in ln units, per site × stratum (zm/L class × day/night). No FFP runs.

## 4. Decision D1 (G, after the results)

Proposed criteria, to discuss: at the measured sites the median |ln ratio| of `gamma_ec` is ≤ those of `clim` and
`gamma_med` and clearly below `default`; no site clearly worse than `default`; γ_best and the γ_EC site median agree
within their intervals at most sites; the bracketing share is small (e.g. < 5 % of rows per site). Options: `gamma_ec`
replaces `gamma_med` as the recommended EC-only setting, and/or becomes the code default. A default change reruns
every step (spec 11 amendment 1, "cost"); **the footprints are not rerun** if the bracketing check holds (G,
2026-11-25); the paper then reports the footprint results for the three sources and states that PBLHmod-γEC lies
inside their range.

## 5. Computation

PC: BE-Bra, FI-Hyy, FR-Gri, CH-Dav (local inputs), all parts except C and E. Cluster: 17a and step 02 at the 39 sites
(two passes), script 11 (A, B), step 04 pair 9, 17b (C, D, E: `scripts/17b_gamma_ec_collect.R`). Runbook
`docs/local/cluster_commands_spec17.md` in the usual style (`~/pbl_vars.sh`; job names not used before; check the
`#SBATCH --array` line before every `sbatch`; `--time` on the command line; no dependencies on ended jobs; checks and
what to bring back). Report the run-time estimate first.

## 6. Tests

θ conversion (hand case). Method B recovers γ on a synthetic encroachment day (θ and h generated with known γ, A,
h0 = 0) within 1 %. Method A recovers γ on a synthetic day generated with `bg_step()` (mechanical terms on, h0 > 0)
within 2 %, and with θ noise of 0.1 K the median bias over 200 synthetic days is < 5 %. Bound optimum → NA. Window
rules (max θ, noon + offset, minimum n, gaps). Selection flags. Monthly table readable by `run_pblh_model()`;
`ec_clim` / `ec_daily` fallbacks; `constant` identity.

## 7. Docs

`docs/STATUS.md` (stand-by reopened for spec 17), spec 02 amendment 15, `CLAUDE.md` (γ sources, step 17a, order),
`docs/project_plan.md` (γ option "EC-derived"), `docs/figure_inventory.md`. Methods note for the paper: both methods,
the selection, the bounds, the A-independence of Method B, the H-bias absorption. References: Batchvarova & Gryning
1991, Tennekes 1973, Kljun et al. 2015 (check the keys in the bib; add nothing unchecked).

## Acceptance criteria

Tests pass; identity with `constant` at four sites; daily and monthly γ_EC tables and parts A, B, D at the four local
sites (labelled "test"); runbook written; run-time and rerun-cost estimate reported before the cluster run.

## Summary to report back

As in `CLAUDE.md`, plus: used days per site and the flag counts; γ_A vs γ_B (median ratio); the local γ_EC climatology
next to the ERA5 one and 0.01; the local part-A metrics by case; anything in the forcing that limits the fit (tower
temperature source, SW_IN_POT availability, early-morning shapes).
