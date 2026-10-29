# Spec 11: free-atmosphere lapse rate γ from ERA5 profiles; Davos MLH; evaluation at four sites

Phase 1 (PBLHmod option) and phase 3. Branch: `spec-11-gamma-era5` (from `main`, after the spec-10 merge).
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/project_plan.md` (phase 1: γ options 1–3; "Scope"; evaluation and
compression rules), `docs/specs/02_pblh_model.md` (B&G), `docs/specs/04_era5_standard_compare.md` (pressure-level
reader), `docs/specs/10_measured_aerosol_mlh.md`. Agreed with G on 2026-10-29. Light work only on this PC.

## Why

At FR-Gri PBLHmod is ≈ 40 % below the SIRTA ceilometer MLH in the daytime window (median ln ratio −0.53), while
ERA5 matches it (+0.002). With γ = 0.0025 K m⁻¹ the deficit disappears (+0.08); start height and gap stops are ruled
out (spec 10). The default γ = 0.01 is about twice the ERA5 850–700 hPa gradient at BE-Bra, FI-Hyy and FR-Gri.
B&G growth scales roughly with γ^(−1/2), so a realistic γ should remove most of the bias. This spec implements γ
options 2 and 3 of the project plan and evaluates them; **the default stays 0.01 until G decides** (end of spec).
The decision is needed before the cluster footprint run (spec 09), which uses PBLHmod.

## Part A — γ from ERA5 pressure levels (`R/era5_gamma.R`)

Inputs: the spec-04 pressure-level reader (6-hourly T and geopotential on 13 levels, tower cell and neighbours),
static orography, site altitude from `config/sites.csv`.

- Per profile: θ = T (1000/p)^(R_d/c_p) (dry θ, as in the model; constants from `pblh_constants()`); height of each
  level above mean sea level = Φ / g; above ground = that − site altitude (also report against the ERA5 cell
  orography).
- γ of a profile = slope of a least-squares fit of θ on height over the levels with height above ground in
  [`gamma_zmin_agl`, `gamma_zmax_agl`] (defaults 500 m and 3000 m). **Levels below the ground are excluded** (Davos is
  at ≈ 840 hPa: 850 hPa is underground). At least 2 levels (report how often only 2 are available); report R² and
  the levels used.
- Profile time: the 6-hourly profile closest to and not after the local sunrise of the day (the start of the B&G
  integration), i.e. the morning profile above the nocturnal layer; option `gamma_profile_hour` to force a UTC hour.
- Plausibility: γ outside [`gamma_min`, `gamma_max`] (defaults 0.001 and 0.02 K m⁻¹) → NA for that day (report).
- Outputs: `{derived_dir}/gamma/{site}_daily.csv` (date, profile time, γ, n levels, R², levels used, flag) and
  `{site}_monthly.csv` (site × calendar month: median, p25, p75 of daily γ, n days) — the climatology.

## Part B — PBLHmod option `gamma_source`

- `gamma_source = "constant"` (default, uses `gamma`), `"era5_clim"` (site × month climatology), `"era5_daily"`
  (the day's morning γ; falls back to the climatology when NA, flagged `flag_gamma_fallback`).
- New output columns: `gamma_used` [K m⁻¹], `gamma_source_used` (`constant`/`era5_daily`/`era5_clim_fallback`).
- Record as an amendment to spec 02. With `gamma_source = "constant"` the output must be identical to now
  (< 1e-9 m) at BE-Bra, FI-Hyy and FR-Gri.
- New cases in `scripts/02_run_pblh_model.R`: `gamma_era5clim`, `gamma_era5daily`; plus a `gamma_era5clim_zmin1000`
  sensitivity (fit layer from 1000 m agl).

## Part C — Davos MLH (MeteoSwiss STRATfinder)

Input: `{measured_dir}/CH-Dav/MeteoSwiss_STRATfinder_MLH_ABLH_DAV_20210701_20261027.nc` (from MeteoSwiss via the CH-Dav
PI, 2026-10-29) and three illustrative PNGs (not used).
Facts checked by Cowork (2026-10-29; verify): CF-1.8, `time` seconds since 1970 **UTC**, ≈ 60 s; 2021-07-01 →
2026-10-27 (2,731,834 profiles); `mlh`, `ablh` [m above ground at the ceilometer; sensor at ground], `mlh_asl`,
`ablh_asl`; `quality_level_mlh` (0 good, 1 suspect, 2 bad, 9 missing; **use 0 only**, as the provider recommends);
`qc_flag_mlh` bitmask (8 = at lowest detectable height ≤ 30 m, informational; 16 = at upper search limit 4000 m;
4 = maintenance/dirty window, incl. 2025-06-14 → 2025-11-27); `precipitation_flag`; `cloud_base_height`;
`homogeneous_segment`; station 46.8130 N, 9.8435 E, 1590 m (CH-Dav tower 0.98 km away, 1637 m); STRATfinder
(`Kotthaus2020`); licence **CC BY-NC 4.0**; acknowledgement: "Please cite MeteoSwiss as the data provider and Kotthaus
et al. (2020)". Night MLH sits on a floor (median 260–320 m all seasons, almost no spread: night search 0–700 m and
incomplete overlap), so night values are not stable-layer depths.

- Reader `read_meteoswiss_stratfinder()` and standard source `meteoswiss_stratfinder_mlh` (`method_class = aerosol`,
  half-hourly, `time_reference = interval_mean`): median of the quality-0 minutes in each half-hour; `n_samples`;
  `qc_flag` 0 if ≥ 75 % of 30 minutes are quality 0, 1 if ≥ 50 %, else NA; `x_ablh_m` (median quality-0 ABLH),
  `x_share_lowest` (share of bit 8), `x_segment` (homogeneous segment). Sidecar: licence, acknowledgement, citation,
  product version, instrument history summary.
- `config/measured_sites.csv`: row CH-Dav | MeteoSwiss | DAV | coordinates | distance | aerosol | CHM15k | CC BY-NC 4.0.
- CH-Dav inputs: G copies from the cluster the L2 zips (FLUXES, FLUXNET HH, METEO), the BIF and ERA5 (three
  families) into `data/raw/`. Run 01, 02 (all cases incl. the new γ cases), 04 for CH-Dav; report site specifics
  (alpine valley, slope, ERA5 cell orography vs 1637 m, zm vs hc).

## Part D — evaluation

Sites and references: FR-Gri (SIRTA CHM15k, window), FI-Hyy (Cloudnet ε MLH, exact + censored tails), CH-Dav
(MeteoSwiss STRATfinder, window), BE-Bra (ERA5 only). Use the generic spec-10 comparison.

- γ diagnostics: daily and monthly γ per site (figure `11_gamma_by_month_{sites}`), against 0.01 and 0.0025; n
  levels and R²; daily vs climatology spread.
- For each site and each γ case (`default`, `gamma0025`, `gamma_era5clim`, `gamma_era5daily`, `zmin1000`): window
  metrics vs the measurement (median ln ratio, within factor 2, FB, nRMSD, Spearman, symmetric slope), the FI-Hyy
  censored-tail test, and PBLHmod vs ERA5 (convective and all regimes). Table `11_gamma_cases_metrics.csv` (long).
- Clear-day growth curves by solar hour (as spec 10, 3b) for each case at FR-Gri, FI-Hyy, CH-Dav.
- Figures: `11_gamma_cases_{site}` (median ln ratio by case with bootstrap intervals), `11_growth_clearday_{site}`.

## Part E — plotting convention (G, 2026-10-29)

In every comparison figure, hours outside the comparable set (distance class, aerosol window, regime) are **drawn but
faded** (alpha ≈ 0.3, or a light grey band behind them); the comparable window at full opacity; a legend entry
"not comparable (shown for context)". Metrics unchanged. Put it in a shared helper in `R/plot_utils.R`; regenerate
the spec-07 and spec-10 figures with it (tables must stay byte-identical; report).

## Docs

- `docs/project_plan.md`: γ options 2–3 implemented; the plotting convention; **data release (phase 6): measured
  PBLH data are never released, not even derived from them as standard files; the release cites DOIs, references and
  acknowledgements case by case (G, 2026-10-29)**; Davos as a measured source.
- `CLAUDE.md`: γ inputs and options; MeteoSwiss file facts; `docs/specs/standard_format.md`: a 1-min → half-hour
  aggregated aerosol example; `docs/STATUS.md`; spec 02 amendment.
- Bib: nothing new expected (Kotthaus2020 present); add any reference you rely on only after checking it.

## Tests

θ and height computation (hand case); levels below ground excluded; fit with 2 and 3 levels; plausibility bounds →
NA; profile selection at sunrise (sites at different longitudes and seasons); `gamma_source` fallback and identity
with `constant`; MeteoSwiss aggregation on a synthetic file (quality levels, bit 8, coverage → qc); faded-plot helper
(comparable flag respected).

## Acceptance criteria

Tests pass; `constant` reproduces current PBLHmod at three sites; γ tables for BE-Bra, FI-Hyy, FR-Gri, CH-Dav;
Davos standard file validated; evaluation tables and figures; regenerated 07/10 figures with identical tables.

## Out of scope

Footprints; cluster runs; γ from sources other than ERA5; changes to B&G other than γ; night use of aerosol heights.

## Summary to report back

As in `CLAUDE.md`, plus: γ climatology per site (vs 0.01), how often the fit has only 2 levels, the case metrics at
each site (which γ option removes the bias, and does it hold across sites and seasons?), CH-Dav specifics, and a
recommendation for the default γ, with the evidence — G decides.
