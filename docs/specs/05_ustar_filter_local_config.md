# Spec 05: u* plausibility filter in the PBLHmod forcing; local config override

Phase 1 (amends spec 02) and housekeeping. Branch: `spec-05-ustar-filter` (from `main`).
Read first: `CLAUDE.md`, `docs/specs/02_pblh_model.md` (section 3, forcing), `docs/project_plan.md` (open point on
implausible USTAR). Agreed with G on 2026-10-21; implement as written.

## Part A — u* plausibility filter (input QC, not output clipping)

### Why

At BE-Bra, 179 PBLHmod half-hours exceed 5 km (70 exceed 10 km), across B1, B3 and B5. They come from implausible
L2 `USTAR` (median 6.9 m s-1 on those rows, max 49 m s-1; 1,726 rows > 2 m s-1). These values are consistent with
L2 `TAU` (u* from TAU ≈ 2.1 m s-1 on rows with USTAR > 2), so they are EddyPro outputs that passed QC, but they are
physically implausible:
- u*/WS: median 0.21, p90 0.30 over all rows; neutral log law at BE-Bra gives u*/U = k / ln((z − d)/z0) ≈ 0.14–0.19
  (z − d ≈ 24 m, z0 ≈ 1–3 m); the suspect rows have WS ≈ 4 m s-1 and u*/WS > 0.5;
- σw/u* on rows with u* > 2 m s-1: 0.87, vs 1.21 overall (≈ 1.25 near neutral, e.g. Panofsky & Dutton 1984;
  Kaimal & Finnigan 1994).

The model output stays raw (no clipping of `pblh_m`); implausible **inputs** are rejected before the model runs.

### Rule (new step in `prepare_pblh_forcing()`, before the short-gap interpolation)

A half-hour's u* is implausible if either holds:
1. **Ratio test**: `ustar / WS > opts$ustar_ws_max` and `WS >= opts$ws_min_ratio_test` (WS = FLUXES `WS`, same sonic,
   [m s-1]); not applied where WS is missing or below `ws_min_ratio_test` (ratio unstable at low wind);
2. **Absolute cap**: `ustar > opts$ustar_max`.

On implausible rows: `ustar` **and `L`** are set to `NA` (L ∝ u*³, so it is unusable too); `flag_ustar_implausible
= TRUE`. The existing short-gap interpolation (≤ `max_gap_interp`, spec 02) then runs as usual; filled rows keep
both `flag_ustar_implausible` and the `flag_*_interp` flags. `H` and the other inputs are not touched.

### Options (`pblh_options()`)

| option | default | unit | meaning |
|---|---|---|---|
| `ustar_ws_max` | 0.5 | – | max u*/WS; `Inf` disables the ratio test |
| `ws_min_ratio_test` | 1 | m s-1 | ratio test only where WS ≥ this |
| `ustar_max` | 3 | m s-1 | absolute cap; `Inf` disables |

Validate (positive numbers or `Inf`). Document the physical reasoning and references in roxygen and in comments next
to the rule. The defaults change PBLHmod results: record this as **amendment 12** in `docs/specs/02_pblh_model.md`.

### Output and standard format

- New output column of `run_pblh_model()`: `flag_ustar_implausible` (logical). Also add `WS` [m s-1] to the forcing
  and model output (needed to audit the rule).
- `standardize_pblhmod()` (spec 04): `qc_flag` meaning unchanged; rows that were rejected and then interpolated are
  already qc 1 through `flag_ustar_interp`; add `flag_ustar_implausible` to the qc 1 condition anyway, and update
  `qc_flag_meaning`.

### Script `scripts/02_run_pblh_model.R`

- Add two sensitivity cases: `ustar04` (`ustar_ws_max = 0.4`) and `ustar06` (`ustar_ws_max = 0.6`); keep the
  others (`default`, `gamma0025`, `max_B2`, `max_auto`), all with the new defaults. Add a case `nofilter`
  (`ustar_ws_max = Inf, ustar_max = Inf`) reproducing the previous default exactly (test: identical to the spec-02
  output on BE-Bra, < 1e-9 m).
- Report per case: rows rejected by the ratio test, by the cap, by both; share of rows; how many were then
  interpolated vs left `NA`; count of `pblh_m > 5 km` and `> 10 km` before (nofilter) and after; quantiles per regime
  (as now).
- Then rerun `scripts/04_era5_compare.R` (it reads the derived files) and report pair 1 (default vs ERA5 tower cell)
  overall and by regime, before and after the filter.

### Tests

- Ratio test: hand cases (u* = 2.5, WS = 4 → rejected; u* = 0.8, WS = 4 → kept; u* = 0.8, WS = 0.5 → kept because
  WS < 1; WS `NA` → only the cap applies).
- Cap: u* = 3.5 with WS = 20 → rejected.
- L set to `NA` on rejected rows; H untouched; flags set; rejected rows inside a short gap get interpolated and keep
  both flags; long runs stay `NA`.
- `ustar_ws_max = Inf, ustar_max = Inf` gives output identical to before.
- Option validation.

## Part B — local config override (no cluster paths in the public repo)

The repository is public. Absolute cluster paths (with user names) must not be in tracked files from now on.

- New file `config.local.yml` (gitignored; add to `.gitignore`), same structure as `config.yml`
  (`default:` / `cluster:` sections, any subset of keys). Create it locally with the current cluster paths moved out
  of `config.yml` (do not commit it).
- New tracked file `config.local.yml.example`: same keys as the cluster profile, placeholder values, comments on how
  to use it.
- `config.yml`: cluster-profile values that are machine-specific become `"LOCAL"` (keep the relative ones, e.g.
  `sites_file`). Header comment: machine-specific paths go in `config.local.yml`.
- `get_cfg(profile)`: read `config.yml` as now; if `config.local.yml` exists, read it with `yaml::read_yaml()`, take
  its `default` section and then its active-profile section, and merge them over the base configuration
  (`utils::modifyList`); report once (message) which keys came from the local file. The active profile is resolved
  as now (`profile` argument, else `R_CONFIG_ACTIVE`, else `default`).
- `resolve_path()`: a value `"LOCAL"` aborts with a message saying to set it in `config.local.yml` (same handling as
  `"TODO"`).
- Remove the absolute cluster paths from tracked docs: `CLAUDE.md` (ERA5 section) and
  `docs/specs/04_era5_standard_compare.md` (section 1 table): replace them with "cluster path: set in
  `config.local.yml` (see `config.local.yml.example`)". Update the "Environment" section of `CLAUDE.md`
  accordingly. Do not rewrite git history (already-published paths stay in history; agreed with G).
- Tests: merge order (local overrides base; profile section overrides local default), missing local file = base
  only, `"LOCAL"` aborts in `resolve_path()`, unknown keys in the local file allowed but reported.
- Check before committing: `git grep -n "/bigdata"` and `git grep -n "nicolini\|trotta"` return nothing outside
  `config.local.yml.example` placeholders (which must not contain real paths).

## Acceptance criteria

- All tests pass; `scripts/02_run_pblh_model.R` and `scripts/04_era5_compare.R` run on BE-Bra; `nofilter`
  reproduces the old default exactly.
- No absolute cluster paths or user names in tracked files; `config.local.yml` not tracked.
- Amendment 12 in spec 02; `docs/STATUS.md` updated (defaults line, spec 05 row); open point on USTAR in
  `docs/project_plan.md` marked resolved with the counts.

## Out of scope

Other input plausibility rules (H, L sign consistency, T, p) — to be discussed separately if needed; changes to the
legacy comparison (spec 03).

## Summary to report back

As in `CLAUDE.md`, plus: rejected-row counts per criterion and case; outliers > 5 / > 10 km before and after;
pair-1 metrics before and after; anything unexpected (e.g. rows rejected in clearly plausible conditions).
