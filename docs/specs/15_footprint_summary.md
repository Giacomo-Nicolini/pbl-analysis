# Spec 15: network footprint summary for the paper (09_collect), and the γ supplementary tables

Branch: `spec-15-footprint-summary` (from `main`). Agreed with G on 2026-11-12. Light work on this PC (collection
logic and figures from tables copied back); the collection runs on the cluster.
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/figure_inventory.md`, `docs/specs/08_footprint_sensitivity.md`
(outputs; Rebmann2018 criteria, lines ~126 and ~215), `docs/specs/09_cluster_deployment.md` (collection),
`docs/specs/11_gamma_era5_davos.md` (decision 3b recorded in spec 09).

## Why

The footprint impact is the core of the paper. All per-site footprint outputs exist (reference year per site,
cases `pblhmod` (γ 0.01), `pblhmod_clim`, `era5`, plus `measured` and `bound` where available). G's scope
(2026-11-12): the five elements below and nothing else. No NEE weighting, no budgets, no maps by regime (maps stay
annual / seasonal / day-night, from the stored groups).

## Part A: network tables (in `09_collect.R` or a helper it calls)

Case pairs: `era5` − `pblhmod`, `pblhmod_clim` − `pblhmod`, `era5` − `pblhmod_clim` (difference = first minus
second). Validity set `etc` as main; `kljun_norsl` as sensitivity. Subsets `valid` and `valid_qc_used`. Row set
`common_pe` (rows valid in every case compared).
1. **Target-area fraction difference [percentage points]:** per site and network-wide (pooled rows, and median of
   site medians); p10 / median / p90; by stability regime (zm/L classes as in spec 08) × day/night, and all.
2. **Fetch:** ln(x_80 ratio) between cases, same strata and statistics; also x_80 in metres per case (median, p90).
3. **Most affected rows:** the 10 % of rows with the largest |target-area difference| per site (and network): their
   share by regime, day/night, season; and what removes them: share failing `valid_qc` and `valid_qc_used`, per
   validity set (expected network figures from the per-site tables: about 71 % stable, 64 % night, 64 % (`etc`) /
   93 % (`kljun`) removed by QC / used; reproduce and report).
4. **ICOS site-design criteria (Rebmann et al. 2018):** per site and case: share of time with ≥ 90 % of the flux from
   the target area (all, and stable only), and whether the ≥ 50 %-of-time criterion is met; share of rows with
   cumulative footprint ≥ 80 %; a flag where the pass/fail outcome differs between cases (expected: none). State the
   criteria and page numbers in the table header.
5. **Exceptions:** list sites where |p10| or |p90| of the target-area difference exceeds 5 pp in any regime ×
   day/night stratum (expected: CH-Dav and IT-Ren, stable / night, p10 −16 to −20 pp); for them, the stratum, n,
   p10 / median / p90, which case gives the smaller target-area fraction, and median h per case in that stratum.
   The 5 pp threshold is a reporting choice; state it.

Also the ICOS labelling criterion (G, 2026-11-12): per site and case, among the half-hours of the "good data"
subset, the share with a cumulative target-area contribution ≥ 70 %; target met if that share is ≥ 70 %. Report it
next to the Rebmann2018 values (columns `icos_lab_share`, `icos_lab_pass`), computed on `valid_qc`, and flag sites
where pass/fail differs between cases. State the rule verbatim in the table header as "ICOS labelling rule (ETC)".

Outputs: `09_network_fp_target_area.csv`, `09_network_fp_fetch.csv`, `09_network_fp_topaffected.csv`,
`09_network_fp_icos_criteria.csv`, `09_network_fp_exceptions.csv`. Each with a header comment (definitions,
sign convention, units).

## Part B: figures (paper candidates; style as the existing figures)

- `09_fp_ta_diff`: target-area difference (p10–p90 bar, median point) per site, ordered by median h or by
  land cover, facets regime × day/night; network summary row at the bottom; exceptions highlighted.
- `09_fp_fetch_ratio`: same layout for ln(x_80 ratio).
- `09_fp_topaffected`: stacked shares (regime, day/night, removed by QC / used) for the top 10 %.
- `09_fp_icos_criteria`: per site, the share of time meeting the 90 % criterion per case (points per case, 50 % line).
- Supplement: `09s_fp_*` versions with `kljun_norsl` and with `valid` instead of `valid_qc_used`.

## Part C: γ supplementary tables

- `09s_gamma_site.csv`: per site, `gamma_site` (median) with IQR, p10 / p90, n days and months used; land cover.
- `09s_gamma_network.csv`: `gamma_med` (median of site medians) with a site-bootstrap 95 % interval; the pooled
  median; the range of site medians.
- Optional figure `09s_gamma_site`: site medians with IQR, `gamma_med` line, γ = 0.01 line.

## Tests

Sign convention of the differences; the 10 % selection on a synthetic set; the Rebmann criteria on a synthetic site
(known pass/fail); the exceptions rule at the threshold; the γ bootstrap reproducible with a seed.

## On this PC / cluster runbook / acceptance

Run Part A–C from the cluster tables already copied back if they hold enough (list what's missing); otherwise tests
only. Runbook `docs/local/cluster_commands_spec15.md` (the collection only, no footprint rerun). Acceptance: tests
pass; the five tables and four figures (plus supplement versions); γ tables; `docs/figure_inventory.md` updated;
numbers in the report checked against the per-site tables.

## Implementation (2026-11-12)

Code: `R/footprint_network.R` (row table, the five tables, ICOS criteria, γ tables, table headers),
`R/plot_footprint_network.R` (figures), `R/plot_utils.R` (style groups `fp_pair`, `fp_share`, elements
`ta90_line`, `lab70_line`), `scripts/09_collect.R` (`--no-footprint` skips the part); tests
`tests/testthat/test-footprint_network.R`.

- **Data source.** The per-site tables of step 08 hold neither the regime × day/night strata nor the
  `pblhmod_clim` − `era5` pair, and pooled network quantiles need rows: `fp_site_rows()` reads each site's
  footprint input table (`<site>_footprint_inputs.rds`, step 08a: zm/L class, day / night, season, PBLHmod regime,
  validity flags) and the per-half-hour chunk metrics of the reference year (`chunks/*_hh.csv.gz`, step 08b: h,
  target-area fraction, x_80, `ffp_ok`). No footprint is recomputed. So the collection reads more than the per-site
  tables for this part (a few minutes on the cluster).
- **`common_pe`.** Rows valid (set × subset) with an FFP result in all three cases, so the three pairs use the same
  rows. The step-08 per-site tables used rows common to each pair (slightly larger sets); differences are small.
- **Strata.** zm/L class (`zl_class()`) × day / night with the margins (`all`); scopes `site`, `network_pooled`
  (quantiles of all rows of all sites), `network_site_median` (median over the sites of the site p10, median, p90).
- **Top 10 %.** Selected on the `common_pe` rows of subset `valid` (so that the share removed by QC / used can be
  measured), per site and pooled; shares by zm/L class and by PBLHmod regime; removed = not in `common_pe` of
  `valid_qc` / `valid_qc_used`. For comparison, the table also aggregates the step-08 per-site tables
  (`selection = "abs_lnx80_era5_pblhmod (step 08)"`, top 10 % by |ln x_80 ratio| ERA5 vs PBLHmod, PBLHmod regime).
  In the figure, the bar of the network median of the site shares need not sum to 100 %.
- **ICOS criteria.** Rebmann on `common_pe` of the given set × subset (all rows; stable = zm/L class); ICOS
  labelling rule on `common_pe` of (set, `valid_qc`); the pass flags (≥ 50 %, ≥ 70 % share) include the boundary.
- **Exceptions.** Strictly above 5 pp (|p10| or |p90|), crossed strata only, main set and subset; no minimum n
  (n reported).
- **Figures.** Sites ordered by land cover (IGBP, `config/sites.csv`), then median PBLHmod h of the main rows;
  exceptions marked with " *"; network rows at the bottom. The top-affected figure shows the ERA5 − PBLHmod pair.
- **γ.** `gamma_site_table()`: quantiles of the valid daily ERA5 γ per site, year-months and calendar months used;
  `gamma_network_table()`: 2000 site resamples, seed 20261112, percentile interval.

### Checks against the per-site tables (cluster collection of 2026-11-06, `data/cluster_results1/`)

- Top 10 % (step-08 selection, `etc`): median of the site shares 71.2 % stable (PBLHmod regime), 64.5 % night,
  64.1 % removed by QC or used, as expected. **Removed by QC or used: 92.6 % is `kljun_norsl` (and `none`); `kljun`
  gives 83.9 % (17 sites; `kljun` keeps no row at forest sites).** The expected values are medians of the site
  shares; pooled (weighted by n) 72.0 % stable.
- Exceptions preview from `08_target_area` (`etc`, `valid_qc_used`, pairwise common rows, strata not crossed):
  CH-Dav and IT-Ren only, ERA5 − PBLHmod: stable p10 −16.4 / −16.0 pp, night p10 −19.9 / −15.6 pp; neutral −5.3 /
  −9.9 pp; ERA5 gives the smaller target-area fraction.
- Rebmann (`etc`, `valid_qc_used`, all rows): 17 of 39 sites pass with PBLHmod and with ERA5 (stable rows: 9); no
  verdict differs in this set × subset; 20 differ elsewhere (other subsets, night or stable rows), all with shares
  within 0.05 of the 50 % limit.
- ERA5 − PBLHmod, all rows: median of the site medians +0.02 pp (site p10 / p90 medians −0.17 / +0.18 pp); PBLHmod
  target-area fraction: median of the site medians 89 %.

## Amendment 1: fix-15 figures (G 2026-11-13, cluster run of 2026-11-13)

- Titles of `09_fp_*` / `09s_fp_*` were clipped on the right: subtitle shortened to the set, subset and reference
  year; the details (site order, network bar, exception mark) moved to a caption; both wrapped (`fp_wrap()`).
- `09_fp_ta_diff`, `09_fp_fetch_ratio` (and `09s_`): 39 sites + network per facet row overlapped: height from
  `fp_fig_height()` (0.42 cm per site row and facet row, + 7 cm), site labels 6.5 pt (`theme_pblh(y_text_size =)`).
- `09_fp_icos_criteria` (and `09s_`): PBLHmod was hidden under PBLHmod-clim and ERA5: shapes from the style table
  (new column `shape`: PBLHmod filled circle, PBLHmod-clim open circle, ERA5 open triangle, the network figure
  convention) and a vertical dodge within each site row (PBLHmod on top); exception sites marked " *" here too.
- Every figure with " *" explains it in a caption line (`fp_exception_caption()`): |p10| or |p90| of the
  target-area fraction difference > 5 pp in at least one zm/L class × day/night stratum, for any case pair (set
  `etc`, subset `valid_qc_used`; `09_network_fp_exceptions`).
- `09_collect.R` "no chunks" message: a case-month absent from `<derived_dir>/footprint/chunk_list.csv` has no valid
  rows and its missing chunks are expected (info line; e.g. SE-Myc 2025-07); a listed case-month without a chunk on
  disk is a real gap (danger line, counted; rerun with `08c --pending`); without the chunk list the cause is reported
  as unknown (`fp_classify_missing()`, `fp_missing_messages()`).
