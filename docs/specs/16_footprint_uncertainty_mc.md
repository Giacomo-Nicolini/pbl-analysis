# Spec 16: footprint uncertainty from the FFP inputs (phase 5b), and where the PBLH source sits in it

Phase 5b. Branch: `spec-16-footprint-uncertainty` (from `main`). Agreed with G on 2026-11-17.
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/specs/08_footprint_sensitivity.md` (inputs, validity sets, FFP
settings, amendments 1–2), `docs/specs/09_cluster_deployment.md` (chunks, hashes, reference years),
`docs/specs/15_footprint_summary.md` (target-area metrics, ICOS criteria, exceptions), `docs/figure_inventory.md`.
Light work only on this PC (a smoke test on the local test month); the run is on the cluster.

## Why

Spec 15 shows that the PBLH source barely changes the footprint: the target-area (TA) fraction differs by a median
+0.02 pp between ERA5 and PBLHmod across 39 sites, except on stable and neutral nights at the two alpine sites
(p10 −15 to −20 pp). A reviewer will ask: small compared with what? The footprint already carries uncertainty from
its other inputs (σv, L, u*, umean, zm through d, wind direction). This spec quantifies that uncertainty on a
stratified sample of half-hours at six sites, ranks the inputs, and places the PBLH-source effect in it. Expected
use in the paper: one paragraph and one figure in Results 3.2, details in the supplement, and the quantitative basis
for the discussion point "why footprints are robust to h".

Not in scope: footprint climatologies by Monte Carlo, network-wide Monte Carlo, flux or budget uncertainty.

## 1. Sites and half-hours

**Sites (6):** FR-Gri (cropland, measured aerosol MLH), FR-Fon (deciduous forest, measured), FI-Hyy (boreal forest,
turbulence-based MLH), BE-Bra (mixed forest, low TA shares), CH-Dav and IT-Ren (alpine, the two spec-15 exceptions).
Each site's approved reference year (`config/reference_years.csv`). Configurable list, `config/footprint_mc_sites.txt`.

**Half-hours:** rows valid under the `etc` set in the base run (all three PBLH cases valid, i.e. spec-15 `common_pe`),
subset `valid`; `valid_qc_used` reported as a sub-selection (no extra runs). Stratified sample per site:
stability (zm/L class: unstable < −0.1 ≤ neutral ≤ 0.1 < stable) × day/night = 6 strata, 170 rows each (all rows if
a stratum has fewer; report n). Seeded (`set.seed` per site, seed in the outputs). Stratum weights = the stratum's
share of the site's valid rows in the reference year, for any site-level or network statistic (the sample is not
representative of time without them).

## 2. Input uncertainties (G, 2026-11-17)

All perturbations are relative to the published L2 value of that row (the base). One standard deviation (1σ):

| Input | Perturbation | 1σ | Note |
|---|---|---|---|
| σv (`V_SIGMA`) | multiplicative, lognormal | 15 % | sensitivity: 10 % and 20 % |
| u* (`USTAR`) | multiplicative, lognormal | 10 % | |
| L (`MO_LENGTH`) | via 1/L: L' = L0 · (u*'/u*0)³ · ε_L, ε_L lognormal | 20 % | keeps the sign; L ∝ u*³ links L to the u* draw (see note) |
| umean (`WS`) | multiplicative, lognormal | 5 % | |
| hc (canopy height) | multiplicative, normal, truncated at ±3σ | 10 % | changes d = 2/3 hc and zm = EC height − d, and z0 = 0.15 hc |
| z0 roughness factor | multiplicative on 0.15, lognormal | 50 % | z0 is used only by the validity tests (FFP runs with umean, z0 = NaN); see section 3 |
| wind direction (`WD`) | additive, normal, wrapped to [0, 360) | 7.5° | affects the TA fraction only |
| h (PBLH) | multiplicative, lognormal | 30 % | applied to PBLHmod; alternative "empirical" (below) |

**L and u*:** in the one-at-a-time design (section 3A) L is perturbed alone (ε_L only, u* fixed). In the joint draws
(3B) L follows the u* draw through (u*'/u*0)³ and then gets its own ε_L. Report the zm/L shift this produces.

**h, two options (both computed):** `h_pct30` = lognormal 30 % around PBLHmod (main); `h_empirical` = lognormal
with σ_ln taken from the spread of ln(PBLHmod/measured) in the daytime window at the aerosol sites (spec 11 cluster
tables: the median over sites of the IQR / 1.349 of the daily values), applied by day only, and 30 % at night
(no night reference). Report the σ_ln used. This keeps the original phase-5b idea of error distributions from the
measured sites.

**Reference: the PBLH-source effect.** For every sampled row, the spec-15 differences ERA5 − PBLHmod and
PBLHmod-clim − PBLHmod of the TA fraction and of ln x80 (from the existing chunk `_hh` files, no new runs).

## 3. Designs

Same FFP code, wrapper and settings as step 08b (Kljun `calc_footprint_FFP.R` v1.3, umean only, z0 = NaN, domain
±1000 m, dx = dy = 1 m, `rslayer = 1`, `smooth_data = 1`; the site's target-area mask as in 08t). Never modify or
commit the Kljun code. Outputs per row and run: TA fraction [%], x_peak, x50, x70, x80, x90 [m], FFP `flag_err`.

**Validity under perturbation:** the base row is valid (section 1). A perturbed run is kept only if the row stays
valid under the `etc` set recomputed with the perturbed inputs (zm, z0, L, u*, h) and the FFP code's own checks.
Count and report the share of perturbed runs that become invalid, by input and stratum: a validity flip is itself an
uncertainty, and the z0 factor enters only here. The ETC z0 multiplier m stays the site's value from 08a (not
recomputed per draw).

**A. One-at-a-time (all sampled rows):** each input at −1σ and +1σ (quantiles 0.1587 / 0.8413 of its distribution;
wind direction ±7.5°), the others at base: 2 × 8 = 16 runs per row, plus the base run (checked against the spec-15
value of the same row, tolerance 0.01 pp). Per input: ΔTA [pp] and Δ ln x80 for −σ and +σ, and the half-range
(|Δ+| + |Δ−|) / 2 as that input's 1σ effect.

**B. Joint Monte Carlo (subsample):** 50 rows per stratum (300 per site), 50 joint draws per row (all inputs drawn
together, L linked to u* as above), h option `h_pct30` and, in a second set, `h_empirical`. Per row: p05, p16, p50,
p84, p95 of the TA fraction and of ln x80; the 1σ-equivalent width (p84 − p16)/2.

**C. Attribution (from B):** per site and stratum, the share of the variance of the TA fraction and of ln x80
explained by each input: squared standardised regression coefficients of the response on the standardised
log-perturbations (wind direction as its additive deviation), with the R² of that linear model reported (low R² =
non-linear response; then say so and rely on A). No Sobol indices (cost).

## 4. Comparison with the PBLH-source effect

Per site × stratum (and network: stratum-weighted median over the six sites):
- the 1σ effect of each input (A) vs the median |ΔTA| and |Δ ln x80| of ERA5 − PBLHmod and PBLHmod-clim − PBLHmod;
- the ratio |PBLH-source effect| / joint 1σ width (B), and the share of rows where the source effect exceeds the
  joint 1σ and 2σ widths. Expected (to test, not to assume): well below 1 almost everywhere, above 1 on alpine nights.
- ICOS criteria under uncertainty (B, stratum-weighted): per site and case, the probability that the share of rows
  with TA ≥ 90 % (Rebmann 2018) and the labelling share (TA ≥ 70 % of good data) cross their thresholds; flag sites
  where the pass/fail verdict is uncertain (probability between 0.1 and 0.9). The labelling share uses `valid_qc`
  rows of the sample.

## 5. Outputs

Tables under `{output_dir}/tables/` (each with a `#` header: definitions, units, sign convention, seed):
- `16_oat_{site}.csv`: row × input × side (−σ, +σ): ΔTA, Δ ln x80, valid flag, stratum, subset flags;
- `16_mc_{site}.csv`: row-level quantiles of B, per h option;
- `16_attribution_{site}.csv`: C per stratum;
- `16_network_uncertainty.csv`: section 4 summaries per site, stratum and network; `16_icos_uncertainty.csv`.

Figures (PNG, shared theme):
- `16_input_ranking`: per stratum (facets regime × day/night), the 1σ effect of each input on the TA fraction
  (median over rows, IQR bar), inputs ordered by size, with the PBLH-source effect (ERA5 − PBLHmod) as a marked
  reference line; sites as small points, network as the large point. Paper candidate (Results 3.2).
- `16_input_ranking_x80`: same for ln x80.
- `16s_mc_width`: joint 1σ width vs |PBLH-source effect| per row (log-log, by stratum, 1:1 line).
- `16s_validity_flips`: share of perturbed runs that become invalid, by input and stratum.

Update `docs/figure_inventory.md`.

## 6. Computation

- Scripts: `scripts/16a_mc_sample.R` (strata, sample, perturbation table; no FFP), `scripts/16b_mc_run.R` (chunks:
  site × design × part, restartable by inputs hash, as 08b), `scripts/16c_mc_collect.R` (tables, figures).
  `slurm_submit.R --what=mc` renders the array; the chunk list goes to `{derived_dir}/footprint_mc/chunk_list.csv`.
- **Step 0, runtime first:** time 100 rows × 17 runs (A) at BE-Bra on the PC (≤ 2 workers) and extrapolate. Expected
  order: A 6 × 1,020 rows × 17 runs; B 6 × 300 × 50 × 2 h options; at about 2 s per 2D run that is about
  60 + 100 CPU-hours. Report the estimate before the cluster run; if above 300 CPU-hours, propose a smaller B
  (fewer draws or rows) and ask G.
- Perturbation draws are generated once in 16a and stored (`{derived_dir}/footprint_mc/{site}_draws.rds`), so the run
  is reproducible and chunks can be rerun alone.

## 7. Tests

Lognormal draws: median 1 and the requested σ (large n); ±1σ quantiles; wind wrap at 0/360; L–u* link (L' / L0 =
(u*'/u*0)³ when ε_L = 1); stratified sampling and weights (weights sum to 1 per site; strata with few rows); validity
recomputation under perturbation (hand cases at the zm/L, u* and roughness boundaries); base run equals the spec-15
row value; attribution on a synthetic linear response recovers the known shares.

## 8. On this PC

Smoke test on the local test month (BE-Bra, FI-Hyy, 2022-09): 20 rows per stratum, design A only, plus 5 rows × 10
draws of B; check the base-run equality with the local 08b values. Runtime step 0 as above. Nothing else locally.

## 9. Cluster runbook

`docs/local/cluster_commands_spec16.md`, in the usual style: `~/pbl_vars.sh` variables; job names not defined there;
check the `#SBATCH --array` line against the expected range before every `sbatch`; no dependencies on jobs that may
have ended; checks and what to bring back.

## 10. Docs

Spec file as is (this file). `docs/STATUS.md` entry; `docs/project_plan.md` phase 5b (no longer postponed; the design
above); `CLAUDE.md` (scripts 16a–c). Methods note for the paper: the input uncertainties of section 2 and their
sources (to be completed by G: EC random-error references for σv and u*, the hc reliability from the BIF).

## Acceptance criteria

Tests pass; smoke test and runtime estimate on the PC; base runs equal the spec-15 values; tables and figures defined
in section 5 produced by the smoke test (test labels); runbook written; the estimate reported to G before the run.

## Amendment 1 (G, 2026-11-18): hc and z0 without the roughness-sublayer test

The smoke test showed that at FI-Hyy zm/z0 is only 2.5 % above the 08a ETC multiplier, so every +1σ hc and z0 run
fails the `etc` roughness-sublayer test (e_rsl, zm ≥ m z0) and the 1σ effects of hc and z0 are `NA` there (no
valid +1σ run).

- Design A also runs hc and z0 at −1σ / +1σ with e_rsl switched off, all other `etc` conditions kept (base checks,
  e_el, e_tower, e_stab, e_ustar; m of 08a unchanged): validity `rsl_off`, runs `<input>_<minus|plus>_rsl_off`, 21
  runs per row. Their FFP inputs equal those of the `etc` runs (z0: the base run), so they add an FFP call only where
  the `etc` hc run fails e_rsl.
- Reported next to the main runs: `16_oat_<site>` gains `validity`, `rsl_fail`, `valid_run`; `16_network_uncertainty`
  gains the items `hc|rsl_off`, `z0f|rsl_off` and the column `share_rsl_flip` (share of the input's ±1σ runs that
  fail e_rsl); the 16c log prints both effects per site.
- `16_input_ranking` (and `16_input_ranking_x80`): hc and z0 effects from `rsl_off`, labelled "rsl test off", with the
  network share of rsl-flipped runs written beside them ("rsl flips n %"), so a flip never reads as a zero effect.
  Stated in the figure caption and in the methods note. `16s_validity_flips` stays on the `etc` runs.
  **Superseded (G 2026-11-19, close-out):** the flip shares are no longer written beside the bars; the figure caption
  replaces that wording (it names the site × input cells whose runs fail the test and points to
  `16s_validity_flips`, which shows the share per site). See "Close-out figure fixes" below.

## Close-out figure fixes (G, 2026-11-19; re-rendered from the cluster tables, no cluster run)

- `16_input_ranking`, `16_input_ranking_x80`: one legend below the panels (inputs by colour, the PBLH-source
  references as lines; reference lines in neutral greys so they do not share a colour with an input); title and
  subtitle say "reference year per site", the caption lists the sites with their years.
- The network median flip share ("rsl flips n %") is dropped from the main figure (it read 0 % while FI-Hyy and CH-Dav
  lose half of their hc / z0 runs); the caption names the site × input cells with failing runs and points to
  `16s_validity_flips`, which now shows the share per site (one panel per site, runs at −1σ and +1σ together,
  stratum-weighted with the 16a row weights, the e_rsl share in brackets). This replaces the "beside the bars" wording
  of amendment 1.
- `16_oat_<site>` gains `w_row` (16a row weight). `16c --from-tables=DIR [--fp-table=FILE]` redraws the three figures
  from collected tables; for tables without `w_row` the weights come from the eligible rows per stratum of the spec-15
  `09_network_fp_target_area` (set `etc`, subset `valid` = the 16a eligible set; checked: every stratum sampled below
  170 rows has n equal to the spec-15 count).
- Cluster result (reference years): rsl failures under `etc` at FI-Hyy (hc and z0, every +1σ run: 50 % of the runs)
  and CH-Dav (z0, every +1σ run: 50 %); elsewhere < 0.5 %; u* 1–4 % (e_ustar) at every site.

## Summary to report back

As in `CLAUDE.md`, plus: the runtime estimate and the proposed chunking; the smoke-test input ranking at BE-Bra and
FI-Hyy (labelled "test, one month"); the share of validity flips by input; anything in the FFP code that limits a
perturbation (e.g. inputs outside its range).
