# Spec 19: why the PBLH source changes the footprint at CH-Dav and IT-Ren, and not at other mountain sites

Extension after `analysis-v1`. Branch: `spec-19-alpine-exception` (from `main`). Agreed with G on 2026-11-25.
Read first: `CLAUDE.md`, `docs/STATUS.md` (stand-by section), `docs/specs/15_footprint_summary.md` (exceptions, min
n = 50), `docs/specs/16_footprint_uncertainty_mc.md`, `docs/specs/11_gamma_era5_davos.md` (Part C: Davos facts).
Analysis and figures only: **no PBLHmod or FFP runs**, no change to existing outputs.

## Why

Spec 15 flags as exceptions (|p10| or |p90| of the target-area (TA) difference > 5 pp, n ≥ 50) the neutral and stable
nights and the stable days at CH-Dav and IT-Ren (p10 −15 to −20 pp at night) and, mildly, the stable nights at DE-Tha
(−7.5 pp). Reviewers will ask why these sites and not other mountain sites (CZ-BK1, IT-MBo, FR-Lqu). A first look at
the frozen tables (Cowork, 2026-11-25) suggests three ingredients together:
1. relief ERA5 does not resolve: sub-grid orography SD ≈ 370–390 m at CH-Dav and IT-Ren (≤ 124 m at the other forest
   sites); the towers sit 360 m below (CH-Dav) and 380 m above (IT-Ren) the ERA5 cell surface;
2. a large night-time disagreement of the sources: stable-regime median PBLHmod / ERA5 ≈ 3.7 (CH-Dav) and 4.8
   (IT-Ren), ≤ 1.5 elsewhere, where ERA5 is usually the higher one;
3. a tall tower (zm a sizeable fraction of a shallow h) and enough valid stable and neutral nights (CZ-BK1: 32
   stable-night rows in its reference year, below the minimum; IT-MBo: alpine relief but a 2.4 m grassland tower).
This spec tests that explanation on the rows actually used for the footprints and documents the physical causes.

## 1. Site descriptors (all 39 sites) — `scripts/19a_site_descriptors.R`

One row per site, on the reference-year rows of step 08 (`common_pe`, set `etc`, subset `valid_qc_used`; also
`valid` as a column suffix), by day/night and zm/L class where marked [s]:
- terrain (from `09_network_era5_cells.csv`): site altitude, ERA5 cell orography, difference, sub-grid orography SD
  (`sdor`, max over the cells), orography range of the 9 cells;
- tower: EC height, hc, zm (reference-year median), land cover (BIF/IGBP);
- PBLH [s]: median h of ERA5, PBLHmod, PBLHmod-clim; median ln(PBLHmod / ERA5) on the same rows; zm / h_ERA5 and
  zm / h_PBLHmod (median, p90);
- turbulence [s]: u*, σw / u*, L (median, p10, p90), share of rows with u* below the site's L2 u* threshold if
  available in the FLUXNET files (VUT reference; report availability);
- night-time flows: circular mean wind direction by night and by day, their difference (reversal), directional
  steadiness (|vector mean| / scalar mean speed) by night and by day, share of night rows within ±45° of the dominant
  night direction, median night wind speed;
- footprint [s]: n, p10, median, p90 of the TA difference ERA5 − PBLHmod (spec-15 table), the spec-15 status
  (`exception`, `below_min_n`, none); valid stable / neutral night rows.
Output `{output_dir}/tables/19_site_descriptors.csv` (long or wide, with a `#` header: definitions, units, rows used).

## 2. Row level (selected sites) — same script

Sites: CH-Dav, IT-Ren, DE-Tha, CZ-BK1, FR-Fon, BE-Bra, IT-MBo (configurable). Night rows (neutral and stable) of the
reference year: TA difference (ERA5 − PBLHmod) against ln(h_ERA5 / h_PBLHmod) and against zm / h_ERA5; and the x80
ln ratio likewise. Report, per site, how much of the TA-difference variance the two predictors explain (simple linear
model on the night rows, R²; no further modelling). Output `19_row_level.csv` (site × stratum summaries, not rows).

## 3. Which source is closer at night? (existing tables only)

- FI-Hyy (Cloudnet ε, valid at night): ERA5 and PBLHmod against the lidar in the night / stable strata (from the
  spec-10 / spec-11 tables; if a night stratum is missing, say so; do not rerun).
- CH-Dav: daytime ERA5 vs MeteoSwiss STRATfinder bias compared with the other aerosol sites (spec-10 tables), as
  evidence of how ERA5 handles the valley; the night STRATfinder values sit on a floor and are not usable (spec 11 C).
- State plainly that no measured night-time PBLH exists at the alpine sites.

## 4. Figures (supplement)

- `19s_alpine_drivers`: x = night (neutral + stable) median ln(PBLHmod / ERA5), y = p10 of the TA difference on
  stable nights (second panel: neutral nights); point colour = EC height (or zm), size = sub-grid orography SD;
  labels for all sites with |p10| > 2 pp and for CZ-BK1, IT-MBo, FR-Lqu; the −5 pp line; `below_min_n` sites as open
  symbols.
- `19s_row_level`: the night-row scatter of §2 for the selected sites (facets), coloured by zm / h_ERA5.
- `19s_night_flows`: day vs night wind-direction distributions (polar histograms or direction × hour heat maps) for
  CH-Dav, IT-Ren, DE-Tha, CZ-BK1, IT-MBo, BE-Bra.
Shared theme; update `docs/figure_inventory.md`.

## 5. What to report (the causes)

A short numbered list of findings with numbers, one per ingredient (terrain, source disagreement and its turbulence
side, tower height, data availability, night flows), and for each mountain site that is **not** an exception, which
ingredient is missing. No interpretation beyond the numbers: Cowork and G write the paper text and supply the
literature (ERA5 boundary-layer height in complex terrain and at night, valley and drainage flows); add no reference
that is not checked in the bib.

## 6. Computation

The reference-year rows with h per case, u*, L, σw, WD, WS, zm, hc and the TA results are on the cluster (step-08 /
chunk outputs and the L2 files): one light job (`docs/local/cluster_commands_spec19.md`, usual style). Develop and
test locally on the sites with local inputs (BE-Bra, FI-Hyy, FR-Gri, CH-Dav); the figures can be redrawn from the
collected tables (`--from-tables=DIR`).

## 7. Tests

Circular statistics (mean direction, steadiness, reversal) on hand cases including the 0/360 wrap; row selection
equals the spec-15 `common_pe` rows of a test site; descriptor table layout; the plot helpers on a synthetic table.

## Acceptance criteria

Tests pass; local descriptors for the four local sites; runbook written; tables and figures from the cluster
collection; the findings list of §5.

## Summary to report back

As in `CLAUDE.md`, plus the findings list of §5 and anything that contradicts the three-ingredient explanation.
