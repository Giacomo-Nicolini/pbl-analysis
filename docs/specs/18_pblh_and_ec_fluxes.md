# Spec 18: where the PBLH source matters for EC fluxes — surface layer, u* filter, storage, σv

Extension after `analysis-v1` (Results 3.5 of the outline). Branch: `spec-18-pblh-ec-fluxes` (from `main`). Agreed with
G on 2026-11-25 (plan) and 2026-11-27 (go). Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/specs/15_footprint_summary.md`,
`docs/specs/19_alpine_exception.md` (descriptors, row basis, u* thresholds from the BIF), `docs/specs/02_pblh_model.md`
(night-time PBLHmod formula), `docs/specs/standard_format.md`. Analysis only: **no PBLHmod or FFP runs**, no change to
existing outputs. The γEC cases of spec 17 are not used (D1 rejected).

## Why

The footprint barely depends on h (specs 15–16), but h decides other things in EC processing. Monin–Obukhov similarity,
the u*-filter logic and the footprint model assume the sensor is in the surface layer, roughly zm < 0.1 h. On stable
nights h can fall to 100–200 m, so a tall-forest tower (zm 15–30 m) can sit above the surface layer, where flux
divergence and decoupling bias night-time NEE. There the PBLH source matters: at night ERA5 and PBLHmod differ by
factors of 2–5 (spec 19: at CH-Dav and IT-Ren zm/h_ERA5 is 0.15–0.20 on neutral and stable nights). This spec asks how
often towers are above the surface layer, by source; whether the u* filter already removes those half-hours; whether
night-time NEE that passes the filter is biased there; how storage relates to h; and which h best explains σv.

**Circularity, stated once and kept in every result.** PBLHmod's night-time h is computed from u* (and L) (spec 02,
B1), so with PBLHmod "zm/h > 0.1" is close to "low u*", and its overlap with the u* filter is partly built in. ERA5 is
independent of the station's turbulence and is the **main source for §1–3**; PBLHmod and PBLHmod-clim are reported
next to it; at FI-Hyy the Cloudnet ε MLH (valid at night) is a third, measured source.

## 0. Rows and variables (all sites, `config/run_sites.txt`)

- **Period:** every year from the station's labelling date with L2 data (the full record of the standard files); the
  reference year of spec 15 as a flag column (`ref_year`), so results can be checked against the footprint rows.
- **h:** standard files `era5_center`, `pblhmod_default`, `pblhmod_clim`; FI-Hyy also `cloudnet_eps_mlh` (exact values;
  censored ones as their own class, never imputed).
- **Fluxes and turbulence:** FLUXNET HH (`NEE_VUT_REF`, `NEE_VUT_REF_QC`, `RECO_NT_VUT_REF`, `TA_F`, `SW_IN_POT`,
  `NIGHT`, `USTAR`) and L2 FLUXES (`FC`, `NEE`, `V_SIGMA`, `W_SIGMA`, `MO_LENGTH`, `WS`, `H_SCF_STAT`, `LE_SCF_STAT`,
  `FC_SCF_STAT`, `TAU_SCF_STAT`, the FC/NEE QC flags). Report per site which variables exist (a site table); a missing
  variable skips only the part that needs it.
- **Storage:** ICOS L2 NEE = FC + SC where a CO2 profile exists, so SC = NEE − FC on rows with both (report per site
  whether a profile exists: `CO2_DRY_1..n` columns, `MISSING_LEV_*`); sites without a profile have no storage results.
  If a direct SC column exists, use it and check SC ≈ NEE − FC.
- **u* filter:** the ONEFlux VUT reference threshold per year (BIF, as spec 19); u*-filtered = `NEE_VUT_REF_QC > 0`
  or `USTAR` below the threshold (report both and their agreement).
- **Heights:** zm = EC height − 2/3 hc from `config/sites_timevarying.csv` (time-varying, as everywhere).
- **Strata:** night (`SW_IN_POT = 0`) / day; zm/L class as spec 15; season (DJF, MAM, JJA, SON).

## 1. How often is the tower above the surface layer?

Per site × source × day/night × zm/L class × season: share of half-hours with zm/h > 0.1, and with zm/h > 0.05 and
> 0.2 (sensitivity), and the median and p90 of zm/h. Rows: all half-hours with h of all three sources (common rows),
and separately each source on its own rows. Network summary: median over sites, and the list of sites where the share
on stable nights exceeds 25 % for any source. Output `18_sl_share.csv`.

## 2. Does the u* filter already remove them?

On night rows with all three h: a 2 × 2 table per site and source, u*-filtered (yes/no) × above the surface layer
(zm/h > 0.1, yes/no), with counts and shares; the share of above-surface-layer half-hours that the u* filter removes,
and the share of accepted (u*-passed) half-hours that are above the surface layer. ERA5 first. Output
`18_ustar_cross.csv`.

## 3. Is accepted night-time NEE biased above the surface layer?

Rows: night, u*-passed, measured (not gap-filled: `NEE_VUT_REF_QC == 0`), with `RECO_NT_VUT_REF`. Residual
r = NEE_VUT_REF − RECO_NT_VUT_REF [µmol m⁻² s⁻¹] (night NEE minus the night-time respiration model of ONEFlux) and its
relative form r / RECO_NT_VUT_REF. Compare the residual between zm/h ≤ 0.1 and > 0.1 (ERA5 h; then PBLHmod), within
u* bins (u* / u*_thr: 1–1.5, 1.5–2, > 2) and TA_F terciles, so the comparison is not just a u* or temperature effect:
median difference per site with a night-block bootstrap 90 % interval (resampling whole nights), n per cell; a cell
needs ≥ 30 rows on each side. Repeat with FC instead of NEE (FC − RECO_NT) at sites with storage, to see whether the
difference is in the turbulent flux or in the storage term. State clearly that RECO_NT is fitted on the same night
data, so a bias would be diluted, not created. Output `18_nee_residual.csv`; figure `18_nee_residual` (per site:
median difference with interval, by u* bin, sites ordered by their stable-night share of §1).

## 4. Storage against h (sites with a CO2 profile)

- Night accumulation: sum of positive SC over each night (sunset to sunrise, in µmol m⁻²), against the night median h
  (ERA5, PBLHmod); daily points, Spearman ρ per site, and the binned median curve.
- Morning release: sum of negative SC from sunrise to sunrise + 4 h, against the morning growth of h (PBLHmod: the B&G
  growth over the same hours; ERA5: hourly difference), Spearman ρ per site.
- Share of night-time NEE that is storage: median of SC / NEE on night rows with |NEE| > 0.5 µmol m⁻² s⁻¹, by zm/h class.
Output `18_storage_h.csv`; figure `18_storage_h` (accumulation vs night h, facets by site, log x).

## 5. Which h best explains σv? (Panofsky test)

In unstable conditions in the surface layer, σv / u* = (12 − 0.5 h / L)^(1/3) (Panofsky et al. 1977). Rows: day,
−zm/L between 0.05 and 2, zm/h < 0.1 for the source tested, valid V_SIGMA and USTAR. Per site and source (ERA5,
PBLHmod, PBLHmod-clim, and the measured MLH at the measured sites by day): median ln(observed / predicted), the median
absolute ln ratio, and the slope of ln(σv/u*) on ln(12 − 0.5 h/L) (expected 1/3). A site-specific offset (terrain,
meandering) is expected, so rank the sources **within** each site by the median absolute ln ratio after removing the
site's median offset, and report how often each source ranks first. Output `18_panofsky.csv`; figure `18_panofsky`
(binned σv/u* against −h/L per source for 6 sites: FR-Gri, CH-Dav, FI-Hyy, BE-Bra, DE-Tha, IT-Ren; the Panofsky curve).

## 6. Spectral correction factors (descriptive only)

`H_SCF_STAT`, `LE_SCF_STAT`, `FC_SCF_STAT`, `TAU_SCF_STAT` (total factors, EddyPro) against zm/h and −h/L in unstable
conditions: binned medians per site and network, no model. Caption and summary must say that the low-frequency part of
these factors is analytic and contains no h, so any relation reflects their shared dependence on wind speed, zm and L,
not a measured loss. Output `18_scf_h.csv`; figure `18s_scf_h` (supplement).

## 7. Outputs and figures

Tables under `{output_dir}/tables/` with a `#` header (definitions, units, rows, sources); figures as named above plus
`18_sl_share` (stable-night share above the surface layer per site, three sources side by side, sites ordered by zm).
Shared theme; update `docs/figure_inventory.md`. Every figure caption names the h source and the circularity note
where PBLHmod is shown.

## 8. Computation

Scripts: `scripts/18a_pblh_ec_site.R` (per site: rows, §1–6 site tables; restartable) and
`scripts/18b_pblh_ec_collect.R` (network tables, figures, `--from-tables=DIR`). The `R/` changes of specs 17–19 change
every run_site hash, so **call the scripts directly, never through `run_site.R`**. PC: develop and test on the local
sites (BE-Bra, FI-Hyy, FR-Gri, CH-Dav), with timing. Cluster: one array over the 39 sites + one collection job;
runbook `docs/local/cluster_commands_spec18.md` in the usual style (`~/pbl_vars.sh`; new job names; check the
`#SBATCH --array` line; `--time` on the command line; no dependency on ended jobs; checks; bring back). Report the
cost estimate before the cluster run.

## 9. Tests

zm/h shares on a synthetic series; the 2 × 2 table; night grouping (sunset to sunrise across midnight, polar day and
night handled); night-block bootstrap; SC = NEE − FC and the profile detection; the Panofsky prediction on hand values
and the within-site ranking; u*-filter agreement (QC flag vs threshold) on a hand case; skipping a part when a variable
is missing.

## 10. Docs

`docs/STATUS.md`, `CLAUDE.md` (scripts 18a–b), `docs/project_plan.md` (Results 3.5), `docs/figure_inventory.md`.
Methods note: the surface-layer criterion and its sensitivity, the u* filter definition, the residual test and its
dilution caveat, the storage definition, Panofsky et al. (1977) (check the key in the bib; add nothing unchecked).

## Acceptance criteria

Tests pass; local tables and figures for the four local sites (labelled "test"); runbook and cost estimate; nothing
existing overwritten.

## Summary to report back

As in `CLAUDE.md`, plus: variable availability per site (storage profile, SCF, RECO_NT, V_SIGMA); the local §1–5
results in one line each; anything that limits a part (e.g. too few u*-passed night rows above the surface layer).
