# Methods note: where the PBLH source matters for EC fluxes (spec 18; Results 3.5)

Draft for the methods / supplement text (2026-11-27). Code: `R/pblh_ec_fluxes.R`, `R/pblh_ec_collect.R`,
`scripts/18a_pblh_ec_site.R`, `scripts/18b_pblh_ec_collect.R`. Analysis only: no PBLHmod or footprint run.

## Rows and sources

Every half-hour of a station's L2 record from its ICOS labelling date (the first valid PBLHmod value, as
`09_network_record_start`) to the end of the FLUXES file, with the exclusion periods of `config/exclusion_periods.csv`
removed (counted per site in `18_variables`). The reference year of spec 15 is kept as a flag. Heights: zm = EC height
− 2/3 canopy height from the time-varying site metadata. Night = `SW_IN_POT` = 0; stability class from zm/L as spec 15
(unstable < −0.1 ≤ neutral ≤ 0.1 < stable, L = L2 `MO_LENGTH`); season from the UTC month.

h: ERA5 BLH of the tower cell (instantaneous, interpolated linearly to the half-hour mid, as the footprint cases),
PBLHmod (constant γ), PBLHmod-clim. Measured: the site's primary source. The Cloudnet ε MLH (FI-Hyy) is valid at night
and enters the surface-layer shares with its exact values; censored values (below the lowest gate, signal-limited
lower bounds) are counted as their own classes and never imputed. Aerosol MLH (SIRTA, MeteoSwiss) sits on a night-time
floor that is not a stable-layer depth: it enters only the daytime Panofsky test.

**Circularity.** PBLHmod's night-time h is computed from u* and L (spec 02, B1). With PBLHmod, "zm/h > 0.1" is close
to "low u*", and its overlap with the u* filter is partly built in. ERA5 is independent of the station's turbulence and
is the main source of sections 1–3; PBLHmod and PBLHmod-clim are reported next to it.

## Surface-layer criterion

The tower is taken to be above the surface layer when zm/h > 0.1 (the surface layer as the lowest ≈ 10 % of the
boundary layer); 0.05 and 0.2 as sensitivity. Shares and the median and p90 of zm/h per source × day/night × zm/L class
× season, on the rows with h of all three model sources (common) and on each source's own rows.

## u* filter

The ONEFlux VUT reference u* threshold of each (local) year, from the BIF file (`GRP_FLUX_REF`, as spec 19).
u*-filtered = FLUXNET `USTAR` below that threshold (main definition); `NEE_VUT_REF_QC` > 0 is reported alongside (it
also contains the other NEE removals and gaps), with the agreement of the two on night rows. The 2 × 2 table u*-filtered
× above the surface layer (night rows with all three h) gives the share of above-surface-layer half-hours removed by
the filter and the share of accepted half-hours above the surface layer.

## Residual test

Night, u*-passed, measured (`NEE_VUT_REF_QC` = 0) half-hours with `RECO_NT_VUT_REF` and `TA_F`. Residual r =
`NEE_VUT_REF` − `RECO_NT_VUT_REF` [µmol m⁻² s⁻¹] and r / RECO_NT (rows with RECO_NT > 0.5 µmol m⁻² s⁻¹). In each cell
u* / u*_thr bin (1–1.5, 1.5–2, > 2) × `TA_F` tercile, the difference of the medians above (zm/h > 0.1) and within the
surface layer, valid with at least 30 rows on each side; pooled per u* bin and overall as the mean of the valid cell
differences weighted by the smaller side's n. 90 % interval: night-block bootstrap (whole nights resampled with
replacement, 1000 replicates). The same with the L2 FC instead of NEE at sites with a profile storage term, to separate
the turbulent flux from the storage term. **Dilution caveat:** RECO_NT is fitted on the same night-time data, so a bias
of the above-surface-layer half-hours enters the fitted respiration too; the test can dilute a bias, not create one.

## Storage

ICOS L2 NEE = FC + SC (checked: |NEE − FC − SC| ≤ 0.011 µmol m⁻² s⁻¹ on all rows, values rounded to two decimals).
`SC_FLAG` 0 = profile storage (rows with flag 1 have SC equal to the single-point estimate `SC_1P`; sites without
`CO2_DRY_<n>` profile levels have flag 1 only). Only profile rows are used; sites without a profile have no storage
results. Night accumulation = sum of positive SC × 1800 s over a night (sunset to sunrise, at least 80 % of the rows
valid) [µmol m⁻²], against the night median h; morning release = sum of negative SC × 1800 s over the 4 h from the
first day row (at least 6 of 8 valid), against the growth of h over the same rows (PBLHmod: its B&G growth; ERA5: the
interpolated hourly values). Spearman ρ per site and binned medians. Share of night-time NEE that is storage: median
SC / NEE (L2 NEE with QC 0) on night rows with |NEE| > 0.5 µmol m⁻² s⁻¹, by zm/h class.

## Panofsky test

In the unstable surface layer σv / u* = (12 − 0.5 h / L)^(1/3) (Panofsky et al. 1977, Boundary-Layer Meteorol. 11,
355–361; **not yet in `docs/references/*.bib`: G to check and add**). Rows: day, 0.05 ≤ −zm/L ≤ 2, valid L2
`V_SIGMA` and `USTAR`, zm/h < 0.1 for the source tested. Per site and source: median ln(observed / predicted), median
|ln ratio|, and the OLS slope of ln(σv/u*) on ln(12 − 0.5 h/L) (1/3 expected). A site-specific offset (terrain,
meandering) is expected, so the sources are ranked within each site on common rows by the median absolute deviation
of the ln ratio from its own median (offset removed); how often each source ranks first is counted over the sites.

## Spectral correction factors

Descriptive only: EddyPro total factors `H_SCF_STAT`, `LE_SCF_STAT`, `FC_SCF_STAT` (`TAU_SCF_STAT` is not in the L2
files checked) against zm/h and −h/L on unstable rows, binned medians per site and over sites. The low-frequency part
of these factors is analytic and contains no h: any relation reflects their shared dependence on wind speed, zm and L,
not a measured loss.
