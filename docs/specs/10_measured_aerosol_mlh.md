# Spec 10: aerosol-based measured MLH — SIRTA ReOBS at FR-Gri (Lanzhot CZ-Lnz when the data arrive)

Phases 2–3. Branch: `spec-10-measured-aerosol` (from `main`, after the spec-08 merge).
Read first: `CLAUDE.md`, `docs/project_plan.md` ("Scope", "Comparing different method classes", distance rule,
censoring convention), `docs/specs/standard_format.md`, `docs/specs/07_measured_pblh_cloudnet.md` (comparison code
and amendment 1), `docs/references/literature_review_pblh.md` (section 3). Agreed with G on 2026-10-27.

## Goal

Add the second measurement class, **aerosol-based** mixing-layer heights from ceilometers, using products that are
already computed (no retrieval of our own: G, 2026-10-27), and compare PBLHmod and ERA5 with them.

- Part A (now): SIRTA ReOBS hourly MLH at Palaiseau, ≈ 23 km from FR-Gri → daytime convective comparisons only.
- Part B (when the data arrive): Lanzhot BL-View L3 MLH at CZ-Lnz.
- Part C (optional, only with G's go): Cloudnet ε MLH at Palaiseau (2025-01 →), the spec-07 chain, giving a
  turbulence-based and an aerosol-based height at the same place.

## Part A — SIRTA ReOBS (FR-Gri)

### A1. Inputs

- Measured: `data/raw/measured/FR-Gri/Palaiseau_SIRTA/MLH4CK.zip` (from SIRTA via the FR-Gri PI, 2026-10). Contains
  `reobs_sirta_1D_request_2026_08_1_20101214_000000_15y.nc` and four overview PNGs. Read the NetCDF by extracting it
  to a temporary directory at read time (`ncdf4` cannot read inside a zip); never unzip permanently into the repo.
- ICOS L2 FR-Gri (FLUXES, FLUXNET HH, METEO, FLUXNET BIF) and ERA5 FR-Gri (hourly single levels, pressure levels,
  static fields): G copies them from the cluster into the usual `data/raw/` folders. If anything is missing, stop
  and list it.

### A2. Facts already checked (Cowork, 2026-10-27; verify and report differences)

- SIRTA 1D ReOBS dataset (AERIS; `Chiriaco2018`), product version 2026.08.1, ReOBS 4.0.0; hourly values
  2010-12-14 00:00 → 2026-06-30 23:00 (136,272 hours). `time`: minutes since 1970-01-01, **UTC**, start of the hour;
  `time_bnds` = [t, t + 1 h). Station 48.718° N, 2.207° E, 156 m.
- Variables used (heights above ground, m; `_FillValue` NaN for values, −999 for counts):
  - `mlh_chm15k` (+ `_std`, `_nb`): MLH from the Lufft CHM15k with **STRATfinder** (high-SNR ceilometers;
    `Kotthaus2020`); 2015 →; `_nb` up to 60 (1-min values per hour);
  - `mlh_cl31` (+ `_std`, `_nb`): MLH from the Vaisala CL31 with **CABAM** v2019.3 (`Kotthaus2018`,
    `Kotthaus2020`); 2011 →; `_nb` up to 4 (15-min values per hour);
  - `cbh_chm15k`, `cbh_cl31` (cloud base), `tcf_cl31` (total cloud fraction, oktas): ancillary.
- `ablh_chm15k` ("atmosphere boundary layer thickness") is the top of the aerosol layer (≈ 1.4 km on summer nights,
  i.e. residual layer), **not** a mixing height: not used as a reference (may be kept as `x_` column).
- Licence: "Inspired by CC-BY" (file attribute); acknowledgement text in the global attribute `acknowledgement`.
- Timing check: JJA median MLH minimum ≈ 04 UTC, maximum ≈ 15 UTC (solar noon at SIRTA ≈ 11:50 UTC); CHM15k vs CL31
  Spearman 0.80 (2020 →); JJA 15 UTC median 1,590 m (CHM15k) vs 1,400 m (CL31); DJF midday within ≈ 100 m.

### A3. Reader and standard files (`R/measured_reobs.R`)

- `read_reobs(path, vars)`: tidy hourly table with UTC start/end from `time_bnds`, values, std, nb; global
  attributes (product version, uuid, acknowledgement, licence, references) returned as attributes.
- Two standard files, both `method_class = aerosol`, `time_reference = interval_mean`, hourly:
  - `sirta_mlh_chm15k` (primary reference): `n_samples` = `_nb`, `uncertainty_m` = `_std`;
  - `sirta_mlh_cl31` (second aerosol method).
- `qc_flag`: 0 if `_nb` ≥ 75 % of nominal (≥ 45 of 60; ≥ 3 of 4), 1 if ≥ 50 %, otherwise `pblh_m = NA`.
  Write `qc_flag_meaning`. Ancillary `x_cbh_m`, `x_tcf_oktas`, `x_ablh_m` (CHM15k file only).
- Sidecar: instrument, `retrieval_algorithm` (STRATfinder / CABAM v2019.3), `licence` as in the file,
  `citation` = `Chiriaco2018`, `Kotthaus2020` (+ `Kotthaus2018` for CABAM) and the file's acknowledgement text,
  product version and uuid in `processing_history`; `height_reference = agl`.
- `config/measured_sites.csv`: new row FR-Gri | SIRTA-ReOBS | provider site "SIRTA Palaiseau" | 48.718, 2.207 |
  distance from `haversine_km()` | aerosol | CHM15k, CL31 | licence | notes ("hourly means; daytime convective
  only by distance rule"). Keep the existing Cloudnet `palaiseau` row.

### A4. FR-Gri runs

`scripts/01_check_inputs.R FR-Gri`, `02_run_pblh_model.R FR-Gri` (all current cases), `04_era5_compare.R FR-Gri`
(no `--with-oper`). Report site-specific points: crop site with seasonal canopy (time-varying hc, d, z0 from
`config/sites_timevarying.csv`), zm ≈ 2 m (very low: neutral B2/B3 and stable B1 behaviour, values below zm), the
reference period (FLUXES range), u* filter and L check counts, ERA5 cell distance and orography.

### A5. Comparison (generalise the spec-07 code)

Refactor the spec-07 comparison into functions that take any measured standard source (its time reference, qc
semantics and distance class), so FI-Hyy (half-hourly, censored) and FR-Gri (hourly means, not censored) use the
same code; FI-Hyy results must stay identical (regression check). Script `scripts/10_measured_compare.R SITE`.

Matching to the hourly measured means [t, t + 1 h):
- PBLHmod: mean of the two half-hours inside the hour (both required); the hour is convective if both half-hours
  are in the convective regime;
- ERA5 (instantaneous at the hour): mean of the values at t and t + 1 h (both required), as an approximation of
  the hourly mean; document it.

Strata:
- primary: daytime convective = PBLHmod convective **and** local solar time within `opts$aerosol_window`
  (default 10–16 h; aerosol-based heights track the mixed layer late morning to mid-afternoon: `Schween2014`,
  `Kotthaus2018`). Sensitivity: all convective hours;
- distance class from `measured_sites.csv` (23 km → daytime convective only): night and transition hours are shown
  in the diurnal-cycle figure for context only, labelled, and never enter metrics.

Pairs (metrics as `pblh_metrics()`; by season and hour; ERA5 tower cell; bilinear in tables only):
1. PBLHmod vs `sirta_mlh_chm15k`; 2. ERA5 vs `sirta_mlh_chm15k` (cross-class: Ri vs aerosol);
3. PBLHmod vs `sirta_mlh_cl31`; 4. ERA5 vs `sirta_mlh_cl31`;
5. `sirta_mlh_cl31` vs `sirta_mlh_chm15k` — the inter-method spread, reported as a "definitional uncertainty" band
   (project plan, "Comparing different method classes");
6. PBLHmod vs ERA5 at FR-Gri (for the record, as spec 04).

New metric for all pairs (and added to FI-Hyy's exact-value tables): the **range-compression slope**, Theil–Sen
slope of ln(model) on ln(measured) (slope < 1: the model compresses the range, i.e. too deep when shallow, too
shallow when deep — the candidate finding from spec 07), with a bootstrap 95 % interval; plus the median ln ratio
by tercile of the measured MLH.

Figures (PNG, ERA5 tower cell only): `10_diurnal_by_season_FR-Gri` (CHM15k, CL31, PBLHmod, ERA5; aerosol window
shaded; night labelled "not comparable"), `10_scatter_FR-Gri` (window, log axes, with the Theil–Sen line),
`10_logratio_by_measured_tercile_FR-Gri`. Tables: `10_metrics_FR-Gri.csv` (long format), `10_compression_{site}.csv`
(FR-Gri and FI-Hyy).

## Part B — Lanzhot BL-View L3 (CZ-Lnz), when the data arrive

If `data/raw/measured/CZ-Lnz/` holds more than the two sample days, implement; otherwise skip and report.
Facts from the two sample files (Cowork, 2026-10-26; verify on the full set): CL61, BL-View 2.3.2; 1-min records;
up to 3 candidates `aerosol_layer_heights` with `quality_indices` 1–3 (3 best), `detected_layers`,
`precipitation_detection`, `fog_detection`, `cloud_base_heights`; settings constant (method 0, sensitivity 10,
30–4000 m, variable averaging ≈ 3120 s); `mean_mixing_layer_height` is a trailing hourly mean (not used);
latitude/longitude/elevation = 0 (coordinates from the PI); **time base: probably local time** (`time_zone` = 2,
daylight signal centred ≈ 12:30 on 7 Oct) despite "seconds since 1970": convert with the per-record `time_zone`, check
DST transitions (gaps/duplicates); timestamps mark the **end** of the averaging window: centre at t − T/2 with
T = `time_averaging_period`.
Selection rule (documented, options): lowest candidate with quality 3 (sensitivity ≥ 2), no precipitation or fog,
below the lowest cloud base, continuity within ±30 min (|ΔMLH| ≤ `opts$max_jump`, default 300 m); aggregate to
half-hours (median, ≥ 50 % of minutes valid). Standard source `lanzhot_blview_mlh`, `method_class = aerosol`.
Comparison as A5 (distance class from the PI's coordinates).

## Part C — optional: Cloudnet ε at Palaiseau

Only after G's go: the spec-07 chain for `palaiseau` (2025-01 →; inventory: two lidars, ≈ 9.2 GB; report the size
first). Gives a turbulence-based and an aerosol-based MLH at the same site; compare them (pair 5 analogue).

## Bibliography

`Chiriaco2018` and `Kotthaus2020` added to `docs/references/pblh_literature_2026-10.bib` by Cowork (2026-10-27).

## Tests

ReOBS reader on a synthetic NetCDF fixture (time bounds, fill values, nb → qc); hourly matching (PBLHmod pair
required, ERA5 t and t + 1 h); solar-time window; Theil–Sen slope and bootstrap on synthetic data with a known slope;
FI-Hyy regression after the refactor (identical tables); Part B selection rule on synthetic candidate sets
(quality, precipitation, cloud base, continuity; time-zone and window-centre conversion).

## Docs

`CLAUDE.md` (ReOBS facts; BL-View facts if Part B runs), `docs/specs/standard_format.md` (aerosol hourly
interval-mean example), `docs/project_plan.md` (measured sources: SIRTA added; PI data status: FR-Gri SIRTA
received, Trappes/Melun and FI-Lom pending, CZ-Lnz pending; compression slope as a metric), `docs/STATUS.md`.

## Acceptance criteria

- Tests pass; FI-Hyy comparison identical after the refactor; scripts 01, 02, 04 and 10 run for FR-Gri.
- Standard files for both SIRTA sources written and validated; nothing under `data/` or `output/` committed.
- Metrics restricted to the daytime convective window; compression slopes for FR-Gri and FI-Hyy.

## Out of scope

Own MLH retrieval from backscatter (G, 2026-10-27); `ablh_chm15k` as a reference; footprint runs at FR-Gri (cluster,
spec 09); E-PROFILE data.

## Summary to report back

As in `CLAUDE.md`, plus: verification of the A2 facts; FR-Gri model/ERA5 results vs BE-Bra and FI-Hyy; window
metrics per pair; CL31 vs CHM15k spread; compression slopes (FR-Gri, FI-Hyy); anything unexpected.
