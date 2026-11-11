# Spec 13: measured MLH from MeteoSwiss / E-PROFILE (five stations) and CZ-Lnz preparation

G, 2026-11-06. Branch `spec-13-eprofile`. Light work on the PC; standard files and QC report are light (≤ 20 s per
site).

## Inputs

`data/raw/measured/` (cluster: `measured_dir`), one folder per ICOS station plus `_shared/`:

| ICOS | File | Station (code, WIGOS) | Instrument, method, step | Period |
|---|---|---|---|---|
| DE-RuS | `MeteoSwiss_STRATfinder_MLH_ABLH_JUE_20210318_20261101.nc` | Jülich (JUE, 0-276-0-10508) | CHM15k, STRATfinder, 60 s | 2021-03-18 → 2026-11-01 |
| DE-HoH | `MeteoSwiss_STRATfinder_MLH_ABLH_MAG_20210101_20261101.nc` | Magdeburg (MAG, 0-20000-0-10361) | CHM15k, STRATfinder, 60 s | 2021-01-01 → 2026-11-01 |
| FI-Lom | `MeteoSwiss_CABAM_MLH_KEN_20210114_20240516.nc` | Kenttärova (KEN, 0-246-0-101987) | CL31, CABAM, 15 min | 2021-01-14 → 2024-05-16 |
| FR-Fon | `MeteoSwiss_CABAM_MLH_MEL_20220823_20261101.nc` | Melun (MEL, 0-20000-0-07153) | CL31, CABAM, 15 min | 2022-08-23 → 2026-11-01 |
| FR-Lus | `MeteoSwiss_CABAM_MLH_POI_20210101_20261101.nc` | Poitiers (POI, 0-20000-0-07335) | CL31, CABAM, 15 min | 2021-01-01 → 2026-11-01 |

`_shared/meteoswiss_eprofile_2026-11-05/`: `README.md`, `qc_summary.json`, `qc_events_stratfinder.csv`, PNGs.
Product 1.0 (2026-11-05), same layout as the Davos file (spec 11 C). Data owners: University of Cologne (JUE), DWD
(MAG), FMI (KEN), Météo-France (MEL, POI); processing MeteoSwiss.

**Licence: to be confirmed. Internal use only; contact the provider before sharing results** (MeteoSwiss
condition, 2026-11-06). Sidecar `licence` says so; nothing derived from these files leaves the project.

## 1. Readers and standard sources (done)

- STRATfinder (`read_meteoswiss_stratfinder()`, unchanged aggregation of spec 11 C): source
  `meteoswiss_stratfinder_mlh`, half-hourly median of the quality-0 profiles, qc from the distinct minutes (0 ≥ 75 %,
  1 ≥ 50 %). Differences from Davos, handled:
  - bit 8192 `no_temperature_overlap_model` (Jülich TUB120017, 2021-03 → 2026-01): informational, quality stays 0;
    kept, share per half-hour in `x_share_no_temp_overlap`;
  - no maintenance log: bit 4 only from the L1 window transmission (daily median < 60 %), quality 1, so excluded by
    the quality-0 rule. The window-transmission days are in the file attribute `instrument_events` and in
    `qc_events_stratfinder.csv`; 11b and 13a check that the file flags them (all listed days flagged on ≥ 99.9 % of
    the profiles; no flagged day missing from the list);
  - every qc bit must be in `meteoswiss_qc_bits()` (declared `flag_masks` and set bits); a new bit aborts.
- CABAM (`read_meteoswiss_cabam()`, `cabam_halfhour()`): source `meteoswiss_cabam_mlh` (`method_class` aerosol).
  `time` = interval centre, interval from `time_bnds` [start, end); checked on all three files: every step is 15 min,
  starts on the quarter hour, time = start + 7.5 min. Half-hour [t, t + 30 min) = the steps starting at t and
  t + 15 min. `pblh_m` = mean of the quality-0 steps; qc 0 both good, 1 one good, else NA. Heights already
  tilt-corrected, above ground. Bit 16384 `incomplete_interval` → quality 1 (excluded). `dlh` kept as `x_dlh_m`
  (mean of the steps with a value; only the dirty-window removal applies to it); `alh` (9 unordered slots) not read.
  Other columns: `x_coverage`, `x_share_lowest`, `x_segment`, `x_n_steps`, `x_n_profiles`, `x_tilt_deg`.
- Sidecar extra fields (`standard_extra_fields()`, not counted in the completeness): `data_owner` (operator),
  `processor` (MeteoSwiss, method and version, client product version), `known_limitations` (file attribute).
  Citation: `Kotthaus2020` (STRATfinder); `Kotthaus2018` (Part 1, qj.3299) and `Kotthaus2018b` (Part 2, qj.3298,
  added to the bib, checked 2026-11-06) for CABAM.
- `scripts/11b_meteoswiss_mlh.R` handles both file types (station code and distance checked against
  `measured_sites.csv`).

## 2. Pairing table (done)

`config/measured_sites.csv` new columns: `station_code`, `wigos_id`, `provider_alt_m`, `elevation_diff_m` (tower
minus ceilometer), `distance_class` (checked against `distance_class()` on read), `status` (`used`, `inventory`,
`pending`; a primary source only on a `used` row, at most one per site).

| ICOS | Station | km | Elevation diff. [m] | Class | Note |
|---|---|---|---|---|---|
| DE-RuS | JUE | 5.30 | −5 | daytime_convective | JOYCE (same site as the Cloudnet row) |
| DE-HoH | MAG | 24.48 | +132 | daytime_convective | file coordinates rounded to 0.01° (±0.5 km): borderline |
| FI-Lom | KEN | 1.80 | −100 | all_regimes | Kenttärova hill 100 m above the mire |
| FR-Fon | MEL | 16.64 | +12 | daytime_convective | |
| FR-Lus | POI | 24.87 | +31 | daytime_convective | borderline |

Borderline pairings (G 2026-11-06): DE-HoH (24.48 km) and FR-Lus (24.87 km) keep the daytime-convective class;
column `borderline` (logical, `FALSE` when missing) in `measured_sites.csv`. Every per-site comparison table of
`scripts/07_measured_compare.R` and `scripts/10_measured_compare.R` carries `site_id`, `measured_provider_site`,
`measured_distance_km` and `borderline` (`tag_measured_pairing()`), so network results can be shown with and
without them. (DE-RuS JUE, 5.30 km, is 0.3 km above the 5 km all-regimes limit; not flagged, G to say.)

Primary source: the MeteoSwiss row at all five. At FI-Lom and DE-RuS the Cloudnet rows are now `inventory` (their ε
files cover 2022-08 → 2023-12 / from 2023-04, not the reference years; never downloaded).

## 3. Additional QC (done; internal, `scripts/13a_meteoswiss_qc.R`, outputs not committed)

Tables and figures per site under `output/` (`13_qc_*`). Findings (2026-11-06):
- **KEN (FI-Lom)**: 69 % good; 212 of 1,219 days missing: 115 polar night (Dec–mid Jan, not processed), 30 polar day,
  67 other; May–Aug 56 missing (MeteoSwiss: 49 midsummer days stopped in the CABAM noon search). Window flag 0.4 %.
  Good values ≤ 100 m 7 %, minimum 39 m. Share good drops from 0.75 after sunrise to 0.64 around sunset. Record ends
  2024-05-16: **no overlap with the reference year 2025** (no measured footprints). The FI-Lom ICOS record starts
  2024-01, so the comparison covers only 2024-01 → 2024-05-16 (winter and spring, few convective hours).
  **Correction (2026-11-11): no overlap at all.** The FLUXES rows start 2024-01, but the EC values (and PBLHmod) start
  2024-11-04 22:00 UTC; see spec 11, amendment 1, Part B correction.
- **MEL (FR-Fon)**: window flag 10.8 %; ≥ 30 % in 2023-05, 2023-06, 2025-09, 2026-08, 2026-09, 2026-10. In those
  months the remaining midday median is 1,414 m vs 1,012 m elsewhere (night 512 vs 384 m): consistent with spurious
  deep layers. Breakpoint 2025-12-16 (firmware 321, new corrections): no step in the 60-day medians.
- **POI (FR-Lus)**: window flag ≥ 30 % in 2021-03, 2021-09, 2021-10 (the README's 2021-02 → 2021-10 period, lower
  shares in the other months); breakpoints 2021-10-26 and 2025-11-05 (new correction periods) coincide with a drop of
  the 60-day midday median (893 → 509 m, 835 → 494 m): mostly seasonal (Sep–Oct vs Nov–Dec), not separable here.
- **JUE (DE-RuS)**: 96 % good; no data 2025-11-19 → 2026-06-17. Bit 8192 (no temperature overlap model, until the
  optical unit change) vs TUB160062, same months (Jun–Nov): night median 352 vs 382 m, mode 320 vs 290 m, night
  p90 1,429 vs 561 m, share at the lowest height 1.3 % vs 0.3 %: the night distribution differs (the README warns that
  low nocturnal layers are less reliable without the model). Daytime: midday median 1,010 vs 1,160 m (different
  years). Night MLH is outside the comparable set anyway (daytime convective class).
- **MAG (DE-HoH)**: 98 % good; breakpoint 2023-09-22 (firmware 1.12, new overlap model): night median 322 → 262 m and
  share at the lowest height 0 → 3 % in the 60 days after (the README: do not compare night values across it).

**Exclusions (approved G 2026-11-06, applied to the measured source only):** FR-Fon the six MEL months above
(2023-05, 2023-06, 2025-09, 2026-08, 2026-09, 2026-10; dirty window); FR-Lus the whole README period 2021-02-01 →
2021-10-31 (dirty window). None at KEN, JUE, MAG (night issues are outside the daytime comparable set).

Mechanism: `config/measured_exclusions.csv` (`measured_exclusions_file` in `config.yml`), columns `site_id`,
`source_id`, `date_start`, `date_end` (UTC calendar dates, end inclusive), `reason`. Read by
`read_measured_exclusions()` (all fields required, `source_id` must be a source of `pblh_sources_style`), applied by
`apply_measured_exclusions()` in every measured standard-file builder (`standardize_meteoswiss_mlh()`,
`standardize_reobs_mlh()`, `standardize_cloudnet_mlh()`): rows whose `timestamp_start_utc` falls in a period are
dropped, counted in the log and in the sidecar `processing_history` (with the periods and reasons). Unlike
`config/exclusion_periods.csv`, the ICOS records stay. Local run: FR-Fon 8,784 half-hours removed (183 days), FR-Lus
13,104 (273 days). `meas_std` and `meas_cmp` hash `measured_sites.csv` and `measured_exclusions.csv` (input group
`measured_config`), so an edit of either reruns them (and, through the standard files, `08a`).

## 4. Comparison and footprints

- `config/measured_run_sites.txt`: 8 sites (FI-Hyy, FR-Gri, CH-Dav, DE-RuS, DE-HoH, FI-Lom, FR-Fon, FR-Lus).
- MeteoSwiss comparison = the generic spec-10 script (`measured_site_scripts()`, also at CH-Dav, which keeps
  `scripts/11_gamma_evaluation.R`). `scripts/10_measured_compare.R` now takes the reference from the site's primary
  source: one half-hourly aerosol source gives pairs 1, 2, 2b and 6 (PBLHmod on the half-hours, ERA5 hours vs the mean
  of the two adjacent half-hours), primary set window (10–16 h AST) and PBLHmod convective; FR-Gri unchanged.
- Footprint `measured` case (aerosol: comparable set only) in the approved reference years: DE-RuS 2022, DE-HoH 2023,
  FR-Fon 2023, FR-Lus 2023 (12 months with window values each, FR-Fon 10 after the exclusion of 2023-05 and
  2023-06); FI-Lom 2025 does not overlap (record ends 2024-05; reference year unchanged). Extra chunks (upper bound,
  before the PBLHmod-convective condition): 46 new `measured` chunks, and the `pblhmod`, `pblhmod_clim`, `era5`
  chunks of the same months (138) change hash because their `common_pme` groups now exist; ≈ 85 CPU-hours nominal. Exact list: `scripts/08c_make_chunks.R --pending`
  (`chunk_list_pending.csv`, new / changed by the hash of 08b, `footprint_task_hash()`).

## 4b. Fix of step 08 at the non-extra sites (G 2026-11-06)

`scripts/08_footprint_compare.R` failed on the cluster at the 33 sites without the supplementary cases (after
10–15 min, at the join with `kljun2015_table5_h20()`): the ±20 % anchor list was empty. Now the anchor comparison
is skipped when `pblhmod_x0.8` / `pblhmod_x1.2` are absent: one log line, and `08_anchor_<site><period>.csv` holds a
single `note` row saying why. The other steps only use the supplementary cases through `intersect()` / `setdiff()`
(`oper`, `measured`, `bound`, `era5_bilinear`, `gamma001` / `0025` enter the generic loops only when present; ERA5
blocks are guarded by `"era5" %in% cases`); empty tables print "(no rows)". Checked with a fake non-extra site built
from the BE-Bra 2022-09 test chunks with only `pblhmod`, `pblhmod_clim` (a copy of `pblhmod`) and `era5`: the old
script stops with the cluster error, the new one writes every table and figure. Step 08 has to be rerun at all 39
sites (`docs/local/cluster_commands_spec13.md`).

## 5. CZ-Lnz (pending)

Row in `measured_sites.csv` with `status = pending`: Vaisala CL61 with BL-View at the Lanzhot grassland tower, ~3 km
from CZ-Lnz, ~3 m elevation difference; L3 `mean_mixing_layer_height` and `mean_mixing_layer_quality` (hourly,
1–3); usable from 2024-05 (BL-View 3.1.1), reprocessed with 3.2.2 since 2025-09, firmware 2.3.2 from 2025-12-11.
Reader prototype `read_blview_l3(path, time_offset_h)` (reads only `time_aerosol` and the two variables): the hourly
value is repeated on every ~1-min profile of the hour and changes at hh:00:50; one row per hour. Sample facts: 1,436
profiles per day, stored time 00:00:50 → 23:55:50, `sw_version` 2.3.2, latitude / longitude / elevation 0, a
`time_zone` variable = 2. To confirm with the site: time zone of the stored time (UTC or CEST), whether the value
of the profiles hh:00–hh:59 belongs to [hh, hh + 1 h) or to the hour ending at hh, the meaning of the quality index.
No standard file yet.
