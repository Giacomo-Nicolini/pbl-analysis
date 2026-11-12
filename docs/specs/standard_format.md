# Standard PBLH format (as implemented, spec 04)

Code: `R/standard_format.R` (`standard_meta()`, `validate_standard_pblh()`, `write_standard_pblh()`,
`read_standard_pblh()`). Design: `docs/project_plan.md`, "Standard PBLH format".

## Files

One pair per site × source in `<standard_dir>/<site>/`:

- `<site>__<source_id>.csv.gz`: data table, comma-separated, missing values `NA`;
- `<site>__<source_id>.yml`: sidecar.

## Data table

| column | type | unit | notes |
|---|---|---|---|
| `timestamp_start_utc` | ISO 8601 UTC `YYYY-MM-DDTHH:MM:SSZ` | – | required, never missing |
| `timestamp_end_utc` | idem | – | required; ≥ start |
| `pblh_m` | numeric | m above ground | required; finite or `NA` |
| `source_id` | character | – | required; equal to the sidecar `source_id` |
| `qc_flag` | integer | – | recommended (may be all `NA`); meaning in sidecar `qc_flag_meaning` |
| `uncertainty_m` | numeric | m | recommended (may be all `NA`) |
| `n_samples` | integer | – | recommended (may be all `NA`) |
| `x_*` | any | documented by the writer | source-specific extras; any other column name is rejected |

Columns are written in this order (`x_*` sorted); missing recommended columns are added as `NA`.
Start times are unique.

Time reference (sidecar `time_reference`):
- `instantaneous` (e.g. ERA5 BLH): `timestamp_start_utc == timestamp_end_utc` on every row;
- `interval_mean` (e.g. half-hourly PBLHmod): `timestamp_start_utc < timestamp_end_utc` on every row.

## Sidecar (YAML)

Required (never `unknown`): `site_id`, `source_id`, `method_class` ∈ {`turbulence`, `aerosol`, `thermodynamic`,
`model`, `parameterisation`, `unknown`}.

Optional (value or `unknown`), 12 fields: `instrument`, `retrieval_algorithm`, `height_reference`,
`instrument_altitude_m`, `time_reference`, `original_resolution`, `original_qc`, `qc_flag_meaning` (map
flag → meaning), `processing_history`, `contact`, `licence`, `citation`.

Additional fields (spec 13, `standard_extra_fields()`; written only when given, not counted in the completeness):
`data_owner` (instrument operator / owner of the data), `processor` (who derived the product, method and version),
`known_limitations` (provider text).

Added automatically on write:
- `metadata_completeness` = fraction of the **12 optional fields** that are not `unknown` (a field is unknown if
  absent, empty, `NA` or the string `unknown`, case-insensitive). `created_utc` and `code_version` are excluded
  from numerator and denominator (spec 04, amendment 4);
- `created_utc`: write time, ISO 8601 UTC;
- `code_version`: `git describe --always --dirty --abbrev=7` of the repository (`-dirty` = uncommitted changes),
  or `unknown`.

## Validation

`validate_standard_pblh()` runs on write and on read and aborts listing all problems: missing required columns,
non-prefixed extra columns, wrong types, missing timestamps, start > end, inconsistency with `time_reference`,
duplicated start times, infinite `pblh_m`, `source_id` mismatch, sidecar required fields missing/unknown,
invalid `method_class`.

## Example of a measured source: `cloudnet_eps_mlh` (spec 07)

Turbulence-based mixing-layer height from the Cloudnet Doppler-lidar dissipation rate (surface-connected
ε ≥ 1e-4 m² s⁻³ from the lowest reliable gate, Manninen et al. 2018), half-hourly (`time_reference = interval_mean`,
UTC half-hours [start, start + 30 min)), `method_class = turbulence`, `height_reference = agl` (height a.s.l. minus
the altitude in the Cloudnet file; ±31 m ground-reference uncertainty at FI-Hyy noted in `processing_history`).
Written by `standardize_cloudnet_mlh()` (`scripts/07c_cloudnet_mlh.R`); the file uuids and PIDs of the inputs are
listed in `<site>__cloudnet_eps_mlh_files.csv` next to the standard file, referenced by the sidecar `citation`.

**Censoring convention** (`qc_flag`):

| `qc_flag` | `pblh_m` | Meaning | Use in comparisons |
|---|---|---|---|
| 0 | value | all valid profiles measured | exact |
| 2 | value | partly censored: median of the measured profiles | exact |
| 3 | value | signal-limited: more than half of the measured profiles ended at a gate without signal | **lower bound**: true MLH ≥ value |
| 1 | `NA` | more than half of the valid profiles censored below z_min | **censored**: true MLH in (0, `x_censor_below_m`] |
| `NA` | `NA` | coverage below `min_coverage`, or no valid profile | not usable |

Main `x_` columns: `x_censor_below_m` (z_min [m agl]), `x_coverage`, `x_n_expected`, `x_n_measured`,
`x_n_censored_weak`, `x_n_lowest_gate_missing`, `x_n_no_signal`, `x_n_precipitation`, `x_n_gap_in_column`,
`x_frac_signal_limited`, `x_precip_checked`, `x_stare_seconds`, `x_first_profile_utc`, `x_last_profile_utc`,
`x_cloud_base_agl_m`, `x_cloud_fraction`, and sensitivity cases as `x_pblh_m_<case>` / `x_qc_flag_<case>`
(`precip_column`, `precip_strict`, `bridged`). Comparisons code the three kinds as
`survival::Surv(left, right, type = "interval2")` intervals (`measured_interval()`).

## Example of an aerosol-based hourly source: `sirta_mlh_chm15k`, `sirta_mlh_cl31` (spec 10)

Mixing-layer heights from the SIRTA ceilometers at Palaiseau, as distributed in the SIRTA 1D ReOBS product (AERIS;
Chiriaco et al. 2018): CHM15k with STRATfinder (primary) and CL31 with CABAM v2019.3 (Kotthaus et al. 2018, 2020).
Hourly interval means (`time_reference = interval_mean`, UTC hours [start, start + 1 h)), `method_class = aerosol`,
`height_reference = agl`. Written by `standardize_reobs_mlh()` (`scripts/10a_sirta_reobs.R`); product version,
ReOBS version and file uuid in `processing_history`, licence and acknowledgement text from the file.

| `qc_flag` | `pblh_m` | Meaning |
|---|---|---|
| 0 | value | `n_samples` ≥ 75 % of nominal (60 one-minute values for CHM15k, 4 fifteen-minute values for CL31) |
| 1 | value | `n_samples` ≥ 50 % and < 75 % of nominal |
| `NA` | `NA` | no value, or `n_samples` < 50 % of nominal |

Not censored: every non-missing `pblh_m` is an exact value in comparisons. `uncertainty_m` = within-hour standard
deviation of the native values (not a retrieval uncertainty). `x_` columns: `x_cbh_m` (cloud base of the same
instrument), `x_tcf_oktas` (CL31 total cloud fraction), `x_coverage` (`n_samples` / nominal), and for CHM15k only
`x_ablh_m` (aerosol-layer top, not a mixing height).

Comparisons (`measured_time_base()`): PBLHmod = mean of the two half-hours inside the hour (both required; regime
common to both, else `mixed`); ERA5 (instantaneous) = mean of the values at t and t + 1 h (both required).

## Example of an aerosol-based source aggregated from 1-min values: `meteoswiss_stratfinder_mlh` (spec 11)

STRATfinder mixed-layer height of the MeteoSwiss CHM15k ceilometer at Davos (DAV), delivered at the native ~60 s
resolution with a quality level per profile. Aggregated by `standardize_meteoswiss_mlh()` (`scripts/11b_meteoswiss_mlh.R`)
to half-hourly interval values (`time_reference = interval_mean`, UTC [start, start + 30 min)), `method_class = aerosol`,
`height_reference = agl` (above ground at the ceilometer; `instrument_altitude_m` = 1590 m). `pblh_m` = median of the
quality-level-0 profiles in the half-hour (provider recommendation: level 0 only); `n_samples` = number of distinct
minutes with a level-0 value (profile times jitter around 60 s); `uncertainty_m` = `NA`.

| `qc_flag` | `pblh_m` | Meaning |
|---|---|---|
| 0 | value | level-0 values in ≥ 75 % of the 30 minutes |
| 1 | value | level-0 values in ≥ 50 % and < 75 % of the minutes |
| `NA` | `NA` | no level-0 value, or in < 50 % of the minutes |

Not censored. `x_` columns: `x_coverage` (`n_samples` / 30), `x_ablh_m` (median level-0 ABLH, aerosol-layer top),
`x_share_lowest` (share of level-0 profiles at the lowest detectable height, bit 8), `x_segment` (homogeneous
instrument segment), `x_n_profiles` (all profiles). Comparisons (`measured_time_base()`, half-hourly base without
censoring): PBLHmod on the same half-hours; ERA5 hours against the mean of the two adjacent half-hours (both
required); qc is a coverage flag, not a censoring code.

Spec 13: the same source at the E-PROFILE stations Jülich (DE-RuS) and Magdeburg (DE-HoH); `x_share_no_temp_overlap`
(share of level-0 profiles with bit 8192, no temperature-dependent overlap model; informational) added for all.

## Example of an aerosol-based source from 15-min means: `meteoswiss_cabam_mlh` (spec 13)

CABAM mixed-layer height of the Vaisala CL31 at the E-PROFILE stations Kenttärova (FI-Lom), Melun (FR-Fon) and
Poitiers (FR-Lus), delivered by MeteoSwiss as 15-min means (`time` = interval centre, `time_bnds` = [start, end) on
the quarter hour), heights tilt-corrected, above ground. Aggregated by `cabam_halfhour()` (`scripts/11b_meteoswiss_mlh.R`)
to half-hourly interval values (`interval_mean`, UTC [start, start + 30 min) = the steps starting at start and
start + 15 min), `method_class = aerosol`, `height_reference = agl`. `pblh_m` = mean of the quality-level-0 steps;
`n_samples` = number of level-0 steps (0–2); `uncertainty_m` = `NA`.

| `qc_flag` | `pblh_m` | Meaning |
|---|---|---|
| 0 | value | both 15-min steps level 0 (their mean) |
| 1 | value | one step level 0 (its value) |
| `NA` | `NA` | no level-0 step |

`x_` columns: `x_coverage` (`n_samples` / 2), `x_dlh_m` (mean CABAM disconnected-layer height of the steps with one),
`x_share_lowest` (bit 8 among the level-0 steps), `x_segment`, `x_n_steps` (steps in the file), `x_n_profiles` (raw
profiles of the level-0 steps), `x_tilt_deg`. Sidecar `licence`: to be confirmed, internal use only; `data_owner`,
`processor`, `known_limitations` filled. Comparisons as for `meteoswiss_stratfinder_mlh`.

## Column types of `x_` columns (spec 14)

`read_standard_pblh()` reads every column as text and types it explicitly (no guessing from the first rows; before
spec 14 a column empty in its first rows, e.g. the CABAM `x_dlh_m`, was guessed as logical and lost). Required and
recommended columns as in the table above. Documented non-numeric `x_` columns (`standard_x_col_types()`):

| column | type |
|---|---|
| `x_regime`, `x_method`, `x_h_source`, `x_version` | character |
| `x_precip_checked`, `x_consistent` | logical (`TRUE` / `FALSE`; any other value aborts) |
| `x_first_profile_utc`, `x_last_profile_utc` | ISO 8601 UTC time |

Every other `x_` column is double when all its values parse as numbers (an empty column too), otherwise character.
A writer that adds a non-numeric `x_` column adds it to `standard_x_col_types()` and to this table.

## Secondary measured sources (spec 14)

Manufacturer algorithms near ICOS towers: shown and tabulated, never pooled with the E-PROFILE sites, never a primary
source (`status = secondary` in `config/measured_sites.csv`). Internal use only.

**`tud_skyvue_mlh`** (DE-Tha, TU Dresden Campbell SkyVUE PRO, built-in MLH; `scripts/14a_skyvue_mlh.R`,
`standardize_skyvue_mlh()`): half-hourly interval values [start, start + 30 min) (a record belongs to the half-hour
of its time stamp minus 1 s; the stamp is the end of the ~8.6 s record), `method_class = aerosol`, heights above the
instrument (0.5 m above ground). **No layer chosen: `pblh_m` and `qc_flag` are `NA` on every row** until G decides a
layer rule. `n_samples` = records with a valid CRC-16. `x_` columns per layer k = 1–3: `x_layer<k>_m` (median height
[m]), `x_layer<k>_conf` (median quality 1–3, 3 = highest confidence), `x_layer<k>_n` (records with the layer),
`x_layer<k>_share_conf3`; `x_n_records`, `x_cloud_share` (records with detection status 1–4), `x_window_tr_pct`.

**`vaisala_blview_mlh`** (CZ-Lnz, CzechGlobe CL61 with Vaisala BL-View L3; `scripts/14b_blview_mlh.R`,
`standardize_blview_mlh()`): hourly interval values [hour end − 1 h, hour end), `method_class = aerosol`. `pblh_m` =
the BL-View hourly mean MLH with **all quality values kept**; `qc_flag` `NA` (no filtering until G decides);
`n_samples` = profiles of the hour. `x_quality` (1–3; meaning not documented, indicative 3 = best), `x_version`
(3.1.1, 3.1.2, 3.2.2, from the file variables), `x_consistent` (`FALSE`: value changing within the hour, then
`pblh_m` `NA`), `x_n_values`, `x_time_zone_h` (3.2.2 file variable), `x_time_shift_h` [h, corrected minus stored
hour end]. Time corrected by default (spec 14 amendment 1 item 11, G 2026-11-12): 3.2.2 stamps read as Europe/Prague
local time (summer time) and converted to UTC, then every version −1 h (empirical alignment from the lag test,
mechanism unknown); details in the sidecar.
