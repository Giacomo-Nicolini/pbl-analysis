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
