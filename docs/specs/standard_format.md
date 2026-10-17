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
