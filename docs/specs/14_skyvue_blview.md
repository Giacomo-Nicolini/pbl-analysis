# Spec 14: readers for two secondary measured sources: DE-Tha SkyVUE PRO (built-in MLH) and CZ-Lnz CL61 BL-View

Branch: `spec-14-skyvue-blview` (from `main`). Agreed with G on 2026-11-11. Light work only on the PC (samples
only; the full data are on the cluster).
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/specs/13_measured_eprofile.md` (CZ-Lnz BL-View prototype
`read_blview_l3`, `config/measured_sites.csv` columns, standard source conventions), `docs/specs/standard_format.md`,
`docs/specs/10_measured_aerosol_mlh.md` (aerosol window), `docs/figure_inventory.md`.

## Why

Two more ceilometers near ICOS towers. Both are manufacturer algorithms, so they are **secondary sources**: shown
and tabulated, never pooled with the E-PROFILE (CABAM / STRATfinder) sites. Measured data are never released or
committed; nothing from the files goes into tracked files except variable names and formats.

## Part A: DE-Tha, Campbell SkyVUE PRO (TU Dresden)

Facts from the site (Markus Hehn, 2026-11-11) and from Cowork's check of `samples/07_2026/E1227_150726_*` (verify):
- 113 m from the flux tower at bearing 275°, base 0.5 m above ground, not tilted; profiles from 2025-02-06;
  the built-in MLH (KNMI-based algorithm) from 2025-06-25. Time in UTC, timestamp = end of the interval.
- Daily files `E1227_<ddmmyy>_Midnight.{csv,txt,nc}`. **The `.nc` files are empty** (header only, 0 records): ignore.
- `.csv` (about 195 MB/day): one row per record (about 8.6 s, 10,063 rows on 2026-07-15), 26 columns; header
  `UTC Timestamp` (`dd.mm.yyyy HH:MM:SS`), cloud bases 1–4, vertical visibility, `Ceil Status`, `Backscatter
  Profile( (100000*srad*km)^-1 )` (2048 values in one quoted field), sky condition 1–5 and heights,
  `Mixing Layer 1..3( Meters )`, `Confidence 1..3( Confidence 1-3 )` (values 0–3), `Highest Signal`,
  `Ceil Resolution` (05). **Decimal comma** throughout; `/////` = missing. The profile field mixes decimal commas
  with comma separators (`3,65E-06,9,7E-06,...,0,`), so it is ambiguous: **don't parse the profile from the CSV.**
- `.txt` (about 105 MB/day): the instrument's message, 7 lines per record, same record count and order as the CSV,
  **no timestamp**. Line 1 `CS0016006` (after a control character), then status, sky condition, the profile header
  (`00100 05 2048 100 +39 00 0096 0020 30 131`: scale, resolution, n gates, laser energy, temperature, ...; take the
  field meanings from the SkyVUE PRO / CS135 manual, message 006, and cite the section), the MLH line
  (`00375 00002 00650 00001 01009 00001` = 3 heights with confidence), the profile as 2048 × 5 hex digits, and a
  checksum line.

Tasks:
1. **Reader `read_skyvue_day(csv, txt)`:** timestamps and MLH layers 1–3 with confidence from the CSV; the profile
   header fields and the decoded profile from the TXT. Pair the records by order, and check the pairing: equal
   counts, and the CSV cloud base and MLH equal to the TXT values for every record. Fail loudly on a mismatch.
   Decode the hex profile with the scale factor as the manual says, and check it against the CSV profile values
   where the CSV parse is unambiguous (e.g. the first 20 gates). Verify the checksum if the manual documents it.
2. **Standard source `tud_skyvue_mlh`** (`method_class = aerosol`, secondary): half-hourly, in the same way as the
   other aerosol sources; keep all three layers (`x_layer1_m`, `x_layer2_m`, `x_layer3_m` and confidences), and
   **don't choose a layer yet.** `pblh_m` = NA until G decides the rule; diagnostics below. Sidecar: provider, the
   "internal use only" status, instrument facts.
3. **Profiles:** don't write profile files yet. Just provide the decode function and its test: the profiles are for a
   later retrieval (possibly MeteoSwiss, pending).
4. `config/measured_sites.csv`: a DE-Tha row (SkyVUE PRO, position from the facts above, distance class < 5 km,
   status `secondary`).

## Part B: CZ-Lnz, Vaisala CL61 BL-View L3 (all versions)

Facts from the site contact (2026-11-11) and Cowork's check of the three samples in `data/raw/measured/CZ-Lnz/samples/`
(verify): the cluster folder holds his NCO extracts (`mlh_*`), 1,370 files in folders `BL-View v3.1.1`,
`v3.1.2`, `v3.2.2` (one year of 3.2.2; many 3.1.1 files without MLH).
- 3.1.1 / 3.1.2: `L3_OFFLINE__<YYYYMMDDhhmm>_<settings>.nc`, variables `time` (s since 1970), `date_stamp` (empty),
  `Mean_Layer_Height`, `Mean_Layer_QualityIndex` (both n×1), 3.1.2 also `location_utc_offset`, `Mean_Layer_Calculation_Time`.
- 3.2.2: `l3_<dd.mm.yyyy>.nc` (also the full L3 sample `l3_07.10.2026.nc`), `time_aerosol` ("end of the averaging
  period"), `mean_mixing_layer_height`, `mean_mixing_layer_quality` ("for the period of one hour"), fill −999.
- **1-minute records, but the value changes only on the full hour**: an hourly product repeated. Quality 1–3; which
  end is good is **unknown** (find it in the Vaisala BL-View documentation if you can, and cite it; otherwise leave it
  open).
- **Time:** the stamps look like UTC (3.1.x file names; offset 0), but the 3.2.2 day files run 23:54 → 22:59 UTC
  (cut at local midnight UTC+1). Not proven.

Tasks:
1. Extend the prototype `read_blview_l3()` to all three versions (detect by variables, not by folder name); collapse
   to **one value per hour** (check the value is constant within the hour; report the hours where it isn't); the
   hour's timestamp = end of the averaging hour; carry version, quality and a `-999` → NA mapping.
2. Standard source `vaisala_blview_mlh` (aerosol, secondary, hourly like SIRTA), all quality values kept, quality as
   a column; no filtering until G decides.
3. Report per version: date range, days with any valid MLH, share of valid hours, quality distribution, and the
   3.1.1 files with no MLH.

## Part C: diagnostics (run on the cluster; G runs, with a runbook)

For DE-Tha (test month July 2026, and later the full record) and CZ-Lnz (all versions):
- **Time convention test:** median diurnal cycle of the measured height (and, for DE-Tha, each layer) vs PBLHmod and
  ERA5 at the site, for clear days (kt ≥ 0.65, as in spec 10); the lag (−3 … +3 h, 1 h steps) that maximises the
  correlation of the morning growth (sunrise → solar noon), per season and per BL-View version. A lag of 0 supports
  UTC; ±1 or ±2 h means local time (CET / CEST). Report; don't correct anything automatically.
- **DE-Tha layers:** per hour of day, median of layers 1–3 and their confidence, vs PBLHmod and ERA5, so that G can
  decide a layer rule (or none).
- **CZ-Lnz quality:** window metrics vs PBLHmod and ERA5 by quality class (1, 2, 3).
- Output under `{output_dir}/tables/internal/14_*` and `figures/internal/14_*` (internal, not for the paper).

## Tests

SkyVUE: the CSV decimal-comma parser on the sample header and a record; the TXT record splitter; the pairing check
(a deliberately shifted record must fail); hex decode vs known CSV values; checksum if documented. BL-View: each
version's sample read; the hourly collapse (a synthetic file with a change inside an hour is reported); −999 → NA.
Add `docs/figure_inventory.md` entries for the new tables.

## Acceptance criteria

Tests pass on the PC samples; both standard sources defined; `measured_sites.csv` rows; diagnostics script and
runbook `docs/local/cluster_commands_spec14.md` (style of the earlier runbooks: `~/pbl_vars.sh` variables, job names
not used in it, array range checked before `sbatch`, no dependencies on ended jobs).

## Further tasks

- The standard-file reader guesses column types, so the CABAM `x_dlh_m` column is read as logical and lost (vroom
  "parsing issues" warnings at FR-Fon and FR-Lus in the script-11 cluster run). Read standard files with explicit
  column types from `docs/specs/standard_format.md` (unknown `x_*` columns as character or double, never
  logical-guessed); test on a file whose first rows of an `x_` column are empty.
- Docs and config only: (1) FI-Lom: set the CABAM (Kenttärova) row in `config/measured_sites.csv` to status
  inventory, reason "no overlap: EC record starts at the ICOS labelling date 2024-11-05, Kenttärova CL31 ends
  2024-05-16"; check in the Cloudnet inventory (spec 07) whether the Kenttärova ε MLH record continues after
  2024-11-05 and report its date range (no download). (2) Record in `docs/project_plan.md` (methods notes) and
  `docs/STATUS.md` (G, 2026-11-12): each site's EC record starts at its ICOS labelling date, the first valid L2 EC
  value; rows before it are empty grid rows. Add a per-site table of these dates to the network collection
  (`09_network_record_start.csv`, from the first valid `pblhmod_default` value: site, labelling date, last valid
  value, n_valid), for the supplement.

## Implementation (2026-11-12)

Code: `R/measured_skyvue.R`, `R/measured_blview.R` (rewritten), `R/measured_secondary_diag.R`,
`R/standard_format.R` (explicit types), `R/network_collect.R` (`network_record_start()`), `R/measured_cloudnet.R`
(status `secondary`), `R/plot_utils.R` (style rows); scripts `14a_skyvue_mlh.R`, `14b_blview_mlh.R`,
`14c_secondary_diagnostics.R`, `09_collect.R`; tests `test-measured_secondary.R` (and updates in
`test-measured_eprofile.R`, `test-run_status.R`, `test-gamma_evaluation.R`, `test-network_collect.R`).

### SkyVUE PRO message 006 (`CampbellSkyVUEPRO2026`, Revision 04/2026)

| line | content (Sect. 6.3.6, pp. 87–91) |
|---|---|
| 1 | SOH `CS` ID(1) OS(3) N(3) STX — `CS0016006`: ID 0, OS 016, message 006 |
| 2 | S(1) detection status (0 none, 1–4 cloud bases, 5 obscuration, 6 transparent, `/` suspect); WA(1) `0` / `W` / `A`; tr(3) window transmission %; h1–h4(5) heights (cloud bases; status 5: vertical visibility, highest signal); flags(12) alarm / warning bits (message 001) |
| 3 | sky condition: amount (oktas) and height (10 m units) of up to 5 layers |
| 4 | scale(5) Attenuated_SCALE %, default 100; res(2) gate [m]; n(4) gates; energy(3) laser pulse energy %; lt(3) laser temperature °C; ti(2) total tilt °; bl(4) background light mV; pulse(4) pulses × 1000; rate(2) sample rate MHz; sum(3) sum of normalised backscatter |
| 5 | mh1 q1 mh2 q2 mh3 q3 (5 characters each): MLH heights and quality 1–3, **3 = highest confidence** (Sect. 2.1.4); `/////` = no layer; all `/////` if the MLH option is not active |
| 6 | profile: n × 5 hex characters, signed 20-bit two's complement (values > 2^19 − 1 minus 2^20), × 1e-8 = attenuated backscatter [m⁻¹ sr⁻¹], scaled by Attenuated_SCALE (Sect. 2.1.3, p. 6); last two gates always 0; not tilt-corrected |
| 7 | ETX CRC-16(4) EOT: CRC-16 (polynomial 0x1021, initial 0xFFFF, output XOR 0xFFFF; = CRC-16/GENIBUS) over the characters after SOH up to and including ETX (Sect. 6.2, p. 68, with C code) |

Not in the manual: the height of a gate (centre or bottom). MLH algorithm: gradient method after the KNMI operational
algorithm (de Haij et al. 2007), quality from the backscatter contrast below / above the MLH (Sect. 2.1.4).

Checks on the sample day 2026-07-15 (PC): 10,063 records in both files (10,430 bytes each); every CRC valid; layers,
confidences and cloud base 1 identical in CSV and TXT on every record; decoded profile vs CSV on all 10,063 records,
maximum relative difference 3.5e-16 (exact). The CSV parse of the profile is ambiguous at a zero gate (`0,6E-08` =
`0` then `6E-08`): the comparison stops at the first zero. Scale 100 % (the scale direction for other values is our
reading, flagged), 2048 gates × 5 m, tilt field 0 and 1° (site: not tilted), record steps 1–13 s (median 11 s, mean
8.6 s). 24 s per day on the PC.

`tud_skyvue_mlh`: half-hours [t, t + 30 min) of the record time stamp minus 1 s; per layer the median height,
median confidence, number of records and share of confidence 3 (`x_layer<k>_m`, `x_layer<k>_conf`,
`x_layer<k>_n`, `x_layer<k>_share_conf3`), `x_n_records`, `x_cloud_share`, `x_window_tr_pct`; `n_samples` =
records with a valid CRC; `pblh_m`, `qc_flag` `NA`. Heights above the instrument (0.5 m above ground), used as is.

### BL-View versions (detected from the variables)

| version | variables | samples (2026-11-12) |
|---|---|---|
| 3.1.1 | `time` (int, attribute `Units`), `date_stamp` (filled, equal to `time` in UTC), `Mean_Layer_Height`, `Mean_Layer_QualityIndex` | 2024-10-01: 23 / 24 valid hours |
| 3.1.2 | as 3.1.1 + `location_utc_offset` (0), `Mean_Layer_Calculation_Time` (= first profile of the hour, hh:00:55) | 2025-03-23: 24 / 24 |
| 3.2.2 | `time_aerosol` ("end of the averaging period"), `mean_mixing_layer_height`, `mean_mixing_layer_quality` (fill −999); full files also `time_zone` (2) | 2025-10-13 (extract, firmware 1.2.7, 23:54 → 22:59), 2026-10-07 and 10-10 (full, firmware 2.3.2, 00:00:50 → 23:55:50) |

In every sample the value changes on the first profile after the full hour. Hour convention (implemented): the
profiles hh:00–hh:59 form the hour ending at hh:00 ([hh − 1 h, hh)); supported by the 3.1.2 calculation time, an
assumption for 3.1.1 and 3.2.2 (lag test). Hours with a changing value: `NA`, `x_consistent = FALSE`. Quality index:
the BL-View User Guide (M211185EN-E, NetCDF L3 File Format) gives no meaning of the values nor the hour convention;
Münkel et al. (ISARS 2012) show the gradient-minimum index 1–3 as red / yellow / green, higher for a stronger
gradient and a larger distance to other minima, lower with clouds (indicative: 3 = best; open for the L3 hourly
mean). Overlaps between files (reprocessing, day files): newest version kept (3.2.2 > 3.1.2 > 3.1.1), counted in
`14_blview_overlap_<site>.csv` (rule flagged for G). `vaisala_blview_mlh`: `pblh_m` = hourly value (all quality),
`qc_flag` `NA`, `x_quality`, `x_version`, `x_consistent`, `x_n_values`, `x_time_zone_h`; stored time read as UTC
(`--time-offset=0`).

### Further tasks

- `read_standard_pblh()`: all columns read as character, then typed: required / recommended as before; documented
  non-numeric `x_` columns from `standard_x_col_types()` (character, logical, ISO datetime); every other `x_` column
  double when all values parse, else character. The CABAM `x_dlh_m` is no longer lost.
- FI-Lom: the MeteoSwiss (KEN) row is `inventory` (no primary source), FI-Lom left `config/measured_run_sites.txt`.
  Cloudnet inventory (spec 07, 2026-10-22): Kenttärova ε-lidar FMI HALO 146 2022-08-23 → 2023-12-15 (426 files), FMI
  HALO 34 2026-04-03 → 2026-10-22 (203 files, volatile, continuing): after the labelling date ε only from
  2026-04-03. The Cloudnet row had the two HALO names swapped (corrected).
- `09_network_record_start.csv` (`network_record_start()`): labelling date (local date of the first valid
  `pblhmod_default` value), first and last valid value (UTC), n_valid; PC: BE-Bra 2020-10-26, FI-Hyy 2018-11-06.
