# Project status (handover note)

Last update: 2026-11-26. Read this first in a new session, then `docs/project_plan.md` and `CLAUDE.md`.

## Analysis on stand-by from 2026-11-19 (reopened 2026-11-25 for specs 19 and 17, G)

Everything up to spec 16 and the close-out fixes is merged on `main` (spec 16 figures and spec 15 exception minimum:
`57a9460`; SkyVUE reader fixes: `e5e52b7`). Pending on the cluster (not blocking): 14a + 14c for DE-Tha with the new
SkyVUE reader (`docs/local/cluster_commands_spec14.md`, section "Close-out"; the 2026-11-20 rerun left 4 June 2026 days
failing, fixed by spec 14 amendment 3 on `main`: rerun 14a + 14c for DE-Tha when the analysis reopens, no hurry, G).
Nothing else is running or scheduled.

**Spec 17, decision D1: rejected (pre-test 2026-11-26, local, BE-Bra, FI-Hyy, FR-Gri, CH-Dav; no spec-17 cluster
run).** Methods / supplement note:
- *Question.* The first local run gave γ_EC 2–3 × the ERA5 site γ and above the measured benchmark γ_best
  (constant γ where the median ln(PBLHmod / measured) on the primary set crosses 0; spec 17 B): FR-Gri 0.0030,
  CH-Dav 0.0037 K m⁻¹ (FI-Hyy: no crossing on the grid 0.002–0.012; BE-Bra: no measured reference). Suspected cause:
  the early-morning θ rise at the tower (erosion of the nocturnal inversion while the model grows from the night-time
  h0; superadiabatic surface layer).
- *Test.* Option `gamma_ec_start_offset_h` (config, default 0; `17a --start-offset-h=`): the fit starts 2 or 3 h after
  the first positive heat flux, the end unchanged; h is still integrated from the first positive flux (h at the fit
  start model-consistent for each candidate γ), only the earlier θ residuals are left out. Offset 0 reproduces the
  existing tables (checked). Criterion fixed beforehand (G): at FR-Gri and CH-Dav the site median γ_A within a factor
  1.5 of γ_best.
- *Results* (used days; share of the fitted days with the optimum at the upper bound 0.02; median γ_A [K m⁻¹];
  median γ_B / γ_A; part A median ln(PBLHmod-γEC / measured), default −0.53 FR-Gri, −0.43 CH-Dav, +0.44 FI-Hyy):

  | site | offset 0 | +2 h | +3 h |
  |---|---|---|---|
  | FR-Gri | 112; 0.78; 0.0141; 0.92; −0.64 | 212; 0.63; 0.0120; 0.94; −0.63 | 256; 0.53; 0.0115; 0.94; −0.60 |
  | CH-Dav | 518; 0.44; 0.0074; 0.87; −0.32 | 586; 0.21; 0.0047; 0.90; −0.14 | 466; 0.16; 0.0052; 0.90; −0.14 |
  | FI-Hyy | 423; 0.50; 0.0094; 0.74; +0.43 | 516; 0.35; 0.0067; 0.80; +0.49 | 532; 0.27; 0.0063; 0.81; +0.53 |
  | BE-Bra | 346; 0.45; 0.0102; 0.78; – | 452; 0.24; 0.0066; 0.83; – | 422; 0.18; 0.0060; 0.86; – |

  Median γ_A / γ_best: FR-Gri 4.7 → 4.0 (+2 h) → 3.9 (+3 h); CH-Dav 2.0 → 1.26 → 1.38. CH-Dav passes with either
  offset, FR-Gri fails both: **D1 = reject** (G's rule: stop, no cluster run).
- *What we learned.* (1) The early-morning rise explains part of the excess: a later start halves the share of
  upper-bound optima and lowers γ_A by 15–35 %, bringing CH-Dav close to γ_best and to the ERA5 γ (0.0044). (2) It
  does not explain FR-Gri, where γ_A stays ≈ 4 × γ_best with a later start: the θ rise is still too fast for the
  heat put in by the measured H. Consistent with H too small for the θ budget (the fitted γ absorbs any H bias, see
  `docs/methods_gamma_ec.md`), and with advection and radiative warming that the slab model attributes to the surface
  flux; FR-Gri's daytime EBR is 0.97, so non-closure alone does not account for it. (3) Even where γ_A approaches
  γ_best, PBLHmod-γEC gains little over PBLHmod-clim (CH-Dav −0.14 vs −0.08), and at FR-Gri it is worse than the
  default (−0.60 vs −0.53). (4) γ_B / γ_A ≈ 0.8–0.94: the closed form agrees within 20 %, so the excess is in the θ
  rise itself, not in the integration. A tower θ-rise γ is therefore not recommended as an operational γ; PBLHmod-clim
  (ERA5) and the constant `gamma_med` stay the options. Tables of the pre-test: `output/tables/17d1s2_*`,
  `17d1s3_*` (script 11 `--prefix=17d1s2` / `17d1s3`, local, not collected); 17a summaries and daily tables of the
  variants in the session scratchpad only (the baseline gamma_EC tables and outputs were restored).
- *Supplement network run* (G 2026-11-26; branch `spec-17-supp-network`): 17a at the 39 sites with offsets 0 and
  +2 h (`<derived_dir>/gamma_ec_off2h/`) and 17b parts C, D, with the new height table `17_gamma_ratio_vs_height`
  and figure `17s_gamma_ratio_vs_height`; no PBLHmod-γEC case (step 02 runs the γEC cases only when named), no
  script 11, no pair 9. Scripts called directly (run_site would rerun 01/11a/02): ≈ 4.5 CPU-hours.
  Runbook: `docs/local/cluster_commands_spec17.md` section 0. Cluster run done (17b log `gec17b_19467.out`,
  results in `data/cluster_spec17/`): 38 sites with a ratio (FI-Lom none); Spearman rho(ratio, zm) −0.28 (+0 h,
  p 0.09) and −0.30 (+2 h, p 0.07); land-cover medians +0 h / +2 h: crop / grass 3.11 / 2.69, DBF 2.91 / 2.22,
  ENF / MF / EBF 2.26 / 1.40, wetland 1.89 / 1.60. Figure redrawn as one panel (G 2026-11-27,
  `scripts/17r_gamma_ratio_height_figure.R`, ggrepel added to renv). IT-BCi (+2 h ratio 1.00 vs 2.36 at +0 h;
  checked 2026-11-27): no fit fault (r² 0.93, no bound hits, flags as elsewhere), but the window ends at the day's
  θ maximum at ≈ 12 UTC (other crop sites ≈ 14 UTC; 46 % of the used days end before 12 UTC), so the +2 h fit has
  9 rows and a θ rise of 2.5 K (6.5 K over the full window): θ levels off around midday, probably the sea breeze
  (coast ≈ 10 km; not checked with the wind direction). Per-day γ_A(+2 h) / γ_A(+0 h) median 0.35 (FR-Gri 0.69).

The analysis reopens for one of the following; steps listed, nothing run.

**1. Ummendorf ceilometer for DE-HoH** (about 9 km instead of Magdeburg MAG 24.48 km), if MeteoSwiss provides it.
- Data: the E-PROFILE client file under `<measured_dir>/DE-HoH/`. 11b reads the *single*
  `MeteoSwiss_(STRATfinder|CABAM)_*.nc` of the site folder: move the Magdeburg file out of `measured_dir` (e.g. an
  archive folder beside it) first. New qc bits abort the reader (`meteoswiss_qc_bits()`): review them if it stops.
- `config/measured_sites.csv`: new DE-HoH row for Ummendorf (provider MeteoSwiss, `provider_site` / `station_code` =
  the file's station code, `wigos_id`, `provider_lat/lon/alt_m`, `distance_km` with `haversine_km()`,
  `elevation_diff_m`, `distance_class` (9 km stays `daytime_convective`, 5–25 km), `borderline` FALSE, method class,
  instruments, licence, notes), `primary_source_id` (`meteoswiss_stratfinder_mlh` for a CHM15k, `meteoswiss_cabam_mlh`
  for a CL31), `status` `used`; the MAG row: `primary_source_id` empty, `status` `inventory`.
  `config/measured_exclusions.csv`: only if the QC report shows source periods to drop (G decides).
  `config/measured_run_sites.txt`: unchanged (DE-HoH already listed).
- Reruns, DE-HoH only: `scripts/13a_meteoswiss_qc.R` (QC report) → `run_site.R DE-HoH --steps=meas_std,meas_cmp,08a`
  (11b standard file, 10 comparison; 08a hashes the measured standard file, so check that the DE-HoH reference year
  of the proposal is still the approved one in `config/reference_years.csv`, else G re-approves) → `08c --pending` → the pending chunk array (`measured`
  case and the `pblhmod`, `pblhmod_clim`, `era5` chunks of the months with measured values, `common_pme`) →
  `run_site.R DE-HoH --steps=08` → `09_collect.R` (network tables, spec-15 summary). Then `11_gamma_evaluation.R`
  (one job; DE-HoH is one of its sites; decision 3b stays unless G reopens it).
- Spec 16: `h_empirical` is the median σ_ln over the primary aerosol sources, DE-HoH included. `16a
  --recompute-h-empirical` prints the new σ_ln; only if it moves does design B need a rerun (B chunks hold both h
  options: ≈ 95 CPU-hours at 3.3 s) and 16c. DE-HoH is not one of the six Monte Carlo sites.
- Paper and docs: distance table, `CLAUDE.md` E-PROFILE section, `docs/figure_inventory.md` only if a table changes.

**2. BE-Maa: RMI Diepenbeek CL51 data** (about 16 km from the tower, approximate station coordinates; check from the
files). New reader for the secondary set (spec 14 conventions):
- Before coding: file format, time reference, the MLH algorithm (Vaisala BL-View / manufacturer = `secondary`; an
  RMI retrieval of its own could be a primary candidate: G decides), licence and owner (internal use until cleared).
- Code (new branch, spec 14 amendment or a new spec): `R/measured_rmi.R` (reader with unit and qc-bit checks, tests
  on a synthetic fixture), a standard source id (e.g. `rmi_cl51_mlh`), a script `scripts/14d_rmi_mlh.R` (standard
  file, exclusions applied and counted), documented non-numeric `x_` columns in `standard_x_col_types()` and
  `docs/specs/standard_format.md`, a style row in `pblh_sources_style`, the `config/measured_sites.csv` row
  (`status` `secondary`), the BE-Maa case in `scripts/14c_secondary_diagnostics.R` (window metrics, lag test, diurnal
  cycle vs PBLHmod, PBLHmod-clim, ERA5), `docs/figure_inventory.md` (`14_*` tables and figures).
- Cluster: the new 14d script and 14c for BE-Maa only; secondary sources are never pooled and enter no footprint case.

**3. Different 1σ values for the spec-16 inputs after the co-author review.**
- Change the CVs in `mc_input_spec()` (`R/footprint_mc.R`; σv also as an argument) and the header line of
  `scripts/16c_mc_collect.R` that lists them; the methods note in `docs/project_plan.md` (5b) and spec 16.
- Reruns: `16a` (same seed per site, so the sample and the standard deviates are unchanged; the new spec is stored
  in the draws file; with `--recompute-h-empirical` only if the h option changes) → `16b` array (every chunk: the
  chunk hash includes the perturbed inputs) → `16c`. The spec-15 footprints and the base runs are unchanged.
- CPU estimate from the last cluster run (collected tables, reference years, six sites): 277,997 runs, 192,268 FFP
  calls after removing invalid and duplicated runs (A 88,757, B 103,511) ≈ 176 CPU-hours at 3.3 s per call (PC
  rate; ≈ 107 at 2 s). Larger σ make more runs invalid (fewer calls), smaller σ fewer; if every run were computed:
  ≈ 255 CPU-hours at 3.3 s (below the 300 CPU-hour limit of spec 16). The 16b wall times are in `sacct` of the 16b array (not copied back). Re-render only (no rerun) with
  `16c --from-tables=DIR`.

## PBLH sources (names used in the code and in the standard files)

| Source | `source_id` | What it is | Class |
|---|---|---|---|
| PBLHmod v2 ("PBLHmod") | `pblhmod_default` | Our EC-based model (Kljun et al. 2015 App. B): stable B1 Nieuwstadt, neutral B2/B3 auto, convective B5 Batchvarova–Gryning slab daily cycle. `run_pblh_model()`, `R/pblh_*.R`; EC data only, constant γ 0.01 K/m ("PBLHmod (literature γ)"; kept by decision 3b, G 2026-11-11) | parameterisation |
| PBLHmod-clim | `pblhmod_clim` | PBLHmod with the site × month ERA5 γ climatology (`gamma_source = "era5_clim"`, spec 11; G 2026-10-30); constant γ where a month has no value (flagged) | parameterisation |
| PBLHmod, reference γ values (supplement) | `pblhmod_gamma001`, `pblhmod_gamma0025` | γ = 0.01 K/m (Batchvarova & Gryning 1991 as cited by Kljun et al. 2015; = `pblhmod_default` for now) and 0.0025 K/m (legacy ETC code); tables and supplementary figures only | parameterisation |
| ERA5 BLH | `era5_blh_center`, `era5_blh_bilinear` | ECMWF reanalysis diagnostic (bulk Ri = 0.25, IFS Cy41r2 §3.10.1); tower cell, or bilinear to the tower from 4 cells | model |
| Operational ICOS PBLH | `icos_oper` | `PBLH` column of L2 FLUXES, produced by the ETC pipeline with `legacy/calc_PBL_height_operational.R` + `legacy/proc_PBL_height.R`; `qc_flag` 1 = 50 m floor, 2 = run-mean fill. Used in the phase-1 regression only (spec 03), not in the evaluation (project plan "Scope", G 2026-10-23); script 04 compares it only with `--with-oper` | parameterisation |
| Legacy (spec 03) | `legacy_<variant>` | Port of the operational legacy code with switches, for the regression (in `data/derived`, not a standard file) | parameterisation |
| Measured MLH, Cloudnet (spec 07) | `cloudnet_eps_mlh` | Surface-connected mixing from the Doppler-lidar dissipation rate ε (ε ≥ 1e-4 m² s⁻³ upward from the lowest reliable gate; Manninen et al. 2018); half-hourly; `qc_flag` 0/2 exact, 1 censored below z_min (105 m at FI-Hyy), 3 signal-limited = lower bound; precipitation = Cloudnet `rain_detected` | turbulence |
| Measured MLH, SIRTA ReOBS (spec 10) | `sirta_mlh_chm15k` (primary), `sirta_mlh_cl31` | Ceilometer MLH at Palaiseau (CHM15k STRATfinder; CL31 CABAM v2019.3), hourly means from the SIRTA 1D ReOBS product; qc 0/1 from the sample count; FR-Gri, 23.37 km: daytime convective only | aerosol |
| Measured MLH, MeteoSwiss STRATfinder (specs 11, 13) | `meteoswiss_stratfinder_mlh` | CHM15k STRATfinder MLH (MeteoSwiss client files), half-hourly median of the quality-0 1-min values; CH-Dav (Davos, 0.98 km), DE-RuS (Jülich, 5.30 km), DE-HoH (Magdeburg, 24.48 km); E-PROFILE files: licence to be confirmed, internal use only | aerosol |
| Measured MLH, MeteoSwiss CABAM (spec 13) | `meteoswiss_cabam_mlh` | CL31 CABAM MLH (MeteoSwiss E-PROFILE client files), half-hourly mean of the quality-0 15-min steps (qc 0 both, 1 one); FR-Fon (Melun, 16.64 km), FR-Lus (Poitiers, 24.87 km); FI-Lom (Kenttärova, 1.80 km) inventory only since 2026-11-12 (no overlap with the EC record); licence to be confirmed, internal use only | aerosol |
| Secondary measured MLH, SkyVUE PRO (spec 14) | `tud_skyvue_mlh` | Campbell SkyVUE PRO built-in MLH (KNMI-based gradient method, 3 layers, confidence 1–3; TU Dresden), DE-Tha (113 m); half-hourly medians of each layer, `pblh_m` NA; **inventory since G 2026-11-12** (no daytime growth in any layer on the July 2026 test month, confidence 1 almost always); raw profiles kept for a later retrieval; internal use only | aerosol (inventory) |
| Secondary measured MLH, BL-View (spec 14) | `vaisala_blview_mlh` | Vaisala CL61 BL-View L3 hourly MLH (versions 3.1.1, 3.1.2, 3.2.2), CZ-Lnz (Lanzhot grassland, ~3 km); all quality values kept (`x_quality`); time zone not proven; internal use only; never pooled | aerosol (secondary) |
| Other measured PBLH (later) | – | PI data (SE-Htm, SE-Nor, SE-Svb, FI-Hyy, CZ-BK1 requested), ceilometer / radiosonde (E-PROFILE, IGRA) | aerosol / thermodynamic / turbulence |

## Specs

| Spec | Content | Phase | State |
|---|---|---|---|
| 01 | Site metadata, time-varying heights, exclusions, ICOS L2 readers, unit checks | 0 | merged (`main`) |
| 02 | PBLHmod v2 | 1 | merged (`main`, 58d45ad) |
| 04 | ERA5 readers (3×3 cells, pressure levels, static), standard PBLH format, first PBLHmod vs ERA5 comparison at BE-Bra (amendments 1–8) | 0, 2, 3 (first step) | merged (`main`, aba3394) |
| 05 | u* plausibility filter in the PBLHmod forcing (spec 02 amendment 12); local config override (`config.local.yml`, no cluster paths in tracked files) | 1 | merged (`main`, 25ae55c) |
| 06 | Consistency check of L2 `MO_LENGTH` against the forcing (option C, G 2026-10-21): L2 L kept where consistent (same sign, within 3×), otherwise replaced by L from the forcing and flagged (spec 02 amendment 13); checked only where the heat flux is positive (spec 06 amendment 1) and in daytime (amendment 2) | 1 | merged (`main`, 2026-10-21; amendment 2 merged 2026-10-21) |
| 03 | Legacy regression: reproduce operational PBLH from L2 (operational code `calc_PBL_height_operational.R` + `proc_PBL_height.R`), port it, attribute differences to PBLHmod v2 (amendments 1–2) | 1 | merged (`main`, e224c4b) |
| 07 | Cloudnet inventory (FI-Hyy, FI-Lom, DE-RuS, SE-Nor, FR-Gri); turbulence-based MLH from Doppler-lidar ε (Manninen et al. 2018) in the standard format; FI-Hyy runs (PBLHmod, ERA5) and first comparison against measured MLH (amendment 1: no operational PBLH, ERA5 tower cell in figures); reference period = FLUXES range; METEO TA levels from the header | 2, 3 | merged (`main`, b22f7e2) |
| 08 | Footprint sensitivity to the PBLH source (FFP, Kljun et al. 2015) at BE-Bra and FI-Hyy; one-month PC test (FI-Hyy 2022-09); ONEFlux BIF reader and u* thresholds | 5a | merged (`main`, 320f4ec); production run in spec 09; amendment 2 (DE-Tha m override 4.1, reference year 2025) on `spec08-m-override` (2026-11-10), to review |
| 10 | Aerosol-based measured MLH: Part A SIRTA ReOBS at FR-Gri (reader, two standard sources, FR-Gri runs 01/02/04, generalised spec-07 comparison, range-compression statistics also for FI-Hyy; review additions 2026-10-29: FB/nRMSD, symmetric slope + Bland–Altman, censored-tail null, PBLHmod-deficit diagnostics, NEE column check); Part B Lanzhot not started (two sample days only); Part C (Cloudnet ε at Palaiseau) only with G's go | 2, 3 | Part A merged (`main`, 63497d8) |
| 11 | γ from ERA5 profiles (`gamma_source` option, spec 02 Amendment 14); MeteoSwiss STRATfinder MLH at CH-Dav; CH-Dav runs 01/02/04; four-site evaluation of the γ cases; faded-context plotting convention | 1, 3 | merged (`main`, ccc6c08); decided: default unchanged (EC-only), ERA5 climatology = variant PBLHmod-clim |
| 09 | Cluster deployment. Part A (PC): config (BIF folder with glob, FLUXNET year spans, `pblh_gamma_default`), per-site driver with run status, SLURM templates and submit script, footprint chunks (c5 storage, chunk list, chunk mode, new cases, measured primary source), network collection and reference-year proposal, one figure style (`pblh_sources_style`, `theme_pblh()`). Part B (cluster runs) by G | 1–5 | Part A merged (`main`, 9656cd1); Part B: BE-Bra/FI-Hyy cluster regression passed (2026-11-02); first network run (job 11449) failed at FR-Bil, FR-Lus, UK-AMo, GL-ZaF (ERA5 0-360 longitudes) and SE-Sto (polar months, cli message): amendment 1 merged (`main`, 6328d4e). Network run on 39 stations: IT-SR2 not run (no L2 data in the L2I 2026 release; `config/sites_not_run.csv`). First collection (job 11546): amendment 2 merged (`main`, ebf4f71): network figure with the PBLHmod IQR bar, GL-ZaF NEE QC: `NEE_CUT_REF_QC` fallback where `NEE_VUT_REF_QC` is entirely missing, and reference-year exception (`config/reference_year_exceptions.csv`), measured steps `meas_std`, `meas_cmp`. Amendment 3 merged (`main`, c9eb5f3): reference year prefers a measured-overlap year at ≥ 80 % of the best year (FI-Hyy → 2022), γ strategy recorded, side-result constants `gamma_med` / `gamma_site` in step 02 (two-pass order), supplementary network figure. Reference years approved (G 2026-11-04, cluster collection 2026-11-04): `config/reference_years.csv` (39 sites, column `approved`), merged (`main`, 268bd34). Footprint run (jobs 11769/11796, 2026-11-05): 27 chunks failed (24 CZ-Lnz target-area mask from the first KML layer only, 1 mask read race, 2 worker ports); fixes merged (`main`, c654af9): mask from all KML layers (union), atomic mask write, worker port with retries, `scripts/08t_target_area_masks.R` (rebuild, invalidation, rerun list) |
| 13 | Measured MLH from MeteoSwiss / E-PROFILE at DE-RuS, DE-HoH (STRATfinder) and FI-Lom, FR-Fon, FR-Lus (CABAM, new source `meteoswiss_cabam_mlh`); pairing table with station code, WIGOS id, elevation difference, distance class, status; internal QC report (`scripts/13a_meteoswiss_qc.R`); generic spec-10 comparison for one half-hourly aerosol source; footprint pending list (`08c --pending`); measured-source exclusions (`config/measured_exclusions.csv`) and borderline flag; fix of step 08 at the sites without the ±20 % anchor cases; CZ-Lnz BL-View L3 reader prototype (pending data). Licence of the E-PROFILE files to be confirmed: internal use only | 2, 3, 5a | merged (`main`, 90c1da8); FI-Lom KEN row inventory since 2026-11-12 (no overlap) |
| 14 | Secondary measured sources: DE-Tha SkyVUE PRO built-in MLH (reader with CSV/TXT pairing, CRC-16 and profile-decode checks; `tud_skyvue_mlh`, 3 layers, no layer chosen) and CZ-Lnz BL-View L3 (all versions; `vaisala_blview_mlh`); Part C diagnostics (time-convention lag test, SkyVUE layers, BL-View quality); status `secondary`; explicit column types of standard files; FI-Lom KEN inventory; `09_network_record_start.csv` (ICOS labelling dates) | 2, 3 | merged (`main`, e88dda8); cluster runs by G (`docs/local/cluster_commands_spec14.md`) |
| 15 | Network footprint summary for the paper (target-area difference, fetch, top 10 % rows, Rebmann criteria and ICOS labelling rule, exceptions; figures and supplement) and γ supplementary tables, in `09_collect.R` | 5 | merged (`main`, 3e1625f; figure fixes fix-15, c672f8a) |
| 16 | Footprint uncertainty from the FFP inputs (phase 5b): one-at-a-time (A) and joint Monte Carlo (B) perturbation of σv, u*, L, umean, hc, z0, WD, h at six sites, attribution (C), the PBLH-source effect of spec 15 placed in it; `16a`–`16c` | 5b | merged (`main`, 573e9e2; amendment 1 rsl_off); cluster run done (six sites, reference years; collection job 18905); close-out figure fixes merged (`main`, 57a9460) |
| 19 | Why the PBLH source changes the footprint at CH-Dav and IT-Ren and not at other mountain sites: site descriptors on the spec-15 rows (terrain, tower, PBLH, turbulence, night-time flows, footprint, spec-15 status), night-row R2, existing night-time measured evidence, supplementary figures `19s_*`; `scripts/19a_site_descriptors.R`, `R/alpine_drivers.R`; analysis only (no PBLHmod / FFP run) | 5 | `spec-19-alpine-exception` (2026-11-25); local test on BE-Bra, FI-Hyy, FR-Gri, CH-Dav (row basis flags_only, no local chunks); cluster job pending (`docs/local/cluster_commands_spec19.md`) |
| 17 | γ from the station's own data (EC-derived γ, PBLHmod-γEC): daily γ from the tower θ rise in the B&G window (Method A, model-consistent; Method B, closed-form encroachment), day selection by flags, site × month climatology; `gamma_source` `ec_clim` / `ec_daily` (spec 02 amendment 15); cases `gamma_ec` (standard source `pblhmod_gamma_ec`), `gamma_ec_daily`, constant grid `gamma_g002`–`gamma_g012`; step 04 pair 9; evaluation `11 --prefix=17` (A, B), `17b` (C, D, E); `R/gamma_ec.R`, `scripts/17a_gamma_ec.R`, `scripts/17b_gamma_ec_collect.R` | 1, 3 | `spec-17-gamma-ec` (2026-11-25); local test on BE-Bra, FI-Hyy, FR-Gri, CH-Dav; D1 pre-test (2026-11-26, `spec-17-d1-pretest`, option `gamma_ec_start_offset_h`): **D1 rejected**, no cluster run (note at the top) |

Phases (see `project_plan.md`): 0 done; 1 done (model, input QC, legacy regression); 2 standard format, three
sources (spec 04) and the first measured source (spec 07, FI-Hyy) done; 3 BE-Bra (vs ERA5) and FI-Hyy (vs ERA5 and
measured MLH) done; 4–6 not started.

- Tests: 1,483 expectations passing (`spec-09-cluster`, 2026-10-30). Environment: R 4.4.1 + renv (Posit Package
  Manager), same as the cluster. Machine-specific paths in the untracked `config.local.yml` (template
  `config.local.yml.example`).
- PBLHmod defaults: `stable_method = "nieuwstadt"` (B1), `neutral_method = "auto"` (B2 if 10 < N/|f| < 70, else
  B3), `combine = "bg_only"`, `gamma = 0.01` K/m, `h_source = "H_F_MDS"` with `h_qc_max = 1`; u* input filter
  `ustar_ws_max = 0.5` (where WS ≥ `ws_min_ratio_test = 1` m s-1), `ustar_max = 3` m s-1 (spec 02 amendment 12;
  u* and L set to NA, short gaps interpolated); L check `l_source = "L2_checked"` (`l_consistency_factor = 3`,
  daytime rows with wtheta0 > `l_check_wtheta_min = 0.005` K m s-1 only, `l_check_daytime_only = TRUE`;
  inconsistent L2 L replaced by the forcing L; spec 02 amendment 13, spec 06 amendments 1–2); raw output (no smoothing, clipping or filling). Forcing T and p from
  measured meteo (BE-Bra: FLUXNET `TA_F`/`PA_F` 2020-01 → 2026-08-22, METEO `TA`/`PA` after; before 2026-10-26
  FLUXNET 2020–2025 and METEO for all of 2026).
- BE-Bra PBLHmod (default, FLUXNET 2020–2026, 2026-10-26): convective median 707 m (p5–p95 75–1581), stable 257 m,
  all rows 414 m; γ = 0.0025 raises the convective median to 1122 m. With FLUXNET 2020–2025 (spec 06 amendment 2):
  690 m (74–1575), 259 m, 412 m; 1095 m. Spec 06 amendment 1: 690 m (74–1575), 260 m, 413 m;
  1095 m. Spec 05: 698 m (72–1586), 261 m, 413 m; 1104 m. Spec 02: 704 m (72–1609), 271 m, 423 m; 1109 m.
- Spec 10 Part A at FR-Gri (2026-10-29; PC RAM test passed 2026-10-29, runs single-process):
  - Inputs: reference period (FLUXES range) 2021-01-01 → 2026-08-31 local, but valid L2 data from 2021-11-17;
    FLUXES has NEE / NEE_DATA_FLAG / SH (an earlier note said they were missing: a misread of the script-01 table,
    corrected 2026-10-29; 128 rows with NEE_DATA_FLAG = 0 and NEE missing); zm 2–3 m, hc 0.06–0.97 m (0.2 m after the
    harvest used as is, G); no exclusions. ERA5 tower cell 11.05 km, orography 127.5 m (tower 125 m), sdor 25 m.
  - PBLHmod default: 81 % valid; convective median 415 m (p5–p95 55–1134; BE-Bra 707, FI-Hyy 700), stable 104 m;
    u* rejected 0.14 %; L replaced 1,834; B&G stopped at a gap on 19.9 % of days; 56 rows below zm. Operational
    PBLH: run-mean fills 781.37 m (2021–2025) and 808.35 m (2026) on 33.4 % of the non-missing values (day; with the
    night fills 214.42 / 185.67 m 39.3 %; ETC QC report items 12–13).
  - PBLHmod vs ERA5 (spec 04, 40,212 h): median log ratio −0.84 (convective −0.54, stable −1.24), 39 % within a
    factor 2, Spearman 0.76 (BE-Bra −0.14 / 72 % / 0.78); EC H vs ERA5 H bias +1.5 W m-2 (not a forcing problem).
  - SIRTA ReOBS facts of spec 10 A2 verified (details in `CLAUDE.md`). Comparison period 2021-11-16 23:00 →
    2026-06-30 19:00 UTC; window = PBLHmod convective (both half-hours) and apparent solar time at the interval
    midpoint 10–16 h. Window metrics (n, median ln ratio, within a factor 2, Spearman):
    PBLHmod vs CHM15k 6,531, −0.53, 66 %, 0.72; ERA5 vs CHM15k +0.002, 93 %, 0.80; PBLHmod vs CL31 4,691, −0.41,
    71 %, 0.73; ERA5 vs CL31 +0.06, 89 %, 0.80; CL31 vs CHM15k 4,295, −0.018 (MAD 0.10), 94 %, 0.86 (definitional
    spread small in the window; CL31 lower in the evening decay); PBLHmod vs ERA5 in the window −0.54, 63 %.
    PBLHmod low in every season (−0.47 JJA to −0.60 MAM vs CHM15k); JJA diurnal medians at 14 UTC: CHM15k 1,575,
    ERA5 1,431, PBLHmod 931 m.
  - FB / nRMSD (`Chang2004`, `Sinclair2022`; FB > 0 = observation deeper), window: PBLHmod vs CHM15k +0.53 / 0.73,
    ERA5 vs CHM15k −0.01 / 0.31, PBLHmod vs CL31 +0.41 / 0.64, ERA5 vs CL31 −0.09 / 0.35, CL31 vs CHM15k +0.08 / 0.29.
  - Compression (G 2026-10-30): **supported at FI-Hyy** by the censored-tail test (main evidence: beyond the
    scatter-only bound, one-sided, conservative); **mild at FR-Gri** (symmetric slope 0.89 in the window). The
    symmetric slope and Bland–Altman trend (day-block bootstrap 95 %) are supporting evidence only: source and
    measurement differ in definition, not only in noise (equation error, `Warton2006`); Theil–Sen diagnostics only.
    Window: PBLHmod vs CHM15k 0.89 [0.86, 0.93], trend −0.14
    [−0.18, −0.09]; ERA5 vs CHM15k 0.94 [0.91, 0.98], −0.07; PBLHmod vs CL31 0.74 [0.71, 0.78]; ERA5 vs CL31 0.84
    [0.80, 0.88]; CL31 vs CHM15k 1.16 [1.12, 1.21]. All convective hours: > 1 for every pair (PBLHmod vs CHM15k
    1.17): the sign depends on the stratum (CHM15k does not resolve the shallow morning and evening layers below
    its overlap height of a few hundred metres). FI-Hyy (exact values): symmetric slopes > 1 (convective PBLHmod 1.29
    [1.19, 1.40], ERA5 1.26); censored-tail null (`10_censored_null_FI-Hyy.csv`): convective, MLH < 105 m: PBLHmod
    ≤ 105 m on 14 % vs ≥ 25 % under scatter alone (ERA5 13 % vs ≥ 33 %); MLH a lower bound: PBLHmod below it on 40 %
    vs ≤ 30 % (stable 43 % vs ≤ 29 %), ERA5 33 % vs ≤ 26 % from the stable regime (37 % vs ≤ 22 %; convective 31 %
    vs ≤ 35 %). Details in `project_plan.md` open points.
  - PBLHmod deficit (`scripts/10b_pblhmod_deficit.R`, existing outputs only): γ = 0.0025 case window median ln ratio
    +0.08 (90 % within ×2), max_B2 −0.09, all other cases −0.51 to −0.54; h0 median 88 m (JJA 57, DJF 135; 85 % from
    B1), ERA5 at the same time 231 m; clear days (341, kt ≥ 0.65): PBLHmod close to CHM15k/ERA5 until 8 h solar time,
    then grows too slowly; 1.2 % of window hours on days where B&G stopped (51 of 1,541 B&G days; 79 of 1,687 days).
    **Decided (G 2026-10-30): γ explains the deficit (γ = 0.0025: −0.53 → +0.08); start height and gap stops ruled
    out; clear-day definition (kt ≥ 0.65, SW_IN coverage ≥ 90 %) confirmed. Next: spec 11 (γ from ERA5 profiles).**
  - CL31 < 150 m with CHM15k > 300 m: 95 of 4,295 window convective hours (2.2 %).
  - Refactor check: FI-Hyy script 07 outputs (6 tables, 4 figures, console report) byte-identical after moving the
    matching into `measured_time_base()`. After adding `fb`/`nrmsd` (2026-10-29): 5 tables and 4 figures
    byte-identical; `07_metrics_FI-Hyy.csv` keeps its 1,441 rows unchanged and in order, plus 262 `fb`/`nrmsd` rows.
- Spec 11 (2026-10-30, `spec-11-gamma-era5`):
  - CH-Dav inputs: reference period 2018-12-31 → 2026-08-31 (134,400 half-hours), no exclusions; zm 35 m, hc 19–20 m;
    ERA5 tower cell 10.85 km, orography 1999 m (tower 1637 m, −362 m), sdor 366 m; EC H vs ERA5 H bias −18.6 W m-2.
    PBLHmod default 86 % valid, convective median 483 m (γ 0.0025: 810 m), stable 168 m; B&G stopped at a gap on
    14.6 % of days. PBLHmod vs ERA5 (all regimes) median ln ratio +0.84 (convective +0.09). 850 hPa lies below the
    ground at CH-Dav, so the spec-04 850–700 hPa γ table is not meaningful there.
  - MeteoSwiss STRATfinder (`meteoswiss_stratfinder_mlh`): 93,344 half-hours, 72,488 with a value (qc 0 71,309);
    79–89 % valid per year, 2025 50 % (window contamination); night floor 262 m; midday medians DJF 741 → JJA
    1698 m. Facts and three corrections to the spec (first profile 08:20, 59–61 s steps, bit-4 spells to 2026-02)
    in `CLAUDE.md`.
  - ERA5 γ (500–3500 m above the site): medians BE-Bra 0.0047, FI-Hyy 0.0046, FR-Gri 0.0048, CH-Dav 0.0044 K m-1
    (≈ 0.47 × default); 2 levels on 0.7 % (BE-Bra), 28 % (FI-Hyy), 4 % (FR-Gri), 100 % (CH-Dav, 700/600 hPa) of days;
    implausible < 0.2 %; `zmin1000`: 1 level only on 71 FI-Hyy days (no γ), identical to the default layer at CH-Dav.
  - `gamma_source = "constant"` reproduces every existing case at BE-Bra, FI-Hyy, FR-Gri (all 37 old columns of the
    written files identical; the automatic check flagged only FI-Hyy `nofilter` on L, a CSV round-off at
    |L| = 2.3 × 10^10 m). 07, 10, 10b tables byte-identical (10b gains rows for the new cases); 07 FI-Hyy figures
    byte-identical; 10 FR-Gri figures regenerated with the faded context.
  - Evaluation (`scripts/11_gamma_evaluation.R`, 11.6 min): see `project_plan.md` open points (spec 11). Window bias
    FR-Gri −0.53 → −0.20 (ERA5 clim), CH-Dav −0.43 → −0.08; γ 0.0025 +0.08 / +0.17. γ in N affects only the neutral
    B3 heights (no B2 rows at any site), not the window.
- NEE columns for spec 09 (`scripts/10c_fluxes_nee_columns.R`, header check + value counts; run on the cluster for
  all 40 stations of `config/sites.csv`): BE-Bra, FI-Hyy, FR-Gri have FLUXES NEE + NEE_DATA_FLAG and FLUXNET
  NEE_VUT_REF_QC.
- `config/sites.csv` (fixed 2026-10-30): the last line (UK-AMo) had 12 fields instead of 13 and no final newline;
  after lines ending with an empty `notes` field, readr dropped it without reporting a problem, so `read_sites()`
  returned 39 stations since the file was created (2026-10-08). Fixed by completing the line (format only, no value
  changed); `read_sites()` now aborts when the number of rows read differs from the number of data lines. 40
  stations = 40 target-area KMLs. No all-station run had been made, so only the NEE check above was affected.
- Cloudnet re-check, Hyytiälä (2026-10-29, by instrument PID, `showLegacy=true`, no legacy files returned): HALO 46
  (21.12132/3.421272f219be4f97): `lidar` 2014-02 → 2021-08 (gaps 2014-04 → 2016-11, 2017-11 → 2018-01),
  `doppler-lidar` / `epsilon-lidar` / `doppler-lidar-wind` 2022-01 → 2024-06-13 (epsilon to 2024-04-11), then at
  Sammaltunturi from 2026-09-21. UH CL61 (uuid 241bda14-…, PID 21.12132/3.241bda142975460b): `lidar` 2021-08-23 →
  2026-09-08 at Hyytiälä, Sammaltunturi from 2026-09-23.
- Spec 08 at BE-Bra and FI-Hyy (2026-10-26; **test, one month, September 2022**; no definitive conclusions):
  - Verification (`scripts/08v_ffp_verify.R`): A (paper's equations) = `calc_footprint_FFP()` on 1,000 half-hours
    per site (max rel. diff 5e-16); C (summed B fields) = `calc_footprint_FFP_climatology()` (identical, smoothing
    off; EBImage not installed, so `ffp_smooth()` is verified against an FFT circular convolution only).
  - Coverage (full period, 08a): `kljun` set 0 valid rows at both sites (20 z0 < zm never holds); ETC multiplier m
    = 5 (BE-Bra), 3.7 (FI-Hyy). Test month valid (`none`): BE-Bra 1,202 of 1,440 half-hours, FI-Hyy 1,003.
  - ERA5 vs PBLHmod, ln(x_80 ratio), subset `valid`, set `none`: median −0.020 (BE-Bra) / −0.008 (FI-Hyy); rows
    with |x_80 change| > 10 %: 30 % / 9 %; > 25 %: 14 % / 2 %. Mostly stable nights (BE-Bra stable: > 10 % on 43 %,
    convective 12 %); large where h approaches zm. In `valid_qc_used` the differences almost vanish (BE-Bra stable
    > 10 % on 13 %): 97–100 % of the top-10 % affected half-hours are removed by the NEE QC or the u*/gap-fill
    subset. Target-area fraction: median +0.4 pp (BE-Bra), p90 |diff| 9.9 / 2.6 pp.
  - Monthly climatologies (common rows): Jaccard of the 80 % regions ERA5 vs PBLHmod 0.77 (BE-Bra; night 0.70) and
    0.93 (FI-Hyy); 0.97 / 0.99 in `valid_qc_used`. The BE-Bra night `valid` region reaches the ±1000 m domain edge.
  - Rebmann et al. (2018) ≥ 90 % from the target area on ≥ 50 % of half-hours: fails at both sites with every h
    (≤ 1.5 % of half-hours; target areas 0.115 / 0.105 km²); no verdict changes with the PBLH source.
  - Anchor (h ± 20 %): median |Δx_peak| unstable/neutral/stable 0.6/0.7/3.4 % (FI-Hyy, ×0.8), Kljun et al. (2015)
    Table 5: 0.3/0.7/3.7 %; BE-Bra 1.0/1.3/6.8 %.
  - FI-Hyy measured MLH (297 common half-hours): ln(x_80 PBLHmod / measured) median −0.008, > 10 % on 13 %.
  - Runtime: 1.7–2.1 s per half-hour (single), ≈ 4.2 s per worker with 4 PC workers; test month 2.0–2.2 h per site.
    Full period ≈ 460–480 k FFP half-hours per site ≈ 250–270 CPU-hours per site: cluster.
  - Final test-month chunks (hash v3, 2026-10-27, 2 workers) and integrity checks (2026-10-28). The PC crashed four
    times during these runs and once silently flipped data in memory: one float32 cell of one summed group of the
    FI-Hyy `era5` `_sums.rds` (rel. 2.9e-6) differed between two runs; a third computation matched the earlier file,
    which was kept. **Treat this PC as able to corrupt data in memory; no heavy runs on it until the RAM is tested.**
    Checks: every chunk of both sites read back (hash = current task, rows/timestamps/h = task, no non-finite or
    negative cells, group counts consistent); climatologies rebuilt from the float32 sums match the stored ones (max
    rel. diff ≤ 3e-6, Jaccard 80 % = 1; all 288 of BE-Bra; FI-Hyy `bound`, `measured`, `pblhmod_x0.8`/`x1.2`, 180);
    FI-Hyy `pblhmod`, `era5`, `era5_bilinear`, `bound` byte-identical across two independent runs (hh, clim, sums). Script 08 run twice per site: all 22 outputs
    (tables and figures) byte-identical (SHA-256). The headline numbers above are unchanged by the final chunks.
  - Storage options for the cluster sums (`scripts/dev/08_storage_options.R`, for spec 09; tables
    `08_storage_*.csv`): current format (float32 + gzip) 54 GB (FI-Hyy, 7 cases) / 33 GB (BE-Bra, 5 cases) per
    site-year with all sets, 18 / 11 GB with the `etc` set only; 5 m block sums (`c5`) 2.2 / 1.3 GB (all sets), within
    every tolerance (|Δ area_80| ≤ 0.16 %, |Δ xmax_80| ≤ 1.1 m, |Δ f_target_area| ≤ 0.008 pp, Jaccard vs full
    ≥ 0.991, |Δ Jaccard ERA5 vs PBLHmod| ≤ 3e-4); sparse 1 − 1e-3 / 1 − 1e-4 pass but save little (42–49 / 30–33
    GB); `c5_sp4` 2.0 / 1.3 GB. Decision on the storage design for spec 09: G.
- BE-Bra FLUXNET 2020–2026 replacement (G, 2026-10-26; old outputs in `*/pre_bebra_fluxnet2026/`):
  - Inputs, old vs new file on 2020–2025 (105,216 common rows): only gap-filled values differ. `TA_F` 14,318 rows (all
    QC 2, 99 % in 2020; max 0.033 K), `PA_F` 15,530 rows (all QC 2; max 28 Pa), `H_F_MDS` 69 rows (all gap-filled,
    2025-12-26 → 12-31; max 63 W m-2), `NEE_VUT_REF_QC` 754 rows (0.8 %); NA pattern unchanged; `H_CORR` differs on
    48 % (not used). 2026: FLUXNET `TA_F`/`PA_F` identical to METEO where both are valid; FLUXNET ends 2026-08-22
    22:30 UTC.
  - PBLHmod (13 cases, scripts 01, 02, 04, 03 rerun): before 2025-12-26 max |diff| ≤ 0.16 m (gap-filled TA/PA);
    2025-12-26 → 12-31 up to 534 m (default; 16 rows > 1 m, all on the days with changed `H_F_MDS`); 2026 changed
    because H now comes from `H_F_MDS` instead of L2 `H` (default: valid H rows 6,318 → 10,553, convective 1,860 →
    4,337, B&G stopped at gaps 8,959 → 2,148, convective median 594 → 783 m). `l_forcing` gains 3,387 valid 2026 rows.
  - Pair 1 vs ERA5 (49,400 h): RMSE 337.5 → 333.5 m, bias −49 → −42 m, median log ratio −0.153 → −0.146, within a
    factor 2 72.6 → 72.8 %, Spearman 0.797 → 0.801; convective n 15,242 → 16,442.
  - Spec 03: the reproduction of the operational PBLH is unchanged (all and per-year strata identical; regime strata
    change because regimes come from PBLHmod); variants with our forcing (`legacy_forcing`, `cum_06`–`cum_12`) change
    as expected. The `cleanset` variant is no longer produced (cleanset moved to `data/archive/`; old result kept in
    the backup).
  - ONEFlux u* check now runs at BE-Bra: VUT thresholds 0.18–0.28 m s-1 (CUT 0.2345); 100 % of u*-filtered rows
    gap-filled; u* explains 80.7 % of the QC-good NEE gap-filled in FLUXNET (day 71.7 %, night 84.2 %; FI-Hyy 74 %).
- Spec 06 amendment 2 at BE-Bra (2026-10-21): L check in daytime only. Rows checked 26,198; replaced 2,191 (sign
  1,170, magnitude 1,021; convective 1,720, unstable_diag 471), night replacements 285 → 0 (all in Nov 2022–Jan
  2023, positive gap-filled H_F_MDS with L2 H QC-rejected); `flag_sign_inconsistent` 335; values > 5 km 1.
  Regression chain: `nofilter_l2` = spec 02, `l_L2` = spec 05, `l_check_allday` = amendment 1 (all < 1e-9 m).
  Pair 1 vs ERA5 practically unchanged vs amendment 1 (RMSE 337.6 → 337.5 m; unstable_diag RMSE 449 → 454 m).
- Spec 06 at BE-Bra (2026-10-21, amendment 1 = V1): 26,652 rows checked, 2,476 replaced (sign 1,349, magnitude
  1,127; 69 % convective, 31 % unstable_diag; 88 % daytime; on 956 days; winter clusters, e.g. 2023-01-19 → 01-28); L2 L used
  on 96.5 % of rows with L; `flag_sign_inconsistent` 1,228 → 335. PBLHmod > 5 km: 11 → 1 (2024-01-05 09:30, B1,
  u* 2.67 m s-1: spec 05 filter limitation). Regression chain: `nofilter_l2` = spec 02 and `l_L2` = spec 05
  (max |diff| ≤ 4e-12 m). Pair 1 vs ERA5, spec 05 → 06: RMSE 349 → 338 m, Pearson 0.725 → 0.737, within a
  factor 2 72.4 → 72.6 %, convective RMSE 385 → 353 m, stable unchanged (310 → 309 m).
- Spec 05 at BE-Bra (2026-10-21): the filter rejects 2,149 of 100,186 u* values (2.1 %: ratio 1,328, cap 46,
  both 775); 618 then interpolated, 1,531 left NA. PBLHmod > 5 km: 179 → 11; > 10 km: 70 → 0. `nofilter`
  reproduces the spec-02 default (max |diff| 3.6e-12 m). Pair 1 vs ERA5: RMSE 623 → 349 m, Pearson 0.49 → 0.73,
  bias −9 → −46 m, median log ratio −0.136 → −0.151, within a factor 2 71.8 → 72.4 %.
- Spec 07 at FI-Hyy (2026-10-23):
  - Inputs: FLUXES/FLUXNET/METEO 2018-01 → 2026-08 (data from 2018-11-06); unit checks pass; METEO TA_1–TA_9.
    PBLHmod default: 80 % valid, convective median 700 m (stable 315 m), u* rejected 1.8 %, L replaced 2,013, none
    > 5 km; 2,811 values below zm = 27 m (stable nights 2,282, early-morning B&G starts 430). ERA5 tower cell 11.1 km
    from the tower, 42 m lower. PBLHmod vs ERA5 (2020-01 → 2026-08, 51,584 h): median log ratio −0.10 (DJF +0.16,
    JJA −0.35), 73 % within a factor 2, Spearman 0.75 (BE-Bra −0.14 / 72 % / 0.78). Operational fills
    3316.87/2000.85 m (2018-11-06 → 2025-12-31) and 3365.22/1532.14 m (2026-01-01 → 2026-08-26), 28 % of rows.
  - Measured MLH (Cloudnet `epsilon-lidar`, 688 days 2022-01 → 2024-04, 390 MB; classification 682 days): z_min
    105 m agl; profiles: measured 67 %, censored ε < 1e-4 at 105 m 22 %, 105 m gate missing 6 %, no signal 1 %,
    precipitation (`rain_detected`) 3.5 %; half of the measured profiles end at a gate without signal (clean air).
    Half-hours (30,727): qc 0 3,171, qc 2 6,870, qc 3 (lower bound) 11,889, qc 1 (censored) 7,365, NA 1,432. JJA
    daytime: qc 3 on 54–62 % of half-hours, so the censored median is not estimable there.
  - Comparison on exact values (qc 0/2), convective regime: PBLHmod n 2,703, bias +202 m, median log ratio +0.44,
    56 % within a factor 2, Spearman 0.47; ERA5 tower cell n 880, +111 m, +0.32, 69 %, 0.53 (bilinear practically
    identical). Exact values over-represent shallow layers. Where the lidar MLH is below 105 m, PBLHmod is ≤ 105 m on
    26 % (ERA5 24 %); where it is only a lower bound, PBLHmod is below the bound on 40 % (JJA 62 %; ERA5 33 %, JJA
    45 %): "shallow too deep / deep too shallow", candidate finding. Outputs with the operational PBLH kept in
    `output/*/reference_spec07_with_oper/` (and `output/tables/reference_spec04_with_oper/`).
  - Sensible-heat timing (`scripts/07d_h_timing.R`, 2020–2025, values at interval midpoints): ERA5 H earlier than
    EC H, best-fit lag FI-Hyy 25–37 min (JJA 37), BE-Bra 16–25 min; FI-Hyy JJA afternoon crossings ≈ 63 min earlier,
    morning ≈ 15–18 min (asymmetric). Adding the L2 storage term `SH` (air below the EC system, ±12 W m-2 at
    FI-Hyy) reduces the lag (FI-Hyy 12–27 min, BE-Bra −1 to 14 min) but leaves the FI-Hyy afternoon difference at
    51–64 min. Diagnostic only (open point in `project_plan.md`).
- Spec 04 at BE-Bra (2026-10-17): overlap 2020-10-25 → 2026-08-25; PBLHmod default vs ERA5 tower cell 50,269 matched
  hours, bias −9 m, median log ratio −0.14, 72 % within a factor 2, Spearman 0.78; stable −0.23, convective −0.07;
  JJA evenings (19–22 UTC) PBLHmod ≈ half of ERA5. Operational PBLH vs ERA5: median log ratio +0.57 (qc 0).
  Open points from spec 04 in `project_plan.md` ("Open points").

## Next steps

1. ~~Review and merge spec 04 (`spec-04-era5-compare`).~~ Done (2026-10-21, merged at aba3394; `mad_log_ratio`
   switched to the scaled MAD at review, amendment 8).
2. ~~Review and merge spec 05 (`spec-05-ustar-filter`).~~ Done (2026-10-21, merged at 25ae55c).
3. ~~Spec 06 — consistency of L2 `MO_LENGTH` with the forcing.~~ Done (2026-10-21). Why option C (G): at BE-Bra,
   10 of the 11 values > 5 km left after the u* filter came from B5 on days where L2 H is QC-rejected but
   `MO_LENGTH` (from H_UNCLEANED ≈ 0, |L| ≈ 3e6 m) is kept. A (L NA where L2 H is NA) would stop B5 early on 959 of
   1,770 convective days, because L2 H is QC-rejected on 47 % of rows with L; B (inconsistent L set to NA) on 268
   days; D (L always from the forcing) is kept as a sensitivity case (`l_forcing`). Checking stable rows too made
   results worse (spec 06 amendment 1).
4. ~~Review and merge spec 03 (`spec-03-legacy-regression`).~~ Done (2026-10-22, merged at e224c4b). BE-Bra results:
   - Part A: the operational code on inputs rebuilt from L2 (EddyPro engine formulas for AT, AP, rho, cp) reproduces
     the operational column: median |rel. diff| 0.009 %, p95 0.06 %, fill values within 0.006 %, 1,341 of 1,342 fill
     blocks; the local cleanset gives 0.0004 % (remaining differences = rounding of the published inputs). Runs:
     2020-10-26 → 2026-01-01 and 2026-01-01 → 2026-08-26 (local standard time). The development version is out of
     scope (G, 2026-10-22).
   - Part B: port = original (max |diff| 0 m on the full series, < 1e-8 m on the synthetic test).
   - Part C, on common rows (valid in every variant and in PBLHmod v2; spec 03 amendment 2 note 9). One at a
     time (n = 79,263; baseline median 707 m, PBLHmod v2 364 m): neutral formula −211 m (convective −740 m: the
     operational unstable bound −0.003 puts many daytime rows in the neutral class, 0.3 u*/|f| ≈ 1.3 km), smoothing
     −192 m, gap interpolation +176 m, stable formula −38 m, gamma −11 m (convective −101 m); fill, floor, mol_locf,
     forcing ≈ 0. Cumulative (n = 75,689; 724 → 366 m): smoothing −206 m, neutral formula −82 m (convective −970 m),
     stable formula −46 m; the last step equals pblhmod_default (residual 4.6e-13 m). On all valid rows the fill
     (−383 m) and forcing (±160–190 m) effects are row-set changes, not value changes. `at_meteo` (sonic temperature
     offset): negligible.
5. Spec 07 — implemented on `spec-07-measured-cloudnet` (2026-10-23), to review and merge. Scripts: 07a inventory,
   07b download (size report first, `--confirm`; `--match=PRODUCT`), 07c MLH (profiles in `data/derived/cloudnet/`,
   standard file, sensitivity columns), 07 comparison, 07d sensible-heat timing; 04 with optional `--with-oper`.
   Diurnal figures place every value at its interval midpoint (hour bins centred on the full hour). Distance rule and censoring
   convention in `project_plan.md`; Cloudnet facts in `CLAUDE.md`; example in `docs/specs/standard_format.md`.
6. After spec 07: the other Cloudnet sites (FI-Lom, DE-RuS, FR-Gri; on the cluster), PI data into the standard
   format, ETC QC report. Before paper results: re-check the volatile Cloudnet files and create a Cloudnet collection
   (open point in `project_plan.md`).
7. Spec 13 (2026-11-06, merged `main`, 90c1da8): cluster runs (meas_std / meas_cmp at the five E-PROFILE sites,
   08a, `08c --pending`, pending chunks, step 08 at all 39 sites after the anchor fix, collection); measured-source
   exclusions approved and applied (G 2026-11-06, `config/measured_exclusions.csv`: Melun six dirty-window months,
   Poitiers 2021-02 → 2021-10); borderline flag (DE-HoH, FR-Lus) in the comparison tables; ask MeteoSwiss / the owners before sharing
   any result from these files. CZ-Lnz: BL-View data received (spec 14).
8. Spec 08 amendment 2 (2026-11-10, merged `main`, 7fc99a6): ETC z0 multiplier override m = 4.1 at DE-Tha
   (`config/etc_z0_multiplier.csv`; loop value 4.6 rejects every row from 2024-06-11), `m_loop` / `m_used` /
   `m_source` in the 08a table and `09_network_etc_multiplier.csv`; DE-Tha reference year 2025 (overlap with the
   requested TU Dresden SkyVUE ceilometer data). Cluster: 08a at DE-Tha, check 2025 ≥ 9 usable months, then
   `08c --pending`, the DE-Tha 2025 chunks, step 08 at DE-Tha, collection (`docs/local/cluster_commands_detha_m.md`).
   Methods note in `project_plan.md` (phase 5).
9. Spec 11 amendment 1 (2026-11-11, `spec-11-gamma-network`, to review): γ evaluation at the 8 measured sites +
   BE-Bra for decision 3b (cases `default`, `clim`, `gamma_med`, `gamma_site`, `gamma0025`; decision table and
   figure; evidence only, default stays 0.01); FI-Lom `meas_cmp` failure (empty stratum set) fixed: empty sets write
   `n = 0` with `no_data_reason`. Cluster: one job for the evaluation, `meas_cmp` at FI-Lom
   (`docs/local/cluster_commands_gamma_eval.md`); then G and Cowork decide 3b. Merged; the cluster run stopped at
   FI-Lom: **no overlap** (EC values from 2024-11-04 22:00 UTC, KEN CABAM ends 2024-05-16).
10. Follow-up fix (2026-11-11, `fix-11-empty-overlap`, merged `main`, 63450a8): an empty or inverted overlap gives
   n = 0 rows with the reason (scripts 10 and 11); in script 11, one site never stops the job (`status`,
   `no_data_reason` in `11_gamma_decision.csv`; exit 1 if any site is `error` / `error_growth`). Cluster: range check
   of every PBLHmod case file, rerun of script 11 and of FI-Lom `meas_cmp`
   (`docs/local/cluster_commands_filom_overlap.md`).
11. **Decision 3b (G, 2026-11-11; spec 09 Part B 3b):** the EC-only constant stays γ = 0.01 K m⁻¹
   (`pblhmod_default`, "PBLHmod (literature γ)"; no rerun); PBLHmod-clim is the realistic-γ variant in the main
   figures; `gamma_med` = 0.004572 K m⁻¹ is the recommended operational EC-only value (paper and ETC), a side result.
   Evidence: script 11 on the cluster (job 18120): median |ln ratio| at the six aerosol sites 0.43 (default) → 0.20
   (`gamma_med`; `clim`, `gamma_site` within 0.02), overshoot ≈ +0.2 at DE-HoH and FR-Fon; vs ERA5 across the
   network −0.13 → +0.15. Open for later: does the residual site spread follow the energy-balance ratio?
12. **EC record start = ICOS labelling date (G, 2026-11-12):** each site's EC record starts at its ICOS labelling
   date, the first valid L2 EC value; the FLUXES rows before it (from 1 January, local) are empty grid rows. Per-site
   table `09_network_record_start.csv` (09_collect, from the first valid `pblhmod_default` value) for the supplement.
   FI-Lom: labelling date 2024-11-05, after the end of the Kenttärova CL31 record (2024-05-16): the MeteoSwiss KEN
   row is `inventory`, FI-Lom left `config/measured_run_sites.txt` (7 measured sites). Kenttärova ε-lidar (Cloudnet
   inventory 2026-10-22): 2022-08-23 → 2023-12-15 and from 2026-04-03 (FMI HALO 34, volatile).
13. Spec 14 (2026-11-12, `spec-14-skyvue-blview`, to review): secondary sources DE-Tha SkyVUE PRO
   (`tud_skyvue_mlh`) and CZ-Lnz BL-View (`vaisala_blview_mlh`); sample checks pass (SkyVUE day 2026-07-15: 10,063
   records, CRC, layers and profile decode exact; BL-View: 3.1.1, 3.1.2, 3.2.2 samples). Open for G: SkyVUE layer
   rule; BL-View quality rule (meaning of 1–3 not documented, indicative 3 = best), time convention (lag test),
   version overlap rule (newest kept). Explicit column types in `read_standard_pblh()` (the CABAM `x_dlh_m` was lost).
   Cluster: `docs/local/cluster_commands_spec14.md`. Merged (`main`, e88dda8).
14. Spec 15 (2026-11-12, `spec-15-footprint-summary`, to review): network footprint summary for the paper in
   `09_collect.R` (target-area difference, fetch ratio, top 10 % rows, Rebmann criteria and ICOS labelling rule,
   exceptions > 5 pp; figures `09_fp_*`, supplement `09s_fp_*`) from the per-half-hour chunk metrics of the
   reference years (no footprint rerun); γ supplementary tables `09s_gamma_site.csv`, `09s_gamma_network.csv`.
   Checks against the per-site tables: top 10 % (step-08 selection, `etc`) 71 % stable, 64 % night, 64 % removed as
   expected; the expected 93 % removed is `kljun_norsl`, not `kljun` (84 %). Cluster: the collection only
   (`docs/local/cluster_commands_spec15.md`). Merged (`main`, 3e1625f).
15. fix-14c (2026-11-12, `fix-14c-diagnostics`, to review; spec 14 amendment 1): 14c lag table fixed (data masking;
   one row per site × series × reference × season), lags −3 … +6 h, SkyVUE window metrics per layer, lag figure
   grouping; **SkyVUE built-in MLH set to inventory (G 2026-11-12)**, reader / standard source / profile decoder kept
   (raw profiles = DE-Tha asset for a later retrieval); option `--blview-322-tz=Europe/Prague` (diagnostics only,
   suffix `_tz322`) to test the hypothesis that 3.2.2 stamps are local time with summer time; the +1 h baseline is not
   a matching artefact (synthetic chain gives lag 0): it is in the BL-View stamps (CET stamps or a one-hour delay,
   open); `09_fp_topaffected` network bar = pooled shares. Cluster: 14c only (`docs/local/cluster_commands_spec14.md`,
   section "fix-14c"). Merged (`main`, 4338a3a).
16. fix-14a (2026-11-12, `fix-14a-count-mismatch`, to review; spec 14 amendment 1 item 10): a SkyVUE day with
   different CSV / TXT record counts is skipped (`skip_reason` "record count mismatch: CSV n, TXT m", `n_csv` /
   `n_txt` in `14_skyvue_days_<site>.csv`, mismatch days per month in the log), not a failed check; no alignment.
   Item 11 (G 2026-11-12): BL-View time correction in the CZ-Lnz standard file by default (14b,
   `blview_time_correction()`): 3.2.2 Europe/Prague local time → UTC, then all versions −1 h (empirical, mechanism
   unknown); `x_time_shift_h`, lag evidence and DST hours in the sidecar; 14c refuses `--blview-322-tz` on a corrected
   file. Cluster: 14a full record, then the DE-Tha 14c runs; 14b + 14c for CZ-Lnz, expected lag 0 h
   (`docs/local/cluster_commands_spec14.md`, sections "fix-14c" and "fix-14a"). Merged (`main`, 7998608).
17. fix-15 (2026-11-13, `fix-15-figures`, to review): spec 15 figures (subtitle short and wrapped, details in a
   caption explaining the " *"; heights and site-label size for 39 sites + network; ICOS criteria with shapes from
   the style table and a vertical dodge); **CZ-Lnz BL-View main set = quality ≥ 2 (G 2026-11-12)**, `all` and `q3`
   as sensitivity rows (column `quality_min` of `config/measured_sites.csv`; standard file unchanged; spec 14
   amendment 1 item 12); 09_collect.R tells expected missing months (not in `chunk_list.csv`, e.g. SE-Myc 2025-07)
   from real gaps. Cluster: `09_collect.R` and 14c CZ-Lnz (`docs/local/cluster_commands_spec15.md`, section
   "fix-15"). Merged (`main`, c672f8a).
18. Spec 16 (2026-11-17, `spec-16-footprint-uncertainty`, to review): footprint uncertainty from the FFP inputs
   (phase 5b). Step 0 on the PC: 100 rows × 17 runs of design A at BE-Bra (test month) in 2,830 s with 2 workers =
   3.33 s per run per core (1.7–2.1 s single-worker in spec 08); full design ≤ 104,040 A + ≈ 135,000 B runs (the
   h_empirical night runs equal h_pct30 and are computed once) ≈ 239,000 runs = ≈ 133 CPU-hours at 2 s, ≈ 221 at
   3.3 s (below the 300 CPU-hour limit; invalid runs are not computed). Chunks sized at 3.3 s, 2 CPUs, ≤ 2 h. PC smoke
   test on the test month (BE-Bra, FI-Hyy 2022-09, without PBLHmod-clim; 103 / 101 A rows, 6 B rows × 10 draws):
   16b 3.2–3.3 s per FFP call per core (A); base run = spec-15 value at every row (max |diff| 0 pp, x80 0 m). Test,
   one month: TA 1σ effects u* > hc ≈ umean > σv ≈ WD > h > L ≈ 0 (z0 validity only); joint width ≈ 3.3 pp vs
   |ERA5 − PBLHmod| ≈ 0.35 pp. Validity flips: u* 2–3 % (e_ustar); FI-Hyy every +1σ hc and z0 run (e_rsl: zm/z0
   only 2.5 % above the 08a multiplier 3.7). Amendment 1 (G 2026-11-18): hc and z0 also run with the e_rsl test off
   (`rsl_off`, 21 runs per row in A); smoke rerun: base still equal (0 pp), FI-Hyy hc 2.4 pp TA / 0.12 ln x80 with
   the test off (`NA` under `etc`), rsl flips 50 % (BE-Bra 0 %, effects unchanged); 101 extra FFP calls at FI-Hyy, 0
   at BE-Bra; `16_input_ranking` shows hc and z0 from `rsl_off` with "rsl flips n %". Cluster: 16a, then (after G)
   the 16b array and 16c (`docs/local/cluster_commands_spec16.md`). Merged (`main`, 573e9e2); cluster run done.
19. Close-out (G 2026-11-19, before the stand-by): spec 16 figures (one legend below the panels, "reference year per
   site", validity flips per site in `16s_validity_flips`: FI-Hyy hc and z0 and CH-Dav z0 lose every +1σ run to the
   rsl test, 50 % of their runs; no network median in the main figure; `w_row` in `16_oat`; `16c --from-tables`);
   spec 15 exceptions with a minimum of 50 rows per site × stratum (FR-Fon n 3, DE-HoH n 28 and an IT-Ren cell n 27
   now `below_min_n`; exceptions CH-Dav, DE-Tha, IT-Ren; `scripts/09r_fp_figures.R`, `09_network_fp_site_order.csv`),
   re-rendered from the cluster tables; merged (`main`, 57a9460). SkyVUE reader (spec 14 amendment 2): message 004,
   truncated last message, CSV without a header; count-mismatch rule unchanged (G); merged (`main`, e5e52b7). Cluster:
   14a + 14c for DE-Tha.
20. fix-14a June 2026 (G 2026-11-20, spec 14 amendment 3; after the DE-Tha rerun, 4 days still failing): TXT
   fragments anywhere in the file dropped and logged, the CSV row of each found by matching (or paired by order when
   the counts are equal); up to 5 records and 0.1 % of a day failing a pairing check dropped and logged; the run
   summary counts fragments and header-less CSVs on all days. Merged (`main`). Cluster: 14a + 14c for DE-Tha when the
   analysis reopens.


## Data

ERA5 from C. Trotta (Oct 2026), 0.25°, tower cell + 8 neighbours, all stations (details in `CLAUDE.md`):
hourly BLH, cloud base height, sensible heat flux (2020-01-01 → 2026-10-07); 6-hourly temperature and geopotential
on 13 pressure levels (enough for a γ climatology, not for time-varying γ); static orography, land–sea mask, sdor.
SE-Myc uses the corrected coordinates; all stations structured like BE-Bra. BE-Bra tower cell: 6.55 km from the
tower, orography 8.9 m, lsm 0.96, sdor 7.6 m.

Literature: `docs/references/literature_review_pblh.md`, `pblh_literature_2026-10.bib` (~50 verified papers),
`pblh_collection.bib`. Key papers and ECMWF IFS Cy41r2 Part IV downloaded; second-wave papers checked against the PDFs (2026-10-22; Durre2008 from Crossref only, no PDF).
Novelty check: no published network-scale evaluation of EC-derived PBLH against ERA5 and measurements.

FI-Hyy (spec 07): ICOS L2 (FLUXES, FLUXNET 2018–2026, METEO) and ERA5 (three families) under `data/raw/`; Cloudnet
Hyytiälä `epsilon-lidar` and `classification` under `data/raw/cloudnet/FI-Hyy/` (with provenance tables).

L2 data issues found so far (BE-Bra: implausible USTAR, MO_LENGTH kept where H is QC-rejected, positive night-time
H_F_MDS in winter 2022–23, sonic temperature offset, latent pipeline bug; FI-Hyy: USTAR up to 189 m s-1 and WS
median 122 m s-1 on cap rows, T_SONIC 6–220 K in winter, L2 PBLH up to 309 km and 28 % run-mean fills) are
collected for the ETC QC report in `docs/local/etc_qc_report.md` (local working file, not tracked); all open.

## Waiting on others

- Verification of forest canopy-height jumps (FR-Fon, FR-Hes, SE-Svb).

## Open scientific points (details in `docs/project_plan.md`)

- B1 on windy stable nights tends to the B2 limit (p95 ≈ 1.5 km at BE-Bra); candidate: bound by B3. Decide after
  comparison with ERA5/measurements.
- Polar-day convective rule (GL-ZaF, SE-Sto, FI-Lom).
- γ (spec 11, G 2026-10-30): PBLHmod (EC-only, constant γ) and PBLHmod-clim (ERA5 γ climatology). Open: residual
  −0.20 of PBLHmod-clim at FR-Gri (crop site; candidates: mechanical term, entrainment constants of B5); FI-Hyy
  night-time discrepancy = stable formula, not γ (deep end: below the lidar lower bound 43 % vs ≤ 30 %, every γ
  case; the shallow end is too deep in the convective set, 13 % vs ≥ 16 % with PBLHmod-clim); EC-only constant:
  decided (3b, G 2026-11-11: stays 0.01; recommended operational value `gamma_med` = 0.004572 K m⁻¹, side result);
  γ tables of all sites as a deliverable (ERA5 attribution). New: does the residual site spread (crop / grass too
  shallow, forest fine or too deep) follow the energy-balance ratio?

## Local setup notes

- PC rebuild procedure: `docs/local/setup/pc_rebuild_checklist.md` (+ `backup_before_format.ps1`; local only).
- Gitignored but needed: `data/raw/` (ICOS L2 zips for BE-Bra, ERA5 CSVs in `era5/`, `era5_pressure_levels/`,
  `era5_single_levels/`, `target_areas/*.kml`), the local secrets folder, `docs/local/` (local working notes).
- Reference PDFs live in a local bibliography folder outside the repo (not in git; the `file` fields were removed
  from the tracked `.bib` files on 2026-10-26).
- BE-Bra: the replaced FLUXNET 2020–2025 file, the cleanset and the EddyPro outputs are in `data/archive/` (not read
  by the pipeline).
