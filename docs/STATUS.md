# Project status (handover note)

Last update: 2026-10-21. Read this first in a new session, then `docs/project_plan.md` and `CLAUDE.md`.

## PBLH sources (names used in the code and in the standard files)

| Source | `source_id` | What it is | Class |
|---|---|---|---|
| PBLHmod v2 | `pblhmod_default`, `pblhmod_gamma0025` | Our EC-based model (Kljun et al. 2015 App. B): stable B1 Nieuwstadt, neutral B2/B3 auto, convective B5 Batchvarova–Gryning slab daily cycle. `run_pblh_model()`, `R/pblh_*.R`; label = options case (γ = 0.01 default, 0.0025 sensitivity) | parameterisation |
| ERA5 BLH | `era5_blh_center`, `era5_blh_bilinear` | ECMWF reanalysis diagnostic (bulk Ri = 0.25, IFS Cy41r2 §3.10.1); tower cell, or bilinear to the tower from 4 cells | model |
| Operational ICOS PBLH | `icos_oper` | `PBLH` column of L2 FLUXES, produced by the ETC pipeline with `legacy/calc_PBL_height.R` + `legacy/proc_PBL_height.R`; `qc_flag` 1 = 50 m floor, 2 = run-mean fill | parameterisation |
| Legacy (spec 03, to come) | `legacy_<variant>` | Port of the legacy code with switches, for the regression | parameterisation |
| Measured PBLH (later) | – | Ceilometer / lidar / radiosonde (ACTRIS/Cloudnet, E-PROFILE, IGRA) | turbulence / aerosol / thermodynamic |

## Specs

| Spec | Content | Phase | State |
|---|---|---|---|
| 01 | Site metadata, time-varying heights, exclusions, ICOS L2 readers, unit checks | 0 | merged (`main`) |
| 02 | PBLHmod v2 | 1 | merged (`main`, 58d45ad) |
| 04 | ERA5 readers (3×3 cells, pressure levels, static), standard PBLH format, first PBLHmod vs ERA5 comparison at BE-Bra (amendments 1–8) | 0, 2, 3 (first step) | merged (`main`, aba3394) |
| 05 | u* plausibility filter in the PBLHmod forcing (spec 02 amendment 12); local config override (`config.local.yml`, no cluster paths in tracked files) | 1 | merged (`main`, 25ae55c) |
| 06 | L (Obukhov length) computed from the PBLHmod forcing instead of L2 `MO_LENGTH` (open point in `project_plan.md`) | 1 | to be drafted |
| 03 | Legacy regression: reproduce operational PBLH from L2, port legacy code, attribute differences to PBLHmod v2 | 1 | drafted, not started |

Phases (see `project_plan.md`): 0 done; 1 model done, legacy regression pending (spec 03); 2 standard format and
three sources done (spec 04), measured PBLH pending; 3 started at BE-Bra only; 4–6 not started.

- Tests: 645 passing (`main`, after spec 05). Environment: R 4.4.1 + renv (Posit Package
  Manager), same as the cluster. Machine-specific paths in the untracked `config.local.yml` (template
  `config.local.yml.example`).
- PBLHmod defaults: `stable_method = "nieuwstadt"` (B1), `neutral_method = "auto"` (B2 if 10 < N/|f| < 70, else
  B3), `combine = "bg_only"`, `gamma = 0.01` K/m, `h_source = "H_F_MDS"` with `h_qc_max = 1`; u* input filter
  `ustar_ws_max = 0.5` (where WS ≥ `ws_min_ratio_test = 1` m s-1), `ustar_max = 3` m s-1 (spec 02 amendment 12;
  u* and L set to NA, short gaps interpolated); raw output (no smoothing, clipping or filling). Forcing T and p from
  measured meteo (FLUXNET `TA_F`/`PA_F` 2020–2025, METEO `TA`/`PA` 2026).
- BE-Bra PBLHmod (default, with the u* filter): convective median 698 m (p5–p95 72–1586), stable 261 m, all rows
  413 m; γ = 0.0025 raises the convective median to 1104 m. Before the filter (spec 02): 704 m (72–1609), 271 m,
  423 m; 1109 m.
- Spec 05 at BE-Bra (2026-10-21): the filter rejects 2,149 of 100,186 u* values (2.1 %: ratio 1,328, cap 46,
  both 775); 618 then interpolated, 1,531 left NA. PBLHmod > 5 km: 179 → 11; > 10 km: 70 → 0. `nofilter`
  reproduces the spec-02 default (max |diff| 3.6e-12 m). Pair 1 vs ERA5: RMSE 623 → 349 m, Pearson 0.49 → 0.73,
  bias −9 → −46 m, median log ratio −0.136 → −0.151, within a factor 2 71.8 → 72.4 %.
- Spec 04 at BE-Bra (2026-10-17): overlap 2020-10-25 → 2026-08-25; PBLHmod default vs ERA5 tower cell 50,269 matched
  hours, bias −9 m, median log ratio −0.14, 72 % within a factor 2, Spearman 0.78; stable −0.23, convective −0.07;
  JJA evenings (19–22 UTC) PBLHmod ≈ half of ERA5. Operational PBLH vs ERA5: median log ratio +0.57 (qc 0).
  Open points from spec 04 in `project_plan.md` ("Open points").

## Next steps

1. ~~Review and merge spec 04 (`spec-04-era5-compare`).~~ Done (2026-10-21, merged at aba3394; `mad_log_ratio`
   switched to the scaled MAD at review, amendment 8).
2. ~~Review and merge spec 05 (`spec-05-ustar-filter`).~~ Done (2026-10-21, merged at 25ae55c).
3. Spec 06 — L computed from the PBLHmod forcing (H_F_MDS, u*, T, p) instead of L2 `MO_LENGTH`: at BE-Bra, 10 of
   the 11 values > 5 km left after the u* filter come from B5 on days where L2 H is QC-rejected but `MO_LENGTH`
   (from H_UNCLEANED ≈ 0, |L| ≈ 3e6 m) is kept; 1,237 convective rows have L2 L off by > 3× from the forcing L.
   Before spec 03, which compares the legacy forcing with PBLHmod's.
4. Spec 03 — `docs/specs/03_legacy_regression.md`, branch `spec-03-legacy-regression`. Inputs of the ETC call
   confirmed by G (2026-10-17); only published L2 files are used: AT, AP, rho, cp (EddyPro ambient values, not in
   L2) are rebuilt from L2 with the EddyPro formulas; MOL gap-filled by LOCF before the call.
   Already established: `proc_PBL_height.R` floors at 50 m (4,055 rows at BE-Bra) and fills NA with the day/night
   mean of each processing run (1654.58 / 843.07 m in the 2020–2025 run, 1591.68 / 682.86 m in the 2026 run; ≈ 18 %
   of operational rows). PBLH is NA Jan–Oct 2020.

## Data

ERA5 from C. Trotta (Oct 2026), 0.25°, tower cell + 8 neighbours, all stations (details in `CLAUDE.md`):
hourly BLH, cloud base height, sensible heat flux (2020-01-01 → 2026-10-07); 6-hourly temperature and geopotential
on 13 pressure levels (enough for a γ climatology, not for time-varying γ); static orography, land–sea mask, sdor.
SE-Myc uses the corrected coordinates; all stations structured like BE-Bra. BE-Bra tower cell: 6.55 km from the
tower, orography 8.9 m, lsm 0.96, sdor 7.6 m.

Literature: `docs/references/literature_review_pblh.md`, `pblh_literature_2026-10.bib` (~50 verified papers),
`pblh_collection.bib`. Key papers and ECMWF IFS Cy41r2 Part IV downloaded; second-wave papers: G will notify.
Novelty check: no published network-scale evaluation of EC-derived PBLH against ERA5 and measurements.

## Waiting on others

- Second-wave papers (G).
- Relaxed FFP validity thresholds used at ETC (for the phase-5 spec).
- Verification of forest canopy-height jumps (FR-Fon, FR-Hes, SE-Svb).

## Open scientific points (details in `docs/project_plan.md`)

- B1 on windy stable nights tends to the B2 limit (p95 ≈ 1.5 km at BE-Bra); candidate: bound by B3. Decide after
  comparison with ERA5/measurements.
- Polar-day convective rule (GL-ZaF, SE-Sto, FI-Lom).
- γ options 2–3 (ERA5 climatology / time-varying) after phase 1.

## Local setup notes

- PC rebuild procedure: `docs/setup/pc_rebuild_checklist.md` (+ `backup_before_format.ps1`).
- Gitignored but needed: `data/raw/` (ICOS L2 zips for BE-Bra, ERA5 CSVs in `era5/`, `era5_pressure_levels/`,
  `era5_single_levels/`, `target_areas/*.kml`), `C:\WORK\secrets`, `docs/local/` (local working notes).
- Reference PDFs live outside the repo in `C:\WORK\bibliography_new_win11\` (not in git; the `.bib` `file` fields
  point there).
