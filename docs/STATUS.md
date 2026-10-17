# Project status (handover note)

Last update: 2026-10-17. Read this first in a new session (Cowork or Claude Code), then `docs/project_plan.md`
and `CLAUDE.md`.

## Way of working

- **Cowork (Claude Desktop)**: scientific discussion and design; writes `docs/` (plan, specs, notes) and reviews code.
- **Claude Code (VS Code)**: implements one spec at a time on a branch `spec-NN-...`, tests, commits, merges to `main`.
- Each task: discuss → `docs/specs/NN_*.md` → Claude Code in Plan mode → approve → implement → review → merge.

## Done

| Spec | Content | State |
|---|---|---|
| 01 | Site metadata, time-varying heights, exclusions, ICOS L2 and ERA5 readers, unit checks | merged (`main`) |
| 02 | PBLHmod v2: Kljun 2015 App. B (B1 stable, B2/B3 neutral auto, B5 Batchvarova–Gryning daily cycle, RK4) | merged (`main`, 58d45ad) |

- Tests: 375 passing. Environment: R 4.4.1 + renv (Posit Package Manager), same as the cluster.
- Defaults: `stable_method = "nieuwstadt"` (B1), `neutral_method = "auto"` (B2 if 10 < N/|f| < 70, else B3),
  `combine = "bg_only"`, `gamma = 0.01` K/m, `h_source = "H_F_MDS"` with `h_qc_max = 1`, raw output (no smoothing,
  clipping, or filling).
- BE-Bra first results (default): convective median 704 m (p5–p95 72–1609), stable 271 m, all rows 423 m.
  γ = 0.0025 raises the convective median to 1109 m.

## Next steps

Both specs drafted 2026-10-17, ready for Claude Code (Plan mode first). Suggested order: spec 04, then spec 03
(independent; spec 03 reuses `save_figure()` from spec 04 if merged first).

- `docs/specs/04_era5_standard_compare.md` — branch `spec-04-era5-compare`.
- `docs/specs/03_legacy_regression.md` — branch `spec-03-legacy-regression`. ETC call inputs confirmed by G
  (2026-10-17): AT, AP, rho, cp from EddyPro full_output (AP from site altitude, no biomet); MOL gap-filled by
  LOCF. Exact reproduction needs the EddyPro full_output columns (optional, from G); otherwise they are rebuilt from
  L2 with the EddyPro formulas (approximate acceptance thresholds in the spec).
- PBLHmod v2 forcing already uses measured meteo: FLUXNET `TA_F`/`PA_F` 2020–2025 (METEO measurements, gap-filled),
  METEO `TA`/`PA` in 2026. No change needed (G, 2026-10-17: prefer METEO over EddyPro ambient values).

1. **Spec 03 — legacy regression.** Unblocked (2026-10-16): the ETC post-processing is in
   `legacy/proc_PBL_height.R` (negatives → NA; floor at 50 m; remaining NA filled with the mean of all "day" or
   "night" values of the processing run, day = hour 06–18 of the timestamp). The operational artefacts at BE-Bra are
   now **explained and reproduced exactly** (FLUXES L2, day/night from the hour of `TIMESTAMP_END`, local standard
   time; means taken after the 50 m floor):
   - 50 m floor: 4,055 rows;
   - 2020–2025 run: day fill 1654.58 m (15,649 rows) = day mean, night fill 843.07 m (2,739) = night mean;
   - 2026 run: day fill 1591.68 m (2,649), night fill 682.86 m (426; one more row in 2022 is a coincidence).
   - So ≈ 18 % of 2020–2026 operational PBLH rows are run-mean fills and ≈ 4 % are floored. The day fill value is
     inflated by outliers (median of non-filled day values ≈ 1509 m; max 91 km). PBLH is NA Jan–Oct 2020.
   Remaining spec 03 work: run legacy `calc_PBL_height()` + `proc_PBL_height()` on BE-Bra, reproduce the operational
   column, then switch legacy steps/fixes one at a time (γ, smoothing, gap order, units, floor, fill) to quantify each
   effect against PBLHmod v2.
2. **Spec 04 (drafted 2026-10-17)** — new ERA5 readers (format changed, see below), standard PBLH format
   (phase 2) and a first PBLHmod vs ERA5 comparison at BE-Bra (2020–2026). Independent of spec 03.
   Must include: reader for the new ERA5 hourly / pressure-level / static files (cell position from header, not
   file name; sshf accumulated J m-2, sign flip; 9 cells), `config.yml` path update (local + cluster), ERA5 cell
   metadata per site (cell centre, distance, orography height vs site altitude, lsm, sdor), 3×3 spread as
   representativeness metric. Frame the comparison against `literature_review_pblh.md`.

## Data received (2026-10-15)

ERA5 from C. Trotta, 0.25°, tower cell + 8 neighbours, all stations (verified on BE-Bra; details in `CLAUDE.md`):

- Hourly 2020-01-01 → 2026-10-07: `boundary_layer_height`, `cloud_base_height`, `surface_sensible_heat_flux`
  (accumulated J m-2, positive downward). Complete at BE-Bra (59,328 h, no gaps).
- Pressure levels, 6-hourly, 13 levels 50–1000 hPa: temperature, geopotential. Enough for a γ climatology (option 2,
  e.g. θ gradient 850–700 hPa); too coarse for time-varying γ (option 3) or for inversion diagnosis: that would need
  hourly data and the intermediate levels 975–750 hPa. Request only if option 3 is pursued (agreed 2026-10-16).
- Confirmed by G (2026-10-16): SE-Myc extraction uses the corrected coordinates; all stations structured like BE-Bra.
- Static: surface geopotential, land–sea mask, sdor. BE-Bra tower cell: 51.25 N 4.50 E, 6.5 km from the tower,
  orography 8.9 m, lsm 0.96, sdor 7.6 m.
- The ERA5 extraction question (nearest cell? neighbours?) is resolved: nearest 0.25° cell, 3×3 neighbours available.

Literature (2026-10-15): `docs/references/literature_review_pblh.md` + `pblh_literature_2026-10.bib` (~50 verified
papers: ERA5 BLH evaluations, remote-sensing ABLH methods and datasets, stable/convective parameterisations,
PBLH from flux towers, footprint sensitivity to h). Novelty check: no published network-scale evaluation of
EC-derived PBLH against ERA5 and measurements was found. `pblh_collection.bib` fixed (DOIs, typos).

## Waiting on others

- Key papers downloaded (2026-10-16). Still missing: `ECMWF2016` (IFS Cy41r2 Part IV, free from the ECMWF
  e-library). Second-wave papers: G will notify.
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
  `era5_single_levels/`, `target_areas/*.kml`), `C:\WORK\secrets`.
- Reference PDFs live outside the repo in `C:\WORK\bibliography_new_win11\` (not in git; the `.bib` `file` fields
  point there).
