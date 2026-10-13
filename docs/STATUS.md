# Project status (handover note)

Last update: 2026-10-13. Read this first in a new session (Cowork or Claude Code), then `docs/project_plan.md`
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

1. **Spec 03 — legacy regression.** Blocked: needs the ETC operational code that calls `calc_PBL_height()` and
   post-processes its output, to be copied into `legacy/` (e.g. `legacy/etc_pipeline_pblh.R`).
   Findings to explain: operational `PBLH` at BE-Bra has a 50 m floor (4,055 rows) and repeated constants —
   1654.58 m (15,649 rows) and 843.07 m (2,739) in 2020–2025, 1591.68 m (2,649) and 682.86 m (427) in 2026 —
   probably a percentile cap and a fill value per processing run (hypothesis, unconfirmed).
2. **Spec 04 (proposed, not started)** — standard PBLH format (phase 2) and a first PBLHmod vs ERA5 comparison at
   BE-Bra (2025–2026). Independent of spec 03.

## Waiting on others

- ERA5 BLH 2020–2024 and 3×3 neighbours + invariant fields (geopotential, land–sea mask, sdor) at 0.25°
  (C. Trotta). SE-Myc must use the corrected coordinates (58.36503 N, 12.1694 E). ERA5-Land has no BLH.
- Relaxed FFP validity thresholds used at ETC (for the phase-5 spec).
- Verification of forest canopy-height jumps (FR-Fon, FR-Hes, SE-Svb).

## Open scientific points (details in `docs/project_plan.md`)

- B1 on windy stable nights tends to the B2 limit (p95 ≈ 1.5 km at BE-Bra); candidate: bound by B3. Decide after
  comparison with ERA5/measurements.
- Polar-day convective rule (GL-ZaF, SE-Sto, FI-Lom).
- γ options 2–3 (ERA5 climatology / time-varying) after phase 1.

## Local setup notes

- PC rebuild procedure: `docs/setup/pc_rebuild_checklist.md` (+ `backup_before_format.ps1`).
- Gitignored but needed: `data/raw/` (ICOS L2 zips for BE-Bra, ERA5 CSVs, `target_areas/*.kml`), `C:\WORK\secrets`.
