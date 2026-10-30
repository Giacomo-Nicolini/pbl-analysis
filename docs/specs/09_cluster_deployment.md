# Spec 09: cluster deployment — all stations, reference-year footprints, batch infrastructure

Phases 1–5 (network scale). Branch: `spec-09-cluster` (from `main`).
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/project_plan.md` ("Scope"), `docs/specs/08_footprint_sensitivity.md`
(sections 3, 4b, 4c and the storage analysis), `docs/specs/07_measured_pblh_cloudnet.md`,
`docs/specs/11_gamma_era5_davos.md`. Agreed with G on 2026-10-28; γ strategy added 2026-10-30.

## Goal

Run the pipeline at the 40 class 1–2 ICOS ecosystem stations (`config/sites.csv`) on the Unitus cluster (SLURM):
inputs check, PBLHmod, ERA5 comparison, footprint validity and coverage (full period), and the footprint sensitivity
on one reference year per station. Claude Code implements and tests on G's PC (light work only, no heavy runs:
RAM issue on the PC, see STATUS). **G runs everything on the cluster** with Cowork's runbook
(`docs/local/cluster_runbook.md`, untracked). The repo is public: no cluster paths, user names, host names or e-mail
addresses in tracked files (config.local.yml and an untracked `slurm/local.env` hold them).

## Cluster facts (G, 2026-10-28)

SLURM, partition set by G (local setting); R 4.4.1, renv works, packages can be installed; storage sufficient;
code in a home `dev/pbl-analysis` clone, data under a separate big-data project folder.
Partition facts (`sinfo`/`scontrol`, 2026-10-29): 5 nodes × 48 CPUs, ≈ 240 GB memory per node, **no wall-time
limit**, shared with other users. G usually requests 6 CPUs per job; `local.env` sets CPUs and memory per job
(default 6 CPUs, 4 GB per CPU). Chunks are kept short anyway (target ≤ 2 h, option `chunk_target_min`) so that a
failure costs little and the array spreads well across nodes; restartable as in spec 08.

## Part A — infrastructure (implemented on the PC, tested without a cluster)

### A1. Configuration

- New keys in `config.yml`:
  - `oneflux_bif_dir`: default `"{icos_l2_dir}/{site}_output_L2I_2026"` (BIF next to the L2 zips, as on the PC);
    cluster `"LOCAL"`. On the cluster the BIF files live in one ONEFlux run folder per site and processing run, so the
    value may contain a glob (`*`). Resolution: list the matching `ICOSETC_{site}_FLUXNET_BIF_*.csv`; keep those whose
    `PRODUCT_NAME` years match the HH FLUXNET zip in use (`read_oneflux_bif()` check); if several remain, take the
    newest processing timestamp (from the folder name or `PRODUCT_PROCESSING_DATE`) and report the choice. Abort if
    none matches.
  - `measured_dir`: default `"data/raw/measured"` (spec 10 inputs); cluster `"LOCAL"`.
  - cluster profile: `target_areas_dir`, `standard_dir`, `derived_dir`, `output_dir` from `"TODO"` to `"LOCAL"`.
- Update `config.local.yml.example` (placeholders only) and the reader of the BIF to use `oneflux_bif_dir`.
- Readers must find FLUXNET HH zips with any year span (`..._HH_<y1>-<y2>_v*.zip`); several matches → abort with a
  clear message.

### A2. Per-site driver and run status

- `scripts/run_site.R SITE --steps=01,11a,02,04,08a [--force] [--with-oper]`: runs the steps in order for one site, each
  in a fresh R process (`callr` or `Rscript`), logging to `{output_dir}/logs/{site}/{step}_{timestamp}.log`.
- Status file per site × step in `{derived_dir}/run_status/{site}__{step}.json`: start, end, exit code, git commit,
  R version, host name (written only to the untracked derived folder), inputs hash. A step whose status is "ok" with
  the same inputs hash is skipped unless `--force`.
- A failing step stops the site and is reported; other sites continue (one SLURM array task per site).

- **γ strategy (G, 2026-10-30).** Two PBLHmod variants at every site: **PBLHmod** (`pblhmod_default`, EC data
  only, constant γ: the code default) and **PBLHmod-clim** (`pblhmod_clim`, site × month γ climatology from ERA5,
  spec 11). Step `11a` computes the site's γ climatology before step 02; step 02 writes both variants (plus the
  existing sensitivity cases). The EC-only constant is **revised after the network γ climatology** (decision point
  in Part B), so the constant is an option read from `config.yml` (`pblh_gamma_default`, now 0.01), not hard-coded.
- **Reference cases for the supplement (G, 2026-10-30):** the two constants found in the literature and in practice
  are kept as named standard sources at every site, written by step 02/04 and compared like the main variants:
  `pblhmod_gamma001` (γ = 0.01 K m⁻¹, Batchvarova & Gryning 1991 as cited by Kljun et al. 2015) and
  `pblhmod_gamma0025` (γ = 0.0025 K m⁻¹, the value of the legacy ETC code; check whether it also has a literature
  source and cite it if so). Tables always; supplementary figures (same layout as the main ones, labelled "reference
  γ values"); not in the main figures.

### A3. SLURM templates (tracked, no machine data)

- `slurm/sites_array.sbatch`: array over the row index of a site list (`config/run_sites.txt`, default: all 40);
  calls `scripts/run_site.R`. `slurm/footprint_chunks.sbatch`: array over the rows of a chunk list (A5).
- Machine data (partition, mail address, R setup line, memory, wall time) from an untracked `slurm/local.env`,
  template `slurm/local.env.example`; `.gitignore` updated. Templates use `--cpus-per-task` from it (default 6).
- `scripts/slurm_submit.R --what=sites|chunks [--dry-run]`: writes the list, checks `local.env`, prints or runs the
  `sbatch` command. `--dry-run` renders the final script for review.

### A4. Network collection

`scripts/09_collect.R`: reads the per-site outputs and writes network tables `{output_dir}/tables/09_network_*.csv`:
input availability and reference period (01), γ climatology per site × month (11a) and the network summary
(median, p25, p75 per month and overall, sites with only 2 fit levels), PBLHmod and PBLHmod-clim summary per regime (02), u* filter / L check counts, PBLHmod vs
ERA5 metrics (04, tower cell and bilinear), ERA5 cell distance and orography, footprint coverage per case × set ×
subset (08a), reference-year table (spec 08 section 4c), run status. Flags sites with problems (missing inputs,
unit-check failures, > 5 km values, low coverage). Figures: network maps are out of scope; a dot plot of median ln
ratio PBLHmod/ERA5 by site and regime.

### A5. Footprint run on the cluster

Settings as spec 08 (ETC grid ±1000 m at 1 m, umean only, smoothing EBImage-compatible), with:
- **Storage `c5`** (G, 2026-10-28): summed grids stored as 5 m × 5 m block sums (mass-conserving, tower cell at a
  block centre), float32 + gzip, for all validity sets; per-half-hour metrics computed at 1 m and stored as tables.
  Monthly sums kept, annual climatologies built from them (sum, upsample, smooth, normalise).
- **Cases:** `pblhmod` (with the revised EC-only constant), `pblhmod_clim` and `era5` at every site; the reference
  cases `pblhmod_gamma001` and `pblhmod_gamma0025` at the extra sites below (supplement); `measured` where a measured source exists; `era5_bilinear` and the
  ±20 % anchors (`pblhmod_x0.8`, `pblhmod_x1.2`) at the sites listed in `config/footprint_extra_sites.txt` (default:
  BE-Bra, FI-Hyy, FR-Gri, CH-Dav, SE-Svb, IT-Cp2: one per climate/canopy type; G may edit); `oper` only with
  `--with-oper` (internal outputs).
- **Sets/subsets:** `etc`, `kljun`, `kljun_norsl`, `none`; `valid`, `valid_qc`, `valid_qc_used`, `valid_qc_h`.
- **Chunks:** `scripts/08c_make_chunks.R` writes the list site × case × month for each site's reference year
  from the approved reference years only: `09_collect.R` writes the proposal (`output/tables/08_reference_years.csv`);
  after review, the approved list is committed as `config/reference_years.csv` (site, year, reason); 08c aborts if a
  site has no approved year. Each chunk carries an estimated duration from the spec-08 timings;
  months are split into half-months when the estimate exceeds `chunk_target_min`. Each array task runs one chunk with
  the job's cores as workers; chunk outputs keyed by inputs hash (restartable, as spec 08).
- After all chunks of a site: `08_footprint_compare.R SITE` (reference year), then `09_collect.R`.
- Memory per worker from the spec-08 measurements (≈ 2.5 GB with c5 sums); `local.env` sets `--mem`.

### A6. Tests (no cluster needed)

BIF glob resolution (several candidates, product mismatch, newest processing run); FLUXNET HH year-span patterns;
run-status skip logic (same / changed hash, `--force`, failure stops the site); chunk list generation (reference
year approved / not approved, splitting); sbatch rendering in `--dry-run` with a fake `local.env` (no machine data in
the rendered output beyond what `local.env` provides); c5 block sums (mass conservation, tower-centred blocks,
annual = sum of months) — reuse spec-08 tests where they exist.

## Part B — runs on the cluster (G, with the runbook)

1. Setup (runbook): pull, `renv::restore()`, `config.local.yml` and `slurm/local.env` on the cluster, Kljun code and
   target-area KML copied to untracked folders, data folders created, test suite as a small job.
2. **Regression:** BE-Bra and FI-Hyy steps 01/02/04 on the cluster; compare with the PC outputs (copied over):
   PBLHmod max |diff| < 1e-9 m, identical NA/regime/method, identical pair-1 metrics. Same inputs required (check
   the md5 of the L2 zips and BIF on both machines; if they differ, report and stop).
3. All sites, steps 01, 11a, 02, 04, 08a (array). Review `09_network_*.csv` with Cowork: input problems, per-site
   specifics, the reference-year table → G approves the years.
3b. **Decision point (before any footprint run):** from the network γ table, G chooses the recommended EC-only
   constant (e.g. the network median, rounded) and the γ tables to publish (ERA5-derived: Copernicus attribution).
   Set `pblh_gamma_default`, commit `config/gamma_era5_clim.csv` (site × month), rerun step 02 and 04 at all sites
   (cheap; γ = 0.01 and 0.0025 stay as the reference sources `pblhmod_gamma001` / `pblhmod_gamma0025`), collect
   again.
4. Footprint chunks (array) on the reference years, then per-site compare and collection.
5. Measured (after spec 10 is merged, which generalises the comparison code): Cloudnet for FI-Lom, DE-RuS, FR-Gri
   (07b download on the login node, size first; 07c; comparison), FI-Hyy Cloudnet re-downloaded on the cluster
   (checksums) rather than copied; SIRTA (FR-Gri) and later PI data copied into `measured_dir`.

## Docs

`CLAUDE.md`: cluster workflow (driver, status files, SLURM templates, local.env, chunks), the γ step and the two
PBLHmod variants, the new config keys, the
"no machine data in tracked files" rule extended to `slurm/`. `docs/STATUS.md`: spec 09 row. `docs/project_plan.md`:
execution on the cluster, cases and extra-site list. The runbook stays in `docs/local/` (Cowork).

## Acceptance criteria

- Part A: tests pass on the PC; `--dry-run` renders valid sbatch scripts; `git grep` for `/bigdata`, `/home/`,
  `@unitus`, `g.nicolini`, `c.trotta` returns nothing outside `docs/local/` (untracked) and placeholders.
- Part B (reported by G): BE-Bra/FI-Hyy regression on the cluster passes; network tables for all sites; reference
  years approved; footprint outputs for every site's reference year; collection tables.

## Out of scope

Network maps; Monte Carlo (phase 5b); own MLH retrievals; Python port; changes to PBLHmod other than the γ constant.

## Summary to report back

Part A as in `CLAUDE.md`, plus: the dry-run sbatch scripts (rendered with the example env), estimated chunk count and
CPU-hours for all sites (from spec-08 timings), and anything that needs a cluster fact to be confirmed.
