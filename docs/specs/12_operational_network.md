# Spec 12: operational ICOS PBLH vs PBLHmod and ERA5 across the network (internal)

Phase 3, internal ETC output. Branch: `spec-12-oper-network` (from `main`). Agreed with G on 2026-11-02.
Read first: `CLAUDE.md`, `docs/STATUS.md`, `docs/project_plan.md` ("Scope", operational PBLH open points),
`docs/specs/04_era5_standard_compare.md` (pair 5, amendment 7: fill detection), `docs/specs/09_cluster_deployment.md`
(A4, collection). Light work only on this PC; the network run happens on the cluster.

## Why

The operational ICOS PBLH (column `PBLH` in L2 FLUXES, standard source `icos_oper`) is used today for the ICOS
footprints. Scope (G, 2026-10-23, confirmed 2026-11-02): it stays **out of the paper's three-way comparison**, but an
internal network table is wanted, for ETC: how far is the operational product from PBLHmod v2 and from ERA5, and how
much of it is floor (50 m) or run-mean fill. A summary figure may follow (supplement at most; decided later).

## What exists (check first, report)

Step 04 already writes `{site}__icos_oper` as a standard file at every site and computes pair 5 (operational vs
PBLHmod, on half-hours). Check what step 04 writes for the operational source today (pairs, regimes, fill flags) and
reuse it. If step 04 must change, keep the change small: re-running 04 at all sites on the cluster is cheap (minutes
per site) but must be stated in the report.

## Part A — per-site comparison (in step 04 or a helper it calls)

- Pairs: `icos_oper` vs `pblhmod_default`, vs `pblhmod_clim` (half-hours); `icos_oper` vs `era5_blh_center` and
  `era5_blh_bilinear` (ERA5 hours, as the other ERA5 pairs). Reference period: the site's FLUXES range, as everywhere.
- Operational row classes (per half-hour): `floor` (exactly 50 m), `fill` (run-mean fill, spec 04 amendment 7),
  `implausible` (> 5 km; also count > 10 km), `regular` (the rest), `missing`.
- Fill detection is per calendar year at present; processing runs may not align with calendar years at other sites.
  Report per site and year the candidate fill values (day/night) and their counts, and flag sites where the rule looks
  unreliable (e.g. no single dominant value, or several per year). Don't change the rule silently.
- Metrics (same set and code as the other pairs): n, median ln ratio, within factor 2, FB, nRMSD, Spearman, scaled MAD;
  by regime (from PBLHmod: stable / neutral / convective, plus all) and by day/night. Computed twice: on `regular`
  rows only, and on all non-missing rows (to show the effect of floor and fills).
- Output: `{output_dir}/tables/internal/12_oper_{site}.csv` (long format) and `12_oper_classes_{site}.csv` (counts and
  shares per class, per year).

## Part B — network collection (`09_collect.R`)

- `{output_dir}/tables/internal/09_network_oper.csv`: one row per site × pair × regime × row set; plus
  `09_network_oper_classes.csv` (share of floor, fill, implausible, missing per site and year) and the fill-rule flags.
- Header comment and file names make clear the tables are **internal** (not for the paper). They live under
  `output/` (never committed) like all outputs.
- No figure in this spec (G will decide on a summary figure after seeing the table).

## Tests

Row classification on a synthetic series (floor, a fill value repeated per year/day-night, > 5 km, missing);
fill-rule flag when a year has two dominant values; metrics identical to the existing pair-5 values at BE-Bra on
`all` rows; collection with one site missing the operational source (row reported, not an error).

## Acceptance criteria

Tests pass; BE-Bra and FI-Hyy per-site tables on the PC; pair-5 values unchanged where they existed; collection runs
on the two sites; report says whether step 04 changed (→ G re-runs 04 at all sites on the cluster) or not.

## Out of scope

Operational PBLH in the paper's figures or main tables; the `oper` footprint case (spec 09, `--with-oper`, step 4);
the measured comparison of the operational source (FI-Hyy spec-07 reference test stays as is).

## Summary to report back

As in `CLAUDE.md`, plus: what step 04 already provided, whether 04 changed, the class shares and the metrics at BE-Bra
and FI-Hyy, and any site/year where the fill rule looks unreliable.
