# Spec 06: consistency check of L2 Obukhov length against the PBLHmod forcing (option C)

Phase 1 (amends spec 02). Branch: `spec-06-obukhov-consistency` (from `main`).
Read first: `CLAUDE.md`, `docs/specs/02_pblh_model.md` (section 3 and amendment 12), `docs/specs/05_ustar_filter_local_config.md`,
`docs/project_plan.md` (open point (b) on the values > 5 km). Agreed with G on 2026-10-21 (option C); implement as written.

## Why

After spec 05, 11 PBLHmod values > 5 km remain at BE-Bra. Ten come from B5 on 2021-08-18 and 2022-03-11: at those
half-hours the L2 QC rejected `H` (NA), but `MO_LENGTH` was kept. EddyPro computed it from `H_UNCLEANED` ≈ 0 W m-2,
so |L| ≈ 3 × 10⁶ m, while the forcing uses gap-filled `H_F_MDS` ≈ 60 W m-2. That pairing makes the B&G growth rate
explode in one step, and the height is carried for the rest of the day.

More generally, the model takes u*, H and T from the forcing but L from L2, and the two can disagree. Over the
23,621 convective half-hours at BE-Bra, L2 L vs the L implied by the forcing has a median ratio of 0.97, but 1,237
rows differ by more than 3× and 152 by more than 10×. L2 `H` is QC-rejected on 47 % of the rows that have
`MO_LENGTH` (35–56 % per year), so rejecting L wherever L2 H is missing (option A) would stop B5 early on 959 of
1,770 convective days. Rejecting inconsistent L without replacing it (option B) would stop B5 early on 268 days.

Option C (G): use the published L2 `MO_LENGTH` wherever it is consistent with the forcing; where it is not, replace
it with the L computed from the forcing, and flag those rows. All inputs remain published (DOI) data.

## Rule (new step in `prepare_pblh_forcing()`)

Order of the forcing steps after this spec:
1. excluded rows → NA (as now);
2. u* plausibility filter (spec 05; rejected rows get `ustar` and `L` NA);
3. **L consistency check (this spec)**, on the values *before* short-gap interpolation;
4. short-gap interpolation (as now; L via 1/L);
5. `rho`, `wtheta` (as now).

For step 3, compute a provisional kinematic heat flux from the not-yet-interpolated H, T and p of the same row:
`wtheta0 = H / (rho0 * cp)`, `rho0 = p_Pa / (r_d * T_K)`, and the Obukhov length implied by the forcing:

`L_forcing = - ustar^3 * T_K / (k * g * wtheta0)`  [m]

(Obukhov 1946; e.g. Stull 1988 — cite the section/equation after checking the PDF; T in place of θv: the moisture contribution to buoyancy is
neglected. At BE-Bra this gives a median L2/L_forcing of 0.97, i.e. about 3 %, well inside the tolerance below.
Constants `k`, `g`, `r_d`, `cp` from `pblh_constants()`.)

A row is **checked** if L2 `L`, `ustar`, `wtheta0` and `T_K` are all available and `|wtheta0| > opts$l_check_wtheta_min`.
Near-neutral rows (`|wtheta0|` ≤ the minimum) are not checked: there |L| is large and ill-determined in both
estimates, and L2 L is kept.

A checked row is **inconsistent** if either holds:
- sign: `sign(L) != sign(L_forcing)`;
- magnitude: `max(L / L_forcing, L_forcing / L) > opts$l_consistency_factor` (same sign).

With `opts$l_source = "L2_checked"` (default), inconsistent rows get `L = L_forcing`, `flag_L_replaced = TRUE`,
`L_source = "forcing"`. All other rows keep L2 L (`L_source = "L2"`). Rows where L2 L is missing stay missing and go
through the usual short-gap interpolation (`L_source = "interpolated"` where filled). L2 L is **not** filled from
the forcing where it is missing (not part of option C).

Other values of `opts$l_source`:
- `"L2"`: no check (behaviour of spec 05); used for the identity check;
- `"forcing"`: L computed from the forcing on every row where `ustar`, `wtheta0`, `T_K` are available (option D), as
  a sensitivity case; rows with `wtheta0 = 0` get `NA`.

## Options (`pblh_options()`)

| option | default | unit | meaning |
|---|---|---|---|
| `l_source` | `"L2_checked"` | – | `"L2_checked"`, `"L2"`, `"forcing"` |
| `l_consistency_factor` | 3 | – | max ratio between L2 L and L_forcing (> 1) |
| `l_check_wtheta_min` | 0.005 | K m s-1 | rows with \|wtheta0\| at or below this are not checked |

Validate (`l_consistency_factor > 1`, `l_check_wtheta_min >= 0`, `l_source` one of the three). Record as
**amendment 13** in `docs/specs/02_pblh_model.md`.

## Output

- Forcing and `run_pblh_model()` output gain: `L_l2` (L2 `MO_LENGTH` as read, after the u* filter) [m],
  `L_forcing` [m], `L_source` (`"L2"`, `"forcing"`, `"interpolated"`, `NA`), `flag_L_replaced` (logical).
  `L` remains the value used by the model.
- `standardize_pblhmod()` (spec 04): `flag_L_replaced` joins the qc 1 condition; `qc_flag_meaning` for 1 becomes
  "forcing interpolated or replaced (…)".

## Script `scripts/02_run_pblh_model.R`

- Backup first: copy the current derived PBLHmod files and `output/tables/04_metrics_BE-Bra.csv` to
  `data/derived/pblhmod/pre_spec06/` (gitignored).
- Cases: all existing cases now use `l_source = "L2_checked"`. Add `l_L2` (`l_source = "L2"`), `l_forcing`
  (`l_source = "forcing"`), `lcons2` (`l_consistency_factor = 2`) and `lcons5` (`l_consistency_factor = 5`).
- Identity check: `l_L2` equals the spec-05 default from the backup (max difference < 1e-9 m, same NA, regime and
  method).
- Report per case: rows checked; rows replaced (by sign, by magnitude), overall and by regime; number of days with at
  least one replacement; share of rows with `L_source = "L2"`; `flag_sign_inconsistent` count; count of
  `pblh_m > 5 km` and `> 10 km`; quantiles per regime (as now).
- Then rerun `scripts/04_era5_compare.R` and report pair 1 (default vs ERA5 tower cell) overall and by regime,
  spec 05 vs spec 06.

## Tests

- `L_forcing` hand case (show the calculation): e.g. u* = 0.5 m s-1, H = 200 W m-2, T = 293.15 K, p = 100000 Pa
  → rho0, wtheta0, L_forcing.
- Sign mismatch → replaced, flagged, `L_source = "forcing"`.
- Magnitude: ratio 2.9 kept, 3.1 replaced (factor 3), in both directions (L2 larger and smaller).
- `|wtheta0|` ≤ minimum → not checked, L2 kept.
- L2 L missing → stays missing, interpolated if the gap is short (`L_source = "interpolated"`), never filled from
  the forcing.
- Rows rejected by the u* filter → L NA, no check.
- `l_source = "L2"` gives output identical to spec 05; `"forcing"` uses `L_forcing` everywhere it exists.
- The BE-Bra jump case as a synthetic fixture: H_F_MDS ≈ 60 W m-2 with L2 L = −3 × 10⁶ m inside a convective
  window → replaced; B5 growth over the step stays below a few hundred metres.
- Option validation.

## Acceptance criteria

- All tests pass; scripts 02 and 04 run on BE-Bra; `l_L2` reproduces spec 05 exactly.
- No PBLHmod value from the two jump days above 5 km in the default case; report any value > 5 km that remains,
  with its rows (as in the spec 05 review).
- Amendment 13 in spec 02; `docs/STATUS.md` (spec 06 row, defaults line, before/after numbers); open point (b) in
  `docs/project_plan.md` marked resolved with the counts.

## Out of scope

Filling missing L2 L from the forcing; moisture correction of L_forcing (θv, LE); changes to the u* filter; spec 03.

## Summary to report back

As in `CLAUDE.md`, plus: rows checked/replaced per case and by regime; share of L2 L used; remaining values > 5 km;
`flag_sign_inconsistent` before/after; pair-1 metrics spec 05 vs spec 06; anything unexpected (e.g. many
replacements at night, or clusters of replacements on specific days or periods).
