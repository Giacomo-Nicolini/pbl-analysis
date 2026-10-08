# Legacy code review (`legacy/`)

Reviewed 2026-10-08. Reference only: new code goes in `R/` and `scripts/`.

## What's there

| File | Purpose | Reuse |
|---|---|---|
| `calc_PBL_height.R` | Core function. PBLH from EC data (L, u*, H, T) following Kljun et al. 2015 App. B: diagnostic formulas for neutral/stable, Batchvarova & Gryning (1991) slab model for convective conditions. Stability classes from 1/L (Essa et al. 2006). | **Port and refactor.** Core science is here. |
| `PBL_validation.R` | Validation at FI-Hyy (SMEAR II), 2016: Doppler lidar MLH (.mat), ceilometer cloud base, SMEAR fluxes and meteo. Diurnal cycles (year/season/month), Bland–Altman, error boxplots. | Reuse the analysis design. Rewrite the I/O. |
| `PBL_ERA_comparison_test.R` | BE-Bra 2024: PBLH from ETC L2 footprint output vs ERA5 `boundary_layer_height`. Unfinished, does not run. | Reuse the idea (ERA5 as a second reference). |
| `FFP_sens_to_PBL.R` | Sensitivity of the FFP (Kljun) 80% fetch to PBLH. Depends on external `calc_footprint_FFP*.R`. | Possible later module. |
| `ERA-PBL-stations.txt` | 85 stations (ICOS + NEON) with lat/lon for ERA5 extraction. | Keep as site metadata. |
| `STabilityClasses_on_MOL.*` | Alternative stability classes from L (a–h). | Reference. |

## Issues found in `calc_PBL_height.R` and the validation script

1. **Temperature units (likely major).** The function expects `AT` in K, but `PBL_validation.R` passes `AT.24` in °C. In the B&G denominator the mechanical term `C·T·u*²` is then about 19× too small, so convective growth is overestimated. This may explain the TODO "PBL sovrastima molto".
2. **Order of gap filling and B&G.** `na.approx(maxgap = 2)` is applied to `PBL.h.0` *before* the B&G step. Short convective runs (≤ 2 steps) get linearly interpolated instead of modelled.
3. **Inconsistent convective trigger.** Stability classes use 1/L, but B&G is triggered by `H > 0`. The two criteria can disagree.
4. **Neutral formula `0.3 u*/|f|`.** At 61°N (f ≈ 1.28e-4 s⁻¹) and u* = 0.5 m/s it gives ≈ 1170 m. A likely source of the "higher in winter" TODO (windy, near-neutral).
5. **First time step.** At `i = 1`, `PBL.h.0[i-1]` is empty, so `if()` errors if `H[1] > 0`.
6. **`smooth.PBL = FALSE` breaks the function.** `PBL.h` is never assigned. The comment says 9-point smoothing while the code uses 5 points. `filter()` is masked if dplyr is loaded (use `stats::filter`).
7. **Pressure units.** Documented as Pa, passed as kPa (only harmless because `rho` is supplied).
8. **`gamma = 0.0025 K/m`.** The comment says ≈ 0.01. This needs a justified choice, ideally a parameter.
9. **Scripts don't run as-is.** `MLH.TS` is undefined (`MLH.df.TS` is the one created). dplyr and lubridate are used without being loaded. There are hardcoded `C:/WORK/DATA/...` paths, `install.packages()` inside a function, and `windows()` devices.
10. **Time conventions not handled.** GMT is assumed everywhere. The timestamp convention (start vs end of averaging period) and the 30-min vs 60-min ERA5 alignment are open TODOs.

## Open TODOs carried over from the legacy notes

- QC: exclude rain, low cloud and fog periods. MLH < 60 m is not resolvable by the lidar.
- Uncertainty of the method, of the lidar MLH, and of temporal aggregation.
- Bland–Altman per season. Error boxplots at fixed hours (Haefele 2014). 4-h means (Hellén et al. 2019).
- Is PBLH (from EC) comparable with lidar MLH at all? Address this in the design.
