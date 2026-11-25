# Methods note: γ from the station's own data (EC-derived γ, PBLHmod-γEC; spec 17)

Draft for the paper's methods section (G and co-authors write the final text). Code: `R/gamma_ec.R`,
`scripts/17a_gamma_ec.R`; options `gamma_source = "ec_clim"` / `"ec_daily"` of `run_pblh_model()` (spec 02
amendment 15). References by BibTeX key where the key is in `docs/references/`.

## Principle

The convective growth of PBLHmod is the zero-order slab model of Batchvarova & Gryning (`Batchvarova1991`), as given in
Kljun et al. (`Kljun2015`, Eq. B5). In that model the mixed layer warms as

  dθ/dt = (1 + A) (w'θ')_s / h,

with the entrainment heat flux equal to −A times the surface flux (A = 0.2), while it grows into a layer with
potential-temperature gradient γ. Given the heat put in by the surface, the rise of θ at the tower therefore carries γ:
a steeper γ keeps the layer shallower, and the same heat warms a shallower layer faster. γ is fitted from the
station's own sensible heat flux and tower air temperature, without ERA5.

## Data and fit window

θ = T (100000 Pa / p)^(R_d / c_p) (dry potential temperature; R_d = 287.05 J kg⁻¹ K⁻¹, c_p = 1005 J kg⁻¹ K⁻¹), with T
and p of the model forcing (gap-filled FLUXNET `TA_F`, `PA_F`; METEO after the end of FLUXNET). Half-hours whose
temperature comes from the sonic temperature are not used. The fit window of a day starts at the first half-hour with
a positive kinematic heat flux in the B&G window of the default run and ends at the end of that window, at the
half-hour of the day's maximum θ if earlier, and at local solar noon + 3 h if earlier. A window needs at least 6 valid
half-hours and no gap in the forcing.

## Method A (primary, model-consistent)

For a candidate γ, h(t) is integrated from the day's initial height h0 with the B&G growth equation exactly as in
PBLHmod (4th-order Runge–Kutta, 30 sub-steps per half-hour, h held during short negative-flux spells and where the
flux is positive under stable L), and the mixed-layer temperature θ_mod(t) = θ0 + ∫ (1 + A) (w'θ')_s / h dt on the
same sub-steps. Half-hourly values are interval means. θ0 is the mean residual (analytic least squares); γ minimises
Σ (θ_obs − θ_mod)² by Brent's method on log γ within [0.001, 0.02] K m⁻¹ (the bounds of the ERA5 γ, spec 11). An
optimum at a bound gives no value (flagged). Reported per day: γ_A, θ0, RMSE, R², n.

## Method B (closed form, sensitivity)

Encroachment from h0 = 0 without the mechanical terms (Tennekes 1973 — not in the project bib; to be checked and added
by G — and `Batchvarova1991`): dh/dt = (1 + 2A) (w'θ')_s / (γ h) gives h² = 2 (1 + 2A) Q / γ with
Q = ∫ (w'θ')_s dt [K m] from the window start, and dθ/dt = (1 + A) (w'θ')_s / h = γ (1 + A)/(1 + 2A) dh/dt gives
θ − θ_s = γ (1 + A)/(1 + 2A) h. Hence θ_obs = a + s √Q with s = √γ (1 + A) √(2 / (1 + 2A)), and

  γ_B = s² (1 + 2A) / (2 (1 + A)²),

fitted by linear least squares. The factor (1 + 2A)/(1 + A)² is 0.97 for A = 0.2 relative to A = 0: γ_B hardly
depends on the entrainment ratio. Method B ignores h0 and the mechanical terms, so it is a sensitivity estimate.

## Day selection (flags; nothing deleted)

A day is used when it is clear (mean SW_IN / SW_IN_POT over the window ≥ 0.6, FLUXNET `SW_IN_F`, `SW_IN_POT`; days
without these are flagged `clear_unknown` and not used), has an initial height, no forcing gap in the window or the
model's gap flag, a positive heat flux, a θ rise ≥ 1 K over the window, R²_A ≥ 0.8 and γ_A inside the bounds. The
climatology is the site × calendar-month median (and p25, p75) of γ_A over the used days; a month with fewer than 5
used days has no value, and PBLHmod then uses the constant γ (flagged `constant_fallback`).

## H bias

A γ fitted with the same H that drives the model also absorbs a systematic error of H, e.g. energy-balance
non-closure: with H too small, the observed θ rise needs a shallower layer, i.e. a larger γ. This is reported
(γ_EC / γ_ERA5 against the daytime energy-balance ratio, spec 17 D) and not corrected (no H post-correction,
G 2026-11-12).

## Limitation seen in the first local runs (2026-11-25, four sites, test)

On most fitted days the optimum lies at the upper bound (0.02 K m⁻¹): the morning θ rise at the tower is faster
than a mixed layer that grows from the night-time height into the free-atmosphere γ allows. Erosion of the nocturnal
inversion and surface-layer warming early in the window are the likely reasons; the used days are therefore a
selection, and their γ_A (medians 0.007–0.014 K m⁻¹) describe the morning inversion as much as the free atmosphere.
