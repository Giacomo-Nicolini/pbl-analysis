# Notes: Kljun et al. (2015), GMD 8, 3695–3713 — Appendix B and Sect. 6.2

Source: https://doi.org/10.5194/gmd-8-3695-2015 (open access, CC BY 3.0). Extracted 2026-10-12 for spec 02.
Symbols: L Obukhov length [m], u* friction velocity [m s-1], f = 2 Ω sin(φ) Coriolis parameter [s-1],
k = 0.4 von Kármán constant, g gravity, (w'θ')0 surface kinematic heat flux [K m s-1], T mean temperature [K],
γ potential-temperature gradient above the CBL [K m-1], N Brunt–Väisälä frequency above the BL [s-1].

## Neutral and stable conditions (diagnostic)

- **(B1)** Nieuwstadt (1981) interpolation formula, neutral to stable:
  `h = L/3.8 * ( -1 + sqrt(1 + 2.28 * u* / (f L)) )`
- **(B2)** Neutral limit (large L): `h_n = c_n * u* / |f|`, `c_n = 0.3` (Hanna & Chang 1993; Tennekes 1973).
  Strictly valid only if the free atmosphere is near-neutral, `10 < N/|f| < 70`.
- **(B3)** Otherwise, for very large L: `h_nn = c_nn * u* / sqrt(|f| N)`, `c_nn = 1.36` (Zilitinkevich et al. 2012).
- **(B4)** Strongly stable limit: `h_s = c_s * sqrt(u* L / |f|)` (Zilitinkevich 1972), applicability
  `4 <= u* L/|f| <= 100` (Zilitinkevich & Mironov 1996). `c_s ≈ 0.4` implied by (B1);
  `c_s ≈ 0.63` (Zilitinkevich & Mironov 1996); `c_s ≈ 0.3` (Zilitinkevich et al. 2007, LES).

## Convective conditions (prognostic)

- h cannot be diagnosed; integrate from **sunrise** (before the surface heat flux first becomes positive),
  with the initial height diagnosed by one of the expressions above.
- **(B5)** Batchvarova & Gryning (1991) slab model:
  `dh/dt = ((w'θ')0 / γ) * [ h^2 / ((1+2A) h - 2 B k L) + C u*^2 T / (γ g ((1+A) h - B k L)) ]^(-1)`
  with `A = 0.2`, `B = 2.5`, `C = 8`. May be solved iteratively, or using h(t_i) to compute the rate giving h(t_i+1).
- γ often unavailable: "can be approximated by a constant (e.g. γ ≈ 0.01 K m-1, a typical mid-latitude value adopted
  from Batchvarova and Gryning, 1991)". dh/dt is **quite sensitive** to γ.

## FFP sensitivity to h (Table 5, Sect. 6.1)

±20 % change in h changes footprint peak location by ≤ ~4 % (stable scenarios) and ≤ ~0.1 % (convective).
Overestimated h shifts the stable-case peak closer to the receptor.

## FFP validity limits (Sect. 6.2, Eq. 27)

- `20 z0 < zm < h_e`, with `h_e ≈ 0.8 h` (entrainment-layer height; Holtslag & Nieuwstadt 1986).
  `20 z0` is of the order of the roughness-sublayer height z*; may be adjusted for other z*.
- `-15.5 <= zm / L`.
- (zm here is the measurement height above the displacement height, z - d.)
- Source-area fraction R <= 0.9 recommended (R = 0.8 usually sufficient).
- The official FFP code also uses a u* threshold (check `calc_footprint_FFP.R`: typically u* > 0.1 m s-1)
  — verify in the code before writing the phase-5 spec.
