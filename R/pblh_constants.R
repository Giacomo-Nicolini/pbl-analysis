# Physical constants and model parameters of PBLHmod (single source of truth).
# Spec: docs/specs/02_pblh_model.md section 1. Kljun et al. (2015, GMD 8, 3695-3713), App. B.

#' Physical constants and PBLH model parameters
#'
#' @return Named list:
#'   - `k` von Karman constant [-]: 0.4, as in Kljun et al. (2015) App. B.
#'   - `g` standard gravity [m s-2]: 9.80665 (CGPM 1901 conventional value).
#'   - `omega` Earth's rotation rate [rad s-1]: 7.2921159e-5, i.e. 2 pi per
#'     sidereal day (86164.1 s). The legacy code used 2 pi / 86400 (solar day,
#'     7.2722e-5), 0.27 % smaller.
#'   - `r_d` gas constant of dry air [J kg-1 K-1]: 287.05.
#'   - `cp` specific heat of dry air at constant pressure [J kg-1 K-1]: 1005.
#'   - `bg_A`, `bg_B`, `bg_C` [-]: 0.2, 2.5, 8; Batchvarova & Gryning (1991),
#'     as in Kljun et al. (2015) Eq. B5.
#'   - `nieuwstadt_a`, `nieuwstadt_b` [-]: 3.8, 2.28; Nieuwstadt (1981), Kljun
#'     et al. (2015) Eq. B1.
#'   - `c_n` [-]: 0.3; neutral-limit coefficient, Kljun et al. (2015) Eq. B2
#'     (Hanna & Chang 1993).
#'   - `c_nn` [-]: 1.36; neutral coefficient with a stratified free atmosphere,
#'     Kljun et al. (2015) Eq. B3 (Zilitinkevich et al. 2012).
#'   - `n_over_f_min`, `n_over_f_max` [-]: 10, 70; range of N/|f| in which the
#'     neutral formula B2 is valid (Kljun et al. 2015, App. B); outside it B3.
#'   - `p0` reference sea-level pressure [Pa]: 101325 (standard atmosphere).
#'   - `scale_height` isothermal scale height [m]: 8434 = r_d T / g with
#'     T ~ 288 K (e.g. Wallace & Hobbs 2006, Sect. 3.2); barometric fallback.
#'   - `h0_search_halfhours` [half-hours]: 4 (2 h) window before sunrise in
#'     which the initial height h0 is searched (spec 02, section 6).
#'   - `neutral_limit_eps` [-]: 1e-6; below this |u* / (f L)| the Nieuwstadt
#'     formula is replaced by its neutral limit B2 (loss of precision).
pblh_constants <- function() {
  list(
    k = 0.4,                 # von Karman constant; Kljun et al. 2015 App. B
    g = 9.80665,             # standard gravity [m s-2]
    omega = 7.2921159e-5,    # Earth's rotation rate, sidereal day [rad s-1] (legacy: 2*pi/86400)
    r_d = 287.05,            # gas constant of dry air [J kg-1 K-1]
    cp = 1005,               # specific heat of dry air, constant pressure [J kg-1 K-1]
    bg_A = 0.2,              # Batchvarova & Gryning 1991; Kljun et al. 2015 Eq. B5
    bg_B = 2.5,              # idem
    bg_C = 8,                # idem
    nieuwstadt_a = 3.8,      # Nieuwstadt 1981; Kljun et al. 2015 Eq. B1
    nieuwstadt_b = 2.28,     # idem
    c_n = 0.3,               # Kljun et al. 2015 Eq. B2 (Hanna & Chang 1993)
    c_nn = 1.36,             # Kljun et al. 2015 Eq. B3 (Zilitinkevich et al. 2012)
    n_over_f_min = 10,       # B2 valid for 10 < N/|f| < 70; Kljun et al. 2015 App. B
    n_over_f_max = 70,       # idem
    p0 = 101325,             # standard sea-level pressure [Pa]
    scale_height = 8434,     # isothermal scale height r_d * 288 K / g [m]; Wallace & Hobbs 2006
    h0_search_halfhours = 4, # h0 search window before sunrise (2 h); spec 02 section 6
    neutral_limit_eps = 1e-6 # switch B1 -> B2 when |u*/(f L)| < eps; spec 02 section 4
  )
}
