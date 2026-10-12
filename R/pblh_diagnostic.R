# Diagnostic PBL height for neutral and stable conditions.
# Kljun et al. (2015, GMD 8, 3695-3713) App. B, Eqs. B1, B2, B4
# (notes: docs/references/kljun2015_appB_notes.md). Spec 02 section 4.

#' Neutral PBL height (Kljun et al. 2015, Eq. B2)
#'
#' `h = c_n u* / |f|` (Hanna & Chang 1993; Tennekes 1973).
#'
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (sign ignored).
#' @param c_n Coefficient [-]; default `pblh_constants()$c_n` (0.3).
#' @return PBL height [m].
h_neutral <- function(ustar, f, c_n = pblh_constants()$c_n) {
  c_n * ustar / abs(f) # Kljun et al. 2015 Eq. B2
}

#' Brunt-Vaisala frequency above the boundary layer
#'
#' `N = sqrt(g gamma / T)` (dry static stability, gamma = potential-temperature
#' gradient), as used with Kljun et al. (2015) Eq. B3.
#'
#' @param gamma Potential-temperature gradient [K m-1].
#' @param T_K Air temperature [K].
#' @return Brunt-Vaisala frequency [s-1] (`NA` for gamma < 0 or missing input).
brunt_vaisala <- function(gamma, T_K) {
  x <- pblh_constants()$g * gamma / T_K
  x[!is.na(x) & x < 0] <- NA_real_
  sqrt(x)
}

#' Neutral PBL height with a stratified free atmosphere (Kljun et al. 2015, Eq. B3)
#'
#' `h = c_nn u* / sqrt(|f| N)`, `c_nn = 1.36` (Zilitinkevich et al. 2012).
#'
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (sign ignored).
#' @param N Brunt-Vaisala frequency above the boundary layer [s-1].
#' @param c_nn Coefficient [-]; default `pblh_constants()$c_nn` (1.36).
#' @return PBL height [m].
h_neutral_nn <- function(ustar, f, N, c_nn = pblh_constants()$c_nn) {
  c_nn * ustar / sqrt(abs(f) * N) # Kljun et al. 2015 Eq. B3
}

#' Neutral PBL height by the selected method
#'
#' `"B2"`: Eq. B2 always. `"B3"`: Eq. B3 always. `"auto"`: Eq. B2 where
#' `10 < N/|f| < 70` (its validity range, Kljun et al. 2015 App. B), Eq. B3
#' otherwise, with `N = sqrt(g gamma / T)`.
#'
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (sign ignored).
#' @param neutral_method `"auto"`, `"B2"` or `"B3"`.
#' @param T_K Air temperature [K] (needed for `"auto"` and `"B3"`).
#' @param gamma Potential-temperature gradient above the BL [K m-1] (idem).
#' @return Tibble: `h` [m], `method` (`"B2"`, `"B3"`, `NA` if undecidable).
h_neutral_select <- function(ustar, f, neutral_method = "auto", T_K = NULL, gamma = NULL) {
  const <- pblh_constants()
  n <- length(ustar)
  f <- rep_len(f, n)
  if (neutral_method == "B2") {
    return(dplyr::tibble(h = h_neutral(ustar, f), method = "B2"))
  }
  if (is.null(T_K) || is.null(gamma)) {
    cli::cli_abort("{.arg T_K} and {.arg gamma} are needed for neutral_method = {.val {neutral_method}}.")
  }
  N <- brunt_vaisala(rep_len(gamma, n), rep_len(T_K, n))
  use_b2 <- if (neutral_method == "auto") {
    nf <- N / abs(f)
    nf > const$n_over_f_min & nf < const$n_over_f_max
  } else {
    rep(FALSE, n)
  }
  h <- ifelse(use_b2, h_neutral(ustar, f), h_neutral_nn(ustar, f, N))
  method <- ifelse(use_b2, "B2", "B3")
  method[is.na(use_b2) | is.na(h)] <- NA_character_
  h[is.na(method)] <- NA_real_
  dplyr::tibble(h = h, method = method)
}

#' Stable PBL height, Nieuwstadt (1981) (Kljun et al. 2015, Eq. B1)
#'
#' `h = L/a (-1 + sqrt(1 + b u* / (|f| L)))`, `a = 3.8`, `b = 2.28`, for
#' `L > 0`. Evaluated in the algebraically identical form
#' `h = (b u*/|f|) / (a (1 + sqrt(1 + b x)))`, `x = u* / (|f| L)`, which avoids
#' cancellation; for `|x| < 1e-6` the neutral limit (Eq. B2, `c_n = b / (2a) = 0.3`)
#' is used (spec 02).
#'
#' @param L Obukhov length [m]; `L <= 0` gives `NA`; `Inf` gives the neutral limit.
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (sign ignored).
#' @return PBL height [m].
h_nieuwstadt <- function(L, ustar, f) {
  const <- pblh_constants()
  a <- const$nieuwstadt_a
  b <- const$nieuwstadt_b
  n <- max(length(L), length(ustar), length(f))
  L <- rep_len(L, n)
  ustar <- rep_len(ustar, n)
  af <- rep_len(abs(f), n)
  h <- rep(NA_real_, n)
  ok <- !is.na(L) & !is.na(ustar) & !is.na(af) & L > 0
  x <- ustar[ok] / (af[ok] * L[ok])
  # Kljun et al. 2015 Eq. B1: L/a * (sqrt(1 + b x) - 1) = (b u*/|f|) / (a (1 + sqrt(1 + b x)))
  hk <- (b * ustar[ok] / af[ok]) / (a * (1 + sqrt(1 + b * x)))
  neutral <- abs(x) < const$neutral_limit_eps
  hk[neutral] <- h_neutral(ustar[ok][neutral], af[ok][neutral]) # neutral limit, Eq. B2
  h[ok] <- hk
  h
}

#' Strongly stable PBL height, Zilitinkevich (1972) (Kljun et al. 2015, Eq. B4)
#'
#' `h = c_s sqrt(u* L / |f|)` for `L > 0` (finite), else `NA`.
#'
#' @param L Obukhov length [m].
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (sign ignored).
#' @param c_s Coefficient [-]; default 0.4 (implied by Eq. B1; 0.63 Zilitinkevich
#'   & Mironov 1996; 0.3 Zilitinkevich et al. 2007).
#' @return PBL height [m].
h_zilitinkevich <- function(L, ustar, f, c_s = 0.4) {
  n <- max(length(L), length(ustar), length(f))
  L <- rep_len(L, n)
  ustar <- rep_len(ustar, n)
  af <- rep_len(abs(f), n)
  h <- rep(NA_real_, n)
  ok <- !is.na(L) & !is.na(ustar) & is.finite(L) & L > 0
  h[ok] <- c_s * sqrt(ustar[ok] * L[ok] / af[ok]) # Kljun et al. 2015 Eq. B4
  h
}

#' Diagnostic PBL height with method label
#'
#' - `L > 0`: `opts$stable_method` `"nieuwstadt"` (B1), `"zilitinkevich"` (B4,
#'   `opts$c_s`) or `"neutral"` (neutral formula by `opts$neutral_method`).
#' - `L < 0` (unstable outside the convective window): neutral formula by
#'   `opts$neutral_method` (B2, B3 or auto switch, [h_neutral_select()]).
#' - `L` or `ustar` missing, `L == 0`: `NA`.
#'
#' @param L Obukhov length [m].
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (scalar or vector; sign ignored).
#' @param opts Output of [pblh_options()].
#' @param T_K Air temperature [K]; needed unless `opts$neutral_method == "B2"`.
#' @param gamma Potential-temperature gradient above the BL [K m-1], scalar or
#'   vector; default `opts$gamma`.
#' @return Tibble: `h_diag_m` [m], `method` (`"B1"`, `"B2"`, `"B3"`, `"B4"`, `NA`).
h_diagnostic <- function(L, ustar, f, opts = pblh_options(), T_K = NULL, gamma = opts$gamma) {
  n <- max(length(L), length(ustar))
  L <- rep_len(L, n)
  ustar <- rep_len(ustar, n)
  f <- rep_len(f, n)
  if (!is.null(T_K)) T_K <- rep_len(T_K, n)
  gamma <- rep_len(gamma, n)
  neutral <- function(idx) {
    h_neutral_select(ustar[idx], f[idx], opts$neutral_method, T_K[idx], gamma[idx])
  }
  h <- rep(NA_real_, n)
  method <- rep(NA_character_, n)
  ok <- !is.na(L) & !is.na(ustar) & L != 0
  st <- ok & L > 0
  un <- ok & L < 0
  if (opts$stable_method == "nieuwstadt") {
    h[st] <- h_nieuwstadt(L[st], ustar[st], f[st])
    method[st] <- "B1"
  } else if (opts$stable_method == "zilitinkevich") {
    h[st] <- h_zilitinkevich(L[st], ustar[st], f[st], c_s = opts$c_s)
    method[st] <- "B4"
  } else {
    nt <- neutral(st)
    h[st] <- nt$h
    method[st] <- nt$method
  }
  nt <- neutral(un)
  h[un] <- nt$h
  method[un] <- nt$method
  method[is.na(h)] <- NA_character_
  dplyr::tibble(h_diag_m = h, method = method)
}
