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
#'   `opts$c_s`) or `"neutral"` (B2).
#' - `L < 0` (unstable outside the convective window): B2.
#' - `L` or `ustar` missing, `L == 0`: `NA`.
#'
#' @param L Obukhov length [m].
#' @param ustar Friction velocity [m s-1].
#' @param f Coriolis parameter [s-1] (scalar or vector; sign ignored).
#' @param opts Output of [pblh_options()].
#' @return Tibble: `h_diag_m` [m], `method` (`"B1"`, `"B4"`, `"B2"`, `NA`).
h_diagnostic <- function(L, ustar, f, opts = pblh_options()) {
  n <- max(length(L), length(ustar))
  f <- rep_len(f, n)
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
    h[st] <- h_neutral(ustar[st], f[st])
    method[st] <- "B2"
  }
  h[un] <- h_neutral(ustar[un], f[un])
  method[un] <- "B2"
  method[is.na(h)] <- NA_character_
  dplyr::tibble(h_diag_m = h, method = method)
}
