# Convective PBL growth: Batchvarova & Gryning (1991) slab model,
# as given in Kljun et al. (2015, GMD 8, 3695-3713) Eq. B5. Spec 02 section 5.

#' Growth rate of the convective boundary layer (Kljun et al. 2015, Eq. B5)
#'
#' `dh/dt = (wtheta / gamma) / ( h^2 / ((1+2A) h - 2 B k L) +
#'   C u*^2 T / (gamma g ((1+A) h - B k L)) )`, with `A = 0.2`, `B = 2.5`,
#' `C = 8` (Batchvarova & Gryning 1991). Returns `NA` where either denominator
#' `(1+2A) h - 2 B k L` or `(1+A) h - B k L` is <= 0, or an input is missing.
#'
#' @param h Boundary-layer height [m].
#' @param wtheta Surface kinematic heat flux [K m s-1].
#' @param ustar Friction velocity [m s-1].
#' @param L Obukhov length [m].
#' @param T_K Air temperature [K].
#' @param gamma Potential-temperature gradient above the boundary layer [K m-1].
#' @param const Output of [pblh_constants()].
#' @return Growth rate dh/dt [m s-1].
bg_rate <- function(h, wtheta, ustar, L, T_K, gamma, const = pblh_constants()) {
  A <- const$bg_A
  B <- const$bg_B
  C <- const$bg_C
  k <- const$k
  d1 <- (1 + 2 * A) * h - 2 * B * k * L # Kljun et al. 2015 Eq. B5, first denominator
  d2 <- (1 + A) * h - B * k * L         # second denominator
  rate <- (wtheta / gamma) / (h^2 / d1 + C * ustar^2 * T_K / (gamma * const$g * d2))
  rate[is.na(d1) | is.na(d2) | d1 <= 0 | d2 <= 0] <- NA_real_
  rate
}

#' Integrate Eq. B5 over one averaging interval
#'
#' `method = "rk4"` (default, PBLHmod v2): classical 4th-order Runge-Kutta with
#' `n_sub` explicit sub-steps, forcing held constant over the interval.
#' `method = "euler"` (legacy `calc_PBL_height()`, used only by the spec 03
#' attribution; spec 03 Amendment 1): one explicit Euler step over the whole
#' interval, `h_end = h0 + dt * rate(h0)`, and `h_mean = h_end` (the legacy
#' reports the end-of-interval height). If [bg_rate()] returns `NA` at any
#' stage, `h` is held constant from that (sub-)step on and `ok = FALSE` (the
#' caller flags it).
#'
#' @param h0 Height at the start of the interval [m].
#' @param wtheta Surface kinematic heat flux [K m s-1].
#' @param ustar Friction velocity [m s-1].
#' @param L Obukhov length [m].
#' @param T_K Air temperature [K].
#' @param gamma Potential-temperature gradient above the boundary layer [K m-1] (scalar).
#' @param dt Interval length [s]; default 1800.
#' @param n_sub Number of sub-steps [-]; default 30 (60 s); ignored by `"euler"`.
#' @param const Output of [pblh_constants()].
#' @param method `"rk4"` (default) or `"euler"`.
#' @return List: `h_end` height at the end of the interval [m], `h_mean`
#'   interval mean height (trapezoidal mean of the `n_sub + 1` nodes; `"euler"`:
#'   `h_end`) [m], `ok` [logical, `FALSE` if the rate was undefined].
bg_step <- function(h0, wtheta, ustar, L, T_K, gamma, dt = 1800, n_sub = 30,
                    const = pblh_constants(), method = "rk4") {
  # Precompute the constant parts of Eq. B5 for speed (scalar forcing).
  A <- const$bg_A
  B <- const$bg_B
  bkl <- B * const$k * L
  num <- wtheta / gamma
  c2 <- const$bg_C * ustar^2 * T_K / (gamma * const$g)
  rate <- function(h) {
    d1 <- (1 + 2 * A) * h - 2 * bkl
    d2 <- (1 + A) * h - bkl
    if (is.na(d1) || is.na(d2) || d1 <= 0 || d2 <= 0) return(NA_real_)
    num / (h^2 / d1 + c2 / d2) # Kljun et al. 2015 Eq. B5
  }
  if (identical(method, "euler")) {
    r <- rate(h0)
    if (is.na(r)) return(list(h_end = h0, h_mean = h0, ok = FALSE))
    h1 <- h0 + dt * r # explicit Euler, one step per interval (legacy calc_PBL_height())
    return(list(h_end = h1, h_mean = h1, ok = TRUE))
  }
  if (!identical(method, "rk4")) cli::cli_abort("{.arg method} must be {.val rk4} or {.val euler}.")
  hs <- dt / n_sub
  nodes <- numeric(n_sub + 1)
  nodes[1] <- h0
  h <- h0
  ok <- TRUE
  for (i in seq_len(n_sub)) {
    if (ok) {
      k1 <- rate(h)
      k2 <- if (!is.na(k1)) rate(h + hs / 2 * k1) else NA_real_
      k3 <- if (!is.na(k2)) rate(h + hs / 2 * k2) else NA_real_
      k4 <- if (!is.na(k3)) rate(h + hs * k3) else NA_real_
      if (is.na(k4)) {
        ok <- FALSE
      } else {
        h <- h + hs / 6 * (k1 + 2 * k2 + 2 * k3 + k4) # RK4
      }
    }
    nodes[i + 1] <- h
  }
  list(
    h_end = h,
    h_mean = (sum(nodes) - (nodes[1] + nodes[n_sub + 1]) / 2) / n_sub, # trapezoidal mean
    ok = ok
  )
}
