# Analysis of the footprint sensitivity to the PBLH source (spec 08, section 4). Pure functions; the script is
# scripts/08_footprint_compare.R. Footprint model: Kljun et al. (2015), `Kljun2015`.

#' Log ratio of footprint distances
#'
#' @param x_case,x_ref Distances [m] (e.g. x_80 with the case h and with PBLHmod h).
#' @return ln(x_case / x_ref) [-]; `NA` where either is missing or not positive.
fp_log_ratio <- function(x_case, x_ref) {
  ok <- !is.na(x_case) & !is.na(x_ref) & x_case > 0 & x_ref > 0
  out <- rep(NA_real_, length(x_case))
  out[ok] <- log(x_case[ok] / x_ref[ok])
  out
}

#' Elasticity of a footprint distance with respect to h
#'
#' Finite-difference elasticity d ln x / d ln h = ln(x_case / x_ref) / ln(h_case / h_ref)
#' (spec 08, section 4).
#'
#' @inheritParams fp_log_ratio
#' @param h_case,h_ref Boundary layer heights [m].
#' @param min_dlnh Smallest |ln(h_case / h_ref)| for which the elasticity is computed [-].
#' @return Elasticity [-]; `NA` where |ln h ratio| < `min_dlnh` or inputs are missing.
fp_elasticity <- function(x_case, x_ref, h_case, h_ref, min_dlnh = 1e-6) {
  lx <- fp_log_ratio(x_case, x_ref)
  lh <- fp_log_ratio(h_case, h_ref)
  out <- lx / lh
  out[is.na(lh) | abs(lh) < min_dlnh] <- NA_real_
  out
}

#' Summary of log ratios (and elasticities) in one stratum
#'
#' @param lr Log ratios [-].
#' @param el Elasticities [-] (optional).
#' @param dta Target-area fraction differences [percentage points] (optional).
#' @return One-row tibble: `n`, `median_lr`, `p25_lr`, `p75_lr` [-], `share_gt10`,
#'   `share_gt25` [0-1] (rows with |x_case / x_ref - 1| > 10 % / 25 %), `median_elasticity`
#'   [-], `median_dta_pp`, `p90_abs_dta_pp` [percentage points].
fp_summary <- function(lr, el = NULL, dta = NULL) {
  ok <- !is.na(lr)
  q <- function(v, p) if (any(!is.na(v))) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
  dev <- abs(exp(lr[ok]) - 1)
  dplyr::tibble(
    n = sum(ok), median_lr = q(lr, 0.5), p25_lr = q(lr, 0.25), p75_lr = q(lr, 0.75),
    share_gt10 = if (sum(ok)) mean(dev > 0.10) else NA_real_,
    share_gt25 = if (sum(ok)) mean(dev > 0.25) else NA_real_,
    median_elasticity = if (is.null(el)) NA_real_ else q(el, 0.5),
    median_dta_pp = if (is.null(dta)) NA_real_ else q(dta, 0.5),
    p90_abs_dta_pp = if (is.null(dta)) NA_real_ else q(abs(dta), 0.9)
  )
}

#' Stability class from zm/L
#'
#' Reporting classes (spec 08 analysis; boundaries chosen here, flagged for G):
#' unstable zm/L < -0.1, neutral -0.1 <= zm/L <= 0.1, stable zm/L > 0.1.
#'
#' @param zm Measurement height above d [m].
#' @param L Obukhov length [m].
#' @return Factor `unstable`, `neutral`, `stable` (NA where missing).
zl_class <- function(zm, L) {
  z <- zm / L
  cl <- ifelse(z < -0.1, "unstable", ifelse(z > 0.1, "stable", "neutral"))
  factor(cl, levels = c("unstable", "neutral", "stable"))
}

#' zm/h class
#'
#' Reporting classes (boundaries chosen here, flagged for G): < 0.02, 0.02-0.05,
#' 0.05-0.1, 0.1-0.2, 0.2-0.8, >= 0.8 (the last fails Kljun et al. 2015 Eq. 27).
#'
#' @param zm Measurement height above d [m].
#' @param h Boundary layer height [m].
#' @return Factor of classes (NA where missing).
zmh_class <- function(zm, h) {
  cut(zm / h, c(0, 0.02, 0.05, 0.1, 0.2, 0.8, Inf),
      labels = c("<0.02", "0.02-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.8", ">=0.8"), right = FALSE)
}

#' Share of a source region inside the target area
#'
#' @param region Integer cell indices of the region (e.g. the 80 % source region).
#' @param mask Target-area mask (logical matrix on the footprint grid).
#' @return Share of the region's cells inside the target area [0-1].
region_share_inside <- function(region, mask) {
  if (length(region) == 0) return(NA_real_)
  mean(c(mask)[region])
}

#' ICOS target-area criteria (Rebmann et al. 2018)
#'
#' A site-design criterion. Rebmann et al. (2018, `Rebmann2018`), p. 475: a flux
#' contribution from the area of interest of at least 90 % (Goeckede et al. 2004 "high
#' quality"), also under stable conditions, for at least 50 % of the time; p. 487: a
#' cumulative footprint (CIF) of at least 80 %, integrated from the peak along levels of
#' equal footprint value, for an acceptable measurement. Operationalised here: per half-hour,
#' `ta90` = target-area fraction >= 0.9 and `cif80` = target-area fraction >= 0.8; the
#' 90 % criterion passes for a set of half-hours if `ta90` holds on >= 50 % of them
#' (checked on all, night and stable rows separately).
#'
#' @param f_ta Footprint fraction inside the target area per half-hour [0-1].
#' @return One-row tibble: `n`, `share_ta90`, `share_cif80` [0-1], `pass_ta90_50pct` [logical].
rebmann_criteria <- function(f_ta) {
  ok <- !is.na(f_ta)
  n <- sum(ok)
  s90 <- if (n) mean(f_ta[ok] >= 0.9) else NA_real_
  dplyr::tibble(n = n, share_ta90 = s90, share_cif80 = if (n) mean(f_ta[ok] >= 0.8) else NA_real_,
                pass_ta90_50pct = s90 >= 0.5)
}

#' Kljun et al. (2015) Table 5: sensitivity to h +-20 %
#'
#' Absolute changes of the peak location x_max and of the crosswind-integrated peak
#' value f^y_max for h changed by 20 %, per stability group (scenarios 10-12
#' convective, 13-15 neutral, 16-18 stable). The signs are not legible in the PDF text;
#' magnitudes only.
#'
#' @return Tibble: `zl_class`, `kljun_abs_dxmax_pct`, `kljun_abs_dfymax_pct` [%].
kljun2015_table5_h20 <- function() {
  dplyr::tibble(zl_class = c("unstable", "neutral", "stable"),
                kljun_abs_dxmax_pct = c(0.3, 0.7, 3.7), kljun_abs_dfymax_pct = c(0.3, 0.7, 3.6))
}

#' Anchor: footprint peak change for h scaled by a factor
#'
#' With umean (Kljun et al. 2015, Eq. 6) x scales with zm / (1 - zm/h) umean / u*, and the
#' crosswind-integrated footprint with its inverse, so the peak value changes by
#' x_peak,ref / x_peak,case - 1.
#'
#' @param x_peak_case,x_peak_ref Peak distances [m] with scaled and original h.
#' @return Tibble: `dxmax_pct`, `dfymax_pct` [%].
fp_peak_change <- function(x_peak_case, x_peak_ref) {
  dplyr::tibble(dxmax_pct = 100 * (x_peak_case / x_peak_ref - 1), dfymax_pct = 100 * (x_peak_ref / x_peak_case - 1))
}
