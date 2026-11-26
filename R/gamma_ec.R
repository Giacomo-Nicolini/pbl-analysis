# Spec 17: gamma from the station's own data (EC-derived gamma, "PBLHmod-gammaEC"). No ERA5.
#
# Principle (zero-order slab model behind Batchvarova & Gryning 1991; Kljun et al. 2015 Eq. B5): the mixed layer warms
# as d theta / dt = (1 + A) w'theta'_s / h (entrainment flux = -A x surface flux) while it grows into a layer of
# gradient gamma, so the rise of theta at the tower, given the heat put in, carries gamma. A gamma fitted with the
# same H that drives the model also absorbs a systematic error of H (e.g. energy-balance non-closure): reported, not
# corrected (no H post-correction, G 2026-11-12).
#
# Method A (primary, model-consistent): h(t) integrated with bg_step() from the day's h0 exactly as run_pblh_model()
#   in the convective window (same sub-steps, same holds), theta_mod(t) = theta0 + int (1 + A) wtheta / h dt on the
#   same sub-steps; least squares over theta0 (analytic) and gamma (Brent on log gamma within the spec-11 bounds).
# Method B (closed form, sensitivity): encroachment with h0 = 0 (Tennekes 1973, not in the project bib; Batchvarova &
#   Gryning 1991 without the mechanical terms): dh/dt = (1 + 2A) wtheta / (gamma h) gives h^2 = 2 (1 + 2A) Q / gamma
#   with Q = int wtheta dt [K m]; d theta / dt = (1 + A) wtheta / h = gamma (1 + A) / (1 + 2A) dh/dt gives
#   theta - theta_s = gamma (1 + A) / (1 + 2A) h. Hence theta = a + s sqrt(Q) with s = sqrt(gamma) (1 + A)
#   sqrt(2 / (1 + 2A)), i.e. gamma_B = s^2 (1 + 2A) / (2 (1 + A)^2); the factor (1 + 2A) / (1 + A)^2 = 0.97 for
#   A = 0.2 (almost independent of A).

#' Options of the gamma_EC fit (spec 17)
#'
#' @param end_offset_h Fit window ends at local solar noon + this [h] (default 3).
#' @param start_offset_h The fit starts this long after the first positive heat flux [h, multiple of 0.5; default 0]
#'   (D1 pre-test, G 2026-11-26). h is still integrated from the first positive heat flux (model-consistent h at the
#'   fit start for each candidate gamma); only the theta residuals of the earlier rows are left out.
#' @param min_n Smallest number of valid half-hours in the fit window [-] (default 6).
#' @param clear_min Smallest mean SW_IN / SW_IN_POT over the window (clear sky) [-] (default 0.6).
#' @param rise_min Smallest theta rise over the window [K] (default 1).
#' @param r2_min Smallest R2 of Method A [-] (default 0.8).
#' @param gamma_min,gamma_max Bounds of gamma [K m-1] (spec 11: 0.001, 0.02; [gamma_options()]).
#' @param min_days Smallest number of used days of a calendar month in the climatology [-] (default 5).
#' @return Named list.
gamma_ec_options <- function(end_offset_h = 3, min_n = 6L, clear_min = 0.6, rise_min = 1, r2_min = 0.8,
                             gamma_min = gamma_options()$gamma_min, gamma_max = gamma_options()$gamma_max,
                             min_days = 5L, start_offset_h = 0) {
  o <- list(end_offset_h = end_offset_h, start_offset_h = start_offset_h, min_n = as.integer(min_n), clear_min = clear_min, rise_min = rise_min,
            r2_min = r2_min, gamma_min = gamma_min, gamma_max = gamma_max, min_days = as.integer(min_days))
  num_ok <- vapply(o, function(v) is.numeric(v) && length(v) == 1 && is.finite(v), logical(1))
  if (!all(num_ok)) cli::cli_abort("gamma_EC options must be single finite numbers: {.val {names(o)[!num_ok]}}.")
  if (!(o$gamma_min > 0 && o$gamma_max > o$gamma_min)) cli::cli_abort("Need 0 < gamma_min < gamma_max.")
  if (o$start_offset_h < 0 || abs(o$start_offset_h * 2 - round(o$start_offset_h * 2)) > 1e-9)
    cli::cli_abort("{.field start_offset_h} must be >= 0 and a multiple of 0.5 h.")
  if (o$min_n < 3 || o$min_days < 1) cli::cli_abort("{.field min_n} must be >= 3 and {.field min_days} >= 1.")
  o
}

#' gamma_EC options from the configuration (`gamma_ec_*` keys of `config.yml`)
#'
#' @param cfg Configuration list ([get_cfg()]); missing keys take the defaults of [gamma_ec_options()].
#' @return Named list ([gamma_ec_options()]).
gamma_ec_options_cfg <- function(cfg) {
  map <- c(end_offset_h = "gamma_ec_end_offset_h", start_offset_h = "gamma_ec_start_offset_h", min_n = "gamma_ec_min_n", clear_min = "gamma_ec_clear_min",
           rise_min = "gamma_ec_rise_min", r2_min = "gamma_ec_r2_min", min_days = "gamma_ec_min_days")
  args <- Filter(Negate(is.null), lapply(map, function(k) cfg[[k]]))
  do.call(gamma_ec_options, args)
}

#' Dry potential temperature
#'
#' theta = T (p0 / p)^(r_d / cp), p0 = 100000 Pa (as spec 11 A; constants of [pblh_constants()]).
#'
#' @param T_K Air temperature [K]; @param p_Pa Air pressure [Pa].
#' @param const Output of [pblh_constants()].
#' @return Potential temperature [K].
theta_dry <- function(T_K, p_Pa, const = pblh_constants()) {
  T_K * (100000 / p_Pa)^(const$r_d / const$cp)
}

#' Mixed-layer height and temperature path over a fit window (Method A)
#'
#' Same rules as [run_pblh_model()] after the first positive heat flux: [bg_step()] where `wtheta > 0` and `L <= 0`;
#' h held where `wtheta <= 0` (short negative spell) or `L > 0` (sign inconsistent), and from an undefined rate on.
#' The mixed-layer temperature increment C(t) = int (1 + A) wtheta / h dt is integrated on the same sub-step nodes
#' (trapezoid); both are reported as interval means (trapezoidal mean of the nodes, as `h_mean`).
#'
#' @param h0 Height at the start of the first row [m].
#' @param wtheta,ustar,L,T_K Forcing of the window rows [K m s-1, m s-1, m, K] (no `NA`).
#' @param gamma Potential-temperature gradient [K m-1].
#' @param n_sub RK4 sub-steps per half-hour [-]; @param const Output of [pblh_constants()].
#' @return List: `h_mean` [m], `c_mean` [K, increment from the window start], `h_end` [m].
bg_theta_path <- function(h0, wtheta, ustar, L, T_K, gamma, n_sub = 30, const = pblh_constants()) {
  n <- length(wtheta)
  hm <- cm <- numeric(n)
  h <- h0
  cum <- 0
  dts <- 1800 / n_sub
  a1 <- 1 + const$bg_A
  m <- n_sub + 1
  for (i in seq_len(n)) {
    if (wtheta[i] > 0 && L[i] <= 0) {
      st <- bg_step(h, wtheta[i], ustar[i], L[i], T_K[i], gamma, dt = 1800, n_sub = n_sub, const = const,
                    nodes = TRUE)
      nd <- st$nodes
      h <- st$h_end
    } else {
      nd <- rep(h, m) # held (run_pblh_model(): negative spell or wtheta > 0 with L > 0)
    }
    inv <- 1 / nd
    cn <- cum + c(0, cumsum(a1 * wtheta[i] * dts * (inv[-1] + inv[-m]) / 2)) # d theta/dt = (1 + A) wtheta / h
    cm[i] <- (sum(cn) - (cn[1] + cn[m]) / 2) / n_sub
    hm[i] <- (sum(nd) - (nd[1] + nd[m]) / 2) / n_sub
    cum <- cn[m]
  }
  list(h_mean = hm, c_mean = cm, h_end = h)
}

#' Method A: gamma from the theta rise with the B&G model
#'
#' @param theta_obs Observed potential temperature of the window rows [K] (`NA` = not used in the fit).
#' @inheritParams bg_theta_path
#' @param opts Output of [gamma_ec_options()].
#' @return List: `gamma` [K m-1; `NA` at a bound], `gamma_opt` [K m-1, optimum also at a bound], `theta0` [K],
#'   `rmse` [K], `r2` [-], `n` [-], `at_bound` [logical], `bound` (`low`, `high`, `NA`).
gamma_ec_fit_a <- function(theta_obs, h0, wtheta, ustar, L, T_K, opts = gamma_ec_options(), n_sub = 30,
                           const = pblh_constants()) {
  ok <- is.finite(theta_obs)
  y <- theta_obs[ok]
  sse <- function(lg) {
    cm <- bg_theta_path(h0, wtheta, ustar, L, T_K, exp(lg), n_sub, const)$c_mean[ok]
    r <- y - cm
    sum((r - mean(r))^2) # theta0 = mean residual (analytic)
  }
  lims <- log(c(opts$gamma_min, opts$gamma_max))
  op <- stats::optimize(sse, lims, tol = 1e-5) # Brent (golden section + parabolic interpolation)
  side <- if ((op$minimum - lims[1]) < 1e-3) "low" else if ((lims[2] - op$minimum) < 1e-3) "high" else NA_character_
  at_bound <- !is.na(side)
  g <- exp(op$minimum)
  cm <- bg_theta_path(h0, wtheta, ustar, L, T_K, g, n_sub, const)$c_mean[ok]
  th0 <- mean(y - cm)
  sst <- sum((y - mean(y))^2)
  list(gamma = if (at_bound) NA_real_ else g, gamma_opt = g, theta0 = th0, rmse = sqrt(op$objective / length(y)),
       r2 = if (sst > 0) 1 - op$objective / sst else NA_real_, n = length(y), at_bound = at_bound, bound = side)
}

#' Method B: closed-form encroachment gamma
#'
#' theta = a + s sqrt(Q), gamma_B = s^2 (1 + 2A) / (2 (1 + A)^2) (derivation in the file header).
#'
#' @param theta_obs Potential temperature of the window rows [K] (`NA` = not used).
#' @param wtheta Kinematic heat flux of the window rows [K m s-1].
#' @param const Output of [pblh_constants()].
#' @return List: `gamma` [K m-1; `NA` if s <= 0 or too few rows], `slope` [K m-1/2], `r2` [-], `q_total` [K m].
gamma_ec_fit_b <- function(theta_obs, wtheta, const = pblh_constants()) {
  A <- const$bg_A
  q_mid <- c(0, cumsum(wtheta * 1800))[seq_along(wtheta)] + wtheta * 900 # Q at the interval midpoint [K m]
  ok <- is.finite(theta_obs) & q_mid > 0
  out <- list(gamma = NA_real_, slope = NA_real_, r2 = NA_real_, q_total = sum(wtheta) * 1800)
  if (sum(ok) < 3) return(out)
  x <- sqrt(q_mid[ok])
  y <- theta_obs[ok]
  if (stats::var(x) == 0) return(out)
  s <- stats::cov(x, y) / stats::var(x)
  res <- y - (mean(y) - s * mean(x)) - s * x
  sst <- sum((y - mean(y))^2)
  out$slope <- s
  out$r2 <- if (sst > 0) 1 - sum(res^2) / sst else NA_real_
  if (s > 0) out$gamma <- s^2 * (1 + 2 * A) / (2 * (1 + A)^2)
  out
}

#' Flags of the gamma_EC day selection that exclude a day
#' @return Character vector.
gamma_ec_excluding_flags <- function() {
  c("polar_day", "no_window", "no_h0", "bg_stopped_gap", "no_positive_h", "gap", "few_rows", "clear_unknown", "cloudy",
    "rise_lt_min", "r2_low", "bound")
}

#' Daily gamma_EC of one site (spec 17 section 1)
#'
#' @param pm Output of [run_pblh_model()] of the `default` case (or its csv: `timestamp_start_utc`, `day_id`,
#'   `is_day`, `regime`, `h_diag_m`, `wtheta`, `ustar`, `L`, `T_K`, `p_Pa`, `T_source`, `flag_no_h0`,
#'   `flag_bg_stopped_gap`, `flag_no_positive_h`, `flag_polar_day`), complete sorted 30-min grid.
#' @param sw_in,sw_in_pot FLUXNET `SW_IN_F` and `SW_IN_POT` on the rows of `pm` [W m-2] (`NA` where missing).
#' @param site_meta One row of [read_sites()] (`lat`, `lon`).
#' @param opts Output of [gamma_ec_options()]; @param n_sub RK4 sub-steps (as the `default` run).
#' @return Tibble per day with daylight: `date` [Date, `day_id`], `window_start_utc`, `window_end_utc` [POSIXct UTC,
#'   start of the first / last fit half-hour; h is integrated from the first positive heat flux, `start_offset_h`
#'   before the fit start], `n` [-, valid theta rows], `n_tsonic` [-, rows with T from T_SONIC, not
#'   used], `gamma_A`, `gamma_B` [K m-1], `r2_A`, `r2_B` [-], `rmse_A` [K], `theta0_A` [K], `theta_rise` [K], `q_total`
#'   [K m], `h0` [m], `clear_index` [-, sum SW_IN / sum SW_IN_POT over the window], `gamma_A_bound` (`low` / `high`:
#'   Method A optimum at that bound, `gamma_A` NA), `flags` (`;`-separated), `used`
#'   [logical], `gamma` [K m-1, = `gamma_A` where `used`, else `NA`: the daily table read by [run_pblh_model()]].
gamma_ec_daily <- function(pm, sw_in, sw_in_pot, site_meta, opts = gamma_ec_options(), n_sub = 30) {
  const <- pblh_constants()
  n <- nrow(pm)
  t <- as_utc_time(pm$timestamp_start_utc)
  theta <- theta_dry(pm$T_K, pm$p_Pa, const)
  tsonic <- pm$T_source %in% "T_SONIC"
  theta[tsonic] <- NA_real_
  day <- as.Date(pm$day_id)
  isd <- pm$is_day %in% TRUE
  days <- sort(unique(day[isd & !is.na(day)]))
  noon <- sun_times(days, site_meta$lat, site_meta$lon)$solar_noon_utc
  fl <- function(k, i) k %in% names(pm) && any(pm[[k]][i] %in% TRUE)
  rows_by_day <- split(seq_len(n), factor(day, levels = days))
  out <- vector("list", length(days))
  for (j in seq_along(days)) {
    rd <- rows_by_day[[j]]
    f <- character()
    if (fl("flag_polar_day", rd)) f <- c(f, "polar_day")
    if (fl("flag_no_h0", rd)) f <- c(f, "no_h0")
    if (fl("flag_bg_stopped_gap", rd)) f <- c(f, "bg_stopped_gap")
    if (fl("flag_no_positive_h", rd)) f <- c(f, "no_positive_h")
    r <- list(date = days[j], window_start_utc = as.POSIXct(NA, tz = "UTC"), window_end_utc = as.POSIXct(NA, tz = "UTC"),
              n = 0L, n_tsonic = 0L, gamma_A = NA_real_, gamma_B = NA_real_, r2_A = NA_real_, r2_B = NA_real_,
              rmse_A = NA_real_, theta0_A = NA_real_, theta_rise = NA_real_, q_total = NA_real_, h0 = NA_real_,
              clear_index = NA_real_, gamma_A_bound = NA_character_)
    conv <- rd[pm$regime[rd] %in% "convective"]
    s <- rd[isd[rd]][1]
    p <- conv[!is.na(pm$wtheta[conv]) & pm$wtheta[conv] > 0][1]
    if (length(conv) == 0 || is.na(p)) {
      out[[j]] <- c(r, flags = paste(unique(c(f, "no_window")), collapse = ";"))
      next
    }
    # h0 as run_pblh_model(): half-hour before sunrise, else the nearest valid diagnostic height within the 2 h before
    cand <- seq.int(s - 1, max(1, s - const$h0_search_halfhours))
    cand <- cand[cand >= 1 & cand < s & !is.na(pm$h_diag_m[cand])]
    r$h0 <- if (length(cand)) pm$h_diag_m[cand[1]] else NA_real_
    # fit window: first wtheta > 0 to the window end, cut at the day's maximum theta and at solar noon + offset
    q <- max(conv)
    lim <- noon[j] + opts$end_offset_h * 3600
    q <- min(q, max(c(p - 1, rd[t[rd] + 900 <= lim])))
    th_day <- theta[rd[isd[rd]]]
    if (any(is.finite(th_day))) q <- min(q, rd[isd[rd]][which.max(th_day)])
    if (q < p) {
      out[[j]] <- c(r, flags = paste(unique(c(f, "few_rows")), collapse = ";"))
      next
    }
    # fit start pf: start_offset_h after the first positive heat flux; h is integrated from p, rows before pf not fitted
    pf <- p + as.integer(round(opts$start_offset_h * 2))
    if (q < pf) {
      out[[j]] <- c(r, flags = paste(unique(c(f, "few_rows")), collapse = ";"))
      next
    }
    w <- p:q
    wf <- pf:q
    th_fit <- theta[w]
    th_fit[w < pf] <- NA_real_
    r$window_start_utc <- t[pf]
    r$window_end_utc <- t[q]
    r$n_tsonic <- sum(tsonic[wf])
    r$n <- sum(is.finite(theta[wf]))
    r$q_total <- sum(pm$wtheta[w]) * 1800
    ok_sw <- is.finite(sw_in[wf]) & is.finite(sw_in_pot[wf]) & sw_in_pot[wf] > 0
    if (any(ok_sw)) r$clear_index <- sum(sw_in[wf][ok_sw]) / sum(sw_in_pot[wf][ok_sw])
    vt <- wf[is.finite(theta[wf])]
    if (length(vt) >= 2) r$theta_rise <- theta[vt[length(vt)]] - theta[vt[1]]
    gap <- anyNA(c(pm$wtheta[w], pm$ustar[w], pm$L[w], pm$T_K[w]))
    if (gap) f <- c(f, "gap")
    if (r$n < opts$min_n) f <- c(f, "few_rows")
    if (!any(ok_sw)) f <- c(f, "clear_unknown") else if (r$clear_index < opts$clear_min) f <- c(f, "cloudy")
    if (is.na(r$theta_rise) || r$theta_rise < opts$rise_min) f <- c(f, "rise_lt_min")
    if (!gap && r$n >= opts$min_n && !is.na(r$h0)) {
      fa <- gamma_ec_fit_a(th_fit, r$h0, pm$wtheta[w], pm$ustar[w], pm$L[w], pm$T_K[w], opts, n_sub, const)
      fb <- gamma_ec_fit_b(th_fit, pm$wtheta[w], const)
      r$gamma_A <- fa$gamma
      r$r2_A <- fa$r2
      r$rmse_A <- fa$rmse
      r$theta0_A <- fa$theta0
      r$gamma_B <- fb$gamma
      r$gamma_A_bound <- fa$bound
      r$r2_B <- fb$r2
      if (fa$at_bound) f <- c(f, "bound")
      if (is.na(fa$r2) || fa$r2 < opts$r2_min) f <- c(f, "r2_low")
    }
    out[[j]] <- c(r, flags = paste(unique(f), collapse = ";"))
  }
  x <- dplyr::bind_rows(lapply(out, function(r) {
    r$window_start_utc <- as.numeric(r$window_start_utc)
    r$window_end_utc <- as.numeric(r$window_end_utc)
    dplyr::as_tibble(r)
  }))
  x$window_start_utc <- as.POSIXct(x$window_start_utc, origin = "1970-01-01", tz = "UTC")
  x$window_end_utc <- as.POSIXct(x$window_end_utc, origin = "1970-01-01", tz = "UTC")
  ex <- gamma_ec_excluding_flags()
  x$used <- vapply(strsplit(x$flags, ";", fixed = TRUE), function(z) !any(z %in% ex), logical(1)) & !is.na(x$gamma_A)
  x$gamma <- ifelse(x$used, x$gamma_A, NA_real_)
  x
}

#' Monthly gamma_EC climatology (spec 17 section 1)
#'
#' Same layout as the ERA5 monthly table of script 11a ([gamma_climatology()]), so [run_pblh_model()] reads either.
#'
#' @param daily Output of [gamma_ec_daily()].
#' @param min_days Smallest number of used days of a month [-]; fewer: `median`, `p25`, `p75` `NA`.
#' @return Tibble per calendar month 1-12: `month`, `median`, `p25`, `p75` [K m-1], `n_days` [-, used days],
#'   `n_fitted` [-, days with a gamma_A], `n_flagged` [-, days with a gamma_A not used].
gamma_ec_climatology <- function(daily, min_days = 5L) {
  m <- as.integer(format(as.Date(daily$date), "%m"))
  dplyr::bind_rows(lapply(1:12, function(k) {
    g <- daily$gamma_A[m == k & daily$used %in% TRUE]
    q <- if (length(g) >= min_days) stats::quantile(g, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
    nf <- sum(m == k & !is.na(daily$gamma_A))
    dplyr::tibble(month = k, median = q[2], p25 = q[1], p75 = q[3], n_days = length(g), n_fitted = nf,
                  n_flagged = nf - length(g))
  }))
}

#' Flag counts of the gamma_EC days
#'
#' @param daily Output of [gamma_ec_daily()].
#' @return Tibble: `flag`, `n_days` (a day can carry several flags), plus rows `days` and `used`.
gamma_ec_flag_counts <- function(daily) {
  fl <- unlist(strsplit(daily$flags[nzchar(daily$flags)], ";", fixed = TRUE))
  lv <- c(gamma_ec_excluding_flags())
  dplyr::bind_rows(dplyr::tibble(flag = c("days", "used"), n_days = c(nrow(daily), sum(daily$used))),
                   dplyr::tibble(flag = lv, n_days = vapply(lv, function(k) sum(fl == k), integer(1), USE.NAMES = FALSE)))
}

#' Daytime energy-balance ratio (spec 17 section 3 D; diagnostic only, no correction)
#'
#' EBR = sum(H + LE) / sum(NETRAD - G) over daytime (SW_IN_POT > 0) half-hours of the selected months with H and LE
#' measured (`H_F_MDS_QC == 0`, `LE_F_MDS_QC == 0`) and NETRAD and G present; also without G where G is missing.
#'
#' @param d FLUXNET table ([read_icos_fluxnet()], prefixed `fn_`) restricted to the period of interest.
#' @param months Calendar months [integer] (default June-August).
#' @return One-row tibble: `n_rows` [-], `ebr` [-, with G], `n_rows_no_g` [-], `ebr_no_g` [-, (H + LE) / NETRAD],
#'   `g_available` [logical], `months`.
energy_balance_ratio <- function(d, months = 6:8) {
  col <- function(k) if (k %in% names(d)) d[[k]] else rep(NA_real_, nrow(d))
  mo <- as.integer(format(d$timestamp_start_utc, "%m", tz = "UTC"))
  base <- mo %in% months & col("fn_SW_IN_POT") > 0 & col("fn_H_F_MDS_QC") %in% 0 & col("fn_LE_F_MDS_QC") %in% 0 &
    is.finite(col("fn_H_F_MDS")) & is.finite(col("fn_LE_F_MDS")) & is.finite(col("fn_NETRAD"))
  base[is.na(base)] <- FALSE
  wg <- base & is.finite(col("fn_G_F_MDS"))
  hle <- col("fn_H_F_MDS") + col("fn_LE_F_MDS")
  dplyr::tibble(n_rows = sum(wg), ebr = if (any(wg)) sum(hle[wg]) / sum(col("fn_NETRAD")[wg] - col("fn_G_F_MDS")[wg]) else NA_real_,
                n_rows_no_g = sum(base), ebr_no_g = if (any(base)) sum(hle[base]) / sum(col("fn_NETRAD")[base]) else NA_real_,
                g_available = any(wg), months = paste(range(months), collapse = "-"))
}

# --- Spec 17 B: gamma benchmark (what the measurement asks for) ----------------------------------------------------

#' Zero crossing of a curve by linear interpolation
#'
#' @param x Abscissa, increasing [any unit]; @param y Ordinate (e.g. median ln ratio) [-].
#' @return x where y first changes sign (linear interpolation; an exact zero at a node is returned as is); `NA` if y
#'   does not cross 0 within the range.
zero_crossing <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) == 0) return(NA_real_)
  o <- order(x)
  x <- x[o]
  y <- y[o]
  if (any(y == 0)) return(x[which(y == 0)[1]])
  k <- which(sign(y[-1]) != sign(y[-length(y)]))
  if (length(k) == 0) return(NA_real_)
  k <- k[1]
  x[k] + (0 - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k])
}

#' gamma benchmark per site (spec 17 B)
#'
#' gamma_best = the constant gamma where the median ln(PBLHmod / measured) on the site's primary set crosses 0
#' (linear interpolation over the constant grid); interval = the crossings of the day-block bootstrap bounds
#' (`mlr_lo`, `mlr_hi`), sorted.
#'
#' @param wide Metrics of script 11 (one row per site x case x reference x stratum set: `site_id`, `case`,
#'   `reference`, `stratum_set`, `median_log_ratio`, `mlr_lo`, `mlr_hi`, `n_days`).
#' @param primary_set Named character: primary stratum set per site (names = site ids).
#' @param grid Named numeric: gamma of each grid case [K m-1] (names = case labels).
#' @return List: `curve` (site x grid gamma: `site_id`, `stratum_set`, `gamma` [K m-1], `median_log_ratio`, `mlr_lo`,
#'   `mlr_hi` [-], `n_days`) and `best` (one row per site: `site_id`, `stratum_set`, `gamma_best`, `gamma_best_lo`,
#'   `gamma_best_hi` [K m-1], `n_days`).
gamma_benchmark_table <- function(wide, primary_set, grid) {
  cur <- wide[wide$case %in% names(grid) & wide$reference != "era5_blh_center", ]
  cur <- cur[cur$stratum_set == unname(primary_set[cur$site_id]) & !is.na(primary_set[cur$site_id]), ]
  cur <- dplyr::tibble(site_id = cur$site_id, stratum_set = cur$stratum_set, gamma = unname(grid[cur$case]),
                       median_log_ratio = cur$median_log_ratio, mlr_lo = cur$mlr_lo, mlr_hi = cur$mlr_hi,
                       n_days = cur$n_days)
  cur <- cur[order(cur$site_id, cur$gamma), ]
  best <- dplyr::bind_rows(lapply(split(cur, cur$site_id), function(d) {
    b <- sort(c(zero_crossing(d$gamma, d$mlr_lo), zero_crossing(d$gamma, d$mlr_hi)), na.last = TRUE)
    dplyr::tibble(site_id = d$site_id[1], stratum_set = d$stratum_set[1], gamma_best = zero_crossing(d$gamma, d$median_log_ratio),
                  gamma_best_lo = if (length(b) == 2 && all(is.finite(b))) b[1] else NA_real_,
                  gamma_best_hi = if (length(b) == 2 && all(is.finite(b))) b[2] else NA_real_,
                  n_days = suppressWarnings(max(d$n_days, na.rm = TRUE)))
  }))
  list(curve = cur, best = best)
}

#' Figure 17_gamma_benchmark
#'
#' @param tab One row per site: `site_id`, `gamma_best`, `gamma_best_lo`, `gamma_best_hi`, `gamma_ec_median`,
#'   `gamma_ec_p25`, `gamma_ec_p75`, `gamma_site` [K m-1].
#' @param gamma_med Network gamma [K m-1] (`NA` = not drawn); @param gamma_default Constant of PBLHmod [K m-1].
#' @param caption Caption.
#' @return ggplot.
plot_gamma_benchmark <- function(tab, gamma_med, gamma_default = 0.01, caption = NULL) {
  tab$site <- factor(tab$site_id, levels = rev(sort(unique(tab$site_id))))
  ids <- c("gamma_best", "pblhmod_gamma_ec", "pblhmod_gamma_site")
  long <- dplyr::bind_rows(
    dplyr::tibble(site = tab$site, id = "gamma_best", x = tab$gamma_best, lo = tab$gamma_best_lo, hi = tab$gamma_best_hi, dy = 0.2),
    dplyr::tibble(site = tab$site, id = "pblhmod_gamma_ec", x = tab$gamma_ec_median, lo = tab$gamma_ec_p25, hi = tab$gamma_ec_p75, dy = 0),
    dplyr::tibble(site = tab$site, id = "pblhmod_gamma_site", x = tab$gamma_site, lo = NA_real_, hi = NA_real_, dy = -0.2))
  long$y <- as.numeric(long$site) + long$dy
  cols <- c(gamma_best = unname(style_of("gamma_best", "colour", "element")),
            pblhmod_gamma_ec = unname(style_of("pblhmod_gamma_ec")), pblhmod_gamma_site = unname(style_of("pblhmod_gamma_site")))
  labs_v <- c(gamma_best = unname(style_of("gamma_best", "label", "element")),
              pblhmod_gamma_ec = "γEC site median (p25-p75 of the used days)",
              pblhmod_gamma_site = "ERA5 γ site median (gamma_site)")
  vl <- dplyr::tibble(x = c(gamma_default, gamma_med), id = c("pblhmod_gamma001", "pblhmod_gamma_med"))
  vl <- vl[is.finite(vl$x), ]
  ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, colour = id)) +
    ggplot2::geom_vline(data = vl, ggplot2::aes(xintercept = x), colour = style_of(vl$id),
                        linetype = style_of(c("reference_line"), "linetype", "element")) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y", width = 0, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = cols, labels = labs_v, breaks = ids, name = NULL) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(tab$site)), labels = levels(tab$site)) +
    ggplot2::labs(x = "γ [K m-1]", y = NULL, caption = caption) +
    theme_pblh() + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
}

#' Write a csv only when its content changes
#'
#' Keeps the file (and its modification time) when the new content is identical, so that the inputs hash of a
#' dependent step does not change (spec 17: step 02 hashes the gamma_EC tables of step 17a, which reads the default
#' output of step 02).
#'
#' @param x Table; @param path Output path [character].
#' @return `TRUE` if written, `FALSE` if unchanged (invisibly).
write_csv_if_changed <- function(x, path) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_csv(x, tmp)
  if (file.exists(path) && unname(tools::md5sum(tmp)) == unname(tools::md5sum(path))) return(invisible(FALSE))
  if (!file.copy(tmp, path, overwrite = TRUE)) cli::cli_abort("Could not write {.path {path}}.")
  invisible(TRUE)
}

# --- Spec 17 C, D, E: network collection (scripts/17b_gamma_ec_collect.R) ------------------------------------------

#' Monthly gamma_EC next to the ERA5 climatology of one site (spec 17 C)
#'
#' @param ec Monthly gamma_EC table ([gamma_ec_climatology()], with or without `site_id`).
#' @param era5 Monthly ERA5 table of script 11a ([gamma_climatology()]).
#' @param site_id Station code.
#' @return Tibble per month: `site_id`, `month`, `gamma_ec`, `gamma_ec_p25`, `gamma_ec_p75` [K m-1], `n_days_ec` [-],
#'   `gamma_era5` [K m-1], `n_days_era5` [-], `ratio_ec_era5` [-].
gamma_ec_vs_era5_table <- function(ec, era5, site_id) {
  m <- 1:12
  g <- function(x, k) if (is.null(x) || !k %in% names(x)) rep(NA_real_, 12) else x[[k]][match(m, x$month)]
  out <- dplyr::tibble(site_id = site_id, month = m, gamma_ec = g(ec, "median"), gamma_ec_p25 = g(ec, "p25"),
                       gamma_ec_p75 = g(ec, "p75"), n_days_ec = g(ec, "n_days"), gamma_era5 = g(era5, "median"),
                       n_days_era5 = g(era5, "n_days"))
  out$ratio_ec_era5 <- out$gamma_ec / out$gamma_era5
  out
}

#' Per-site summary of gamma_EC (spec 17 C, D)
#'
#' @param mon Output of [gamma_ec_vs_era5_table()] (one site); @param daily [gamma_ec_daily()] of the site.
#' @param era5_daily Daily ERA5 gamma table of the site (11a) or `NULL`; @param ebr [energy_balance_ratio()] row.
#' @return One-row tibble: `site_id`, `n_used_days`, `months_with_value` [-], `share_months` [0-1], `gamma_ec_site`
#'   [K m-1, median of the used days], `gamma_era5_site` [K m-1, [gamma_site_constant()]], `ratio_ec_era5` [-, median of
#'   the monthly ratios], `ratio_B_A` [-, median gamma_B / gamma_A on the used days], `ebr`, `ebr_no_g` [-],
#'   `ebr_n_rows` [-], `g_available` [logical].
gamma_ec_site_summary <- function(mon, daily, era5_daily = NULL, ebr = NULL) {
  u <- daily[daily$used %in% TRUE, ]
  e <- function(k) if (is.null(ebr) || !k %in% names(ebr)) NA else ebr[[k]][1]
  dplyr::tibble(site_id = mon$site_id[1], n_used_days = nrow(u), months_with_value = sum(!is.na(mon$gamma_ec)),
                share_months = mean(!is.na(mon$gamma_ec)), gamma_ec_site = med_fin(u$gamma_A),
                gamma_era5_site = if (is.null(era5_daily)) NA_real_ else gamma_site_constant(era5_daily)$gamma,
                ratio_ec_era5 = med_fin(mon$ratio_ec_era5), ratio_B_A = med_fin(u$gamma_B / u$gamma_A),
                ebr = as.numeric(e("ebr")), ebr_no_g = as.numeric(e("ebr_no_g")), ebr_n_rows = as.numeric(e("n_rows")),
                g_available = as.logical(e("g_available")))
}

#' Footprint bracketing check (spec 17 E, for decision D1)
#'
#' On the footprint rows (spec-15 `common_pe` of the reference year), the share of half-hours where h of PBLHmod-gammaEC
#' lies outside [min, max] of (ERA5, PBLHmod, PBLHmod-clim), and the median distance outside in ln units.
#'
#' @param r Rows of [fp_site_rows()] (one site) with `h__pblhmod`, `h__pblhmod_clim`, `h__era5` [m].
#' @param h_ec h of PBLHmod-gammaEC on the rows of `r` [m].
#' @param set,subset Validity set and subset of the row basis (default `etc`, `valid_qc_used`).
#' @return Tibble per regime x day/night (with margins `all`): `site_id`, `regime`, `daynight`, `n` (rows with every h),
#'   `share_outside`, `share_below`, `share_above` [0-1], `median_dist_ln` [-, median ln distance of the outside rows:
#'   ln(min / h) below, ln(h / max) above], `p90_dist_ln` [-].
gamma_ec_bracketing <- function(r, h_ec, set = "etc", subset = "valid_qc_used") {
  cm <- fp_common(r, set, subset)
  hs <- cbind(r$h__era5, r$h__pblhmod, r$h__pblhmod_clim)
  ok <- cm & is.finite(h_ec) & h_ec > 0 & stats::complete.cases(hs) & apply(hs > 0, 1, all)
  lo <- apply(hs[ok, , drop = FALSE], 1, min)
  hi <- apply(hs[ok, , drop = FALSE], 1, max)
  h <- h_ec[ok]
  below <- h < lo
  above <- h > hi
  dist <- ifelse(below, log(lo / h), ifelse(above, log(h / hi), 0))
  long <- fp_strata_long(r$site_id[ok], r$zl_class[ok], r$daynight[ok], below = below, above = above, dist = dist)
  dplyr::ungroup(dplyr::summarise(dplyr::group_by(long, site_id, regime, daynight), n = dplyr::n(),
                                  share_outside = mean(below | above), share_below = mean(below), share_above = mean(above),
                                  median_dist_ln = med_fin(dist[below | above]), p90_dist_ln = q_fin(dist[below | above], 0.9),
                                  .groups = "drop"))
}

#' Figure 17_gamma_ec_vs_era5: monthly gamma_EC against the ERA5 climatology (all sites)
#'
#' @param mon Rows of [gamma_ec_vs_era5_table()] (several sites); @param caption Caption.
#' @return ggplot.
plot_gamma_ec_vs_era5 <- function(mon, caption = NULL) {
  d <- mon[is.finite(mon$gamma_ec) & is.finite(mon$gamma_era5), ]
  ggplot2::ggplot(d, ggplot2::aes(gamma_era5, gamma_ec)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = style_of("mc_one_to_one", "colour", "element"),
                         linetype = style_of("mc_one_to_one", "linetype", "element")) +
    ggplot2::geom_hline(yintercept = 0.01, colour = style_of("pblhmod_gamma001"),
                        linetype = style_of("reference_line", "linetype", "element")) +
    ggplot2::geom_point(colour = style_of("pblhmod_gamma_ec"), size = 1, alpha = 0.6) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "ERA5 γ, monthly median [K m-1]", y = "γEC, monthly median [K m-1]", caption = caption) +
    theme_pblh()
}

#' Figure 17_gamma_ec_vs_era5_seasonal: seasonal cycle of gamma_EC and ERA5 gamma per site
#'
#' @inheritParams plot_gamma_ec_vs_era5
#' @return ggplot.
plot_gamma_ec_seasonal <- function(mon, caption = NULL) {
  ids <- c("pblhmod_gamma_ec", "pblhmod_clim")
  long <- dplyr::bind_rows(dplyr::tibble(site_id = mon$site_id, month = mon$month, id = ids[1], g = mon$gamma_ec,
                                         lo = mon$gamma_ec_p25, hi = mon$gamma_ec_p75),
                           dplyr::tibble(site_id = mon$site_id, month = mon$month, id = ids[2], g = mon$gamma_era5,
                                         lo = NA_real_, hi = NA_real_))
  long <- long[is.finite(long$g), ]
  ggplot2::ggplot(long, ggplot2::aes(month, g, colour = id)) +
    ggplot2::geom_hline(yintercept = 0.01, colour = style_of("pblhmod_gamma001"),
                        linetype = style_of("reference_line", "linetype", "element")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3, na.rm = TRUE) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_colour_manual(values = style_of(ids), labels = c(pblhmod_gamma_ec = "γEC (median, p25-p75)",
                                                                    pblhmod_clim = "ERA5 γ climatology"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(1, 4, 7, 10)) +
    ggplot2::facet_wrap(~site_id) +
    ggplot2::labs(x = "Month", y = "γ [K m-1]", caption = caption) +
    theme_pblh()
}

#' Folder name of the gamma_EC tables for a fit-start offset
#'
#' Offset 0 = the standard tables read by [run_pblh_model()] (`gamma_ec`); other offsets (spec 17 supplement, G
#' 2026-11-26) in their own folder, never read by step 02 (`gamma_ec_off2h` for 2 h).
#'
#' @param offset_h Fit-start offset [h] ([gamma_ec_options()] `start_offset_h`).
#' @return Folder name [character].
gamma_ec_dir_name <- function(offset_h = 0) {
  if (isTRUE(all.equal(offset_h, 0))) "gamma_ec" else paste0("gamma_ec_off", format(offset_h), "h")
}

#' Heights of a site on the used gamma_EC days (spec 17 supplement)
#'
#' EC height and canopy height from [get_site_heights()] at the start of the fit window of every used day; zm =
#' EC height - 2/3 canopy height (footprint convention, d = 2/3 hc).
#'
#' @param daily Daily table of [gamma_ec_daily()]; @param site_id Station code.
#' @param sites,sites_tv Outputs of [read_sites()] and [read_sites_timevarying()].
#' @return One-row tibble: `ec_height_m`, `canopy_height_m`, `zm_m` [m, medians over the used days; `NA` without used
#'   days], `igbp` [-].
gamma_ec_heights <- function(daily, site_id, sites, sites_tv) {
  igbp <- if ("igbp" %in% names(sites)) sites$igbp[sites$site_id == site_id][1] else NA_character_
  t <- as_utc_time(daily$window_start_utc[daily$used %in% TRUE])
  t <- t[!is.na(t)]
  if (!length(t)) {
    return(dplyr::tibble(ec_height_m = NA_real_, canopy_height_m = NA_real_, zm_m = NA_real_, igbp = igbp))
  }
  h <- get_site_heights(site_id, t, sites, sites_tv)
  dplyr::tibble(ec_height_m = med_fin(h$zm_m), canopy_height_m = med_fin(h$canopy_height_m),
                zm_m = med_fin(h$zm_m - 2 / 3 * h$canopy_height_m), igbp = igbp)
}

#' Site gamma_EC / gamma_ERA5 against the measurement height (spec 17 supplement figure)
#'
#' @param x Table `17_gamma_ratio_vs_height` (`site_id`, `offset_h`, `ratio_ec_era5`, `zm_m`, `ec_height_m`, `igbp`).
#' @param caption Caption.
#' @return ggplot: two panels (x = zm, x = EC height, log axes), one point per site and offset, labels site (IGBP).
plot_gamma_ratio_vs_height <- function(x, caption = NULL) {
  x <- x[is.finite(x$ratio_ec_era5), ]
  lab_z <- "zm = EC height - 2/3 hc [m]"
  lab_e <- "EC height [m]"
  long <- dplyr::bind_rows(
    dplyr::mutate(x, height = zm_m, panel = lab_z),
    dplyr::mutate(x, height = ec_height_m, panel = lab_e))
  long <- long[is.finite(long$height) & long$height > 0, ]
  long$panel <- factor(long$panel, levels = c(lab_z, lab_e))
  long$offset <- factor(paste0("fit start +", format(long$offset_h), " h"))
  long$label <- ifelse(long$offset_h == min(long$offset_h), paste0(long$site_id, " (", long$igbp, ")"), NA_character_)
  ggplot2::ggplot(long, ggplot2::aes(height, ratio_ec_era5, shape = offset)) +
    ggplot2::geom_hline(yintercept = 1, colour = style_of("reference_line", "colour", "element")) +
    ggplot2::geom_point(colour = style_of("pblhmod_gamma_ec"), size = 1.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2, vjust = -0.9, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = c(16, 1)[seq_along(levels(long$offset))], name = NULL) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~panel, scales = "free_x", strip.position = "bottom") +
    ggplot2::labs(x = NULL, y = "γEC / γERA5 [-]", caption = caption) +
    theme_pblh() +
    ggplot2::theme(strip.placement = "outside") +
    ggplot2::coord_cartesian(clip = "off")
}
