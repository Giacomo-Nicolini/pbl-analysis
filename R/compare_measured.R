# Comparison of model / reanalysis / operational PBLH with a measured, censored MLH (spec 07 section 5).
# The measured MLH (`cloudnet_eps_mlh`) has three kinds of half-hours (G 2026-10-23):
#   qc 0 / 2: exact value; qc 1: censored below z_min, i.e. in (0, z_min]; qc 3: signal-limited, a lower bound,
#   i.e. in [MLH, Inf). Metrics are computed on qc 0/2; qc 1 and qc 3 are summarised separately.

#' Interval representation of the measured MLH
#'
#' @param pblh_m Measured MLH [m above ground] (`NA` for qc 1).
#' @param qc_flag Measured `qc_flag` (0, 1, 2, 3 or `NA`).
#' @param z_min Censoring height [m above ground] (`x_censor_below_m`).
#' @return Tibble: `left`, `right` [m; `right = NA` = no upper limit],
#'   `kind` (`exact`, `below`, `lower_bound`; `NA` where not usable). Coding as
#'   `survival::Surv(left, right, type = "interval2")`: exact `left = right`,
#'   below z_min `(0, z_min]`, lower bound `[MLH, Inf)`.
measured_interval <- function(pblh_m, qc_flag, z_min) {
  n <- length(pblh_m)
  kind <- rep(NA_character_, n)
  kind[qc_flag %in% c(0L, 2L) & !is.na(pblh_m)] <- "exact"
  kind[qc_flag %in% 1L & !is.na(z_min)] <- "below"
  kind[qc_flag %in% 3L & !is.na(pblh_m)] <- "lower_bound"
  left <- ifelse(kind == "exact" | kind == "lower_bound", pblh_m, ifelse(kind == "below", 0, NA_real_))
  right <- ifelse(kind == "exact", pblh_m, ifelse(kind == "below", z_min, NA_real_))
  dplyr::tibble(left = left, right = right, kind = kind)
}

#' Median of interval-censored values (Turnbull NPMLE)
#'
#' Median of the survival curve fitted with
#' `survival::survfit(Surv(left, right, type = "interval2") ~ 1)` (Turnbull
#' estimator for interval-censored data). Lower bounds `[MLH, Inf)` are passed
#' with a finite upper limit `upper_cap`, a workaround checked on 2026-10-23 with
#' survival 3.8-12: `Surv(l, NA, type = "interval2")` is coded correctly
#' (status 0 = right-censored, as `Surv(time, event)`), but on data with only exact
#' and right-censored values `survfit()` takes its interval-censored path and
#' returns S = 0.787, 0.574, 0.235 (n.risk 6 for 5 values) instead of the
#' Kaplan-Meier 0.8, 0.6, 0.3 (values 100, 200, 300+, 400, 900+); explicit
#' `type = "interval"` coding gives the same. With a finite upper limit the fit
#' equals Kaplan-Meier, and with interval-censored values present both codings
#' agree (tests in `test-compare_measured.R`). The cap only places the unknown
#' mass above all data. `NA` if the median falls above the largest finite limit
#' (more than half the mass in lower bounds above it).
#'
#' @param left,right Interval limits as [measured_interval()] [m].
#' @param upper_cap Upper limit used for lower bounds [m].
#' @return Median [m]; `NA` if not estimable or fewer than 2 values.
censored_median <- function(left, right, upper_cap = 1e7) {
  ok <- !(is.na(left) & is.na(right))
  if (sum(ok) < 2) return(NA_real_)
  l <- left[ok]
  r <- right[ok]
  top <- max(c(l, r), na.rm = TRUE)
  r[is.na(r) | !is.finite(r)] <- upper_cap
  fit <- survival::survfit(survival::Surv(l, r, type = "interval2") ~ 1)
  m <- unname(summary(fit)$table["median"])
  if (is.null(m) || !is.finite(m) || m > top) NA_real_ else as.numeric(m)
}

#' Diurnal cycle of the measured MLH with censored medians
#'
#' @param d Tibble with `left`, `right`, `kind` ([measured_interval()]), `pblh_m`
#'   and the grouping columns.
#' @param by Grouping columns, e.g. `c("season", "hour_utc")`.
#' @return Tibble per group: `median_censored` [m, Turnbull median of all usable
#'   half-hours], `median_exact` [m, plain median of qc 0/2], `n`, `n_exact`,
#'   `n_below`, `n_lower_bound` [-].
measured_diurnal <- function(d, by = c("season", "hour_utc")) {
  d <- d[!is.na(d$kind), ]
  key <- interaction(d[by], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(d)), key), function(i) {
    x <- d[i, ]
    dplyr::bind_cols(x[1, by, drop = FALSE], dplyr::tibble(
      median_censored = censored_median(x$left, x$right),
      median_exact = if (any(x$kind == "exact")) stats::median(x$pblh_m[x$kind == "exact"]) else NA_real_,
      n = nrow(x), n_exact = sum(x$kind == "exact"), n_below = sum(x$kind == "below"),
      n_lower_bound = sum(x$kind == "lower_bound")
    ))
  })
  out <- dplyr::bind_rows(out)
  out[do.call(order, unname(as.list(out[by]))), ]
}

#' Measured MLH on ERA5 hours
#'
#' Mean of the two half-hours adjacent to the hour (spec 04 rule, both required),
#' keeping the censoring: exact if both are exact (mean); lower bound if both are
#' exact or lower bounds and at least one is a lower bound (mean of the two, a
#' lower bound of the hourly mean); below z_min if both are below; otherwise
#' (one below, one not) not usable.
#'
#' @param meas Standard table `cloudnet_eps_mlh` (`timestamp_start_utc`,
#'   `timestamp_end_utc`, `pblh_m`, `qc_flag`, `x_censor_below_m`).
#' @param hours_utc ERA5 times [POSIXct UTC].
#' @return Tibble: `timestamp_utc`, `pblh_m` [m; exact or lower bound], `qc_flag`
#'   (0 exact, 1 below, 3 lower bound, `NA`), `x_censor_below_m` [m].
match_measured_to_hours <- function(meas, hours_utc) {
  hours_utc <- sort(unique(hours_utc))
  i1 <- match(as.numeric(hours_utc), as.numeric(meas$timestamp_end_utc))
  i2 <- match(as.numeric(hours_utc), as.numeric(meas$timestamp_start_utc))
  k1 <- measured_interval(meas$pblh_m[i1], meas$qc_flag[i1], meas$x_censor_below_m[i1])$kind
  k2 <- measured_interval(meas$pblh_m[i2], meas$qc_flag[i2], meas$x_censor_below_m[i2])$kind
  v <- (meas$pblh_m[i1] + meas$pblh_m[i2]) / 2
  up <- c("exact", "lower_bound")
  qc <- ifelse(k1 %in% "exact" & k2 %in% "exact", 0L,
               ifelse(k1 %in% up & k2 %in% up, 3L,
                      ifelse(k1 %in% "below" & k2 %in% "below", 1L, NA_integer_)))
  dplyr::tibble(timestamp_utc = hours_utc, pblh_m = ifelse(qc %in% c(0L, 3L), v, NA_real_), qc_flag = qc,
                x_censor_below_m = pmax(meas$x_censor_below_m[i1], meas$x_censor_below_m[i2]))
}

#' Model values against a lower bound
#'
#' For half-hours (or hours) where the measured MLH is only a lower bound (qc 3).
#'
#' @param df Tibble with `mod` [m], `bound` [m] and the grouping columns.
#' @param by Grouping columns (`NULL` = all).
#' @return Tibble per group: `n`, `frac_mod_below_bound` [0-1], `median_ratio`
#'   [-, median of mod / bound].
lower_bound_summary <- function(df, by = NULL) {
  df <- df[!is.na(df$mod) & !is.na(df$bound) & df$bound > 0, ]
  f <- function(x) dplyr::tibble(n = nrow(x), frac_mod_below_bound = mean(x$mod < x$bound),
                                 median_ratio = stats::median(x$mod / x$bound))
  if (length(by) == 0 || nrow(df) == 0) return(f(df))
  dplyr::bind_rows(lapply(split(df, interaction(df[by], drop = TRUE, lex.order = TRUE)), function(x) {
    dplyr::bind_cols(x[1, by, drop = FALSE], f(x))
  }))
}

#' Model values where the measured MLH is censored below z_min
#'
#' @param df Tibble with `mod` [m], `z_min` [m] and the grouping columns.
#' @param by Grouping columns (`NULL` = all).
#' @return Tibble per group: `n`, `frac_mod_below_zmin` [0-1], `median_mod` [m].
censored_below_summary <- function(df, by = NULL) {
  df <- df[!is.na(df$mod) & !is.na(df$z_min), ]
  f <- function(x) dplyr::tibble(n = nrow(x), frac_mod_below_zmin = mean(x$mod <= x$z_min),
                                 median_mod = stats::median(x$mod))
  if (length(by) == 0 || nrow(df) == 0) return(f(df))
  dplyr::bind_rows(lapply(split(df, interaction(df[by], drop = TRUE, lex.order = TRUE)), function(x) {
    dplyr::bind_cols(x[1, by, drop = FALSE], f(x))
  }))
}

# --- Spec 10 A5: any measured source (censored half-hourly, or hourly interval means) ---

#' Half-hourly standard series on hourly intervals [t, t + 1 h)
#'
#' Spec 10 A5: mean of the two half-hours inside the hour (starting at `t` and `t + 30 min`), both required; regime
#' = the common regime of the two, `"mixed"` if they differ (a convective hour has both half-hours convective);
#' larger `qc_flag`. The rule of [match_halfhours_to_hours()] applied at `t + 30 min`.
#'
#' @param std_hh Half-hourly standard table (`timestamp_start_utc`, `timestamp_end_utc`, `pblh_m` [m]; optional
#'   `qc_flag`, `x_regime`).
#' @param starts_utc Starts of the hourly intervals [POSIXct UTC].
#' @return Tibble: `timestamp_utc` (= interval start), `pblh_m` [m], `n_samples` [0-2], `qc_flag`, `x_regime`.
match_halfhours_to_intervals <- function(std_hh, starts_utc) {
  out <- match_halfhours_to_hours(std_hh, sort(unique(starts_utc)) + 1800, require_both = TRUE)
  out$timestamp_utc <- out$timestamp_utc - 1800
  out
}

#' Instantaneous standard series on hourly intervals [t, t + 1 h)
#'
#' Spec 10 A5: mean of the values at `t` and `t + 1 h`, both required: an approximation of the hourly mean of an
#' instantaneous series such as ERA5 BLH (exact if the series is linear within the hour).
#'
#' @param std_inst Instantaneous standard table (`timestamp_start_utc`, `pblh_m` [m]).
#' @param starts_utc Starts of the hourly intervals [POSIXct UTC].
#' @return Tibble: `timestamp_utc` (= interval start), `pblh_m` [m], `n_samples` [0-2].
match_instants_to_intervals <- function(std_inst, starts_utc) {
  starts_utc <- sort(unique(starts_utc))
  t_in <- as.numeric(std_inst$timestamp_start_utc)
  v1 <- std_inst$pblh_m[match(as.numeric(starts_utc), t_in)]
  v2 <- std_inst$pblh_m[match(as.numeric(starts_utc) + 3600, t_in)]
  n <- (!is.na(v1)) + (!is.na(v2))
  dplyr::tibble(timestamp_utc = starts_utc, pblh_m = ifelse(n == 2, (v1 + v2) / 2, NA_real_),
                n_samples = as.integer(n))
}

#' Measured standard source as usable values
#'
#' Unifies the qc semantics (spec 10 A5): a source with a censoring height (`x_censor_below_m`, e.g.
#' `cloudnet_eps_mlh`) is coded by [measured_interval()] (exact, below z_min, lower bound); a source without
#' censoring (e.g. `sirta_mlh_*`) is exact wherever `pblh_m` is present (the writer already set `pblh_m` to `NA`
#' where its qc rejects the value).
#'
#' @param meas Measured standard table.
#' @return Tibble: `timestamp_utc` (= `timestamp_start_utc`), `pblh_m` [m], `qc_flag`, `z_min` [m], `left`,
#'   `right` [m], `kind`.
measured_values <- function(meas) {
  if ("x_censor_below_m" %in% names(meas)) {
    z <- meas$x_censor_below_m
    iv <- measured_interval(meas$pblh_m, meas$qc_flag, z)
  } else {
    z <- rep(NA_real_, nrow(meas))
    ex <- !is.na(meas$pblh_m)
    iv <- dplyr::tibble(left = ifelse(ex, meas$pblh_m, NA_real_), right = ifelse(ex, meas$pblh_m, NA_real_),
                        kind = ifelse(ex, "exact", NA_character_))
  }
  dplyr::bind_cols(dplyr::tibble(timestamp_utc = meas$timestamp_start_utc, pblh_m = meas$pblh_m,
                                 qc_flag = meas$qc_flag, z_min = z), iv[c("left", "right", "kind")])
}

#' Model and reanalysis series on the time base of a measured source
#'
#' Generalises the matching of spec 07 (spec 10 A5). Half-hourly measured source (spec 07): PBLHmod on the
#' half-hours; ERA5 on its hours, against the measured mean of the two adjacent half-hours
#' ([match_measured_to_hours()]), with the PBLHmod regime of the hour ([match_halfhours_to_hours()]). Hourly
#' interval-mean source [t, t + 1 h) (spec 10): PBLHmod by [match_halfhours_to_intervals()], ERA5 by
#' [match_instants_to_intervals()], both on the measured intervals; ERA5 carries the PBLHmod regime of the
#' interval.
#'
#' @param meas Measured standard table (half-hourly or hourly interval means).
#' @param mod PBLHmod standard table (half-hourly, with `x_regime`).
#' @param era Named list of ERA5 standard tables (instantaneous, hourly, same times).
#' @return List: `base` (`"halfhour"` or `"hour"`); `meas`, `mod`, `mid_s` (measured values, PBLHmod on the
#'   measurement base, offset [s] from the time stamp to the interval midpoint); `era`, `era_meas`, `era_mid_s`
#'   (ERA5 series, measured values on the ERA5 base, offset [s]); `mod_on_era` (PBLHmod on the ERA5 base). Times in
#'   `timestamp_utc`.
measured_time_base <- function(meas, mod, era) {
  step <- stats::median(as.numeric(meas$timestamp_end_utc) - as.numeric(meas$timestamp_start_utc))
  with_regime <- function(e, m) {
    e$x_regime <- m$x_regime[match(as.numeric(e$timestamp_utc), as.numeric(m$timestamp_utc))]
    e
  }
  if (step == 1800) {
    hours <- era[[1]]$timestamp_start_utc
    mod_h <- match_halfhours_to_hours(mod, hours)
    mhr <- match_measured_to_hours(meas, hours)
    mhr <- dplyr::bind_cols(mhr, measured_interval(mhr$pblh_m, mhr$qc_flag, mhr$x_censor_below_m)[c("left", "right", "kind")])
    mhr$z_min <- mhr$x_censor_below_m
    list(base = "halfhour", meas = measured_values(meas), mod = dplyr::rename(mod, timestamp_utc = timestamp_start_utc),
         mid_s = 900, era = lapply(era, function(e) with_regime(dplyr::rename(e, timestamp_utc = timestamp_start_utc), mod_h)),
         era_meas = mhr, era_mid_s = 0, mod_on_era = mod_h)
  } else if (step == 3600) {
    starts <- meas$timestamp_start_utc
    mod_i <- match_halfhours_to_intervals(mod, starts)
    mv <- measured_values(meas)
    list(base = "hour", meas = mv, mod = mod_i, mid_s = 1800,
         era = lapply(era, function(e) with_regime(match_instants_to_intervals(e, starts), mod_i)),
         era_meas = mv, era_mid_s = 1800, mod_on_era = mod_i)
  } else {
    cli::cli_abort("Measured source: time step {step} s not supported (half-hourly or hourly interval means).")
  }
}

#' Keep the regimes allowed by the distance class
#'
#' Distance rule (project plan): `all_regimes` keeps every row; the other classes only the convective regime.
#'
#' @param regime PBLHmod regime per row [character].
#' @param dclass Output of [distance_class()].
#' @return Logical.
keep_regime_by_distance <- function(regime, dclass) {
  if (dclass == "all_regimes") rep(TRUE, length(regime)) else regime %in% "convective"
}
