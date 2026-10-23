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
