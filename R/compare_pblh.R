# Comparison of PBLH sources: matching half-hourly series to ERA5 hours,
# metrics, diurnal summaries and the exploratory extras of spec 04
# (EC H vs ERA5 H, cloud-base classes, free-atmosphere theta gradient).
# Rules: docs/specs/04_era5_standard_compare.md sections 6-7, amendment 2.

#' Match a half-hourly standard table to instantaneous hours
#'
#' For each instant `t` (e.g. ERA5 hour) the mean of the half-hour ending at
#' `t` and the half-hour starting at `t` (decision with G, 2026-10-16).
#' `n_samples` = number of valid (non-missing `pblh_m`) half-hours used. With
#' `require_both = TRUE` the hour is `NA` unless both are valid. Regime and
#' `qc_flag` (amendment 2): if both half-hours are valid, the common regime or
#' `"mixed"` if they differ (a missing regime gives the other one), and the
#' larger `qc_flag`; if only one is valid, its values (only with
#' `require_both = FALSE`; otherwise `NA`).
#'
#' @param std_hh Half-hourly standard table (`timestamp_start_utc`,
#'   `timestamp_end_utc`, `pblh_m` [m]; optional `qc_flag`, `x_regime`).
#' @param hours_utc Instants [POSIXct UTC].
#' @param require_both Require both half-hours [logical].
#' @return Tibble: `timestamp_utc`, `pblh_m` [m], `n_samples` [0-2],
#'   `qc_flag` [integer], `x_regime` [character].
match_halfhours_to_hours <- function(std_hh, hours_utc, require_both = TRUE) {
  hours_utc <- sort(unique(hours_utc))
  i1 <- match(as.numeric(hours_utc), as.numeric(std_hh$timestamp_end_utc))   # half-hour ending at t
  i2 <- match(as.numeric(hours_utc), as.numeric(std_hh$timestamp_start_utc)) # half-hour starting at t
  v1 <- std_hh$pblh_m[i1]
  v2 <- std_hh$pblh_m[i2]
  ok1 <- !is.na(v1)
  ok2 <- !is.na(v2)
  n <- ok1 + ok2
  val <- ifelse(ok1 & ok2, (v1 + v2) / 2, ifelse(ok1, v1, ifelse(ok2, v2, NA_real_)))
  use <- if (require_both) n == 2 else n >= 1
  val[!use] <- NA_real_

  pick <- function(col, combine) {
    if (!col %in% names(std_hh)) return(rep(NA, length(hours_utc)))
    a <- std_hh[[col]][i1]
    b <- std_hh[[col]][i2]
    a[!ok1] <- NA
    b[!ok2] <- NA
    out <- combine(a, b)
    out[!use] <- NA
    out
  }
  qc <- pick("qc_flag", function(a, b) {
    ifelse(is.na(a), b, ifelse(is.na(b), a, pmax(a, b)))
  })
  reg <- pick("x_regime", function(a, b) {
    ifelse(is.na(a), b, ifelse(is.na(b), a, ifelse(a == b, a, "mixed")))
  })
  dplyr::tibble(
    timestamp_utc = hours_utc,
    pblh_m = val,
    n_samples = as.integer(n),
    qc_flag = as.integer(qc),
    x_regime = as.character(reg)
  )
}

#' Match half-hourly values to hourly accumulation intervals
#'
#' For an ERA5 accumulation over (`t - 1 h`, `t`]: the mean of the two
#' half-hours **ending** at `t - 30 min` and `t` (used for EC H vs ERA5 H).
#'
#' @param hh Half-hourly table with `timestamp_end_utc` and `value_col`.
#' @param hours_utc Ends of the accumulation intervals [POSIXct UTC].
#' @param value_col Column to average [character].
#' @param require_both Require both half-hours [logical].
#' @return Tibble: `timestamp_utc` (= interval end), `value` [unit of
#'   `value_col`], `n_samples` [0-2].
match_halfhours_to_accum <- function(hh, hours_utc, value_col, require_both = TRUE) {
  hours_utc <- sort(unique(hours_utc))
  i1 <- match(as.numeric(hours_utc) - 1800, as.numeric(hh$timestamp_end_utc))
  i2 <- match(as.numeric(hours_utc), as.numeric(hh$timestamp_end_utc))
  v1 <- hh[[value_col]][i1]
  v2 <- hh[[value_col]][i2]
  n <- (!is.na(v1)) + (!is.na(v2))
  val <- ifelse(n == 2, (v1 + v2) / 2, ifelse(!is.na(v1), v1, v2))
  val[if (require_both) n < 2 else n < 1] <- NA_real_
  dplyr::tibble(timestamp_utc = hours_utc, value = val, n_samples = as.integer(n))
}

#' Meteorological season of a UTC time
#'
#' @param t Times [POSIXct UTC].
#' @return Factor `DJF`, `MAM`, `JJA`, `SON` (from the UTC month).
season_of <- function(t) {
  m <- as.integer(format(t, "%m", tz = "UTC"))
  s <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")[m]
  factor(s, levels = c("DJF", "MAM", "JJA", "SON"))
}

#' Paired table of two hourly PBLH series
#'
#' Inner join on `timestamp_utc` of rows where both `pblh_m` are valid.
#'
#' @param mod,ref Tibbles with `timestamp_utc` and `pblh_m` [m] (e.g. output
#'   of [match_halfhours_to_hours()] or an instantaneous standard table renamed
#'   to `timestamp_utc`); `mod` may carry `x_regime`, both may carry
#'   `qc_flag`; `ref` may carry `x_blh_cv9` (ERA5 3x3 spread).
#' @param pair Pair name [character].
#' @return Tibble: `pair`, `timestamp_utc`, `mod`, `ref` [m], `regime`,
#'   `qc_mod`, `qc_ref`, `cv9_ref` [-], `season`, `hour_utc` [0-23].
pair_table <- function(mod, ref, pair) {
  col <- function(d, k, na) if (k %in% names(d)) d[[k]] else rep(na, nrow(d))
  # columns computed before tibble(): inside it, the new columns `mod`/`ref` would mask the arguments
  regime <- as.character(col(mod, "x_regime", NA_character_))
  qc_mod <- as.integer(col(mod, "qc_flag", NA_integer_))
  qc_ref <- as.integer(col(ref, "qc_flag", NA_integer_))
  cv9_ref <- as.numeric(col(ref, "x_blh_cv9", NA_real_))
  a <- dplyr::tibble(timestamp_utc = mod$timestamp_utc, mod = mod$pblh_m, regime = regime, qc_mod = qc_mod)
  b <- dplyr::tibble(timestamp_utc = ref$timestamp_utc, ref = ref$pblh_m, qc_ref = qc_ref, cv9_ref = cv9_ref)
  d <- dplyr::inner_join(a[!is.na(a$mod), ], b[!is.na(b$ref), ], by = "timestamp_utc")
  d <- d[order(d$timestamp_utc), ]
  dplyr::bind_cols(
    dplyr::tibble(pair = rep(pair, nrow(d))),
    d,
    dplyr::tibble(season = season_of(d$timestamp_utc),
                  hour_utc = as.integer(format(d$timestamp_utc, "%H", tz = "UTC")))
  )
}

#' PBLH agreement metrics
#'
#' On pairs where both values are finite: `bias = mean(mod - ref)`, `mae`,
#' `rmse`, Pearson and Spearman correlation. In log space, only pairs with
#' `mod > 0` and `ref > 0` (the others counted in `n_log_excluded`):
#' `median_log_ratio = median(log(mod / ref))`, `mad_log_ratio` = median
#' absolute deviation of `log(mod / ref)` from its median, scaled by 1.4826
#' (consistent estimator of the standard deviation for normal data; Rousseeuw
#' & Croux 1993; spec 04 amendment 8), `frac_within_factor2` = fraction with
#' `|log(mod / ref)| <= log 2`.
#' Errors in absolute and log-ratio form as in `docs/project_plan.md`
#' (Evaluation).
#'
#' @param mod,ref PBLH [m].
#' @return One-row tibble: `n` [-], `bias`, `mae`, `rmse` [m], `r_pearson`,
#'   `r_spearman` [-], `n_log` [-], `n_log_excluded` [-],
#'   `median_log_ratio`, `mad_log_ratio` [-, natural log],
#'   `frac_within_factor2` [-].
pblh_metrics <- function(mod, ref) {
  ok <- is.finite(mod) & is.finite(ref)
  m <- mod[ok]
  r <- ref[ok]
  n <- length(m)
  e <- m - r
  cr <- function(method) if (n > 2) suppressWarnings(stats::cor(m, r, method = method)) else NA_real_
  pos <- m > 0 & r > 0
  lr <- log(m[pos] / r[pos])
  nl <- length(lr)
  dplyr::tibble(
    n = n,
    bias = if (n > 0) mean(e) else NA_real_,
    mae = if (n > 0) mean(abs(e)) else NA_real_,
    rmse = if (n > 0) sqrt(mean(e^2)) else NA_real_,
    r_pearson = cr("pearson"),
    r_spearman = cr("spearman"),
    n_log = nl,
    n_log_excluded = n - nl,
    median_log_ratio = if (nl > 0) stats::median(lr) else NA_real_,
    mad_log_ratio = if (nl > 0) stats::mad(lr, constant = 1.4826) else NA_real_,
    frac_within_factor2 = if (nl > 0) mean(abs(lr) <= log(2)) else NA_real_
  )
}

#' PBLH metrics per stratum
#'
#' @param df Data frame with the `mod`, `ref` columns and the stratum columns.
#' @param mod,ref Column names [character].
#' @param by Stratum column names [character]; `NULL` = whole table.
#' @return Tibble: the `by` columns and the columns of [pblh_metrics()].
pblh_metrics_by <- function(df, mod = "mod", ref = "ref", by = NULL) {
  if (length(by) == 0) return(pblh_metrics(df[[mod]], df[[ref]]))
  key <- interaction(df[by], drop = TRUE, lex.order = TRUE)
  parts <- split(seq_len(nrow(df)), key)
  out <- lapply(parts, function(i) {
    dplyr::bind_cols(df[i[1], by, drop = FALSE], pblh_metrics(df[[mod]][i], df[[ref]][i]))
  })
  out <- dplyr::bind_rows(out)
  out[do.call(order, unname(as.list(out[by]))), ]
}

#' Metrics of a pair table in long format
#'
#' Strata: `all`, and each value of `regime`, `season`, `hour_utc`
#' (stratum label `"<column>=<value>"`).
#'
#' @param pt Output of [pair_table()] (one or more pairs).
#' @param strata Stratum columns [character].
#' @return Long tibble: `pair`, `stratum`, `metric`, `value`.
pblh_metrics_long <- function(pt, strata = c("regime", "season", "hour_utc")) {
  res <- list()
  for (p in unique(pt$pair)) {
    d <- pt[pt$pair == p, ]
    res[[length(res) + 1]] <- dplyr::bind_cols(dplyr::tibble(pair = p, stratum = "all"), pblh_metrics(d$mod, d$ref))
    for (s in strata) {
      m <- pblh_metrics_by(d[!is.na(d[[s]]), ], by = s)
      if (nrow(m) == 0) next
      lab <- paste0(s, "=", as.character(m[[s]]))
      m[[s]] <- NULL
      res[[length(res) + 1]] <- dplyr::bind_cols(dplyr::tibble(pair = p, stratum = lab), m)
    }
  }
  wide <- dplyr::bind_rows(res)
  mcols <- setdiff(names(wide), c("pair", "stratum"))
  dplyr::tibble(
    pair = rep(wide$pair, times = length(mcols)),
    stratum = rep(wide$stratum, times = length(mcols)),
    metric = rep(mcols, each = nrow(wide)),
    value = unlist(lapply(mcols, function(k) as.numeric(wide[[k]])), use.names = FALSE)
  )
}

#' Diurnal summary
#'
#' @param df Data frame with `value_cols` and the `by` columns.
#' @param value_cols Columns to summarise [character].
#' @param by Grouping columns [character]; default `c("season", "hour_utc")`.
#' @return Long tibble: the `by` columns, `variable`, `median`, `p25`, `p75`
#'   [unit of the variable], `n` [-].
diurnal_summary <- function(df, value_cols, by = c("season", "hour_utc")) {
  out <- lapply(value_cols, function(v) {
    d <- df[!is.na(df[[v]]), c(by, v)]
    names(d)[names(d) == v] <- "x"
    s <- dplyr::summarise(
      d, .by = dplyr::all_of(by),
      median = stats::median(x),
      p25 = stats::quantile(x, 0.25, names = FALSE),
      p75 = stats::quantile(x, 0.75, names = FALSE),
      n = dplyr::n()
    )
    dplyr::mutate(s, variable = v, .after = dplyr::all_of(by))
  })
  out <- dplyr::bind_rows(out)
  out[do.call(order, unname(as.list(out[c("variable", by)]))), ]
}

#' Sensible heat flux used by PBLHmod
#'
#' Recomputed from the model output as `H = wtheta * rho * cp` with
#' `rho = p / (r_d T)` (ideal gas law for dry air; same constants as spec 02,
#' [pblh_constants()]).
#'
#' @param out Output of [run_pblh_model()] (`wtheta` [K m s-1], `p_Pa` [Pa],
#'   `T_K` [K]).
#' @param const Constants from [pblh_constants()].
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc`, `H_ec_Wm2`
#'   [W m-2, positive upward].
ec_h_from_pblhmod <- function(out, const = pblh_constants()) {
  rho <- out$p_Pa / (const$r_d * out$T_K) # dry-air density [kg m-3]
  dplyr::tibble(
    timestamp_start_utc = out$timestamp_start_utc,
    timestamp_end_utc = out$timestamp_end_utc,
    H_ec_Wm2 = out$wtheta * rho * const$cp
  )
}

#' Cloud-base class
#'
#' @param cbh ERA5 cloud-base height [m above surface; `NA` = no cloud].
#' @return Factor: `no_cloud` (`NA`), `cbh_le_2000` (<= 2000 m), `cbh_gt_2000`.
cloud_class <- function(cbh) {
  factor(ifelse(is.na(cbh), "no_cloud", ifelse(cbh <= 2000, "cbh_le_2000", "cbh_gt_2000")),
         levels = c("no_cloud", "cbh_le_2000", "cbh_gt_2000"))
}

#' Free-atmosphere potential-temperature gradient 850-700 hPa
#'
#' `gamma_850_700 = (theta_700 - theta_850) / (zg_700 - zg_850)` per time
#' from the ERA5 pressure levels (context for the constant gamma of PBLHmod;
#' not a model input in spec 04).
#'
#' @param pl Output of [read_era5_pressure_levels()]; the tower cell is used.
#' @return Tibble: `timestamp_utc`, `gamma_850_700` [K m-1].
gamma_850_700 <- function(pl) {
  pl <- pl[pl$cell_row == 0 & pl$cell_col == 0, ]
  a <- pl[pl$p_hPa == 850, c("timestamp_utc", "theta_K", "zg_m")]
  b <- pl[pl$p_hPa == 700, c("timestamp_utc", "theta_K", "zg_m")]
  if (nrow(a) == 0 || nrow(b) == 0) cli::cli_abort("ERA5 pressure levels: 850 or 700 hPa missing.")
  d <- dplyr::inner_join(a, b, by = "timestamp_utc", suffix = c("_850", "_700"))
  dplyr::tibble(
    timestamp_utc = d$timestamp_utc,
    gamma_850_700 = (d$theta_K_700 - d$theta_K_850) / (d$zg_m_700 - d$zg_m_850)
  )
}

#' Monthly climatology of gamma_850_700
#'
#' @param g Output of [gamma_850_700()].
#' @param hours UTC hours kept [integer]; default `c(0, 12)`.
#' @return Tibble: `hour_utc`, `month` [1-12], `median`, `p25`, `p75`
#'   [K m-1], `n` [-].
gamma_monthly <- function(g, hours = c(0L, 12L)) {
  g$hour_utc <- as.integer(format(g$timestamp_utc, "%H", tz = "UTC"))
  g$month <- as.integer(format(g$timestamp_utc, "%m", tz = "UTC"))
  g <- g[g$hour_utc %in% hours & !is.na(g$gamma_850_700), ]
  out <- dplyr::summarise(
    g, .by = c(hour_utc, month),
    median = stats::median(gamma_850_700),
    p25 = stats::quantile(gamma_850_700, 0.25, names = FALSE),
    p75 = stats::quantile(gamma_850_700, 0.75, names = FALSE),
    n = dplyr::n()
  )
  out[order(out$hour_utc, out$month), ]
}
