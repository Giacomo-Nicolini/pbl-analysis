# PBLHmod v2: options and orchestration (daily cycle).
# Spec: docs/specs/02_pblh_model.md section 6 and Amendments.

#' PBLH model options
#'
#' Returns the model options with defaults; any can be overridden through
#' `...`. Unknown option names or invalid values abort.
#'
#' @param ... Options to override:
#'   - `stable_method` [`"nieuwstadt"` (Eq. B1, default), `"zilitinkevich"`
#'     (Eq. B4), `"neutral"` (Eq. B2)] formula for L > 0;
#'   - `c_s` [-] coefficient of Eq. B4 (default 0.4);
#'   - `gamma` [K m-1] potential-temperature gradient above the CBL, scalar or
#'     one value per half-hour (default 0.01; Batchvarova & Gryning 1991 as cited
#'     in Kljun et al. 2015);
#'   - `h_source` [`"H_F_MDS"` (default), `"H"`, `"H_CORR"`] sensible heat flux;
#'   - `h_qc_max` [-] highest accepted `H_F_MDS_QC` (default 1);
#'   - `h_fallback` [logical] use FLUXES `H` where `H_F_MDS` is missing in the
#'     file (default `TRUE`);
#'   - `max_gap_interp` [half-hours] longest forcing gap filled by linear
#'     interpolation (default 2);
#'   - `n_sub` [-] RK4 sub-steps per half-hour (default 30);
#'   - `n_neg_end` [half-hours] consecutive `wtheta <= 0` ending the convective
#'     window (default 2);
#'   - `neutral_method` [`"auto"` (default), `"B2"`, `"B3"`] neutral formula
#'     (unstable rows outside the window, `combine = "max"`, `stable_method =
#'     "neutral"`): `"auto"` uses Eq. B2 where 10 < N/|f| < 70, else Eq. B3
#'     (Kljun et al. 2015 App. B), N = sqrt(g gamma / T);
#'   - `combine` [`"bg_only"` (default, Kljun et al. 2015 App. B) or `"max"`
#'     (max of B5 and the diagnostic height)] PBLH in the convective window;
#'   - `ustar_ws_max` [-] largest plausible u*/WS (default 0.5; `Inf`
#'     disables the ratio test), [ustar_implausible()], spec 02 Amendment 12;
#'   - `ws_min_ratio_test` [m s-1] ratio test only where WS >= this (default 1);
#'   - `ustar_max` [m s-1] absolute u* cap (default 3; `Inf` disables);
#'   - `label` [character] run label (default `"default"`).
#' @return Named list of options.
pblh_options <- function(...) {
  defaults <- list(
    stable_method = "nieuwstadt",
    c_s = 0.4,
    gamma = 0.01,
    h_source = "H_F_MDS",
    h_qc_max = 1,
    h_fallback = TRUE,
    max_gap_interp = 2,
    neutral_method = "auto",
    n_sub = 30,
    n_neg_end = 2,
    combine = "bg_only",
    ustar_ws_max = 0.5,
    ws_min_ratio_test = 1,
    ustar_max = 3,
    label = "default"
  )
  user <- list(...)
  unknown <- setdiff(names(user), names(defaults))
  if (length(unknown) > 0 || (length(user) > 0 && is.null(names(user)))) {
    cli::cli_abort(c(
      "Unknown PBLH option{?s}: {.val {unknown}}",
      "i" = "Known options: {.val {names(defaults)}}"
    ))
  }
  opts <- utils::modifyList(defaults, user)
  validate_pblh_options(opts)
  opts
}

validate_pblh_options <- function(opts) {
  one_of <- function(name, allowed) {
    v <- opts[[name]]
    if (!(is.character(v) && length(v) == 1 && v %in% allowed)) {
      cli::cli_abort("Option {.field {name}} must be one of {.val {allowed}}, not {.val {v}}.")
    }
  }
  pos_num <- function(name) {
    v <- opts[[name]]
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v) && v > 0)) {
      cli::cli_abort("Option {.field {name}} must be a single positive number.")
    }
  }
  int_ge <- function(name, min) {
    v <- opts[[name]]
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v) && v == round(v) && v >= min)) {
      cli::cli_abort("Option {.field {name}} must be a single integer >= {min}.")
    }
  }
  one_of("stable_method", c("nieuwstadt", "zilitinkevich", "neutral"))
  one_of("h_source", c("H_F_MDS", "H", "H_CORR"))
  one_of("combine", c("max", "bg_only"))
  one_of("neutral_method", c("auto", "B2", "B3"))
  pos_num("c_s")
  g <- opts$gamma
  if (!(is.numeric(g) && length(g) >= 1 && all(is.finite(g)) && all(g > 0))) {
    cli::cli_abort("Option {.field gamma} must be positive and finite (scalar or one value per half-hour) [K m-1].")
  }
  int_ge("h_qc_max", 0)
  if (!(is.logical(opts$h_fallback) && length(opts$h_fallback) == 1 && !is.na(opts$h_fallback))) {
    cli::cli_abort("Option {.field h_fallback} must be TRUE or FALSE.")
  }
  int_ge("max_gap_interp", 0)
  int_ge("n_sub", 1)
  int_ge("n_neg_end", 1)
  for (k in c("ustar_ws_max", "ws_min_ratio_test", "ustar_max")) {
    v <- opts[[k]]
    if (!(is.numeric(v) && length(v) == 1 && !is.na(v) && v > 0)) {
      cli::cli_abort("Option {.field {k}} must be a single positive number or {.code Inf}.")
    }
  }
  if (!(is.character(opts$label) && length(opts$label) == 1 && nzchar(opts$label))) {
    cli::cli_abort("Option {.field label} must be a non-empty string.")
  }
  invisible(TRUE)
}

#' Run the PBLH model (PBLHmod v2)
#'
#' Daily cycle of Kljun et al. (2015) App. B: diagnostic height (Eqs. B1-B4,
#' [h_diagnostic()]) at night and outside the convective window; Batchvarova &
#' Gryning (1991) slab model (Eq. B5, [bg_step()]) integrated from sunrise.
#' Rules: spec 02 section 6 and Amendments 2, 5-8.
#' - A day starts at the first half-hour with `is_day` (solar elevation > 0 at
#'   the midpoint) and lasts until the next one; `day_id` is the local standard
#'   time date of that half-hour (polar day/night dates: the local date).
#' - `h0` = diagnostic height of the half-hour before sunrise, else the nearest
#'   valid value within the 2 h before; none: `flag_no_h0`, diagnostic only.
#' - No `wtheta > 0` between sunrise and sunset: `flag_no_positive_h`,
#'   diagnostic only; if the daylight `wtheta` also has gaps, the day is treated
#'   as a gap day instead (`flag_bg_stopped_gap`). Polar day: `flag_polar_day`,
#'   diagnostic only.
#' - In the window: h held at `h0` until the first `wtheta > 0`; then Eq. B5;
#'   negative spells shorter than `opts$n_neg_end` held; the window ends at the
#'   first run of `opts$n_neg_end` half-hours with `wtheta <= 0` or at sunset.
#'   `wtheta > 0` with `L > 0`, or undefined rate: h held,
#'   `flag_sign_inconsistent`. Missing forcing: integration stops for the rest
#'   of the day (`flag_bg_stopped_gap`).
#' - `pblh_m = h_bg_m` (`combine = "bg_only"`, default) or
#'   `max(h_bg_m, h_diag_m)` (`"max"`) in the window. Output is raw: no smoothing, clipping or filling.
#'
#' @param data Output of [read_icos_site()] (complete 30-min grid, sorted).
#' @param site_meta One row of [read_sites()] (`site_id`, `lat`, `lon`,
#'   `utc_offset_h`; `altitude_m` for the barometric fallback).
#' @param opts Output of [pblh_options()].
#' @return Tibble, one row per input row: `site_id`, `timestamp_start_utc`,
#'   `timestamp_end_utc` [POSIXct UTC], `pblh_m` [m], `regime` (`"stable"`,
#'   `"unstable_diag"`, `"convective"`, `"excluded"`, `NA`), `method` (`"B1"`,
#'   `"B2"`, `"B3"`, `"B4"`, `"B5"`, `NA`), `h_diag_m` [m], `h_bg_m` [m, window only],
#'   `day_id` [Date, local standard time], `is_day` [logical], `ustar` [m s-1],
#'   `L` [m], `wtheta` [K m s-1], `T_K` [K], `p_Pa` [Pa], `H_source`,
#'   `T_source`, `p_source`, `zm_m` [m], `WS` [m s-1], flags (logical):
#'   `flag_ustar_implausible` (spec 02 Amendment 12), `flag_ustar_interp`,
#'   `flag_L_interp`, `flag_H_interp`, `flag_H_qc_masked`, `flag_no_h0`,
#'   `flag_bg_stopped_gap`, `flag_sign_inconsistent`, `flag_no_positive_h`,
#'   `flag_polar_day`, `flag_below_zm` (`pblh_m < zm_m`), and `excl_flag`.
#'   Attribute `pblh_options`: the options used.
run_pblh_model <- function(data, site_meta, opts = pblh_options()) {
  validate_pblh_options(opts)
  if (nrow(site_meta) != 1 || is.na(site_meta$utc_offset_h)) {
    cli::cli_abort("{.arg site_meta} must be one row with non-missing {.field utc_offset_h}.")
  }
  n <- nrow(data)
  if (n > 1 && any(diff(as.numeric(data$timestamp_start_utc)) != 1800)) {
    cli::cli_abort("{.arg data} must be on a complete, sorted 30-min grid (output of {.fn read_icos_site}).")
  }
  gamma <- opts$gamma
  if (length(gamma) == 1) {
    gamma <- rep(gamma, n)
  } else if (length(gamma) != n) {
    cli::cli_abort("Option {.field gamma} must be a scalar or have one value per row ({n}), not {length(gamma)}.")
  }
  const <- pblh_constants()
  f <- prepare_pblh_forcing(data, site_meta, opts)
  dg <- h_diagnostic(f$L, f$ustar, f$f_cor, opts, T_K = f$T_K, gamma = gamma)

  pblh <- dg$h_diag_m
  method <- dg$method
  regime <- ifelse(is.na(f$L), NA_character_, ifelse(f$L > 0, "stable", "unstable_diag"))
  h_bg <- rep(NA_real_, n)
  flag_no_h0 <- flag_stop <- flag_sign <- flag_nopos <- flag_polar <- rep(FALSE, n)

  # --- days: local standard time dates, polar status, day starts ---
  local_date <- as.Date(format(f$timestamp_start_utc + site_meta$utc_offset_h * 3600,
                               "%Y-%m-%d", tz = "UTC"))
  udates <- unique(local_date)
  polar_by_date <- sun_times(udates, site_meta$lat, site_meta$lon)$polar
  polar <- polar_by_date[match(local_date, udates)]
  dd <- dplyr::coalesce(f$is_day, FALSE)
  day_start <- which(dd & c(TRUE, !dd[-n]))
  day_id <- local_date
  seg_end <- integer(0)
  if (length(day_start) > 0) {
    seg_end <- c(day_start[-1] - 1, n)
    for (j in seq_along(day_start)) {
      day_id[day_start[j]:seg_end[j]] <- local_date[day_start[j]]
    }
    # night before the first sunrise of the series: belongs to the previous day
    if (day_start[1] > 1) day_id[seq_len(day_start[1] - 1)] <- local_date[day_start[1]] - 1
  }
  day_id[!is.na(polar)] <- local_date[!is.na(polar)] # polar day/night: local date (Amendment 2)

  # --- convective window per day ---
  for (j in seq_along(day_start)) {
    s <- day_start[j]
    day_rows <- s:seg_end[j]
    e <- s
    while (e < n && dd[e + 1]) e <- e + 1 # last daylight half-hour (sunset)
    light <- s:e

    if (any(!is.na(polar[light]) & polar[light] == "day")) {
      flag_polar[light] <- TRUE # no sunrise: diagnostic only (Amendment 2)
      next
    }
    if (!any(f$wtheta[light] > 0, na.rm = TRUE)) {
      if (anyNA(f$wtheta[light])) {
        flag_stop[day_rows] <- TRUE # no positive value but gaps: undecidable, treated as a gap day
      } else {
        flag_nopos[day_rows] <- TRUE # no positive heat flux: no window (Amendment 5)
      }
      next
    }
    # initial height: half-hour before sunrise, else nearest valid within 2 h (Amendment 8)
    cand <- seq.int(s - 1, max(1, s - const$h0_search_halfhours))
    cand <- cand[cand >= 1 & cand < s]
    cand <- cand[!is.na(dg$h_diag_m[cand])]
    if (length(cand) == 0) {
      flag_no_h0[day_rows] <- TRUE
      next
    }
    h <- dg$h_diag_m[cand[1]]

    seen_pos <- FALSE
    i <- s
    while (i <= e) {
      w <- f$wtheta[i]
      if (anyNA(c(w, f$ustar[i], f$L[i], f$T_K[i], gamma[i]))) {
        flag_stop[day_rows] <- TRUE # forcing gap: stop for the rest of the day
        break
      }
      if (seen_pos && w <= 0) {
        # run of wtheta <= 0 starting here (within daylight, until a gap)
        r <- 0
        while (i + r <= e && !is.na(f$wtheta[i + r]) && f$wtheta[i + r] <= 0) r <- r + 1
        if (r >= opts$n_neg_end) break # window ends before this run
      }
      if (!seen_pos && w > 0) seen_pos <- TRUE
      if (seen_pos && w > 0) {
        if (f$L[i] > 0) {
          flag_sign[i] <- TRUE # positive heat flux with stable L: hold
          h_bg[i] <- h
        } else {
          st <- bg_step(h, w, f$ustar[i], f$L[i], f$T_K[i], gamma[i],
                        dt = 1800, n_sub = opts$n_sub, const = const)
          if (!st$ok) flag_sign[i] <- TRUE
          h_bg[i] <- st$h_mean
          h <- st$h_end
        }
      } else {
        h_bg[i] <- h # before the first positive wtheta, or short negative spell: hold
      }
      regime[i] <- "convective"
      if (opts$combine == "max" && !is.na(dg$h_diag_m[i]) && dg$h_diag_m[i] > h_bg[i]) {
        pblh[i] <- dg$h_diag_m[i] # mechanical (diagnostic) term larger
        method[i] <- dg$method[i]
      } else {
        pblh[i] <- h_bg[i]
        method[i] <- "B5"
      }
      i <- i + 1
    }
  }

  # --- excluded rows ---
  excl <- f$excl_flag
  pblh[excl] <- NA_real_
  regime[excl] <- "excluded"
  method[excl] <- NA_character_
  h_bg[excl] <- NA_real_

  zm <- col_or_na(data, "zm_m")
  site_id <- if ("site_id" %in% names(data)) data$site_id else rep(site_meta$site_id, n)
  out <- dplyr::tibble(
    site_id = site_id,
    timestamp_start_utc = f$timestamp_start_utc,
    timestamp_end_utc = f$timestamp_end_utc,
    pblh_m = pblh,
    regime = regime,
    method = method,
    h_diag_m = dg$h_diag_m,
    h_bg_m = h_bg,
    day_id = day_id,
    is_day = f$is_day,
    ustar = f$ustar,
    WS = f$WS,
    L = f$L,
    wtheta = f$wtheta,
    T_K = f$T_K,
    p_Pa = f$p_Pa,
    H_source = f$H_source,
    T_source = f$T_source,
    p_source = f$p_source,
    zm_m = zm,
    flag_ustar_implausible = f$flag_ustar_implausible,
    flag_ustar_interp = f$flag_ustar_interp,
    flag_L_interp = f$flag_L_interp,
    flag_H_interp = f$flag_H_interp,
    flag_H_qc_masked = f$flag_H_qc_masked,
    flag_no_h0 = flag_no_h0,
    flag_bg_stopped_gap = flag_stop,
    flag_sign_inconsistent = flag_sign,
    flag_no_positive_h = flag_nopos,
    flag_polar_day = flag_polar,
    flag_below_zm = dplyr::coalesce(pblh < zm, FALSE),
    excl_flag = excl
  )
  attr(out, "pblh_options") <- opts
  out
}
