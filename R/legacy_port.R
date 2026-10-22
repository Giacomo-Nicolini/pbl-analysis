# Port of the operational legacy PBLH code (spec 03 section 3, amendments 1-2):
# legacy/calc_PBL_height_operational.R + legacy/proc_PBL_height.R, the versions
# sourced by the ETC pipeline (legacy/FFP_analysis_dev.R). With the default
# switches the port reproduces the original to < 1e-8 m (tests/testthat/
# test-legacy_port.R) including its quirks; each switch moves one legacy choice
# to its PBLHmod v2 value for the attribution (spec 03 section 5).
# The legacy files are never sourced here (CLAUDE.md rule 5).

#' Switches of the legacy port
#'
#' Defaults reproduce the operational legacy code exactly ("legacy values", not
#' in [pblh_constants()]). Use [legacy_switches_pblhmod()] for the PBLHmod v2
#' value of each switch.
#'
#' @param ... Switches to override:
#'   - `gamma` [K m-1] 0.0025 (PBLHmod: 0.01); `k` [-] 0.41 (0.4); `omega`
#'     [rad s-1] 2 pi / 86400 (7.2921159e-5);
#'   - `stability_classes`: `"essa_operational"` (1/L classes of Essa et al.
#'     2006 with the operational unstable bound -0.003, so the near-neutral class
#'     is unreachable) or `"pblhmod"` (regime from L sign + daily window, i.e.
#'     [run_pblh_model()]); must be set together with `convective_trigger`;
#'   - `convective_trigger`: `"h_positive"` (B&G wherever H > 0, row 1 skipped)
#'     or `"daily_window"` (spec 02);
#'   - `stable_formula`: `"zilitinkevich"` (0.4 sqrt(u* L / |f|), Nieuwstadt
#'     fallback in the "stable" class) or `"nieuwstadt"` (Kljun Eq. B1);
#'   - `neutral_formula`: `"hanna_chang"` (0.3 u* / |f|, Nieuwstadt fallback)
#'     or `"auto"` (B2 / B3, [h_neutral_select()]);
#'   - `bg_integration`: `"euler"` (one explicit step per half-hour, legacy
#'     formula) or `"rk4"` ([bg_step()]);
#'   - `gap_interp_before_bg`: `TRUE` (linear interpolation of gaps of at most
#'     2 values in h0 before B&G and in the B&G series, `zoo::na.approx(maxgap =
#'     2)` semantics) or `FALSE`;
#'   - `smoothing`: `TRUE` (5-point centred moving mean, `NA` if any value is
#'     `NA`) or `FALSE`;
#'   - `floor_m` [m]: 50 or `NA` (no floor); `negative_to_na`: `TRUE`;
#'   - `fill`: `"run_mean_day_night"` (NA -> mean of the day / night values of
#'     the run, day = hour 06-18 of the timestamp) or `"none"`;
#'   - `mol_locf`: `TRUE` (L filled by last observation carried forward) or
#'     `FALSE`;
#'   - `forcing`: `"etc"` (ETC inputs: L2 L, u*, H, EddyPro AT / rho / cp) or
#'     `"pblhmod"` ([prepare_pblh_forcing()]: L, u*, H = H_F_MDS, T_K, wtheta).
#'   Fixed legacy constants: `R_d` 287.058 J kg-1 K-1, `g` 9.80665 m s-2, B&G
#'   `A` 0.2, `B` 2.5, `C` 8.
#' @return Named list.
legacy_switches <- function(...) {
  defaults <- list(
    gamma = 0.0025, k = 0.41, omega = 2 * pi / 86400,
    stability_classes = "essa_operational", convective_trigger = "h_positive",
    stable_formula = "zilitinkevich", neutral_formula = "hanna_chang", bg_integration = "euler",
    gap_interp_before_bg = TRUE, smoothing = TRUE, floor_m = 50, negative_to_na = TRUE,
    fill = "run_mean_day_night", mol_locf = TRUE, forcing = "etc",
    R_d = 287.058, g = 9.80665, bg_A = 0.2, bg_B = 2.5, bg_C = 8
  )
  user <- list(...)
  bad <- setdiff(names(user), names(defaults))
  if (length(bad) > 0) cli::cli_abort("Unknown legacy switch{?es}: {.val {bad}}.")
  sw <- utils::modifyList(defaults, user)
  choices <- list(stability_classes = c("essa_operational", "pblhmod"), convective_trigger = c("h_positive", "daily_window"),
                  stable_formula = c("zilitinkevich", "nieuwstadt"), neutral_formula = c("hanna_chang", "auto"),
                  bg_integration = c("euler", "rk4"), fill = c("run_mean_day_night", "none"), forcing = c("etc", "pblhmod"))
  for (k in names(choices)) {
    if (!(length(sw[[k]]) == 1 && sw[[k]] %in% choices[[k]])) {
      cli::cli_abort("Switch {.field {k}} must be one of {.val {choices[[k]]}}.")
    }
  }
  sw
}

#' PBLHmod v2 value of each legacy switch
#'
#' @return Named list (the switches of [legacy_switches()] that differ).
legacy_switches_pblhmod <- function() {
  const <- pblh_constants()
  list(
    fill = "none", floor_m = NA_real_, smoothing = FALSE, gap_interp_before_bg = FALSE, mol_locf = FALSE,
    forcing = "pblhmod", k = const$k, omega = const$omega, gamma = pblh_options()$gamma,
    neutral_formula = "auto", stable_formula = "nieuwstadt",
    stability_classes = "pblhmod", convective_trigger = "daily_window", bg_integration = "rk4"
  )
}

# Pasquill-Gifford classes from 1/L, operational bounds (Essa et al. 2006; the operational
# code moved the unstable bound from -0.016 to -0.003, so "near-neutral" is unreachable).
legacy_classes <- function(inv_l) {
  out <- rep(NA_character_, length(inv_l))
  out[inv_l < -0.056] <- "strongly-unstable"
  out[inv_l >= -0.056 & inv_l < -0.003] <- "unstable"
  out[inv_l >= -0.003 & inv_l < 0.002] <- "neutral" # [-0.004, -0.003) already "unstable"
  out[inv_l >= 0.002 & inv_l < 0.006] <- "stable"
  out[inv_l >= 0.006 & inv_l < 0.022] <- "strongly-stable"
  out[inv_l >= 0.022] <- "nighttime-stable"
  out
}

# Post-processing of proc_PBL_height() on one run: negatives -> NA, floor, run day/night mean fill.
legacy_post <- function(h, hour, sw) {
  if (isTRUE(sw$negative_to_na)) h[!is.na(h) & h < 0] <- NA
  flag_floor <- rep(FALSE, length(h))
  if (!is.na(sw$floor_m)) {
    flag_floor <- !is.na(h) & h < sw$floor_m
    h[flag_floor] <- sw$floor_m
  }
  flag_fill <- rep(FALSE, length(h))
  if (sw$fill == "run_mean_day_night") {
    day <- hour >= 6 & hour <= 18 # clock hours of the timestamp (legacy issue 12)
    flag_fill <- is.na(h)
    m_day <- mean(h[day], na.rm = TRUE)    # computed after the floor (legacy issue 13)
    m_night <- mean(h[!day], na.rm = TRUE)
    h[flag_fill & day] <- m_day
    h[flag_fill & !day] <- m_night
  }
  list(h = h, flag_floor = flag_floor, flag_fill = flag_fill)
}

#' Legacy PBLH, ported (operational ETC version)
#'
#' Reimplements `calc_PBL_height()` (operational version) followed by
#' `proc_PBL_height()`, per processing run, with [legacy_switches()].
#' Legacy algorithm (default switches), per run:
#' 1. `f = 2 omega sin(lat)`; `rho`, `cp` as passed (if `NULL`: `AP / (R_d AT)`,
#'    `1005 + (AT - 273.15 + 23.12)^2 / 3364`; `AP` if `NULL`: `101300 exp(-(altitude
#'    + zm) / 8200)`); `wtheta = H / (rho cp)`.
#' 2. Classes from 1/L ([legacy_classes()]); h0: unstable classes `NA`; neutral
#'    `0.3 u* / |f|` (Hanna & Chang 1993), Nieuwstadt fallback; stable `0.4
#'    sqrt(u* L / |f|)` (Zilitinkevich), Nieuwstadt fallback; strongly / nighttime
#'    stable `0.4 sqrt(u* L / |f|)`. A height is attempted if L or u* is present.
#'    Deviation (documented): where L is missing and u* present, the original
#'    code stops (class test on `NA`; masked by LOCF from the first run row); the
#'    port returns `NA`.
#' 3. Gaps of at most 2 values in h0 interpolated; B&G (Batchvarova & Gryning
#'    1991; Kljun et al. 2015 Eq. B5) where `H > 0` and not on row 1, from the
#'    interpolated h0 of the previous row, else the previous B&G value; one Euler
#'    step `h_prev + dt dh/dt` with the forcing of the row; B&G gaps of at most 2
#'    values interpolated; h0 gaps replaced by B&G.
#' 4. 5-point centred moving mean (`stats::filter`, `NA` near any `NA`).
#' 5. `proc_PBL_height()`: negatives -> `NA`, floor at 50 m, `NA` -> mean of the
#'    run's day (hour 06-18 of `TIMESTAMP`) or night values (after the floor).
#' With `convective_trigger = "daily_window"` (and `stability_classes =
#' "pblhmod"`) the core is [run_pblh_model()] on `icos` (PBLHmod v2 with
#' `gamma` and `bg_integration` from the switches), followed by the enabled
#' post-processing steps; all other switches must have their PBLHmod values.
#'
#' @param inp Output of [legacy_inputs()] (or [legacy_inputs_cleanset()]).
#' @param site_meta One row of [read_sites()] (`lat`, `altitude_m`).
#' @param sw Output of [legacy_switches()].
#' @param forcing Output of [prepare_pblh_forcing()] on the site data
#'   (needed for `forcing = "pblhmod"`), matched to `inp` by time.
#' @param icos Output of [read_icos_site()] (needed for the daily window).
#' @return Tibble, one row per row of `inp`: `pblh_m` [m], `atmos_cond`,
#'   `h0_m` [m, before interpolation], `h_bg_m` [m, B&G before interpolation],
#'   `h_smooth_m` [m, before post-processing], `flag_floor`, `flag_fill`
#'   [logical].
legacy_pblh <- function(inp, site_meta, sw = legacy_switches(), forcing = NULL, icos = NULL) {
  pm <- legacy_switches_pblhmod()
  window <- sw$convective_trigger == "daily_window" || sw$stability_classes == "pblhmod"
  if (window) {
    need <- c("stability_classes", "convective_trigger", "forcing", "mol_locf", "gap_interp_before_bg", "k", "omega",
              "neutral_formula", "stable_formula")
    off <- need[!vapply(need, function(k) identical(sw[[k]], pm[[k]]), logical(1))]
    if (length(off) > 0) cli::cli_abort("The daily window needs the PBLHmod value of: {.field {off}}.")
    if (is.null(icos)) cli::cli_abort("{.arg icos} is needed for the daily window.")
  }
  if (sw$forcing == "pblhmod" && is.null(forcing) && !window) cli::cli_abort("{.arg forcing} is needed for forcing = \"pblhmod\".")
  n <- nrow(inp)
  out <- dplyr::tibble(pblh_m = rep(NA_real_, n), atmos_cond = NA_character_, h0_m = NA_real_, h_bg_m = NA_real_,
                       h_smooth_m = NA_real_, flag_floor = FALSE, flag_fill = FALSE)
  hour <- as.integer(substr(inp$TIMESTAMP, 12, 13))

  if (window) {
    mod <- run_pblh_model(icos, site_meta, pblh_options(gamma = sw$gamma, bg_integration = sw$bg_integration))
    m <- match(as.numeric(inp$timestamp_start_utc), as.numeric(mod$timestamp_start_utc))
    core <- mod$pblh_m[m]
    out$atmos_cond <- mod$regime[m]
    out$h0_m <- mod$h_diag_m[m]
    out$h_bg_m <- mod$h_bg_m[m]
  } else {
    core <- rep(NA_real_, n)
  }

  for (r in sort(unique(inp$run_id))) {
    i <- which(inp$run_id == r)
    h_core <- if (window) core[i] else {
      res <- legacy_core_run(inp[i, ], site_meta, sw, forcing, attr(inp, "pass_null"))
      out$atmos_cond[i] <- res$atmos_cond
      out$h0_m[i] <- res$h0
      out$h_bg_m[i] <- res$bg
      res$h1
    }
    h_s <- if (isTRUE(sw$smoothing)) as.numeric(stats::filter(h_core, rep(1 / 5, 5), sides = 2)) else h_core
    out$h_smooth_m[i] <- h_s
    pp <- legacy_post(h_s, hour[i], sw)
    out$pblh_m[i] <- pp$h
    out$flag_floor[i] <- pp$flag_floor
    out$flag_fill[i] <- pp$flag_fill
  }
  out
}

# Core of the operational calc_PBL_height() on one run (steps 1-3 of legacy_pblh()).
legacy_core_run <- function(x, site_meta, sw, forcing, pass_null) {
  n <- nrow(x)
  f <- 2 * sw$omega * sin(site_meta$lat * pi / 180) # Coriolis parameter (legacy omega = 2 pi / 86400)
  if (sw$forcing == "etc") {
    L <- if (isTRUE(sw$mol_locf)) x$MOL else x$MOL_raw
    ustar <- x$ustar
    H <- x$H
    AT <- x$AT
    AP <- if ("AP" %in% pass_null) 101300 * exp(-(site_meta$altitude_m + x$zm) / 8200) else x$AP # Campbell & Norman 1998
    rho <- if ("rho" %in% pass_null) AP / (sw$R_d * AT) else x$rho
    cp <- if ("cp" %in% pass_null) 1005 + ((AT - 273.15) + 23.12)^2 / 3364 else x$cp
    wtheta <- H / (rho * cp)
  } else {
    m <- match(as.numeric(x$timestamp_start_utc), as.numeric(forcing$timestamp_start_utc))
    L <- forcing$L[m]
    if (isTRUE(sw$mol_locf)) L <- locf(L)
    ustar <- forcing$ustar[m]
    H <- forcing$H[m]
    AT <- forcing$T_K[m]
    wtheta <- forcing$wtheta[m]
  }
  inv_l <- 1 / L
  attempt <- !is.na(inv_l) | !is.na(1 / ustar)
  cls <- legacy_classes(inv_l)
  cls[!attempt] <- NA_character_
  nieu_legacy <- L / 3.8 * (-1 + (1 + 2.28 * ustar / (f * L))^(1 / 2)) # Nieuwstadt (1981), legacy inline form
  zil <- 0.4 * (ustar * L / abs(f))^(1 / 2)                            # Zilitinkevich et al. (2012)
  h0 <- rep(NA_real_, n)
  neu <- cls %in% "neutral"
  if (sw$neutral_formula == "hanna_chang") {
    hn <- 0.3 * ustar / abs(f) # Hanna & Chang (1993)
    hn[is.na(hn)] <- nieu_legacy[is.na(hn)]
  } else {
    hn <- h_neutral_select(ustar, f, "auto", T_K = AT, gamma = sw$gamma)$h
  }
  h0[neu] <- hn[neu]
  st <- cls %in% "stable"
  ss <- cls %in% c("strongly-stable", "nighttime-stable")
  if (sw$stable_formula == "zilitinkevich") {
    hs <- zil
    hs_fb <- ifelse(is.na(hs), nieu_legacy, hs) # Nieuwstadt fallback, "stable" class only
    h0[st] <- hs_fb[st]
    h0[ss] <- zil[ss]
  } else {
    hb1 <- h_nieuwstadt(L, ustar, f)
    h0[st | ss] <- hb1[st | ss]
  }
  h0i <- if (isTRUE(sw$gap_interp_before_bg)) interp_short_gaps(h0, 2)$x else h0
  # B&G (Batchvarova & Gryning 1991; Kljun et al. 2015 Eq. B5), sequential.
  dt <- 1800
  A <- sw$bg_A
  B <- sw$bg_B
  C <- sw$bg_C
  bg <- rep(NA_real_, n)
  const <- utils::modifyList(pblh_constants(), list(k = sw$k, g = sw$g, bg_A = A, bg_B = B, bg_C = C))
  for (i in seq_len(n)) {
    if (i != 1 && !is.na(H[i]) && H[i] > 0) {
      hp <- if (!is.na(h0i[i - 1])) h0i[i - 1] else bg[i - 1]
      if (sw$bg_integration == "euler") {
        dh <- dt * (wtheta[i] / sw$gamma) /
          ((hp^2 / ((1 + 2 * A) * hp - 2 * B * sw$k * L[i])) +
             ((C * AT[i] * ustar[i]^2) / (sw$gamma * sw$g * ((1 + A) * hp - B * sw$k * L[i]))))
        bg[i] <- hp + dh
      } else if (!is.na(hp) && !is.na(wtheta[i]) && !is.na(ustar[i]) && !is.na(L[i]) && !is.na(AT[i])) {
        s <- bg_step(hp, wtheta[i], ustar[i], L[i], AT[i], sw$gamma, dt = dt, const = const, method = "rk4")
        bg[i] <- if (s$ok) s$h_end else NA_real_ # end-of-interval height, as the legacy step
      }
    }
  }
  bgi <- if (isTRUE(sw$gap_interp_before_bg)) interp_short_gaps(bg, 2)$x else bg
  h1 <- ifelse(is.na(h0i), bgi, h0i)
  list(h1 = h1, h0 = h0, bg = bg, atmos_cond = cls)
}
