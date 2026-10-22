# Port of the operational legacy code (R/legacy_port.R) vs the original.
# This file and scripts/03_legacy_regression.R are the only places allowed to source legacy/
# (CLAUDE.md rule 5; spec 03), into a dedicated environment.

utc <- function(x) as.POSIXct(x, tz = "UTC")
site_l <- dplyr::tibble(site_id = "XX-Tst", lat = 51.31, lon = 4.52, altitude_m = 16, utc_offset_h = 1, zm_m = 30)

# Synthetic legacy input table (structure of legacy_inputs()), 3 days of half-hours:
# diurnal H with sign changes, L consistent with H (all 1/L classes), NA runs of 1, 2 and 5 rows in
# H, u*, L and the ambient variables, two processing runs.
make_legacy_inp <- function(n_days = 3) {
  n <- 48 * n_days
  t_loc <- utc("2024-06-01 00:30:00") + 1800 * (seq_len(n) - 1) # TIMESTAMP_END, local standard time
  hod <- as.numeric(format(t_loc - 900, "%H")) + as.numeric(format(t_loc - 900, "%M")) / 60
  H <- 250 * pmax(sin(pi * (hod - 5.5) / 14), 0) - 25 + 10 * sin(seq_len(n))
  ustar <- 0.25 + 0.2 * abs(sin(seq_len(n) / 7))
  AT <- 288 + 6 * sin(pi * (hod - 9) / 12)
  rho <- 1.2 - 0.0002 * (AT - 288)
  cp <- 1010 + 0.1 * (AT - 288)
  L <- -ustar^3 * AT * rho * cp / (0.41 * 9.80665 * H)
  L[c(20, 21, 70)] <- c(2e4, -3e4, 400) # near neutral and a stable spike
  H[c(30, 61, 62, 70:76)] <- NA         # gaps of 1, 2 and 7 rows (the last one by day: B&G chain stops)
  ustar[c(3:9, 103:109)] <- 0.03         # calm nights, very stable: very low h0 (floor at 50 m)
  L[c(3:9, 103:109)] <- 5
  ustar[c(40, 85, 86)] <- NA
  L_raw <- L
  L_raw[c(50, 51, 120:124)] <- NA
  rho[c(90)] <- NA
  dplyr::tibble(
    run_id = ifelse(seq_len(n) <= 96, 1L, 2L),
    timestamp_start_utc = t_loc - 1800 - 3600, timestamp_end_utc = t_loc - 3600,
    TIMESTAMP = format(t_loc, "%Y-%m-%d %H:%M:%S"),
    MOL = unlist(lapply(split(L_raw, ifelse(seq_len(n) <= 96, 1L, 2L)), locf), use.names = FALSE),
    MOL_raw = L_raw, ustar = ustar, H = H, AT = AT, AP = 101100, rho = rho, cp = cp, zm = 30, d = 14
  )
}

# The original pipeline (operational calc + proc), per run.
original_pipeline <- function(inp, lat, altit) {
  skip_if_not_installed("zoo")
  env <- new.env(parent = globalenv())
  sys.source(here::here("legacy", "calc_PBL_height_operational.R"), envir = env)
  sys.source(here::here("legacy", "proc_PBL_height.R"), envir = env)
  out <- rep(NA_real_, nrow(inp))
  for (r in unique(inp$run_id)) {
    i <- which(inp$run_id == r)
    x <- inp[i, ]
    calc <- suppressPackageStartupMessages(env$calc_PBL_height(
      TIMESTAMP = x$TIMESTAMP, LAT = lat, zm = x$zm, d = x$d, MOL = x$MOL, ustar = x$ustar, H = x$H,
      AT = x$AT, AP = x$AP, rho = x$rho, cp = x$cp, altit = altit, plot.atmc = FALSE, plot.PBL = FALSE, save.PBL = FALSE
    ))
    tmp <- withr::local_tempdir()
    ts <- as.POSIXct(x$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
    out[i] <- as.numeric(env$proc_PBL_height(TIMESTAMP = ts, PBL.h = calc, plot.PBL.h.freq = FALSE,
                                             PBL.dir = paste0(tmp, "/")))
  }
  out
}

test_that("port with default switches equals the original operational pipeline (< 1e-8 m, incl. NA)", {
  inp <- make_legacy_inp()
  orig <- original_pipeline(inp, site_l$lat, site_l$altitude_m)
  port <- legacy_pblh(inp, site_l)
  expect_identical(is.na(port$pblh_m), is.na(orig))
  expect_lt(max(abs(port$pblh_m - orig), na.rm = TRUE), 1e-8)
  # the synthetic series exercises the features
  expect_true(any(port$flag_fill) && any(port$flag_floor))
  expect_true(all(c("strongly-unstable", "unstable", "neutral", "stable", "nighttime-stable") %in% port$atmos_cond))
})

test_that("operational classes: unstable bound -0.003, near-neutral unreachable", {
  cl <- legacy_classes(1 / c(-10, -50, -200, -300, -400, 1e9, 300, 100, 20, NA))
  # 1/L: -0.1, -0.02, -0.005 (unstable, not near-neutral), -0.0033, -0.0025, ~0, 0.0033, 0.01, 0.05
  expect_equal(cl, c("strongly-unstable", "unstable", "unstable", "unstable", "neutral", "neutral", "stable",
                     "strongly-stable", "nighttime-stable", NA))
  expect_false("near-neutral" %in% cl)
})

test_that("post-processing switches: negatives, floor, run day/night mean fill", {
  h <- c(-5, 30, 200, NA, 400, NA, 600, 800)
  hour <- c(2, 3, 7, 8, 12, 20, 21, 22) # day = 06-18
  pp <- legacy_post(h, hour, legacy_switches())
  # -5 -> NA -> night fill; 30 -> 50 (floor); day mean = (200 + 400) / 2 = 300; night mean = (50 + 600 + 800) / 3
  expect_equal(pp$h[c(2, 3, 5, 7, 8)], c(50, 200, 400, 600, 800))
  expect_equal(pp$h[4], 300)
  expect_equal(pp$h[c(1, 6)], rep((50 + 600 + 800) / 3, 2))
  expect_equal(pp$flag_floor, c(FALSE, TRUE, rep(FALSE, 6)))
  expect_equal(pp$flag_fill, c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))
  none <- legacy_post(h, hour, legacy_switches(fill = "none", floor_m = NA))
  expect_equal(none$h, c(NA, 30, 200, NA, 400, NA, 600, 800))
  expect_false(any(none$flag_floor | none$flag_fill))
})

test_that("smoothing, gap interpolation and LOCF switches", {
  inp <- make_legacy_inp()
  base <- legacy_pblh(inp, site_l, legacy_switches(fill = "none", floor_m = NA))
  raw <- legacy_pblh(inp, site_l, legacy_switches(fill = "none", floor_m = NA, smoothing = FALSE))
  i <- 60
  expect_equal(base$h_smooth_m[i], mean(raw$h_smooth_m[(i - 2):(i + 2)])) # 5-point centred mean
  expect_true(all(is.na(base$h_smooth_m[1:2])))
  nogap <- legacy_pblh(inp, site_l, legacy_switches(fill = "none", floor_m = NA, smoothing = FALSE,
                                                    gap_interp_before_bg = FALSE))
  expect_gt(sum(is.na(nogap$pblh_m)), sum(is.na(raw$pblh_m)))
  noloc <- legacy_pblh(inp, site_l, legacy_switches(mol_locf = FALSE, fill = "none", floor_m = NA, smoothing = FALSE))
  expect_true(is.na(noloc$atmos_cond[50]) && !is.na(raw$atmos_cond[50])) # L missing there without LOCF
})

test_that("B&G Euler step: hand check of one row", {
  inp <- make_legacy_inp()
  sw <- legacy_switches(fill = "none", floor_m = NA, smoothing = FALSE, gap_interp_before_bg = FALSE)
  res <- legacy_pblh(inp, site_l, sw)
  prev_h0 <- c(NA, res$h0_m[-nrow(inp)])
  i <- which(!is.na(inp$H) & inp$H > 0 & !is.na(prev_h0) & seq_len(nrow(inp)) > 1)[2] # B&G from h0 of the previous row
  hp <- res$h0_m[i - 1]
  w <- inp$H[i] / (inp$rho[i] * inp$cp[i])
  L <- inp$MOL[i]
  rate <- (w / 0.0025) / ((hp^2 / (1.4 * hp - 2 * 2.5 * 0.41 * L)) +
                            (8 * inp$AT[i] * inp$ustar[i]^2 / (0.0025 * 9.80665 * (1.2 * hp - 2.5 * 0.41 * L))))
  expect_equal(res$h_bg_m[i], hp + 1800 * rate)
})

test_that("formula switches use the PBLHmod functions; option checks", {
  inp <- make_legacy_inp()
  sw0 <- legacy_switches(fill = "none", floor_m = NA, smoothing = FALSE, gap_interp_before_bg = FALSE)
  nie <- legacy_pblh(inp, site_l, utils::modifyList(sw0, list(stable_formula = "nieuwstadt")))
  i <- which(nie$atmos_cond %in% "nighttime-stable")[1]
  f <- 2 * (2 * pi / 86400) * sin(site_l$lat * pi / 180)
  expect_equal(nie$h0_m[i], h_nieuwstadt(inp$MOL[i], inp$ustar[i], f))
  aut <- legacy_pblh(inp, site_l, utils::modifyList(sw0, list(neutral_formula = "auto")))
  j <- which(aut$atmos_cond %in% "neutral")[1]
  expect_equal(aut$h0_m[j], h_neutral_select(inp$ustar[j], f, "auto", T_K = inp$AT[j], gamma = 0.0025)$h)
  expect_error(legacy_switches(fill = "sometimes"), "fill")
  expect_error(legacy_switches(foo = 1), "Unknown")
  expect_error(legacy_pblh(inp, site_l, legacy_switches(convective_trigger = "daily_window")), "PBLHmod value")
  expect_error(legacy_pblh(inp, site_l, legacy_switches(forcing = "pblhmod")), "forcing")
})
