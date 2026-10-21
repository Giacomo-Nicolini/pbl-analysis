site <- test_site_meta()
run <- function(d, ...) run_pblh_model(d, site, pblh_options(...))

# start / end index of the convective rows of each day
windows <- function(out) {
  conv <- out$regime %in% "convective"
  split(which(conv), out$day_id[conv])
}

test_that("3-day synthetic series: daily cycle, start at sunrise, no carry-over", {
  d <- make_diurnal(3)
  out <- run(d)
  expect_equal(nrow(out), nrow(d))
  w <- windows(out)
  expect_length(w, 3)

  for (idx in w) {
    s <- idx[1]
    # integration starts at the first daylight half-hour
    expect_true(out$is_day[s])
    expect_false(out$is_day[s - 1])
    # h_bg non-decreasing inside the window (Amendment 7)
    expect_true(all(diff(out$h_bg_m[idx]) >= -1e-9))
    # no carry-over: first window value starts from the diagnostic h0 of the half-hour before
    f <- prepare_pblh_forcing(d, site, pblh_options())
    st <- bg_step(out$h_diag_m[s - 1], f$wtheta[s], f$ustar[s], f$L[s], f$T_K[s], 0.01)
    expect_equal(out$h_bg_m[s], st$h_mean)
    # H > 0 all day: window ends at sunset (last daylight half-hour)
    e <- idx[length(idx)]
    expect_true(out$is_day[e])
    expect_false(out$is_day[e + 1])
    expect_equal(length(idx), e - s + 1)
  }
  # daily maximum well above the morning value; day 2 restarts low
  expect_gt(max(out$h_bg_m[w[[1]]]), 3 * out$h_bg_m[w[[2]][1]])
  # night rows: diagnostic, stable regime (H < 0 -> L > 0)
  night <- !out$is_day
  expect_true(all(out$regime[night] == "stable"))
  expect_true(all(out$method[night] == "B1"))
  expect_equal(attr(out, "pblh_options")$gamma, 0.01)
})

test_that("combine = max vs bg_only, and method labels", {
  d <- make_diurnal(1)
  out <- run(d, combine = "max", neutral_method = "B2")
  conv <- out$regime %in% "convective"
  expect_equal(out$pblh_m[conv], pmax(out$h_bg_m[conv], out$h_diag_m[conv]))
  expect_equal(out$method[conv], ifelse(out$h_diag_m[conv] > out$h_bg_m[conv], "B2", "B5"))
  out_auto <- run(d, combine = "max")
  expect_true(all(out_auto$method[conv] %in% c("B3", "B5"))) # N/|f| = 184 > 70 here
  out2 <- run(d) # default: bg_only
  expect_equal(attr(out2, "pblh_options")$combine, "bg_only")
  expect_equal(out2$pblh_m[conv], out2$h_bg_m[conv])
  expect_true(all(out2$method[conv] == "B5"))
})

test_that("window ends at the first run of n_neg_end non-positive half-hours", {
  d <- make_diurnal(1)
  day <- which(is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon))
  mid <- day[floor(length(day) / 2)]
  # single negative half-hour at midday: short spell, held, window continues
  d1 <- d
  d1$fn_H_F_MDS[mid] <- -10
  d1$MO_LENGTH[mid] <- 50
  out <- run(d1)
  expect_equal(out$regime[mid], "convective")
  expect_equal(out$regime[mid + 1], "convective")
  # held: the held value is the state the next half-hour integrates from
  f <- prepare_pblh_forcing(d1, site, pblh_options())
  expect_gt(out$h_bg_m[mid], out$h_bg_m[mid - 1]) # end-of-interval state > interval mean
  st <- bg_step(out$h_bg_m[mid], f$wtheta[mid + 1], f$ustar[mid + 1], f$L[mid + 1], f$T_K[mid + 1], 0.01)
  expect_equal(out$h_bg_m[mid + 1], st$h_mean)
  # two negative half-hours: window ends at the first of them
  d2 <- d
  d2$fn_H_F_MDS[mid:(mid + 1)] <- -10
  d2$MO_LENGTH[mid:(mid + 1)] <- 50
  out <- run(d2)
  expect_equal(out$regime[mid - 1], "convective")
  expect_equal(out$regime[mid], "stable")
  expect_false(any(out$regime[mid:max(day)] == "convective"))
  # n_neg_end = 3: the same two half-hours are held
  out <- run(d2, n_neg_end = 3)
  expect_equal(out$regime[mid], "convective")
})

test_that("held h before the first positive heat flux", {
  d <- make_diurnal(1)
  day <- which(is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon))
  d$fn_H_F_MDS[day[1:2]] <- -5
  d$MO_LENGTH[day[1:2]] <- 100
  out <- run(d)
  h0 <- out$h_diag_m[day[1] - 1]
  expect_equal(out$h_bg_m[day[1:2]], c(h0, h0))
  expect_equal(out$regime[day[1:2]], c("convective", "convective"))
  expect_false(any(out$flag_sign_inconsistent))
})

test_that("gap rules: short gap interpolated, long gap stops the day", {
  d <- make_diurnal(2)
  day1 <- which(is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon) &
                  seq_len(nrow(d)) <= 48)
  g <- day1[8]
  d_short <- d
  d_short$USTAR[g:(g + 1)] <- NA
  out <- run(d_short)
  expect_true(all(out$flag_ustar_interp[g:(g + 1)]))
  expect_true(all(out$regime[g:(g + 1)] == "convective"))
  expect_false(any(out$flag_bg_stopped_gap))

  d_long <- d
  d_long$USTAR[g:(g + 2)] <- NA
  out <- run(d_long)
  expect_false(any(out$flag_ustar_interp[g:(g + 2)]))
  expect_true(all(out$flag_bg_stopped_gap[out$day_id == out$day_id[g]]))
  expect_equal(out$regime[g - 1], "convective")
  expect_false(any(out$regime[g:max(day1)] %in% "convective"))
  # next day unaffected
  expect_false(any(out$flag_bg_stopped_gap[out$day_id != out$day_id[g]]))
  expect_length(windows(out), 2)
})

test_that("no h0: no integration that day", {
  d <- make_diurnal(2)
  dd <- is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon)
  s2 <- which(dd & c(TRUE, !dd[-length(dd)]))[2]
  d$MO_LENGTH[(s2 - 5):(s2 - 1)] <- NA # 5 > max_gap_interp: not filled
  out <- run(d)
  day2 <- out$day_id == out$day_id[s2]
  expect_true(all(out$flag_no_h0[day2]))
  expect_false(any(out$regime[day2] %in% "convective"))
  expect_false(any(out$flag_no_h0[!day2]))
})

test_that("h0 from the nearest valid value within 2 h", {
  d <- make_diurnal(2)
  dd <- is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon)
  s2 <- which(dd & c(TRUE, !dd[-length(dd)]))[2]
  d$MO_LENGTH[(s2 - 3):(s2 - 1)] <- NA
  out <- run(d, max_gap_interp = 0)
  expect_false(out$flag_no_h0[s2])
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0))
  st <- bg_step(out$h_diag_m[s2 - 4], f$wtheta[s2], f$ustar[s2], f$L[s2], f$T_K[s2], 0.01)
  expect_equal(out$h_bg_m[s2], st$h_mean)
})

test_that("no positive heat flux: no window, flagged", {
  d <- make_diurnal(2)
  day2 <- 49:96
  d$fn_H_F_MDS[day2] <- -10
  d$MO_LENGTH[day2] <- 80
  out <- run(d)
  dd2 <- out$day_id == out$day_id[which(out$is_day & seq_len(nrow(out)) > 48)[1]]
  expect_true(all(out$flag_no_positive_h[dd2]))
  expect_false(any(out$regime[dd2] %in% "convective"))
  expect_length(windows(out), 1)
})

test_that("missing heat flux all day is a gap day, not a no-positive day", {
  d <- make_diurnal(2)
  d$fn_H_F_MDS[49:96] <- NA
  d$H[49:96] <- NA
  out <- run(d)
  day2 <- out$day_id == out$day_id[which(out$is_day & seq_len(nrow(out)) > 48)[1]]
  expect_false(any(out$flag_no_positive_h))
  expect_true(all(out$flag_bg_stopped_gap[day2]))
  expect_false(any(out$regime[day2] %in% "convective"))
})

test_that("sign inconsistency: h held and flagged", {
  d <- make_diurnal(1)
  day <- which(is_day(d$timestamp_start_utc, d$timestamp_end_utc, site$lat, site$lon))
  i <- day[10]
  d$MO_LENGTH[i] <- 30 # H > 0 with L > 0
  out <- run(d)
  expect_true(out$flag_sign_inconsistent[i])
  expect_equal(sum(out$flag_sign_inconsistent), 1)
  expect_equal(out$regime[i], "convective")
  # held: h_bg[i] is the carried state, and interval i + 1 integrates from it
  f <- prepare_pblh_forcing(d, site, pblh_options())
  st <- bg_step(out$h_bg_m[i], f$wtheta[i + 1], f$ustar[i + 1], f$L[i + 1], f$T_K[i + 1], 0.01)
  expect_equal(out$h_bg_m[i + 1], st$h_mean)
  expect_gt(out$h_bg_m[i], out$h_bg_m[i - 1])
})

test_that("excluded rows give NA and regime excluded", {
  d <- make_diurnal(1)
  d$excl_flag[20:22] <- TRUE
  out <- run(d)
  expect_true(all(is.na(out$pblh_m[20:22])))
  expect_true(all(out$regime[20:22] == "excluded"))
  expect_true(all(out$excl_flag[20:22]))
})

test_that("polar day and polar night: diagnostic only, day_id = local date", {
  polar_site <- test_site_meta(lat = 78.2, lon = 15.6, utc_offset_h = 1)
  d <- make_diurnal(2, start = "2024-06-20 23:00:00", site = polar_site)
  out <- run_pblh_model(d, polar_site, pblh_options())
  expect_true(all(out$flag_polar_day))
  expect_false(any(out$regime %in% "convective"))
  expect_equal(out$day_id[1], as.Date("2024-06-21"))

  d <- make_diurnal(2, start = "2024-12-20 23:00:00", site = polar_site)
  out <- run_pblh_model(d, polar_site, pblh_options())
  expect_false(any(out$is_day))
  expect_false(any(out$regime %in% "convective"))
  expect_equal(unique(out$day_id), as.Date(c("2024-12-21", "2024-12-22")))
})

test_that("gamma: vector interface and sensitivity", {
  d <- make_diurnal(1)
  a <- run(d)
  b <- run(d, gamma = rep(0.01, nrow(d)))
  expect_equal(a$pblh_m, b$pblh_m)
  c <- run(d, gamma = 0.0025)
  conv <- a$regime %in% "convective"
  expect_gt(max(c$h_bg_m[conv]), max(a$h_bg_m[conv])) # weaker inversion -> deeper CBL
  expect_error(run(d, gamma = c(0.01, 0.02)), "one value per row")
})

test_that("flag_below_zm and output columns", {
  d <- make_diurnal(1)
  d$zm_m <- 1000
  out <- run(d)
  expect_equal(out$flag_below_zm, dplyr::coalesce(out$pblh_m < 1000, FALSE))
  expect_true(any(out$flag_below_zm))
  expect_true(all(c("site_id", "timestamp_start_utc", "timestamp_end_utc", "pblh_m", "regime",
    "method", "h_diag_m", "h_bg_m", "day_id", "is_day", "ustar", "L", "wtheta", "T_K", "p_Pa",
    "H_source", "T_source", "p_source", "zm_m", "flag_ustar_interp", "flag_L_interp",
    "flag_H_interp", "flag_H_qc_masked", "flag_no_h0", "flag_bg_stopped_gap",
    "flag_sign_inconsistent", "flag_no_positive_h", "flag_polar_day", "flag_below_zm",
    "excl_flag") %in% names(out)))
  expect_s3_class(out$day_id, "Date")
})

test_that("non-regular input aborts", {
  d <- make_diurnal(1)[-5, ]
  expect_error(run(d), "30-min grid")
})

test_that("u* filter: Inf/Inf leaves the forcing raw; rejected rows change only through u* and L", {
  d <- make_diurnal(3)
  d$WS <- 3
  bad <- c(30, 31, 80) # midday and night rows with implausible u*
  d$USTAR[bad] <- 2.5  # u*/WS = 0.83
  off <- run_pblh_model(d, site, pblh_options(ustar_ws_max = Inf, ustar_max = Inf))
  expect_equal(off$ustar, d$USTAR)
  expect_false(any(off$flag_ustar_implausible))
  # without the WS column only the cap applies (2.5 < 3): identical to the disabled filter
  same <- run_pblh_model(dplyr::select(d, -WS), site)
  expect_identical(same$pblh_m, off$pblh_m)

  on <- run_pblh_model(d, site)
  expect_equal(which(on$flag_ustar_implausible), bad)
  expect_equal(on$ustar[80], 0.3) # single row, interpolated
  expect_equal(on$WS, rep(3, nrow(d)))
  expect_true(all(c("WS", "flag_ustar_implausible") %in% names(on)))
  unaffected <- setdiff(seq_len(nrow(d)), c(bad, 29:33, 79:81))
  expect_equal(on$pblh_m[unaffected][1:20], off$pblh_m[unaffected][1:20])
})
