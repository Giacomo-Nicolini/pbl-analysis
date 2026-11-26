test_that("growth display data are unpaired, break at missing hours and keep sparse measured values as points", {
  # 12 clear days in July at lon 0 (apparent solar time ~ UTC); half-hourly rows 00:00-23:30
  t <- as.POSIXct("2023-07-01", tz = "UTC") + rep(0:11, each = 48) * 86400 + rep(0:47, 12) * 1800
  g <- dplyr::tibble(timestamp_utc = t, model = 500, measured = NA_real_)
  h <- as.numeric(format(t, "%H"))
  g$model[h == 14] <- NA_real_                       # missing hour of the model: the line breaks there
  g$measured[h == 10] <- 800                         # 24 values: median
  g$measured[h == 11 & format(t, "%M") == "00" & as.numeric(format(t, "%d")) <= 6] <- 900 # 6 values: points
  kind <- ifelse(is.na(g$measured), NA_character_, "exact")
  kind[h == 12] <- rep(c("below", "lower_bound", "exact", "exact"), length.out = sum(h == 12))
  dd <- growth_display_data(g, c("measured", "model"), "measured", lon = 0, mid_s = 900, n_min = 10, kind = kind)
  ln <- dd$lines[dd$lines$season == "JJA", ]
  expect_equal(nrow(ln), 48) # complete hour grid x 2 variables
  m <- ln[ln$variable == "model", ]
  expect_true(is.na(m$median[m$ast_hour == 14]))
  expect_equal(m$median[m$ast_hour == 13], 500)
  expect_equal(m$n[m$ast_hour == 13], 24) # all clear-day half-hours, independent of the measured values
  me <- ln[ln$variable == "measured", ]
  expect_equal(me$median[me$ast_hour == 10], 800)
  expect_true(is.na(me$median[me$ast_hour == 11]))
  expect_equal(me$n[me$ast_hour == 11], 6)
  expect_equal(nrow(dd$points), 6)
  expect_true(all(dd$points$ast_hour == 11 & dd$points$pblh_m == 900))
  s12 <- dd$strip[dd$strip$season == "JJA" & dd$strip$ast_hour == 12, ]
  expect_equal(s12$share_below, 0.25)
  expect_equal(s12$share_lower_bound, 0.25)
  expect_equal(s12$n_usable, 24)
  expect_null(growth_display_data(g, c("measured", "model"), "measured", 0, 900)$strip)

  pal <- c(measured = unname(style_of("cloudnet_exact")), model = unname(style_of("pblhmod_default")))
  labs <- c(measured = "measured", model = "model")
  p <- plot_growth_clearday(dd, pal, labs, "measured", day_labels = dplyr::tibble(season = "JJA", label = "12 clear days"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  dd0 <- growth_display_data(g, c("measured", "model"), "measured", 0, 900)
  expect_s3_class(ggplot2::ggplot_build(plot_growth_clearday(dd0, pal, labs, "measured", window = c(9, 17))),
                  "ggplot_built")
})
