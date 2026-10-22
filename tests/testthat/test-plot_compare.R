utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("figure functions build; save_figure writes png by default, pdf on request", {
  set.seed(1)
  t <- utc("2020-01-01 00:00:00") + 3600 * 0:(24 * 120 - 1)
  pt <- pair_table(
    dplyr::tibble(timestamp_utc = t, pblh_m = exp(stats::rnorm(length(t), 6, 0.5)),
                  x_regime = sample(c("stable", "convective"), length(t), TRUE)),
    dplyr::tibble(timestamp_utc = t, pblh_m = exp(stats::rnorm(length(t), 6, 0.5)), x_blh_cv9 = 0.1),
    "p"
  )
  pt$cloud <- cloud_class(sample(c(NA, 500, 3000), nrow(pt), TRUE))
  wide <- dplyr::tibble(season = pt$season, hour_utc = pt$hour_utc, pblhmod_default = pt$mod,
                        era5_blh_center = pt$ref, icos_oper_qc0 = pt$mod * 2, icos_oper_all = pt$mod * 2)
  ds <- diurnal_summary(wide, c("pblhmod_default", "era5_blh_center", "icos_oper_qc0", "icos_oper_all"))
  era <- dplyr::tibble(timestamp_start_utc = t, pblh_m = pt$ref, x_blh_cv9 = stats::runif(length(t), 0, 0.3))
  hd <- dplyr::tibble(timestamp_utc = t, H_ec_Wm2 = stats::rnorm(length(t), 50, 50),
                      H_era5_Wm2 = stats::rnorm(length(t), 40, 40), season = season_of(t))
  gm <- dplyr::tibble(hour_utc = rep(c(0L, 12L), each = 12), month = rep(1:12, 2), median = 0.005,
                      p25 = 0.004, p75 = 0.006, n = 30L)

  figs <- list(
    plot_diurnal_by_season(ds), plot_scatter_pair(pt), plot_logratio_by_hour(pt),
    plot_logratio_by_cloud(pt), plot_gamma_monthly(gm)
  )
  for (p in figs) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  two <- plot_era5_spread(era, era)
  expect_length(two, 2)
  expect_length(plot_h_ec_vs_era5(hd), 2)

  dir <- withr::local_tempdir()
  paths <- save_figure(two, "test_fig", 16, 8, dir = dir)
  expect_equal(names(paths), "png") # default: png only
  expect_true(file.exists(paths[["png"]]))
  expect_false(file.exists(file.path(dir, "test_fig.pdf")))
  expect_gt(file.size(paths[["png"]]), 1000)
  both <- save_figure(two, "test_fig2", 16, 8, dir = dir, formats = c("png", "pdf")) # paper figures
  expect_equal(names(both), c("png", "pdf"))
  expect_true(all(file.exists(both)))
  expect_error(save_figure(two, "x", dir = dir, formats = "svg"), "formats")
  expect_equal(unname(names(source_palette())), unname(names(source_labels())))
})
