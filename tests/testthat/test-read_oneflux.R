bif <- function() fixture_path("oneflux/ICOSETC_XX-Tst_FLUXNET_BIF_2020-2021_v1.3_r1.csv")

test_that("read_oneflux_bif returns groups as wide tables and checks the product", {
  fr <- read_oneflux_bif(bif(), "GRP_FLUX_REF")
  expect_equal(nrow(fr), 4)
  expect_equal(fr$REF_VERSION, c("CUT", "VUT", "VUT", "VUT"))
  expect_true(is.na(fr$REF_YEAR[1]))
  expect_equal(attr(fr, "product_name"), "ICOSETC_XX-Tst_FLUXNET_2020-2021_v1.3_r1.zip")
  all <- read_oneflux_bif(bif())
  expect_true(all(c("GRP_FLUX_REF", "GRP_ONEFLUX", "GRP_UST_THR") %in% names(all)))
  expect_equal(all$GRP_UST_THR$UST_THR_VALUE, "0.21")
  # product check against the HH file name
  expect_silent(read_oneflux_bif(bif(), "GRP_FLUX_REF", hh_path = "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2021_v1.3_r1.zip"))
  expect_error(read_oneflux_bif(bif(), "GRP_FLUX_REF", hh_path = "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2020_v1.3_r1.zip"),
               "does not match")
  expect_error(read_oneflux_bif(bif(), "GRP_NOPE"), "no group")
})

test_that("oneflux_ustar_ref picks the reference NEE thresholds", {
  fr <- read_oneflux_bif(bif(), "GRP_FLUX_REF")
  v <- oneflux_ustar_ref(fr, "VUT")
  expect_equal(v$year, c(2020L, 2021L))
  expect_equal(v$ustar_threshold, c(0.2, 0.3)) # the DD threshold (0.9) is not used
  cc <- oneflux_ustar_ref(fr, "CUT")
  expect_equal(cc$ustar_threshold, 0.25)
  expect_true(is.na(cc$year))
  expect_equal(oneflux_ustar_filtered(c(0.1, 0.25, 0.25, NA), c(2020L, 2020L, 2021L, 2021L), v), c(TRUE, FALSE, TRUE, NA))
  expect_equal(oneflux_ustar_filtered(c(0.1, 0.3), c(2020L, 2030L), cc), c(TRUE, FALSE))
  expect_true(is.na(oneflux_ustar_filtered(0.1, 2030L, v)))
})
