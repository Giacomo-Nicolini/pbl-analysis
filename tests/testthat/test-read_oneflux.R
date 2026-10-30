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

test_that("find_oneflux_bif: default folder, glob over run folders, product match, newest run (spec 09 A1)", {
  # default: next to the L2 zips
  tmp <- withr::local_tempdir()
  d0 <- file.path(tmp, "XX-Tst_output_L2I_2026")
  dir.create(d0)
  file.copy(bif(), d0)
  hh <- "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2021_v1.3_r1.zip"
  expect_equal(basename(find_oneflux_bif("XX-Tst", list(icos_l2_dir = tmp), hh)), basename(bif()))
  expect_equal(basename(find_oneflux_bif("XX-Tst", list(icos_l2_dir = tmp, oneflux_bif_dir = "{icos_l2_dir}/{site}_output_L2I_2026"))),
               basename(bif()))
  expect_error(find_oneflux_bif("XX-Tst", list(icos_l2_dir = tmp), "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2022_v1.3_r1.zip"),
               "no ONEFlux BIF file matches")

  # cluster layout: one folder per site and processing run, glob in the configured path
  runs <- withr::local_tempdir()
  b <- readr::read_csv(bif(), col_types = readr::cols(.default = "c"))
  put <- function(folder, product, date) {
    x <- b
    x$DATAVALUE[x$VARIABLE == "PRODUCT_NAME"] <- product
    if (!is.na(date)) x <- rbind(x, data.frame(SITE_ID = "XX-Tst", GROUP_ID = "1", VARIABLE_GROUP = "GRP_ONEFLUX",
                                               VARIABLE = "PRODUCT_PROCESSING_DATE", DATAVALUE = date))
    dir.create(file.path(runs, folder, "99_fluxnet2015"), recursive = TRUE)
    readr::write_csv(x, file.path(runs, folder, "99_fluxnet2015", "ICOSETC_XX-Tst_FLUXNET_BIF_2020-2021_v1.3_r1.csv"))
  }
  put("XX-Tst__ICOS__ex20260101000000_proc202601020000", "ICOSETC_XX-Tst_FLUXNET_2020-2021_v1.3_r1.zip", "20260102")
  put("XX-Tst__ICOS__ex20260901000000_proc202609020000", "ICOSETC_XX-Tst_FLUXNET_2020-2021_v1.3_r1.zip", "20260902")
  put("XX-Tst__ICOS__ex20261001000000_proc202610020000", "ICOSETC_XX-Tst_FLUXNET_2020-2022_v1.3_r1.zip", "20261002")
  cfg <- list(icos_l2_dir = tmp, oneflux_bif_dir = file.path(runs, "{site}__ICOS__ex*_proc*", "99_fluxnet2015"))
  expect_message(p <- find_oneflux_bif("XX-Tst", cfg, hh), "newest processing run")
  expect_match(p, "proc202609020000") # 2026-10 run is another product (2020-2022)
  expect_match(find_oneflux_bif("XX-Tst", cfg, "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2022_v1.3_r1.zip"), "proc202610020000")
  expect_error(find_oneflux_bif("XX-Tst", cfg), "no HH file")
  # same processing date: the folder stamp decides; identical runs abort
  put("XX-Tst__ICOS__ex20260901000000_proc202609021200", "ICOSETC_XX-Tst_FLUXNET_2020-2021_v1.3_r1.zip", "20260902")
  expect_match(suppressMessages(find_oneflux_bif("XX-Tst", cfg, hh)), "proc202609021200")
  put("XX-Tst__ICOS__ex20260901000001_proc202609021200", "ICOSETC_XX-Tst_FLUXNET_2020-2021_v1.3_r1.zip", "20260902")
  expect_error(find_oneflux_bif("XX-Tst", cfg, hh), "same processing run")
  expect_error(find_oneflux_bif("XX-Tst", list(icos_l2_dir = tmp, oneflux_bif_dir = file.path(runs, "nope*"))), "no ONEFlux BIF")
})
