# Spec 08 amendment 2: per-site override of the ETC z0 multiplier (DE-Tha), its reader, the reported values and the
# site-scoped inputs hash.

known <- c("DE-Tha", "BE-Bra", "FI-Hyy")
write_ovr <- function(x) {
  p <- tempfile(fileext = ".csv")
  readr::write_csv(x, p, na = "")
  p
}
ovr_row <- function(...) {
  utils::modifyList(list(site_id = "DE-Tha", z0_multiplier = 4.1, decided_by = "G", date = "2026-11-10",
                         reason = "test"), list(...)) |> dplyr::as_tibble()
}

test_that("read_etc_z0_multiplier: the committed file and the validation", {
  x <- read_etc_z0_multiplier(here::here("config", "etc_z0_multiplier.csv"),
                              read_sites(here::here("config", "sites.csv"))$site_id)
  expect_equal(x$site_id, "DE-Tha")
  expect_equal(x$z0_multiplier, 4.1)
  expect_equal(x$decided_by, "G")
  expect_equal(format(x$date), "2026-11-10")
  expect_equal(nrow(read_etc_z0_multiplier(write_ovr(ovr_row()), known)), 1)
  expect_error(read_etc_z0_multiplier(write_ovr(ovr_row(reason = NA)), known), "required")
  expect_error(read_etc_z0_multiplier(write_ovr(ovr_row(decided_by = NA)), known), "required")
  expect_error(suppressWarnings(read_etc_z0_multiplier(write_ovr(dplyr::select(ovr_row(), -date)), known)),
               "missing[[:space:]]+column")
  expect_error(read_etc_z0_multiplier(write_ovr(ovr_row(site_id = "XX-Nop")), known), "unknown[[:space:]]+site")
  expect_error(read_etc_z0_multiplier(write_ovr(ovr_row(z0_multiplier = 1.9)), known), "between[[:space:]]+2[[:space:]]+and[[:space:]]+5")
  expect_error(read_etc_z0_multiplier(write_ovr(ovr_row(z0_multiplier = 5.1)), known), "between[[:space:]]+2[[:space:]]+and[[:space:]]+5")
  expect_error(read_etc_z0_multiplier(write_ovr(dplyr::bind_rows(ovr_row(), ovr_row(z0_multiplier = 3))), known),
               "more[[:space:]]+than[[:space:]]+one[[:space:]]+row")
})

# Two height periods at EC height 42 m: hc 27.74 m (zm/z0 = 5.65) then 32.59 m (zm/z0 = 4.15)
two_periods <- function(n1, n2) {
  hz <- footprint_heights(rep(42, n1 + n2), c(rep(27.74, n1), rep(32.59, n2)))
  list(zm = hz$zm_m, z0 = hz$z0_m, second = c(rep(FALSE, n1), rep(TRUE, n2)))
}
e_rsl <- function(p, m) {
  v <- footprint_validity(p$zm, p$z0, 42, 1000, 3, -100, 0.5, 0.4, 180, m)
  v$e_rsl
}

test_that("two height periods: the loop rejects the whole second period, the override keeps it", {
  p <- two_periods(70, 30)
  expect_equal(round(p$zm / p$z0, 2)[c(1, 100)], c(5.65, 4.15))
  # second period <= 1/3 of the rows: rejected <= 0.5 x retained already at m = 5
  m <- etc_z0_multiplier_site(p$zm, p$z0, "DE-Tha", NULL)
  expect_equal(m$m_loop, 5)
  expect_equal(m$m_used, 5)
  expect_equal(m$m_source, "loop")
  expect_true(all(!e_rsl(p, m$m_used)[p$second]))
  expect_true(all(e_rsl(p, m$m_used)[!p$second]))
  # override 4.1: no row fails the roughness test; the loop value is still reported
  ov <- read_etc_z0_multiplier(write_ovr(ovr_row()), known)
  mo <- etc_z0_multiplier_site(p$zm, p$z0, "DE-Tha", ov)
  expect_equal(c(mo$m_loop, mo$m_used), c(5, 4.1))
  expect_equal(mo$m_source, "override")
  expect_equal(mo$override$reason, "test")
  expect_true(all(e_rsl(p, mo$m_used)))
  # 60 / 40 % (the split named in the amendment): the loop itself goes down to 4.1 (rejected > 0.5 x retained
  # until no row is rejected), so the legacy rule alone would keep the second period
  expect_equal(etc_z0_multiplier(two_periods(60, 40)$zm, two_periods(60, 40)$z0), 4.1)
})

test_that("a site without an override row: same m and flags as the loop", {
  p <- two_periods(70, 30)
  ov <- read_etc_z0_multiplier(write_ovr(ovr_row()), known)
  m <- etc_z0_multiplier_site(p$zm, p$z0, "BE-Bra", ov)
  expect_equal(m$m_source, "loop")
  expect_null(m$override)
  expect_equal(m$m_used, etc_z0_multiplier(p$zm, p$z0))
  expect_identical(e_rsl(p, m$m_used), e_rsl(p, etc_z0_multiplier(p$zm, p$z0)))
  # regression on the PC's spec-08 input tables (not fixtures; skipped elsewhere): BE-Bra m 5, FI-Hyy m 3.7
  for (s in c("BE-Bra", "FI-Hyy")) {
    f <- here::here("data", "derived", "footprint", s, paste0(s, "_footprint_inputs.rds"))
    if (!file.exists(f)) next
    x <- readRDS(f)
    ms <- etc_z0_multiplier_site(x$tab$zm_m, x$tab$z0_m, s, ov)
    expect_equal(ms$m_used, x$m_etc)
    expect_equal(ms$m_used, c(`BE-Bra` = 5, `FI-Hyy` = 3.7)[[s]])
  }
})

test_that("network_etc_multiplier: new and pre-amendment tables, override visible", {
  new <- dplyr::tibble(site_id = "DE-Tha", m_etc = 4.1, m_loop = 4.6, m_used = 4.1, m_source = "override",
                       n_rows = 10, zm_over_z0_median = 4.9, zm_over_z0_min = 4.15, n_rejected = 0,
                       n_rejected_loop = 3, override_reason = "r")
  old <- dplyr::tibble(site_id = "BE-Bra", m_etc = 5, n_rows = 10, zm_over_z0_median = 12, n_rejected = 0)
  x <- network_etc_multiplier(dplyr::bind_rows(new, old),
                              dplyr::tibble(site_id = c("BE-Bra", "DE-Tha"), year = c(2025L, 2025L)))
  expect_equal(x$m_loop, c(4.6, 5))
  expect_equal(x$m_used, c(4.1, 5))
  expect_equal(x$m_source, c("override", "loop"))
  expect_equal(x$n_rejected_loop, c(3, 0))
  expect_equal(x$reference_year, c(2025L, 2025L))
})

test_that("inputs hash: only the site's own override row enters the 08a hash", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ovr.csv")
  cfg <- list(derived_dir = file.path(tmp, "derived"), output_dir = file.path(tmp, "out"), etc_z0_multiplier_file = f)
  inp <- dplyr::tibble(path = character(), size = numeric(), mtime = character())
  h <- function(site) step_inputs_hash("08a", site, cfg, c(`04` = "x"), input_files = inp)
  readr::write_csv(ovr_row(), f)
  h_bra <- h("BE-Bra")
  h_tha <- h("DE-Tha")
  readr::write_csv(ovr_row(z0_multiplier = 4.2), f)
  expect_equal(h("BE-Bra"), h_bra)
  expect_false(identical(h("DE-Tha"), h_tha))
  readr::write_csv(ovr_row(reason = "other wording"), f)
  expect_equal(h("BE-Bra"), h_bra)
  expect_false(identical(h("DE-Tha"), h_tha))
  # a site without a row: the hash is the one computed without the file (as before the amendment)
  unlink(f)
  expect_equal(h("BE-Bra"), h_bra)
  expect_null(site_config_row(f, "DE-Tha"))
  expect_null(site_config_row(NULL, "DE-Tha"))
  # step 08: the reference year enters only as the site's --year argument
  a <- step_inputs_hash("08", "BE-Bra", cfg, c(`08a` = "x"), args = "--year=2025", input_files = inp)
  expect_false(identical(a, step_inputs_hash("08", "BE-Bra", cfg, c(`08a` = "x"), args = "--year=2023",
                                             input_files = inp)))
})
