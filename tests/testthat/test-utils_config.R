test_that("resolve_path substitutes {site} and keeps absolute paths", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "BE-Bra"))
  p <- resolve_path(file.path(tmp, "{site}"), site = "BE-Bra")
  expect_equal(normalizePath(p), normalizePath(file.path(tmp, "BE-Bra")))
})

test_that("resolve_path resolves relative paths against the project root", {
  expect_equal(resolve_path("config", must_exist = TRUE), here::here("config"))
})

test_that("resolve_path aborts on TODO, missing site, and non-existent paths", {
  expect_error(resolve_path("TODO"), "TODO")
  expect_error(resolve_path("data/{site}", must_exist = FALSE), "site")
  expect_error(resolve_path("no/such/dir_xyz"), "does not exist")
  expect_no_error(resolve_path("no/such/dir_xyz", must_exist = FALSE))
})

test_that("resolve_path substitutes {var} together with {site}", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "BE-Bra", "z"), recursive = TRUE)
  p <- resolve_path(file.path(tmp, "{site}", "{var}"), site = "BE-Bra", var = "z")
  expect_equal(normalizePath(p), normalizePath(file.path(tmp, "BE-Bra", "z")))
  # var ignored when the path has no {var} (local static folder)
  expect_equal(resolve_path(tmp, var = "z"), tmp)
})

test_that("resolve_path aborts on missing var and unresolved placeholders", {
  expect_error(resolve_path("data/{var}", must_exist = FALSE), "var")
  expect_error(resolve_path("data/{site}/{foo}", site = "BE-Bra", must_exist = FALSE), "unresolved")
})

test_that("config.yml has the ERA5 directories in both profiles", {
  for (prof in c("default", "cluster")) {
    cfg <- get_cfg(prof)
    expect_true(all(c("era5_sl_dir", "era5_pl_dir", "era5_static_dir") %in% names(cfg)), info = prof)
    expect_null(cfg$era5_dir)
  }
  expect_match(get_cfg("cluster")$era5_static_dir, "{var}", fixed = TRUE)
})
