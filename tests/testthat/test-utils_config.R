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
