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

test_that("config.yml has the ERA5 directories in both profiles; cluster paths are LOCAL", {
  for (prof in c("default", "cluster")) {
    cfg <- get_cfg(prof, local_file = NULL)
    expect_true(all(c("era5_sl_dir", "era5_pl_dir", "era5_static_dir") %in% names(cfg)), info = prof)
    expect_null(cfg$era5_dir)
  }
  cl <- get_cfg("cluster", local_file = NULL)
  for (k in c("icos_l2_dir", "era5_sl_dir", "era5_pl_dir", "era5_static_dir")) expect_equal(cl[[k]], "LOCAL", info = k)
  expect_equal(cl$sites_file, "config/sites.csv")
})

# Write a local override file in a fresh temporary folder (distinct path = reported again).
local_yml <- function(lines, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  p <- file.path(dir, "config.local.yml")
  writeLines(lines, p)
  p
}

test_that("get_cfg: local file overrides base; profile section overrides local default", {
  p <- local_yml(c(
    "default:",
    "  output_dir: 'out_local'",
    "  test_sites: ['BE-Bra', 'DE-Tha']",
    "cluster:",
    "  output_dir: 'out_cluster'",
    "  icos_l2_dir: '/x/l2'",
    "  era5_static_dir: '/x/static/{var}'"
  ))
  expect_message(d <- get_cfg("default", local_file = p), "2 keys from")
  expect_equal(d$output_dir, "out_local")
  expect_equal(unlist(d$test_sites), c("BE-Bra", "DE-Tha")) # list replaced, not merged
  expect_equal(d$icos_l2_dir, "data/raw/icos_l2")           # cluster section not applied
  c1 <- suppressMessages(get_cfg("cluster", local_file = p))
  expect_equal(c1$output_dir, "out_cluster")                 # profile over local default
  expect_equal(c1$icos_l2_dir, "/x/l2")
  expect_equal(unlist(c1$test_sites), c("BE-Bra", "DE-Tha")) # local default applies to all profiles
  expect_equal(c1$sites_file, "config/sites.csv")            # untouched base value
  expect_match(resolve_path(c1$era5_static_dir, var = "z", must_exist = FALSE), "/x/static/z$")
  # reported once per session for the same file and profile
  expect_no_message(get_cfg("cluster", local_file = p))
})

test_that("get_cfg: profile from R_CONFIG_ACTIVE; missing local file = base only", {
  p <- local_yml(c("cluster:", "  icos_l2_dir: '/y/l2'"))
  withr::local_envvar(R_CONFIG_ACTIVE = "cluster")
  expect_equal(suppressMessages(get_cfg(local_file = p))$icos_l2_dir, "/y/l2")
  expect_equal(get_cfg(local_file = file.path(tempdir(), "no_such_local.yml"))$icos_l2_dir, "LOCAL")
  expect_equal(get_cfg("default", local_file = p)$icos_l2_dir, "data/raw/icos_l2") # no default section
})

test_that("get_cfg: unknown keys in the local file are kept and reported", {
  p <- local_yml(c("default:", "  my_extra_dir: '/z'"))
  expect_message(cfg <- get_cfg("default", local_file = p), "Not in")
  expect_equal(cfg$my_extra_dir, "/z")
})

test_that("resolve_path aborts on LOCAL", {
  expect_error(resolve_path("LOCAL"), "config.local.yml")
})

test_that("config.local.yml.example parses, covers the LOCAL keys, holds placeholders only", {
  ex <- yaml::read_yaml(here::here("config.local.yml.example"))
  base <- get_cfg("cluster", local_file = NULL)
  local_keys <- names(base)[vapply(base, function(v) identical(v, "LOCAL"), logical(1))]
  expect_true(all(local_keys %in% names(ex$cluster)))
  vals <- unlist(ex$cluster)
  expect_true(all(startsWith(vals, "/path/to/")), info = paste(vals, collapse = ", "))
})
