# Spec 09 A2: per-site driver, inputs hash and run-status skip logic (fake runner, temporary folders).

local_pipeline_cfg <- function(env = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = env)
  l2 <- file.path(tmp, "l2", "XX-Tst_output_L2I_2026")
  pl <- file.path(tmp, "era5_pl", "XX-Tst")
  dir.create(l2, recursive = TRUE)
  dir.create(pl, recursive = TRUE)
  writeLines("a", file.path(l2, "ICOSETC_XX-Tst_FLUXES_INTERIM_L2.zip"))
  writeLines("p", file.path(pl, "XX-Tst__reanalysis_era5_pressure_levels-pl__temperature__x.csv"))
  list(icos_l2_dir = file.path(tmp, "l2"), era5_pl_dir = file.path(tmp, "era5_pl", "{site}"),
       derived_dir = file.path(tmp, "derived"), output_dir = file.path(tmp, "output"), pblh_gamma_default = 0.01)
}

fake_runner <- function(fail = character()) {
  calls <- character()
  f <- function(script, args, log) {
    step <- sub("^([0-9a-z]+)_.*$", "\\1", basename(script))
    calls <<- c(calls, paste(basename(script), paste(args, collapse = " ")))
    writeLines("log", log)
    if (any(vapply(fail, function(s) startsWith(basename(script), paste0(s, "_")), logical(1)))) 1L else 0L
  }
  list(run = f, calls = function() calls)
}

test_that("run_site_steps runs, skips current steps, reruns on --force and on changed inputs", {
  cfg <- local_pipeline_cfg()
  steps <- c("01", "11a", "02", "04")
  r1 <- fake_runner()
  res <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = r1$run))
  expect_equal(res$action, rep("ok", 4))
  expect_length(r1$calls(), 4)
  st <- read_run_status(run_status_path("XX-Tst", "02", cfg))
  expect_equal(st[c("site_id", "step", "state", "exit_code")], list(site_id = "XX-Tst", step = "02", state = "ok", exit_code = 0L))
  expect_true(all(c("start", "end", "git_commit", "r_version", "host", "inputs_hash", "log") %in% names(st)))
  expect_true(file.exists(st$log))

  r2 <- fake_runner()
  res2 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = r2$run))
  expect_equal(res2$action, rep("skipped", 4))
  expect_length(r2$calls(), 0)
  expect_equal(res2$inputs_hash, res$inputs_hash)

  r3 <- fake_runner()
  expect_equal(suppressMessages(run_site_steps("XX-Tst", "02", cfg, force = TRUE, runner = r3$run))$action, "ok")
  expect_length(r3$calls(), 1)

  # a changed ERA5 pressure-level file reruns 11a and, through the upstream hash, 02 and 04 (not 01)
  pl <- list.files(resolve_path(cfg$era5_pl_dir, site = "XX-Tst"), full.names = TRUE)
  writeLines(c("p", "q"), pl)
  r4 <- fake_runner()
  res4 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = r4$run))
  expect_equal(res4$action, c("skipped", "ok", "ok", "ok"))
  # a changed argument (--with-oper) reruns the steps that accept it
  r5 <- fake_runner()
  res5 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, with_oper = TRUE, runner = r5$run))
  expect_equal(res5$action, c("skipped", "skipped", "skipped", "ok"))
  expect_match(r5$calls(), "--with-oper")
})

test_that("a failing step stops the site and is rerun next time", {
  cfg <- local_pipeline_cfg()
  steps <- c("01", "11a", "02", "04")
  r1 <- fake_runner(fail = "02")
  res <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = r1$run))
  expect_equal(res$action, c("ok", "ok", "failed", "not_run"))
  expect_equal(res$exit_code[3], 1L)
  expect_equal(read_run_status(run_status_path("XX-Tst", "02", cfg))$state, "failed")
  expect_null(read_run_status(run_status_path("XX-Tst", "04", cfg)))
  r2 <- fake_runner()
  res2 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = r2$run))
  expect_equal(res2$action, c("skipped", "skipped", "ok", "ok"))
  expect_error(run_site_steps("XX-Tst", "99", cfg, runner = r2$run), "Unknown step")
})

test_that("step_is_current and the inputs hash", {
  expect_true(step_is_current(list(state = "ok", inputs_hash = "a"), "a"))
  expect_false(step_is_current(list(state = "ok", inputs_hash = "a"), "a", force = TRUE))
  expect_false(step_is_current(list(state = "failed", inputs_hash = "a"), "a"))
  expect_false(step_is_current(list(state = "ok", inputs_hash = "a"), "b"))
  expect_false(step_is_current(NULL, "a"))
  cfg <- local_pipeline_cfg()
  inp <- site_input_files("XX-Tst", cfg, c("l2", "era5_pl", "bif"))
  expect_equal(sort(unique(inp$group)), c("era5_pl", "l2")) # no BIF: no rows, the step reports it
  h <- step_inputs_hash("02", "XX-Tst", cfg, c(`11a` = "x"))
  expect_equal(h, step_inputs_hash("02", "XX-Tst", cfg, c(`11a` = "x")))
  expect_false(h == step_inputs_hash("02", "XX-Tst", cfg, c(`11a` = "y")))
  expect_false(h == step_inputs_hash("02", "XX-Tst", utils::modifyList(cfg, list(pblh_gamma_default = 0.005)), c(`11a` = "x")))
  expect_false(h == step_inputs_hash("02", "XX-Tst", cfg, c(`11a` = "x"), args = "--with-oper"))
})
