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

test_that("run_site_exit_status: non-zero when a step failed, propagated as the Rscript exit code (spec 09 amendment 1)", {
  cfg <- local_pipeline_cfg()
  steps <- c("01", "11a", "02", "04")
  res <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = fake_runner(fail = "11a")$run))
  expect_equal(res$action, c("ok", "failed", "not_run", "not_run"))
  expect_identical(run_site_exit_status(res), 1L)
  res2 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = fake_runner()$run))
  expect_equal(res2$action, c("skipped", "ok", "ok", "ok"))
  expect_identical(run_site_exit_status(res2), 0L)
  res3 <- suppressMessages(run_site_steps("XX-Tst", steps, cfg, runner = fake_runner()$run))
  expect_identical(run_site_exit_status(res3), 0L) # all skipped
  expect_identical(run_site_exit_status(res[0, ]), 1L)
  # the final lines of scripts/run_site.R in a fresh Rscript process: summary printed, then the exit code
  tail_src <- grep("run_site_exit_status|quit[(]|^if [(]status", readLines(here::here("scripts", "run_site.R")), value = TRUE)
  expect_length(tail_src, 3)
  for (case in list(list(res = res, status = 1L), list(res = res2, status = 0L))) {
    rds <- withr::local_tempfile(fileext = ".rds")
    saveRDS(case$res, rds)
    scr <- withr::local_tempfile(fileext = ".R")
    writeLines(c(sprintf('for (f in list.files(%s, pattern = "[.]R$", full.names = TRUE)) source(f)', deparse(here::here("R"))),
                 sprintf("res <- readRDS(%s); site_id <- 'XX-Tst'", deparse(rds)),
                 'print(as.data.frame(res[c("step", "action", "exit_code")]), row.names = FALSE)', tail_src), scr)
    out <- withr::local_tempfile()
    r <- callr::rscript(scr, stdout = out, stderr = "2>&1", fail_on_status = FALSE, show = FALSE, wd = here::here())
    expect_identical(r$status, case$status)
    expect_true(any(grepl("step +action +exit_code", readLines(out))))
  }
})

test_that("measured steps (spec 09 amendment 2): registered, scripts by provider, sub-scripts and raw files hashed", {
  reg <- pipeline_steps()
  expect_true(all(c("meas_std", "meas_cmp") %in% names(reg)))
  expect_equal(reg$meas_cmp$upstream, c("04", "meas_std"))
  for (s in c(reg$meas_std$script, reg$meas_std$extra_scripts, reg$meas_cmp$extra_scripts)) {
    expect_true(file.exists(here::here(s)), label = s)
  }
  expect_equal(measured_site_scripts("cloudnet")$standard, "scripts/07c_cloudnet_mlh.R")
  expect_equal(measured_site_scripts("SIRTA-ReOBS")$compare, "scripts/10_measured_compare.R")
  expect_equal(measured_site_scripts("MeteoSwiss")$compare, "scripts/10_measured_compare.R") # spec 13
  expect_error(measured_site_scripts("none"), "provider")
  # primary providers of the committed pairings (spec 13: MeteoSwiss E-PROFILE at five more stations)
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  # FI-Lom removed (G 2026-11-12: no overlap of the Kenttarova record with the EC record)
  run <- c("FI-Hyy", "FR-Gri", "CH-Dav", "DE-RuS", "DE-HoH", "FR-Fon", "FR-Lus")
  expect_equal(vapply(run, function(s) measured_primary_provider(ms, s), ""),
               stats::setNames(c("cloudnet", "SIRTA-ReOBS", rep("MeteoSwiss", 5)), run))
  expect_error(measured_primary_provider(ms, "BE-Bra"), "0 primary")
  expect_error(measured_primary_provider(ms, "CZ-Lnz"), "0 primary") # secondary (spec 14): never primary
  expect_error(measured_primary_provider(ms, "FI-Lom"), "0 primary")
  # config/measured_run_sites.txt: stations with a primary measured source
  ids <- trimws(readLines(here::here("config", "measured_run_sites.txt"), warn = FALSE))
  ids <- ids[nzchar(ids) & !startsWith(ids, "#")]
  expect_equal(ids, run)
  expect_true(all(ids %in% ms$icos_site[!is.na(ms$primary_source_id)]))

  cfg <- local_pipeline_cfg()
  cfg$measured_dir <- file.path(dirname(cfg$derived_dir), "measured", "{site}")
  cfg$cloudnet_dir <- file.path(dirname(cfg$derived_dir), "cloudnet", "{site}")
  md <- file.path(dirname(cfg$derived_dir), "measured", "XX-Tst", "sub")
  dir.create(md, recursive = TRUE)
  writeLines("nc", file.path(md, "MeteoSwiss_STRATfinder_x.nc"))
  inp <- site_input_files("XX-Tst", cfg, "measured_raw")
  expect_equal(basename(inp$path), "MeteoSwiss_STRATfinder_x.nc") # recursive; missing Cloudnet folder: no rows
  # the fake runner sees the dispatcher with --part; a changed sub-script changes the hash
  r <- fake_runner()
  res <- suppressMessages(run_site_steps("XX-Tst", c("meas_std"), cfg, runner = r$run))
  expect_equal(res$action, "ok")
  expect_equal(r$calls(), "09m_measured.R XX-Tst --part=standard")
  root <- withr::local_tempdir()
  dir.create(file.path(root, "R"))
  dir.create(file.path(root, "scripts"))
  for (s in c(reg$meas_std$script, reg$meas_std$extra_scripts)) writeLines("x", file.path(root, s))
  h1 <- step_inputs_hash("meas_std", "XX-Tst", cfg, root = root, input_files = inp)
  writeLines("y", file.path(root, "scripts/11b_meteoswiss_mlh.R"))
  expect_false(identical(h1, step_inputs_hash("meas_std", "XX-Tst", cfg, root = root, input_files = inp)))
  # spec 13: the pairing table and the measured-source exclusions are inputs of meas_std and meas_cmp
  expect_true("measured_config" %in% reg$meas_std$inputs)
  expect_true("measured_config" %in% reg$meas_cmp$inputs)
  cf <- file.path(dirname(cfg$derived_dir), c("ms.csv", "mex.csv"))
  for (p in cf) writeLines("a", p)
  cfg$measured_sites_file <- cf[1]
  cfg$measured_exclusions_file <- cf[2]
  mc <- site_input_files("XX-Tst", cfg, "measured_config")
  expect_equal(sort(basename(mc$path)), c("mex.csv", "ms.csv"))
})
