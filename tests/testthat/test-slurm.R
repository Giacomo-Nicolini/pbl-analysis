# Spec 09 A3: SLURM job scripts rendered from the tracked templates and a fake local.env.

fake_env <- function(extra = character()) {
  p <- withr::local_tempfile(.local_envir = parent.frame(), fileext = ".env")
  writeLines(c("# test", "SLURM_PARTITION=part_x", "SLURM_MAIL_USER=", "SLURM_MAIL_TYPE=FAIL  # comment",
               "SITES_CPUS_PER_TASK=1", "SITES_MEM=16G", "CHUNK_CPUS_PER_TASK=6", "CHUNK_MEM_PER_CPU=4G", "SLURM_TIME=",
               "SLURM_ARRAY_MAX=10", "R_SETUP=\"module load R/x; source venv_x\"", "R_CONFIG_ACTIVE=cluster",
               "PROJECT_DIR=/proj_x", "LOG_DIR=/logs_x", extra), p)
  p
}

test_that("read_local_env parses KEY=value with comments and quotes", {
  e <- read_local_env(fake_env())
  expect_equal(e[["SLURM_PARTITION"]], "part_x")
  expect_equal(e[["SLURM_MAIL_USER"]], "")
  expect_equal(e[["SLURM_MAIL_TYPE"]], "FAIL")
  expect_equal(e[["R_SETUP"]], "module load R/x; source venv_x")
  expect_error(read_local_env(file.path(tempdir(), "nope.env")), "local.env.example")
  expect_error(read_local_env(fake_env("not a pair")), "KEY=value")
})

test_that("site and chunk templates render to valid job scripts with only local.env machine data", {
  e <- read_local_env(fake_env())
  s <- render_sbatch(readLines(here::here("slurm", "sites_array.sbatch")),
                     sbatch_values(e, "/lists_x/sites.txt", 40, extra_args = "--with-oper"))
  expect_equal(s[1], "#!/bin/bash")
  expect_false(any(grepl("{{", s, fixed = TRUE)))
  expect_true("#SBATCH --partition=part_x" %in% s)
  expect_true("#SBATCH --array=1-40%10" %in% s)
  expect_true("#SBATCH --cpus-per-task=1" %in% s)
  expect_true("#SBATCH --mail-type=FAIL" %in% s)
  expect_false(any(grepl("--mail-user|--time", s))) # empty values: lines dropped
  expect_true(all(c("module load R/x", "source venv_x") %in% s))
  expect_true(any(grepl("run_site.R \"${SITE}\" --steps=01,11a,02,04,08a --with-oper", s, fixed = TRUE)))
  # SBATCH lines are contiguous after the shebang and comments (sbatch stops reading directives at the first command)
  first_cmd <- which(!grepl("^#|^$", s))[1]
  expect_true(all(grep("^#SBATCH", s) < first_cmd))

  ch <- render_sbatch(readLines(here::here("slurm", "footprint_chunks.sbatch")),
                      sbatch_values(read_local_env(fake_env(c("SLURM_ARRAY_MAX="))), "/lists_x/chunks.csv", 7))
  expect_true("#SBATCH --array=1-7" %in% ch)
  expect_true("#SBATCH --cpus-per-task=6" %in% ch)
  expect_true("#SBATCH --mem-per-cpu=4G" %in% ch)
  expect_true(any(grepl("--chunk-list=\"/lists_x/chunks.csv\"", ch, fixed = TRUE)))
  # only the fake machine values appear
  for (x in list(s, ch)) {
    expect_false(any(grepl("/bigdata|/home/|@unitus|nicolini|trotta", x)))
  }
  expect_error(render_sbatch("x {{NOPE}}", sbatch_values(e, "/l", 1)), "NOPE")
  e2 <- e
  e2[["SLURM_PARTITION"]] <- ""
  expect_error(render_sbatch(readLines(here::here("slurm", "sites_array.sbatch")), sbatch_values(e2, "/l", 1)),
               "SLURM_PARTITION")
  expect_error(sbatch_values(e, "/l", 0), "n_tasks")
})

test_that("tracked SLURM files and the example env hold no machine data", {
  f <- c(list.files(here::here("slurm"), full.names = TRUE), here::here("config.local.yml.example"))
  f <- f[basename(f) != "local.env"]
  txt <- unlist(lapply(f, readLines, warn = FALSE))
  expect_false(any(grepl("/bigdata|/home/|@unitus|g\\.nicolini|c\\.trotta|ONEFLUX", txt)))
})
