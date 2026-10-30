# Spec 09 A2: run pipeline steps for one site, each in a fresh R process, with run-status files.
# Steps (R/run_status.R, pipeline_steps()): 01 input check, 11a ERA5 gamma, 02 PBLHmod (all cases, incl. PBLHmod-clim
# and the reference gamma values), 04 standard files and ERA5 comparison, 08a footprint inputs and coverage (full
# period). A step whose last run is "ok" with the same inputs hash is skipped unless --force; the first failing step
# stops the site (exit code 1). Logs: <output_dir>/logs/<site>/<step>_<timestamp>.log; status:
# <derived_dir>/run_status/<site>__<step>.json.
#
# Usage: Rscript scripts/run_site.R SITE [--steps=01,11a,02,04,08a] [--force] [--with-oper]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster); the child processes inherit it.

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")
steps <- strsplit(opt("steps", "01,11a,02,04,08a"), ",")[[1]]
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
if (!site_id %in% sites$site_id) cli::cli_abort("{.val {site_id}} not in {.file {cfg$sites_file}}.")

cli::cli_h1("{site_id}: steps {paste(steps, collapse = ', ')}")
t0 <- Sys.time()
res <- run_site_steps(site_id, steps, cfg, force = "--force" %in% args, with_oper = "--with-oper" %in% args)
print(as.data.frame(res[c("step", "action", "exit_code")]), row.names = FALSE)
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
if (any(res$action == "failed")) quit(status = 1)
