# Spec 09 A3: write the task list, render the SLURM array script and submit it (or only print it with --dry-run).
#   --what=sites   one task per site: config/run_sites.txt (default: all stations of config/sites.csv) or
#                  --sites-file=PATH; each task runs scripts/run_site.R SITE --steps=... (default 01,11a,02,04,08a)
#   --what=chunks  one task per footprint chunk: the chunk list of scripts/08c_make_chunks.R
#                  (<derived_dir>/footprint/chunk_list.csv, or --chunk-list=PATH)
# Machine settings from slurm/local.env (untracked; template slurm/local.env.example), or --env=PATH.
# The rendered script is written to <output_dir>/slurm/<what>_<timestamp>.sbatch and printed.
#
# Usage: Rscript scripts/slurm_submit.R --what=sites|chunks [--dry-run] [--steps=01,11a,02,04,08a] [--force]
#          [--with-oper] [--sites-file=PATH] [--chunk-list=PATH] [--env=PATH]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
what <- opt("what")
if (!isTRUE(what %in% c("sites", "chunks"))) cli::cli_abort("Use {.code --what=sites} or {.code --what=chunks}.")
dry <- "--dry-run" %in% args
cfg <- get_cfg()
env <- read_local_env(opt("env", here::here("slurm", "local.env")))
out_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "slurm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")

if (what == "sites") {
  src <- opt("sites-file", here::here("config", "run_sites.txt"))
  ids <- trimws(readLines(src, warn = FALSE))
  ids <- ids[nzchar(ids) & !startsWith(ids, "#")]
  known <- read_sites(resolve_path(cfg$sites_file))$site_id
  if (length(setdiff(ids, known))) cli::cli_abort("Unknown site{?s} in {.path {src}}: {.val {setdiff(ids, known)}}.")
  list_file <- file.path(out_dir, paste0("sites_", stamp, ".txt"))
  writeLines(ids, list_file)
  extra <- c(if ("--force" %in% args) "--force", if ("--with-oper" %in% args) "--with-oper")
  vals <- sbatch_values(env, normalizePath(list_file, winslash = "/"), length(ids), opt("steps", "01,11a,02,04,08a"), extra)
  tpl <- readLines(here::here("slurm", "sites_array.sbatch"))
  cli::cli_inform("{length(ids)} site{?s} from {.path {src}}.")
} else {
  src <- opt("chunk-list", file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", "chunk_list.csv"))
  if (!file.exists(src)) cli::cli_abort("Chunk list {.path {src}} not found: run scripts/08c_make_chunks.R first.")
  n <- nrow(readr::read_csv(src, show_col_types = FALSE, progress = FALSE))
  vals <- sbatch_values(env, normalizePath(src, winslash = "/"), n)
  tpl <- readLines(here::here("slurm", "footprint_chunks.sbatch"))
  cli::cli_inform("{n} chunk{?s} from {.path {src}}.")
}

job <- render_sbatch(tpl, vals)
job_file <- file.path(out_dir, paste0(what, "_", stamp, ".sbatch"))
writeLines(job, job_file)
cat(job, sep = "\n")
cli::cli_alert_success("Job script written to {.path {job_file}}.")
if (dry) {
  cli::cli_alert_info("--dry-run: not submitted.")
} else {
  r <- system2("sbatch", shQuote(job_file), stdout = TRUE, stderr = TRUE)
  cat(r, sep = "\n")
  if (!is.null(attr(r, "status")) && attr(r, "status") != 0) cli::cli_abort("sbatch failed.")
}
