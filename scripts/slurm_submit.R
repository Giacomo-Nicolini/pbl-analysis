# Spec 09 A3: write the task list, render the SLURM array script and submit it (or only print it with --dry-run).
#   --what=sites   one task per site: config/run_sites.txt (default: all stations of config/sites.csv) or
#                  --sites-file=PATH; each task runs scripts/run_site.R SITE --steps=... (default 01,11a,02,04,08a)
#   --what=chunks  one task per footprint chunk: the chunk list of scripts/08c_make_chunks.R
#                  (<derived_dir>/footprint/chunk_list.csv, or --chunk-list=PATH); a list longer than
#                  SLURM_MAX_ARRAY_TASKS (local.env, default 1000; cluster MaxArraySize 1001) is split into
#                  independent arrays of at most that many tasks, submitted together (SLURM_ARRAY_MAX applies to
#                  each array)
# Machine settings from slurm/local.env (untracked; template slurm/local.env.example), or --env=PATH.
# The rendered scripts are written to <output_dir>/slurm/<what>_<timestamp>[_partN].sbatch and printed.
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
  tpl <- readLines(here::here("slurm", "sites_array.sbatch"))
  jobs <- list(render_sbatch(tpl, sbatch_values(env, normalizePath(list_file, winslash = "/"), length(ids),
                                                opt("steps", "01,11a,02,04,08a"), extra)))
  names(jobs) <- paste0(what, "_", stamp)
  cli::cli_inform("{length(ids)} site{?s} from {.path {src}}.")
} else {
  src <- opt("chunk-list", file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", "chunk_list.csv"))
  if (!file.exists(src)) cli::cli_abort("Chunk list {.path {src}} not found: run scripts/08c_make_chunks.R first.")
  n <- nrow(readr::read_csv(src, show_col_types = FALSE, progress = FALSE))
  tpl <- readLines(here::here("slurm", "footprint_chunks.sbatch"))
  mx <- if (is.na(env["SLURM_MAX_ARRAY_TASKS"]) || !nzchar(env[["SLURM_MAX_ARRAY_TASKS"]])) 1000L else
    as.integer(env[["SLURM_MAX_ARRAY_TASKS"]])
  parts <- array_parts(n, mx)
  jobs <- lapply(seq_len(nrow(parts)), function(i) {
    render_sbatch(tpl, sbatch_values(env, normalizePath(src, winslash = "/"), parts$n_tasks[i],
                                     row_offset = parts$row_offset[i]))
  })
  names(jobs) <- if (nrow(parts) == 1) paste0(what, "_", stamp) else paste0(what, "_", stamp, "_part", parts$part)
  cli::cli_inform("{n} chunk{?s} from {.path {src}}: {nrow(parts)} array{?s} of at most {mx} tasks.")
  for (i in seq_len(nrow(parts))) {
    cli::cli_inform("part {parts$part[i]}: chunk rows {parts$first_row[i]}-{parts$last_row[i]} ({parts$n_tasks[i]} tasks)")
  }
}

job_files <- file.path(out_dir, paste0(names(jobs), ".sbatch"))
for (i in seq_along(jobs)) {
  writeLines(jobs[[i]], job_files[i])
  cli::cli_h2(basename(job_files[i]))
  cat(jobs[[i]], sep = "\n")
}
cli::cli_alert_success("{length(jobs)} job script{?s} written to {.path {out_dir}}.")
if (dry) {
  cli::cli_alert_info("--dry-run: not submitted.")
} else {
  for (f in job_files) {
    r <- system2("sbatch", shQuote(f), stdout = TRUE, stderr = TRUE)
    cat(r, sep = "\n")
    if (!is.null(attr(r, "status")) && attr(r, "status") != 0) cli::cli_abort("sbatch failed for {.path {f}}.")
  }
}
