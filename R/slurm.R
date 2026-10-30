# SLURM job scripts (spec 09 A3): templates in slurm/*.sbatch with {{KEY}} placeholders, filled from the untracked
# slurm/local.env (machine data: partition, mail, R setup, project and log folders) and from slurm_submit.R (list
# file, number of tasks, steps). Tracked files hold no machine data.

#' Read a local.env file
#'
#' `KEY=value` lines; `#` starts a comment (whole line, or after a space); surrounding quotes are removed; a
#' repeated key takes its last value.
#'
#' @param path File path.
#' @return Named character vector (empty values kept as `""`).
read_local_env <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort(c("{.path {path}} not found.", "i" = "Copy {.file slurm/local.env.example} to {.file slurm/local.env} and fill it in."))
  }
  x <- readLines(path, warn = FALSE)
  x <- sub("(^|\\s)#.*$", "", x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  bad <- !grepl("^[A-Za-z_][A-Za-z0-9_]*=", x)
  if (any(bad)) cli::cli_abort("{.path {path}}: not KEY=value: {.val {x[bad]}}.")
  key <- sub("=.*$", "", x)
  val <- trimws(sub("^[^=]*=", "", x))
  val <- sub("^([\"'])(.*)\\1$", "\\2", val)
  last <- !duplicated(key, fromLast = TRUE) # a repeated key: the last value wins, as in a shell
  stats::setNames(val[last], key[last])
}

#' Render an sbatch template
#'
#' Replaces every `{{KEY}}` by `values[["KEY"]]`. `#SBATCH` lines whose value is empty are dropped (e.g. no mail
#' address, no wall time), so the partition default applies; `R_SETUP` may hold several commands separated by `;`
#' (one line each). Aborts if a placeholder has no value, or if a required key is empty.
#'
#' @param template Template lines [character].
#' @param values Named character vector.
#' @param required Keys that must be non-empty.
#' @return Rendered lines [character].
render_sbatch <- function(template, values, required = c("SLURM_PARTITION", "PROJECT_DIR", "LOG_DIR", "LIST_FILE", "N_TASKS")) {
  empty <- required[!required %in% names(values) | !nzchar(values[required])]
  if (length(empty)) cli::cli_abort("Required SLURM setting{?s} empty or missing: {.field {empty}}.")
  keys <- unique(unlist(regmatches(template, gregexpr("\\{\\{[A-Z0-9_]+\\}\\}", template))))
  keys <- gsub("[{}]", "", keys)
  miss <- setdiff(keys, names(values))
  if (length(miss)) cli::cli_abort("No value for placeholder{?s} {.field {miss}} (slurm/local.env or slurm_submit.R).")
  out <- character()
  for (ln in template) {
    k <- gsub("[{}]", "", unlist(regmatches(ln, gregexpr("\\{\\{[A-Z0-9_]+\\}\\}", ln))))
    if (startsWith(ln, "#SBATCH") && length(k) && any(!nzchar(values[k]) & k != "ARRAY_MAX")) next
    if (identical(k, "R_SETUP") && trimws(ln) == "{{R_SETUP}}") {
      out <- c(out, trimws(strsplit(values[["R_SETUP"]], ";", fixed = TRUE)[[1]]))
      next
    }
    for (kk in k) ln <- gsub(paste0("{{", kk, "}}"), values[[kk]], ln, fixed = TRUE)
    out <- c(out, ln)
  }
  out
}

#' Values for a job script
#'
#' @param env Output of [read_local_env()].
#' @param list_file Absolute path of the site or chunk list (one task per row) [character].
#' @param n_tasks Number of array tasks [integer].
#' @param steps Steps for the site array (`run_site.R --steps`), ignored for chunks.
#' @param extra_args Extra arguments for the R script [character].
#' @return Named character vector for [render_sbatch()].
sbatch_values <- function(env, list_file, n_tasks, steps = "01,11a,02,04,08a", extra_args = character()) {
  if (!(is.numeric(n_tasks) && length(n_tasks) == 1 && n_tasks >= 1)) cli::cli_abort("{.arg n_tasks} must be >= 1.")
  amax <- if (is.na(env["SLURM_ARRAY_MAX"]) || !nzchar(env[["SLURM_ARRAY_MAX"]])) "" else paste0("%", env[["SLURM_ARRAY_MAX"]])
  defaults <- c(SLURM_MAIL_USER = "", SLURM_MAIL_TYPE = "", SLURM_TIME = "", R_SETUP = "", R_CONFIG_ACTIVE = "cluster",
                SITES_CPUS_PER_TASK = "1", SITES_MEM = "16G", CHUNK_CPUS_PER_TASK = "6", CHUNK_MEM_PER_CPU = "4G")
  v <- c(env, defaults[setdiff(names(defaults), names(env))])
  v[c("LIST_FILE", "N_TASKS", "ARRAY_MAX", "STEPS", "EXTRA_ARGS")] <-
    c(list_file, as.character(as.integer(n_tasks)), amax, steps, paste(extra_args, collapse = " "))
  v
}
