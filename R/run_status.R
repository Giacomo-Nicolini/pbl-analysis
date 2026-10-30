# Per-site pipeline driver: step registry, inputs hash and run-status files (spec 09 A2).
# Each step runs one script for one site in a fresh R process (scripts/run_site.R). A status file per site x step
# in <derived_dir>/run_status/<site>__<step>.json records start, end, exit code, git commit, R version, host and the
# inputs hash; a step whose last status is "ok" with the same hash is skipped unless forced. The hash covers the
# step script, all R/ files, the configuration, the site's raw input files (name, size, modification time), the
# extra arguments and the hashes of the upstream steps, so a changed upstream step also reruns its dependants.

#' Pipeline steps
#'
#' @return Named list per step: `script` (path relative to the project root), `upstream` (steps whose outputs it
#'   reads), `inputs` (raw input groups hashed, see [site_input_files()]), `oper` (logical: accepts `--with-oper`).
pipeline_steps <- function() {
  list(
    `01` = list(script = "scripts/01_check_inputs.R", upstream = character(), inputs = c("l2", "era5_sl"), oper = FALSE),
    `11a` = list(script = "scripts/11a_gamma_era5.R", upstream = character(), inputs = c("era5_pl", "era5_static"),
                 oper = FALSE),
    `02` = list(script = "scripts/02_run_pblh_model.R", upstream = "11a", inputs = "l2", oper = FALSE),
    `04` = list(script = "scripts/04_era5_compare.R", upstream = "02", inputs = c("l2", "era5_sl", "era5_pl", "era5_static"),
                oper = TRUE),
    `08a` = list(script = "scripts/08a_footprint_flags.R", upstream = "04", inputs = c("l2", "bif", "measured_std"),
                 oper = TRUE)
  )
}

#' Raw input files of a site
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list ([get_cfg()]).
#' @param groups Input groups: `l2` (ICOS L2 folder), `bif` (ONEFlux BIF), `era5_sl`, `era5_pl` (ERA5 folders),
#'   `era5_static` (static files of the site), `measured_std` (standard files of measured sources).
#' @return Tibble: `group`, `path`, `size` [bytes], `mtime` [character, UTC]; a group whose folder is not configured
#'   or missing gives no rows (the step itself reports missing inputs).
site_input_files <- function(site_id, cfg, groups) {
  safe <- function(expr) tryCatch(expr, error = function(e) character())
  in_dir <- function(key, pattern = NULL) {
    safe({
      d <- resolve_path(cfg[[key]], site = site_id, must_exist = FALSE)
      list.files(d, pattern = pattern, full.names = TRUE)
    })
  }
  site_re <- gsub("([.\\-])", "\\\\\\1", site_id)
  files <- lapply(stats::setNames(groups, groups), function(g) {
    switch(g,
           l2 = safe(list.files(file.path(resolve_path(cfg$icos_l2_dir, site = site_id), paste0(site_id, "_output_L2I_2026")),
                                pattern = "\\.zip$", full.names = TRUE)),
           bif = safe(find_oneflux_bif(site_id, cfg, find_icos_files(site_id, cfg)$fluxnet)),
           era5_sl = in_dir("era5_sl_dir"),
           era5_pl = in_dir("era5_pl_dir"),
           era5_static = safe(unlist(lapply(c("z", "lsm", "sdor"), function(v) {
             d <- resolve_path(cfg$era5_static_dir, site = site_id, var = v, must_exist = FALSE)
             list.files(d, pattern = paste0("^", site_re, "__.*__", v, "__"), full.names = TRUE)
           }))),
           measured_std = in_dir("standard_dir", paste0("^", site_re, "__(cloudnet|sirta|meteoswiss).*\\.csv\\.gz$")),
           cli::cli_abort("Unknown input group {.val {g}}."))
  })
  out <- dplyr::bind_rows(lapply(names(files), function(g) {
    f <- sort(unique(files[[g]]))
    if (length(f) == 0) return(NULL)
    i <- file.info(f)
    dplyr::tibble(group = g, path = f, size = i$size, mtime = format(i$mtime, "%Y-%m-%d %H:%M:%S", tz = "UTC"))
  }))
  if (nrow(out) == 0) out <- dplyr::tibble(group = character(), path = character(), size = numeric(), mtime = character())
  out
}

#' Inputs hash of a step
#'
#' @param step Step name (see [pipeline_steps()]).
#' @param site_id Station code [character].
#' @param cfg Configuration list.
#' @param upstream_hash Named character vector of the hashes recorded for the upstream steps (`NA` if never run).
#' @param args Extra command-line arguments of the step [character].
#' @param root Project root.
#' @param input_files Output of [site_input_files()] (default: computed for the step's groups).
#' @return MD5 hash [character].
step_inputs_hash <- function(step, site_id, cfg, upstream_hash = character(), args = character(), root = here::here(),
                             input_files = NULL) {
  st <- pipeline_steps()[[step]]
  if (is.null(st)) cli::cli_abort("Unknown step {.val {step}}.")
  if (is.null(input_files)) input_files <- site_input_files(site_id, cfg, st$inputs)
  code <- c(file.path(root, st$script), sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)))
  digest::digest(list(
    step = step, site = site_id, code = unname(tools::md5sum(code)),
    cfg = cfg[sort(names(cfg))], inputs = as.data.frame(input_files[c("path", "size", "mtime")]),
    upstream = upstream_hash[sort(names(upstream_hash))], args = sort(args)
  ), algo = "md5")
}

#' Run-status file path
#'
#' @param site_id Station code; @param step Step name; @param cfg Configuration list.
#' @return Path [character].
run_status_path <- function(site_id, step, cfg) {
  file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "run_status", paste0(site_id, "__", step, ".json"))
}

#' Read / write a run-status file
#'
#' @param path Status file path.
#' @param status Named list (`site_id`, `step`, `state` = `running` / `ok` / `failed`, `start`, `end`, `exit_code`,
#'   `inputs_hash`, `git_commit`, `r_version`, `host`, `log`).
#' @return [read_run_status()]: the list, or `NULL` if the file does not exist; [write_run_status()]: `path`,
#'   invisibly.
read_run_status <- function(path) {
  if (!file.exists(path)) return(NULL)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

#' @rdname read_run_status
write_run_status <- function(status, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(status, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  file.rename(tmp, path)
  invisible(path)
}

#' Skip a step?
#'
#' @param status Previous status ([read_run_status()]) or `NULL`.
#' @param hash Current inputs hash.
#' @param force Rerun anyway [logical].
#' @return Logical: `TRUE` if the last run is `"ok"` with the same hash and not forced.
step_is_current <- function(status, hash, force = FALSE) {
  !force && !is.null(status) && identical(status$state, "ok") && identical(status$inputs_hash, hash)
}

# Git commit of the working tree ("<sha>", "<sha>-dirty", or NA without git).
git_commit_id <- function(root = here::here()) {
  sha <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE,
                                           stderr = FALSE)), error = function(e) character())
  if (length(sha) != 1 || !grepl("^[0-9a-f]{40}$", sha)) return(NA_character_)
  dirty <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(root), "status", "--porcelain", "--untracked-files=no"),
                                             stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (length(dirty) > 0) paste0(sha, "-dirty") else sha
}

# Run one script in a fresh R process; output (stdout and stderr) appended to `log`. Returns the exit code.
run_script_rscript <- function(script, args, log) {
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  r <- callr::rscript(script, cmdargs = args, stdout = log, stderr = "2>&1", fail_on_status = FALSE, show = FALSE,
                      wd = here::here())
  r$status
}

#' Run pipeline steps for one site
#'
#' Runs `steps` in the order given; each step is skipped if current ([step_is_current()]), else run with `runner`
#' and its status written before (`running`) and after (`ok` / `failed`). The first failing step stops the site.
#'
#' @param site_id Station code [character].
#' @param steps Step names [character], e.g. `c("01", "11a", "02", "04", "08a")`.
#' @param cfg Configuration list.
#' @param force Rerun current steps [logical].
#' @param with_oper Pass `--with-oper` to the steps that accept it [logical].
#' @param runner Function `(script, args, log)` returning the exit code (default: a fresh `Rscript` process).
#' @param root Project root.
#' @return Tibble per step: `step`, `action` (`skipped`, `ok`, `failed`, `not_run`), `exit_code`, `inputs_hash`, `log`.
run_site_steps <- function(site_id, steps, cfg, force = FALSE, with_oper = FALSE, runner = run_script_rscript,
                           root = here::here()) {
  reg <- pipeline_steps()
  bad <- setdiff(steps, names(reg))
  if (length(bad)) cli::cli_abort(c("Unknown step{?s} {.val {bad}}.", "i" = "Known steps: {.val {names(reg)}}."))
  log_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "logs", site_id)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  res <- dplyr::tibble(step = steps, action = "not_run", exit_code = NA_integer_, inputs_hash = NA_character_,
                       log = NA_character_)
  for (k in seq_along(steps)) {
    s <- steps[k]
    st <- reg[[s]]
    args <- c(site_id, if (with_oper && st$oper) "--with-oper")
    up <- vapply(st$upstream, function(u) {
      x <- read_run_status(run_status_path(site_id, u, cfg))
      if (is.null(x) || !identical(x$state, "ok")) NA_character_ else x$inputs_hash
    }, character(1))
    if (anyNA(up)) {
      cli::cli_alert_warning("{site_id} {s}: upstream step{?s} {.val {names(up)[is.na(up)]}} without a successful run.")
    }
    hash <- step_inputs_hash(s, site_id, cfg, up, args[-1], root)
    res$inputs_hash[k] <- hash
    sp <- run_status_path(site_id, s, cfg)
    if (step_is_current(read_run_status(sp), hash, force)) {
      res$action[k] <- "skipped"
      cli::cli_alert_info("{site_id} {s}: current (same inputs hash), skipped.")
      next
    }
    log <- file.path(log_dir, paste0(s, "_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), ".log"))
    status <- list(site_id = site_id, step = s, state = "running", start = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   end = NA, exit_code = NA, inputs_hash = hash, git_commit = git_commit_id(root),
                   r_version = R.version.string, host = Sys.info()[["nodename"]], args = args, log = log)
    write_run_status(status, sp)
    cli::cli_alert("{site_id} {s}: running {.file {st$script}} (log {.path {log}}).")
    code <- runner(file.path(root, st$script), args, log)
    status$end <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    status$exit_code <- code
    status$state <- if (identical(as.integer(code), 0L)) "ok" else "failed"
    write_run_status(status, sp)
    res$exit_code[k] <- as.integer(code)
    res$log[k] <- log
    res$action[k] <- status$state
    if (status$state == "failed") {
      cli::cli_alert_danger("{site_id} {s}: failed (exit code {code}); later steps not run. Log: {.path {log}}")
      break
    }
    cli::cli_alert_success("{site_id} {s}: ok.")
  }
  res
}
