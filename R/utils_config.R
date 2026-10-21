# Configuration helpers: path profiles from config.yml, with machine-specific
# overrides from the untracked config.local.yml (spec 05 part B).
# Profiles: `default` (local PC) and `cluster`, selected with env var R_CONFIG_ACTIVE.

#' Load the project configuration
#'
#' Reads `config.yml` with `config::get()` for the active profile (`profile`
#' argument, else env var `R_CONFIG_ACTIVE`, else `"default"`). If
#' `local_file` exists (untracked `config.local.yml`, see
#' `config.local.yml.example`), its `default:` section and then its
#' active-profile section are merged over it, key by key at the top level (a
#' local value replaces the base value; lists such as `test_sites` are
#' replaced, not merged). The keys taken from the local file, and those not
#' present in `config.yml`, are reported once per R session.
#'
#' @param profile Optional profile name [character]; `NULL` uses `R_CONFIG_ACTIVE`.
#' @param local_file Path of the local override file [character]; `NULL` or a
#'   missing file = `config.yml` only.
#' @return Named list of configuration values (paths as character, unresolved).
get_cfg <- function(profile = NULL, local_file = here::here("config.local.yml")) {
  active <- if (is.null(profile)) Sys.getenv("R_CONFIG_ACTIVE", "default") else profile
  if (!nzchar(active)) active <- "default"
  base <- config::get(file = here::here("config.yml"), config = active)
  if (is.null(local_file) || !file.exists(local_file)) return(base)

  loc <- yaml::read_yaml(local_file)
  over <- merge_cfg_flat(list(), loc$default)
  if (active != "default") over <- merge_cfg_flat(over, loc[[active]])
  if (length(over) == 0) return(base)
  cfg <- merge_cfg_flat(base, over)

  key <- paste(normalizePath(local_file, winslash = "/"), active)
  done <- getOption("pblh.cfg_local_reported", character(0))
  if (!key %in% done) {
    unknown <- setdiff(names(over), names(base))
    cli::cli_inform(c(
      "Config profile {.val {active}}: {length(over)} key{?s} from {.file {basename(local_file)}}: {.field {names(over)}}.",
      if (length(unknown) > 0) c("!" = "Not in {.file config.yml} (kept): {.field {unknown}}.")
    ))
    options(pblh.cfg_local_reported = c(done, key))
  }
  cfg
}

# Top-level merge: each key of `y` replaces the key of `x` (NULL values ignored).
merge_cfg_flat <- function(x, y) {
  for (k in names(y)) {
    if (!is.null(y[[k]])) x[k] <- list(y[[k]])
  }
  x
}

#' Resolve a configured path
#'
#' Replaces the placeholders `{site}` (station code) and `{var}` (variable name,
#' e.g. ERA5 static fields on the cluster), resolves relative paths against the
#' project root (`here::here()`) and returns absolute paths unchanged. Any
#' `{...}` placeholder left after substitution aborts, and so do the values
#' `"TODO"` (not decided yet) and `"LOCAL"` (machine-specific: set it in
#' `config.local.yml`).
#'
#' @param path Path from `config.yml` [character, length 1].
#' @param site Station code, e.g. `"BE-Bra"` [character]; required if `path`
#'   contains `{site}`.
#' @param var Variable name, e.g. `"z"` [character]; required if `path`
#'   contains `{var}`, ignored otherwise.
#' @param must_exist If `TRUE`, abort when the resolved path does not exist [logical].
#' @return Absolute path [character, length 1].
resolve_path <- function(path, site = NULL, var = NULL, must_exist = TRUE) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    cli::cli_abort("Configured path is missing or empty.")
  }
  if (identical(path, "LOCAL")) {
    cli::cli_abort(c(
      "Configured path is {.val LOCAL}: it is machine-specific.",
      "i" = "Set it in {.file config.local.yml} for the active profile ({.envvar R_CONFIG_ACTIVE} = {.val {Sys.getenv('R_CONFIG_ACTIVE', 'default')}}); see {.file config.local.yml.example}."
    ))
  }
  if (identical(toupper(path), "TODO")) {
    cli::cli_abort(c(
      "Configured path is {.val TODO}.",
      "i" = "Set it in {.file config.yml} for the active profile ({.envvar R_CONFIG_ACTIVE} = {.val {Sys.getenv('R_CONFIG_ACTIVE', 'default')}})."
    ))
  }
  path <- substitute_placeholder(path, "site", site)
  path <- substitute_placeholder(path, "var", var)
  left <- regmatches(path, gregexpr("\\{[^{}]*\\}", path))[[1]]
  if (length(left) > 0) {
    cli::cli_abort("Path {.path {path}} has unresolved placeholder{?s} {.code {left}}.")
  }
  if (!is_absolute_path(path)) {
    path <- here::here(path)
  }
  if (must_exist && !file.exists(path)) {
    cli::cli_abort("Path does not exist: {.path {path}}")
  }
  path
}

# Replace `{name}` in `path` by `value`; abort if the placeholder is present but
# `value` is missing.
substitute_placeholder <- function(path, name, value) {
  ph <- paste0("{", name, "}")
  if (!grepl(ph, path, fixed = TRUE)) {
    return(path)
  }
  if (is.null(value) || length(value) != 1 || is.na(value) || !nzchar(value)) {
    cli::cli_abort("Path {.path {path}} contains {.code {ph}} but no {.arg {name}} was given.")
  }
  gsub(ph, value, path, fixed = TRUE)
}

# TRUE for "/x", "~", "C:/x", "C:\x", "\\server\share".
is_absolute_path <- function(path) {
  grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}
