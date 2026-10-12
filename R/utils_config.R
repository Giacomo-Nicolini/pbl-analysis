# Configuration helpers: path profiles from config.yml.
# Profiles: `default` (local PC) and `cluster`, selected with env var R_CONFIG_ACTIVE.

#' Load the project configuration
#'
#' Wraps `config::get()` on the project `config.yml`. The active profile is
#' chosen by the environment variable `R_CONFIG_ACTIVE` (`default` if unset).
#'
#' @param profile Optional profile name [character]; `NULL` uses `R_CONFIG_ACTIVE`.
#' @return Named list of configuration values (paths as character, unresolved).
get_cfg <- function(profile = NULL) {
  file <- here::here("config.yml")
  if (is.null(profile)) {
    config::get(file = file)
  } else {
    config::get(file = file, config = profile)
  }
}

#' Resolve a configured path
#'
#' Replaces the `{site}` placeholder with the station code, resolves relative
#' paths against the project root (`here::here()`) and returns absolute paths
#' unchanged.
#'
#' @param path Path from `config.yml` [character, length 1].
#' @param site Station code, e.g. `"BE-Bra"` [character]; required if `path`
#'   contains `{site}`.
#' @param must_exist If `TRUE`, abort when the resolved path does not exist [logical].
#' @return Absolute path [character, length 1].
resolve_path <- function(path, site = NULL, must_exist = TRUE) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    cli::cli_abort("Configured path is missing or empty.")
  }
  if (identical(toupper(path), "TODO")) {
    cli::cli_abort(c(
      "Configured path is {.val TODO}.",
      "i" = "Set it in {.file config.yml} for the active profile ({.envvar R_CONFIG_ACTIVE} = {.val {Sys.getenv('R_CONFIG_ACTIVE', 'default')}})."
    ))
  }
  if (grepl("{site}", path, fixed = TRUE)) {
    if (is.null(site) || is.na(site) || !nzchar(site)) {
      cli::cli_abort("Path {.path {path}} contains {.code {{site}}} but no {.arg site} was given.")
    }
    path <- gsub("{site}", site, path, fixed = TRUE)
  }
  if (!is_absolute_path(path)) {
    path <- here::here(path)
  }
  if (must_exist && !file.exists(path)) {
    cli::cli_abort("Path does not exist: {.path {path}}")
  }
  path
}

# TRUE for "/x", "~", "C:/x", "C:\x", "\\server\share".
is_absolute_path <- function(path) {
  grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}
