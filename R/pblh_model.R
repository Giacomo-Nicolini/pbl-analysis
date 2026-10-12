# PBLHmod v2: options and orchestration (daily cycle).
# Spec: docs/specs/02_pblh_model.md section 6 and Amendments.

#' PBLH model options
#'
#' Returns the model options with defaults; any can be overridden through
#' `...`. Unknown option names or invalid values abort.
#'
#' @param ... Options to override:
#'   - `stable_method` [`"nieuwstadt"` (Eq. B1, default), `"zilitinkevich"`
#'     (Eq. B4), `"neutral"` (Eq. B2)] formula for L > 0;
#'   - `c_s` [-] coefficient of Eq. B4 (default 0.4);
#'   - `gamma` [K m-1] potential-temperature gradient above the CBL, scalar or
#'     one value per half-hour (default 0.01; Batchvarova & Gryning 1991 as cited
#'     in Kljun et al. 2015);
#'   - `h_source` [`"H_F_MDS"` (default), `"H"`, `"H_CORR"`] sensible heat flux;
#'   - `h_qc_max` [-] highest accepted `H_F_MDS_QC` (default 1);
#'   - `h_fallback` [logical] use FLUXES `H` where `H_F_MDS` is missing in the
#'     file (default `TRUE`);
#'   - `max_gap_interp` [half-hours] longest forcing gap filled by linear
#'     interpolation (default 2);
#'   - `n_sub` [-] RK4 sub-steps per half-hour (default 30);
#'   - `n_neg_end` [half-hours] consecutive `wtheta <= 0` ending the convective
#'     window (default 2);
#'   - `combine` [`"max"` (default) or `"bg_only"`] PBLH in the convective window;
#'   - `label` [character] run label (default `"default"`).
#' @return Named list of options.
pblh_options <- function(...) {
  defaults <- list(
    stable_method = "nieuwstadt",
    c_s = 0.4,
    gamma = 0.01,
    h_source = "H_F_MDS",
    h_qc_max = 1,
    h_fallback = TRUE,
    max_gap_interp = 2,
    n_sub = 30,
    n_neg_end = 2,
    combine = "max",
    label = "default"
  )
  user <- list(...)
  unknown <- setdiff(names(user), names(defaults))
  if (length(unknown) > 0 || (length(user) > 0 && is.null(names(user)))) {
    cli::cli_abort(c(
      "Unknown PBLH option{?s}: {.val {unknown}}",
      "i" = "Known options: {.val {names(defaults)}}"
    ))
  }
  opts <- utils::modifyList(defaults, user)
  validate_pblh_options(opts)
  opts
}

validate_pblh_options <- function(opts) {
  one_of <- function(name, allowed) {
    v <- opts[[name]]
    if (!(is.character(v) && length(v) == 1 && v %in% allowed)) {
      cli::cli_abort("Option {.field {name}} must be one of {.val {allowed}}, not {.val {v}}.")
    }
  }
  pos_num <- function(name) {
    v <- opts[[name]]
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v) && v > 0)) {
      cli::cli_abort("Option {.field {name}} must be a single positive number.")
    }
  }
  int_ge <- function(name, min) {
    v <- opts[[name]]
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v) && v == round(v) && v >= min)) {
      cli::cli_abort("Option {.field {name}} must be a single integer >= {min}.")
    }
  }
  one_of("stable_method", c("nieuwstadt", "zilitinkevich", "neutral"))
  one_of("h_source", c("H_F_MDS", "H", "H_CORR"))
  one_of("combine", c("max", "bg_only"))
  pos_num("c_s")
  g <- opts$gamma
  if (!(is.numeric(g) && length(g) >= 1 && all(is.finite(g)) && all(g > 0))) {
    cli::cli_abort("Option {.field gamma} must be positive and finite (scalar or one value per half-hour) [K m-1].")
  }
  int_ge("h_qc_max", 0)
  if (!(is.logical(opts$h_fallback) && length(opts$h_fallback) == 1 && !is.na(opts$h_fallback))) {
    cli::cli_abort("Option {.field h_fallback} must be TRUE or FALSE.")
  }
  int_ge("max_gap_interp", 0)
  int_ge("n_sub", 1)
  int_ge("n_neg_end", 1)
  if (!(is.character(opts$label) && length(opts$label) == 1 && nzchar(opts$label))) {
    cli::cli_abort("Option {.field label} must be a non-empty string.")
  }
  invisible(TRUE)
}
