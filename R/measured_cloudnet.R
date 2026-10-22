# ACTRIS/Cloudnet measured mixing-layer height (spec 07).
# Part 1: site pairing, distance rule and Cloudnet data-portal API (inventory).
# API: https://cloudnet.fmi.fi/api (docs: https://docs.cloudnet.fmi.fi/api/data-portal.html).
# Data licence: CC BY 4.0 (ACTRIS/Cloudnet).

#' Read the ICOS site <-> measurement source pairing table
#'
#' Reads `config/measured_sites.csv` (spec 07 section 1): one row per ICOS site x
#' measurement source.
#'
#' @param file Path to `measured_sites.csv` [character].
#' @return Tibble with `icos_site`, `provider`, `provider_site` [-],
#'   `provider_lat` [decimal degrees N], `provider_lon` [decimal degrees E],
#'   `distance_km` [km, tower to measurement site], `method_class`,
#'   `instruments`, `licence`, `url`, `notes` [-].
read_measured_sites <- function(file) {
  ms <- readr::read_csv(
    file,
    na = c("NA", ""),
    col_types = readr::cols(
      provider_lat = readr::col_double(), provider_lon = readr::col_double(),
      distance_km = readr::col_double(), .default = readr::col_character()
    ),
    progress = FALSE
  )
  need <- c("icos_site", "provider", "provider_site", "provider_lat", "provider_lon", "distance_km",
            "method_class", "instruments", "licence", "url", "notes")
  miss <- setdiff(need, names(ms))
  if (length(miss) > 0) cli::cli_abort("{.path {file}}: missing column{?s} {.field {miss}}.")
  key <- paste(ms$icos_site, ms$provider, ms$provider_site)
  if (anyDuplicated(key)) cli::cli_abort("{.path {file}}: duplicated icos_site x provider x provider_site.")
  bad <- is.na(ms$icos_site) | is.na(ms$provider) | is.na(ms$provider_site)
  if (any(bad)) cli::cli_abort("{.path {file}}: icos_site, provider and provider_site are required (rows {which(bad)}).")
  ms
}

#' Distance class of a measured PBLH source
#'
#' Distance rule for measured PBLH (project plan, G 2026-10-22): < 5 km all regimes;
#' 5-25 km (25 included) daytime convective only; 25-50 km (50 included) daytime
#' climatologies only; > 50 km not used.
#'
#' @param distance_km Tower-to-measurement distance [km].
#' @return Character, one of `"all_regimes"`, `"daytime_convective"`,
#'   `"daytime_climatology"`, `"not_used"` (`NA` where `distance_km` is `NA`).
distance_class <- function(distance_km) {
  out <- rep(NA_character_, length(distance_km))
  out[distance_km < 5] <- "all_regimes"
  out[distance_km >= 5 & distance_km <= 25] <- "daytime_convective"
  out[distance_km > 25 & distance_km <= 50] <- "daytime_climatology"
  out[distance_km > 50] <- "not_used"
  out
}

#' GET a Cloudnet API endpoint (one request, polite pause, one retry)
#'
#' @param endpoint Endpoint below the API root, e.g. `"files"`, `"sites"` [character].
#' @param query Named list of query parameters (`NULL` values dropped) [-].
#' @param base API root [character, URL].
#' @param pause_s Pause after each request [s].
#' @param retries Number of retries after a failure [-].
#' @param backoff_s Wait before each retry [s].
#' @return Parsed JSON (data frame or list, as returned by [jsonlite::fromJSON()]).
cloudnet_api_get <- function(endpoint, query = list(), base = "https://cloudnet.fmi.fi/api",
                             pause_s = 0.5, retries = 1, backoff_s = 10) {
  query <- query[!vapply(query, is.null, logical(1))]
  qs <- if (length(query) > 0) {
    paste0("?", paste(names(query), vapply(query, function(v) utils::URLencode(as.character(v), reserved = TRUE),
                                          character(1)), sep = "=", collapse = "&"))
  } else {
    ""
  }
  url <- paste0(base, "/", endpoint, qs)
  for (i in seq_len(retries + 1)) {
    res <- tryCatch(jsonlite::fromJSON(url), error = function(e) e)
    Sys.sleep(pause_s)
    if (!inherits(res, "error")) return(res)
    if (i <= retries) {
      cli::cli_alert_warning("Cloudnet API: {conditionMessage(res)}; retrying in {backoff_s} s ({url}).")
      Sys.sleep(backoff_s)
    }
  }
  cli::cli_abort(c("Cloudnet API request failed after {retries + 1} attempt{?s}.", "x" = conditionMessage(res),
                   "i" = "URL: {url}"))
}

#' Flatten a Cloudnet `/api/files` response
#'
#' Pure function: converts the parsed JSON of `/api/files` (nested `site`,
#' `product`, `instrument`) into one row per file.
#'
#' @param x Parsed JSON as returned by [jsonlite::fromJSON()] (data frame, or an
#'   empty list when no file matches).
#' @return Tibble, one row per file: `uuid`, `filename`, `site`, `product`,
#'   `measurement_date` [Date], `instrument_name`, `instrument_pid`,
#'   `size_bytes` [bytes], `checksum_sha256`, `version`, `pid` (empty while the
#'   file is volatile), `volatile`, `legacy` [logical], `coverage` [0-1, fraction
#'   of the day covered], `error_level`, `download_url`.
cloudnet_flatten_files <- function(x) {
  empty <- dplyr::tibble(
    uuid = character(), filename = character(), site = character(), product = character(),
    measurement_date = as.Date(character()), instrument_name = character(), instrument_pid = character(),
    size_bytes = numeric(), checksum_sha256 = character(), version = character(), pid = character(),
    volatile = logical(), legacy = logical(), coverage = numeric(), error_level = character(),
    download_url = character()
  )
  if (!is.data.frame(x) || nrow(x) == 0) return(empty)
  col <- function(v, sub = NULL, type = as.character) {
    y <- if (is.null(sub)) x[[v]] else x[[v]][[sub]]
    if (is.null(y)) y <- rep(NA, nrow(x))
    type(y)
  }
  inst_name <- if (is.data.frame(x$instrument)) col("instrument", "name") else rep(NA_character_, nrow(x))
  inst_pid <- if (is.data.frame(x$instrument)) col("instrument", "pid") else rep(NA_character_, nrow(x))
  dplyr::tibble(
    uuid = col("uuid"), filename = col("filename"),
    site = if (is.data.frame(x$site)) col("site", "id") else col("site"),
    product = if (is.data.frame(x$product)) col("product", "id") else col("product"),
    measurement_date = as.Date(col("measurementDate")),
    instrument_name = inst_name, instrument_pid = inst_pid,
    size_bytes = col("size", type = as.numeric), checksum_sha256 = col("checksum"),
    version = col("version"), pid = col("pid"),
    volatile = col("volatile", type = as.logical), legacy = col("legacy", type = as.logical),
    coverage = col("coverage", type = as.numeric), error_level = col("errorLevel"),
    download_url = col("downloadUrl")
  )
}

#' List Cloudnet files of one site and product
#'
#' Queries `/api/files` one calendar year at a time (the API documents no paging),
#' one request at a time with a short pause.
#'
#' @param site Cloudnet site id, e.g. `"hyytiala"` [character].
#' @param product Cloudnet product id, e.g. `"epsilon-lidar"` [character].
#' @param date_from,date_to First and last measurement date [Date or "YYYY-MM-DD"].
#' @param show_legacy Also return legacy (pre-CloudnetPy) files [logical].
#' @param ... Passed to [cloudnet_api_get()] (`base`, `pause_s`, `retries`).
#' @return Tibble as [cloudnet_flatten_files()].
cloudnet_files <- function(site, product, date_from, date_to, show_legacy = TRUE, ...) {
  date_from <- as.Date(date_from)
  date_to <- as.Date(date_to)
  years <- seq(as.integer(format(date_from, "%Y")), as.integer(format(date_to, "%Y")))
  res <- lapply(years, function(y) {
    from <- max(date_from, as.Date(sprintf("%d-01-01", y)))
    to <- min(date_to, as.Date(sprintf("%d-12-31", y)))
    q <- list(site = site, product = product, dateFrom = format(from), dateTo = format(to),
              showLegacy = if (show_legacy) "true" else NULL)
    cloudnet_flatten_files(cloudnet_api_get("files", q, ...))
  })
  dplyr::bind_rows(res)
}

#' Monthly inventory of Cloudnet files
#'
#' Pure function: aggregates file metadata to site x product x instrument x month.
#'
#' @param files Tibble as [cloudnet_flatten_files()].
#' @return Tibble with `site`, `product`, `instrument_name`, `instrument_pid`,
#'   `month` ["YYYY-MM"], `n_files` [-, daily files], `size_mb` [MB, 1e6 bytes],
#'   `n_legacy`, `n_volatile` [-], `mean_coverage` [0-1].
cloudnet_inventory_summary <- function(files) {
  if (nrow(files) == 0) {
    return(dplyr::tibble(site = character(), product = character(), instrument_name = character(),
                         instrument_pid = character(), month = character(), n_files = integer(),
                         size_mb = numeric(), n_legacy = integer(), n_volatile = integer(),
                         mean_coverage = numeric()))
  }
  files$month <- format(files$measurement_date, "%Y-%m")
  out <- dplyr::summarise(
    dplyr::group_by(files, .data$site, .data$product, .data$instrument_name, .data$instrument_pid, .data$month),
    n_files = dplyr::n(),
    size_mb = sum(.data$size_bytes, na.rm = TRUE) / 1e6,
    n_legacy = sum(.data$legacy %in% TRUE),
    n_volatile = sum(.data$volatile %in% TRUE),
    mean_coverage = mean(.data$coverage, na.rm = TRUE),
    .groups = "drop"
  )
  dplyr::arrange(out, .data$site, .data$product, .data$instrument_name, .data$month)
}

#' Contiguous coverage periods from an inventory
#'
#' Pure function: collapses the months with files into runs of consecutive months.
#'
#' @param inv Output of [cloudnet_inventory_summary()].
#' @return Tibble with `site`, `product`, `first_month`, `last_month`
#'   ["YYYY-MM"], `n_months` [-], `n_files` [-], `size_mb` [MB].
cloudnet_coverage_periods <- function(inv) {
  if (nrow(inv) == 0) {
    return(dplyr::tibble(site = character(), product = character(), first_month = character(),
                         last_month = character(), n_months = integer(), n_files = integer(),
                         size_mb = numeric()))
  }
  m <- dplyr::summarise(dplyr::group_by(inv, .data$site, .data$product, .data$month),
                        n_files = sum(.data$n_files), size_mb = sum(.data$size_mb), .groups = "drop")
  m <- dplyr::arrange(m, .data$site, .data$product, .data$month)
  mi <- as.integer(substr(m$month, 1, 4)) * 12 + as.integer(substr(m$month, 6, 7))
  grp <- paste(m$site, m$product)
  new_run <- c(TRUE, grp[-1] != grp[-length(grp)] | diff(mi) != 1)
  m$run <- cumsum(new_run)
  dplyr::summarise(
    dplyr::group_by(m, .data$site, .data$product, .data$run),
    first_month = min(.data$month), last_month = max(.data$month), n_months = dplyr::n(),
    n_files = sum(.data$n_files), size_mb = sum(.data$size_mb), .groups = "drop"
  )[, c("site", "product", "first_month", "last_month", "n_months", "n_files", "size_mb")]
}

# ---- Download (spec 07 section 2) ----

#' SHA-256 checksum of a file
#'
#' @param path File path [character].
#' @return Lower-case hex SHA-256 [character]; `NA` where the file does not exist.
cloudnet_sha256 <- function(path) {
  vapply(path, function(p) {
    if (!file.exists(p)) return(NA_character_)
    digest::digest(file = p, algo = "sha256")
  }, character(1), USE.NAMES = FALSE)
}

#' Which Cloudnet files still need downloading
#'
#' A file is skipped when it exists in `dest_dir` with the SHA-256 given by the API.
#'
#' @param files Tibble as [cloudnet_flatten_files()].
#' @param dest_dir Destination folder [character].
#' @return `files` with `local_path` [character], `present` [logical, file exists
#'   with the API checksum] added.
cloudnet_download_plan <- function(files, dest_dir) {
  files$local_path <- file.path(dest_dir, files$filename)
  sha <- cloudnet_sha256(files$local_path)
  files$present <- !is.na(sha) & sha == files$checksum_sha256
  files
}

#' Download Cloudnet files with checksum verification
#'
#' Downloads each file not yet present (see [cloudnet_download_plan()]) to a
#' temporary name, verifies its SHA-256 against the API and renames it. A failed
#' download or checksum mismatch is retried `retries` times, then logged.
#'
#' @param plan Output of [cloudnet_download_plan()].
#' @param retries Retries per file after a failure [-].
#' @param pause_s Pause after each download [s].
#' @return `plan` with `status` (`"present"`, `"downloaded"`, `"failed"`),
#'   `message` [character] and `downloaded_utc` [POSIXct UTC] added.
cloudnet_download <- function(plan, retries = 1, pause_s = 0.5) {
  n <- nrow(plan)
  plan$status <- ifelse(plan$present, "present", NA_character_)
  plan$message <- NA_character_
  plan$downloaded_utc <- as.POSIXct(rep(NA_real_, n), tz = "UTC")
  todo <- which(!plan$present)
  for (j in seq_along(todo)) {
    i <- todo[j]
    tmp <- paste0(plan$local_path[i], ".part")
    msg <- NA_character_
    for (k in seq_len(retries + 1)) {
      res <- tryCatch(utils::download.file(plan$download_url[i], tmp, mode = "wb", quiet = TRUE),
                      error = function(e) e, warning = function(w) w)
      Sys.sleep(pause_s)
      if (inherits(res, "condition")) {
        msg <- conditionMessage(res)
      } else if (!identical(cloudnet_sha256(tmp), plan$checksum_sha256[i])) {
        msg <- "checksum mismatch"
      } else {
        msg <- NA_character_
        break
      }
    }
    if (is.na(msg)) {
      file.rename(tmp, plan$local_path[i])
      plan$status[i] <- "downloaded"
      plan$downloaded_utc[i] <- as.POSIXct(Sys.time(), tz = "UTC")
    } else {
      unlink(tmp)
      plan$status[i] <- "failed"
      plan$message[i] <- msg
      cli::cli_alert_danger("{plan$filename[i]}: {msg}")
    }
    if (j %% 50 == 0) cli::cli_inform("{j} / {length(todo)} files processed.")
  }
  plan
}
