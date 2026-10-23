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
#'   `instruments`, `licence`, `url`, `notes` [-]; optional `z_min_agl_m` [m above ground,
#'   fixed lowest reliable gate of the MLH retrieval; `NA` = per-file default].
read_measured_sites <- function(file) {
  ms <- readr::read_csv(
    file,
    na = c("NA", ""),
    col_types = readr::cols(
      provider_lat = readr::col_double(), provider_lon = readr::col_double(),
      distance_km = readr::col_double(), z_min_agl_m = readr::col_double(),
      .default = readr::col_character()
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

# ---- Doppler-lidar dissipation rate and mixing-layer height (spec 07 section 3) ----
# File structure (Cloudnet `epsilon-lidar`, Hyytiala, checked on all 688 files 2022-01..2024-04, spec 07 report):
# `epsilon` [m2 s-3] (height x time; missing where the signal is too low; no QC / SNR / precipitation field);
# `height` = gate centre [m above mean sea level] = (i + 1/2) * 30 m + `altitude`; `time` [hours since the UTC
# day]; `altitude` [m a.s.l., site altitude, equal to the Cloudnet /api/sites value]; `rolling_window_period`
# 180 s; `ray_accumulation_time` 6 s. Epsilon method: O'Connor et al. (2010), `OConnor2010`.

#' Options of the Doppler-lidar MLH retrieval
#'
#' Thresholds of the surface-connected mixing criterion of Manninen et al. (2018,
#' Sect. 3.2.1; `Manninen2018`) and of the half-hourly aggregation (spec 07).
#'
#' @param eps_surface Threshold on the dissipation rate [m2 s-3]: a gate is
#'   turbulent if epsilon >= `eps_surface` (Manninen et al. 2018: 1e-4).
#' @param z_min_agl Lowest reliable gate, gate centre [m above ground]; `NULL` =
#'   the lowest gate with any valid epsilon in each file (spec 07 default). The
#'   gate whose centre is closest to `z_min_agl` (within half a gate) is used.
#' @param min_coverage Minimum fraction [0-1] of the expected profiles of a
#'   half-hour that must be valid (measured or censored).
#' @param precip_classes Cloudnet `target_classification` values counted as
#'   precipitation (2 = drizzle or rain, 3 = drizzle or rain coexisting with
#'   droplets; Cloudnet classification product).
#' @param precip_top_agl Top of the layer checked for precipitation [m above
#'   ground]: a profile is flagged if any classification gate with centre at or
#'   below this height shows a precipitation class.
#' @param precip_window_s Half-width of the time window around each epsilon
#'   profile in which classification profiles are checked [s] (default: half the
#'   180 s epsilon rolling window).
#' @return Named list of options.
cloudnet_mlh_options <- function(eps_surface = 1e-4, z_min_agl = NULL, min_coverage = 0.5,
                                 precip_classes = c(2L, 3L), precip_top_agl = 300, precip_window_s = 90) {
  list(eps_surface = eps_surface, z_min_agl = z_min_agl, min_coverage = min_coverage,
       precip_classes = precip_classes, precip_top_agl = precip_top_agl, precip_window_s = precip_window_s)
}

#' Read a Cloudnet `epsilon-lidar` file
#'
#' Heights are converted from metres above mean sea level to metres above
#' ground: `height_agl = height_asl - ground_altitude_m`, with the file's own
#' `altitude` by default (spec 07, G 2026-10-23).
#'
#' @param path Path of the daily NetCDF file [character].
#' @param ground_altitude_m Ground altitude used for the conversion [m a.s.l.];
#'   `NULL` = the file's `altitude` variable.
#' @return List: `time` [POSIXct UTC, profile time], `height_asl`, `height_agl`
#'   [m, gate centres], `gate_m` [m, gate spacing], `eps` [m2 s-3, matrix gates x
#'   profiles, `NA` where missing], `altitude_m` [m a.s.l., from the file],
#'   `ground_altitude_m` [m a.s.l., used], `attrs` (file_uuid, instrument_pid,
#'   source_file_uuids, references, software versions, rolling_window_s,
#'   ray_accumulation_s, file name).
read_cloudnet_epsilon <- function(path, ground_altitude_m = NULL) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  miss <- setdiff(c("epsilon", "altitude"), names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {basename(path)}}: variable{?s} {.field {miss}} not found.")
  time <- cloudnet_nc_time(nc, path)
  h_std <- ncdf4::ncatt_get(nc, "height", "standard_name")$value
  if (!identical(h_std, "height_above_mean_sea_level")) {
    cli::cli_abort("{.path {basename(path)}}: height is {.val {h_std}}, expected height_above_mean_sea_level.")
  }
  z <- as.numeric(ncdf4::ncvar_get(nc, "height"))
  eps <- ncdf4::ncvar_get(nc, "epsilon", collapse_degen = FALSE)
  eps <- matrix(as.numeric(eps), nrow = length(z), ncol = length(time))
  eps[!is.finite(eps)] <- NA
  alt <- as.numeric(ncdf4::ncvar_get(nc, "altitude"))
  ground <- if (is.null(ground_altitude_m)) alt else ground_altitude_m
  ga <- ncdf4::ncatt_get(nc, 0)
  get1 <- function(v) if (v %in% names(nc$var)) as.numeric(ncdf4::ncvar_get(nc, v)) else NA_real_
  att <- function(k) if (is.null(ga[[k]])) NA_character_ else as.character(ga[[k]])
  list(
    time = time, height_asl = z, height_agl = z - ground, gate_m = stats::median(diff(z)), eps = eps,
    altitude_m = alt, ground_altitude_m = ground,
    attrs = list(file = basename(path), file_uuid = att("file_uuid"), instrument_pid = att("instrument_pid"),
                 source_file_uuids = att("source_file_uuids"), references = att("references"),
                 doppy_version = att("doppy_version"), cloudnetpy_version = att("cloudnetpy_version"),
                 rolling_window_s = get1("rolling_window_period"), ray_accumulation_s = get1("ray_accumulation_time"))
  )
}

# Profile times of a Cloudnet file: `time` in "hours since YYYY-MM-DD HH:MM:SS +00:00" -> POSIXct UTC.
cloudnet_nc_time <- function(nc, path) {
  t_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  m <- regmatches(t_units, regexec("^hours since ([0-9]{4}-[0-9]{2}-[0-9]{2})[ T]([0-9:]+)", t_units))[[1]]
  if (length(m) != 3) cli::cli_abort("{.path {basename(path)}}: unexpected time units {.val {t_units}}.")
  if (!grepl("(\\+00:00|Z|UTC)\\s*$", t_units) && grepl("[+-][0-9]{2}:[0-9]{2}\\s*$", t_units)) {
    cli::cli_abort("{.path {basename(path)}}: time units {.val {t_units}} not in UTC.")
  }
  t0 <- as.POSIXct(paste(m[2], m[3]), tz = "UTC")
  t0 + as.numeric(ncdf4::ncvar_get(nc, "time")) * 3600
}

#' Surface-connected mixing-layer height per profile
#'
#' Manninen et al. (2018, Sect. 3.2.1; `Manninen2018`), surface-connected
#' turbulent mixing only (the cloud-driven, top-down step is out of scope): from
#' the lowest reliable gate upward while epsilon >= `eps_surface`; MLH = top of
#' the last gate meeting the condition (gate centre + half the gate spacing).
#' Profile status (exclusive; precipitation is applied later, see
#' [cloudnet_precip_flag()]):
#' - `no_signal`: no valid epsilon at any gate (MLH `NA`, not censored);
#' - `lowest_gate_missing`: epsilon missing at the lowest reliable gate but valid
#'   somewhere above (censored below z_min, counted separately);
#' - `censored_weak`: epsilon < `eps_surface` at the lowest reliable gate
#'   (censored below z_min);
#' - `measured`: MLH found; `gap_in_column` = the search stopped at a missing
#'   gate (low signal) rather than at epsilon < `eps_surface`.
#' Gates below the lowest reliable gate are ignored.
#'
#' @param eps Dissipation rate, matrix gates x profiles [m2 s-3].
#' @param z_agl Gate centres [m above ground], increasing.
#' @param gate_m Gate spacing [m].
#' @param opts Output of [cloudnet_mlh_options()].
#' @return Tibble, one row per profile: `status`, `mlh_m` [m above ground, top of
#'   the last turbulent gate; `NA` unless `measured`], `censored` [logical],
#'   `gap_in_column` [logical], `z_min_m` [m above ground, centre of the lowest
#'   reliable gate].
cloudnet_mlh_profiles <- function(eps, z_agl, gate_m, opts) {
  n_t <- ncol(eps)
  n_z <- nrow(eps)
  any_valid <- colSums(!is.na(eps)) > 0
  if (is.null(opts$z_min_agl)) {
    rows_valid <- which(rowSums(!is.na(eps)) > 0)
    i0 <- if (length(rows_valid) > 0) min(rows_valid) else NA_integer_
  } else {
    d <- abs(z_agl - opts$z_min_agl)
    i0 <- which.min(d)
    if (d[i0] > gate_m / 2) {
      cli::cli_abort("No gate within half a gate ({gate_m / 2} m) of z_min_agl = {opts$z_min_agl} m.")
    }
  }
  status <- rep("no_signal", n_t)
  mlh <- rep(NA_real_, n_t)
  gap <- rep(FALSE, n_t)
  if (!is.na(i0) && n_t > 0) {
    e <- eps[i0:n_z, , drop = FALSE]
    ok <- !is.na(e) & e >= opts$eps_surface
    # number of consecutive turbulent gates from the lowest reliable gate upward
    n_ok <- colSums(matrix(apply(ok, 2, cumprod), nrow = nrow(ok)))
    low <- e[1, ]
    status[any_valid & is.na(low)] <- "lowest_gate_missing"
    status[any_valid & !is.na(low) & low < opts$eps_surface] <- "censored_weak"
    meas <- any_valid & !is.na(low) & low >= opts$eps_surface
    status[meas] <- "measured"
    top_idx <- i0 + n_ok - 1
    mlh[meas] <- z_agl[top_idx[meas]] + gate_m / 2
    nxt <- top_idx + 1
    chk <- meas & nxt <= n_z
    gap[chk] <- is.na(eps[cbind(nxt[chk], which(chk))])
  }
  dplyr::tibble(status = status, mlh_m = mlh,
                censored = status %in% c("censored_weak", "lowest_gate_missing"),
                gap_in_column = gap, z_min_m = if (is.na(i0)) NA_real_ else z_agl[i0])
}

#' Read a Cloudnet `classification` file (precipitation and cloud base)
#'
#' @param path Path of the daily NetCDF file [character].
#' @param ground_altitude_m Ground altitude [m a.s.l.]; `NULL` = the file's
#'   `altitude`.
#' @return List: `time` [POSIXct UTC], `height_agl` [m, gate centres above
#'   ground], `target_classification` [matrix gates x profiles, Cloudnet classes
#'   0-10], `cloud_base_agl_m` [m above ground, `NA` = no cloud], `rain_detected`
#'   [logical or `NA` if absent], `altitude_m` [m a.s.l.], `attrs` (file, file_uuid,
#'   pid, cloudnetpy_version, source).
read_cloudnet_classification <- function(path, ground_altitude_m = NULL) {
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  miss <- setdiff(c("target_classification", "altitude"), names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {basename(path)}}: variable{?s} {.field {miss}} not found.")
  time <- cloudnet_nc_time(nc, path)
  h_std <- ncdf4::ncatt_get(nc, "height", "standard_name")$value
  if (!identical(h_std, "height_above_mean_sea_level")) {
    cli::cli_abort("{.path {basename(path)}}: height is {.val {h_std}}, expected height_above_mean_sea_level.")
  }
  z <- as.numeric(ncdf4::ncvar_get(nc, "height"))
  alt <- stats::median(as.numeric(ncdf4::ncvar_get(nc, "altitude")), na.rm = TRUE)
  ground <- if (is.null(ground_altitude_m)) alt else ground_altitude_m
  tc <- ncdf4::ncvar_get(nc, "target_classification", collapse_degen = FALSE)
  tc <- matrix(as.integer(tc), nrow = length(z), ncol = length(time))
  cb <- if ("cloud_base_height_amsl" %in% names(nc$var)) {
    as.numeric(ncdf4::ncvar_get(nc, "cloud_base_height_amsl")) - ground
  } else {
    rep(NA_real_, length(time))
  }
  cb[!is.finite(cb)] <- NA
  rd <- if ("rain_detected" %in% names(nc$var)) {
    as.integer(ncdf4::ncvar_get(nc, "rain_detected")) == 1L
  } else {
    rep(NA, length(time))
  }
  ga <- ncdf4::ncatt_get(nc, 0)
  att <- function(k) if (is.null(ga[[k]])) NA_character_ else as.character(ga[[k]])
  list(time = time, height_agl = z - ground, target_classification = tc, cloud_base_agl_m = cb,
       rain_detected = rd, altitude_m = alt,
       attrs = list(file = basename(path), file_uuid = att("file_uuid"), pid = att("pid"),
                    cloudnetpy_version = att("cloudnetpy_version"), source = att("source")))
}

#' Precipitation flag and cloud base for epsilon profiles
#'
#' A profile is precipitation-flagged if any classification profile within
#' `opts$precip_window_s` of the profile time shows a class in
#' `opts$precip_classes` (drizzle / rain) in a gate with centre at or below
#' `opts$precip_top_agl` (spec 07, G 2026-10-23). The cloud base is that of the
#' nearest classification profile within the window.
#'
#' @param time Epsilon profile times [POSIXct UTC].
#' @param cls Output of [read_cloudnet_classification()], or `NULL` (no file:
#'   no precipitation check).
#' @param opts Output of [cloudnet_mlh_options()].
#' @return Tibble, one row per profile: `precip` [logical, `NA` = not checked],
#'   `rain_detected` [logical, any `rain_detected` in the window; `NA` = not
#'   available], `cloud_base_agl_m` [m above ground, `NA` = no cloud or not
#'   checked].
cloudnet_precip_flag <- function(time, cls, opts) {
  n <- length(time)
  out <- dplyr::tibble(precip = rep(NA, n), rain_detected = rep(NA, n), cloud_base_agl_m = rep(NA_real_, n))
  if (is.null(cls) || length(cls$time) == 0 || n == 0) return(out)
  ord <- order(cls$time)
  ct <- as.numeric(cls$time[ord])
  low <- cls$height_agl <= opts$precip_top_agl
  tc <- cls$target_classification[low, ord, drop = FALSE]
  p_c <- colSums(matrix(tc %in% opts$precip_classes, nrow = nrow(tc))) > 0
  tt <- as.numeric(time)
  lo <- findInterval(tt - opts$precip_window_s, ct, left.open = TRUE) + 1L
  hi <- findInterval(tt + opts$precip_window_s, ct)
  has <- hi >= lo
  cs <- c(0, cumsum(p_c))
  out$precip[has] <- (cs[hi[has] + 1L] - cs[lo[has]]) > 0
  rd <- cls$rain_detected[ord]
  if (!all(is.na(rd))) {
    cr <- c(0, cumsum(rd %in% TRUE))
    out$rain_detected[has] <- (cr[hi[has] + 1L] - cr[lo[has]]) > 0
  }
  # nearest classification profile within the window
  j <- findInterval(tt, ct)
  j_lo <- pmax(j, 1L)
  j_hi <- pmin(j + 1L, length(ct))
  near <- ifelse(abs(ct[j_hi] - tt) < abs(ct[j_lo] - tt), j_hi, j_lo)
  in_win <- abs(ct[near] - tt) <= opts$precip_window_s
  out$cloud_base_agl_m[in_win] <- cls$cloud_base_agl_m[ord][near[in_win]]
  out
}

#' Expected number of profiles per half-hour in a daily file
#'
#' The Doppler lidar alternates vertical stares (about 14-15 min) with other
#' scans every 30 min, so the nominal count per half-hour is the stare count, not
#' 1800 s / profile spacing (spec 07, G 2026-10-23): median number of profiles
#' per UTC half-hour, over the half-hours of the file with at least one profile.
#'
#' @param time Profile times of one file [POSIXct UTC].
#' @return Expected profiles per half-hour [-]; `NA` if `time` is empty.
cloudnet_expected_profiles <- function(time) {
  if (length(time) == 0) return(NA_real_)
  stats::median(as.numeric(table(floor(as.numeric(time) / 1800))))
}

#' Profile-level MLH of one day
#'
#' Reads one `epsilon-lidar` file (and the `classification` file of the same day
#' if given), computes the profile MLH ([cloudnet_mlh_profiles()]) and the
#' precipitation flag ([cloudnet_precip_flag()]). Precipitation overrides the
#' epsilon status: `status = "precipitation"`, MLH `NA`, not censored.
#'
#' @param eps_path Path of the `epsilon-lidar` file [character].
#' @param cls_path Path of the `classification` file, or `NULL`/`NA` (none).
#' @param opts Output of [cloudnet_mlh_options()].
#' @param ground_altitude_m Ground altitude [m a.s.l.]; `NULL` = file altitude.
#' @return Tibble, one row per profile: `time` [POSIXct UTC], `status`, `mlh_m`
#'   [m above ground], `censored`, `gap_in_column` [logical], `z_min_m` [m above
#'   ground], `precip` [logical, `NA` = not checked], `rain_detected` [logical],
#'   `cloud_base_agl_m` [m above ground], `n_expected` [-, expected profiles per
#'   half-hour of the file], `dt_s` [s, median profile spacing of the file],
#'   `file_uuid`. Attributes `eps_attrs`, `cls_attrs`, `altitude_m`, `gate_m`.
cloudnet_mlh_day <- function(eps_path, cls_path = NULL, opts = cloudnet_mlh_options(), ground_altitude_m = NULL) {
  e <- read_cloudnet_epsilon(eps_path, ground_altitude_m)
  p <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, opts)
  has_cls <- !is.null(cls_path) && !is.na(cls_path)
  cls <- if (has_cls) read_cloudnet_classification(cls_path, ground_altitude_m) else NULL
  pr <- cloudnet_precip_flag(e$time, cls, opts)
  out <- dplyr::bind_cols(dplyr::tibble(time = e$time), p, pr)
  rain <- out$precip %in% TRUE
  out$status[rain] <- "precipitation"
  out$mlh_m[rain] <- NA_real_
  out$censored[rain] <- FALSE
  out$gap_in_column[rain] <- FALSE
  out$n_expected <- cloudnet_expected_profiles(e$time)
  out$dt_s <- if (length(e$time) > 1) stats::median(diff(as.numeric(e$time))) else NA_real_
  out$file_uuid <- e$attrs$file_uuid
  attr(out, "eps_attrs") <- e$attrs
  attr(out, "cls_attrs") <- if (has_cls) cls$attrs else NULL
  attr(out, "altitude_m") <- e$altitude_m
  attr(out, "gate_m") <- e$gate_m
  out
}

#' Aggregate profile MLH to half-hours
#'
#' Pure function (spec 07 section 3, G 2026-10-23). Per UTC half-hour
#' [start, start + 30 min): valid profiles = measured + censored (`censored_weak`,
#' `lowest_gate_missing`); `no_signal` and `precipitation` profiles are not valid.
#' Coverage = min(1, valid / expected), expected = `n_expected` of the file
#' ([cloudnet_expected_profiles()]). Rules:
#' - coverage < `opts$min_coverage`: `pblh_m = NA`, `qc_flag = NA`;
#' - more than half of the valid profiles censored: `pblh_m = NA`, `qc_flag = 1`;
#' - otherwise `pblh_m` = median MLH of the measured profiles; `qc_flag = 0` if no
#'   censored profile, `2` if some.
#'
#' @param prof Profile table as [cloudnet_mlh_day()] (several days allowed).
#' @param opts Output of [cloudnet_mlh_options()].
#' @return Tibble, one row per half-hour with at least one profile:
#'   `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `pblh_m` [m above
#'   ground], `qc_flag`, `n_samples` [-, valid profiles], and `x_` columns:
#'   `x_n_expected`, `x_coverage` [0-1], `x_n_measured`, `x_n_censored_weak`,
#'   `x_n_lowest_gate_missing`, `x_n_no_signal`, `x_n_precipitation`,
#'   `x_n_gap_in_column` [-, measured profiles whose search stopped at a missing
#'   gate], `x_censor_below_m` [m above ground, z_min], `x_precip_checked`
#'   [logical, classification available for all profiles], `x_stare_seconds` [s,
#'   valid profiles x profile spacing], `x_first_profile_utc`,
#'   `x_last_profile_utc` [ISO 8601 UTC, first / last valid profile],
#'   `x_cloud_base_agl_m` [m above ground, median over profiles with a cloud base],
#'   `x_cloud_fraction` [0-1, share of profiles with a cloud base, among checked].
cloudnet_mlh_halfhour <- function(prof, opts) {
  hh <- floor(as.numeric(prof$time) / 1800) * 1800
  valid <- prof$status %in% c("measured", "censored_weak", "lowest_gate_missing")
  g <- split(seq_len(nrow(prof)), hh)
  agg <- lapply(g, function(i) {
    st <- prof$status[i]
    v <- valid[i]
    n_valid <- sum(v)
    n_cens <- sum(prof$censored[i] & v)
    n_exp <- stats::median(prof$n_expected[i])
    cov <- min(1, n_valid / n_exp)
    m <- prof$mlh_m[i][st == "measured"]
    if (is.na(cov) || cov < opts$min_coverage) {
      pblh <- NA_real_
      qc <- NA_integer_
    } else if (n_cens > n_valid / 2) {
      pblh <- NA_real_
      qc <- 1L
    } else {
      pblh <- stats::median(m)
      qc <- if (n_cens == 0) 0L else 2L
    }
    tv <- prof$time[i][v]
    cb <- prof$cloud_base_agl_m[i]
    chk <- !is.na(prof$precip[i])
    data.frame(
      pblh_m = pblh, qc_flag = qc, n_samples = as.integer(n_valid), x_n_expected = n_exp, x_coverage = cov,
      x_n_measured = sum(st == "measured"), x_n_censored_weak = sum(st == "censored_weak"),
      x_n_lowest_gate_missing = sum(st == "lowest_gate_missing"), x_n_no_signal = sum(st == "no_signal"),
      x_n_precipitation = sum(st == "precipitation"),
      x_n_gap_in_column = sum(prof$gap_in_column[i] & st == "measured"),
      x_censor_below_m = stats::median(prof$z_min_m[i]), x_precip_checked = all(chk),
      x_stare_seconds = n_valid * stats::median(prof$dt_s[i]),
      x_first_profile_utc = if (n_valid > 0) min(as.numeric(tv)) else NA_real_,
      x_last_profile_utc = if (n_valid > 0) max(as.numeric(tv)) else NA_real_,
      x_cloud_base_agl_m = if (any(!is.na(cb))) stats::median(cb, na.rm = TRUE) else NA_real_,
      x_cloud_fraction = if (any(chk)) mean(!is.na(cb[chk])) else NA_real_
    )
  })
  out <- dplyr::bind_rows(agg)
  start <- as.POSIXct(as.numeric(names(g)), origin = "1970-01-01", tz = "UTC")
  iso <- function(x) ifelse(is.na(x), NA_character_,
                            format(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"))
  out$x_first_profile_utc <- iso(out$x_first_profile_utc)
  out$x_last_profile_utc <- iso(out$x_last_profile_utc)
  dplyr::bind_cols(dplyr::tibble(timestamp_start_utc = start, timestamp_end_utc = start + 1800), out)
}

#' Standardize the Cloudnet epsilon MLH (source `cloudnet_eps_mlh`)
#'
#' @param hh Output of [cloudnet_mlh_halfhour()].
#' @param site_id ICOS station code [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h Site UTC offset [h].
#' @param info List with `instrument` (name and PID), `citation`,
#'   `processing_history`, `instrument_altitude_m` [m a.s.l.] [character/numeric].
#' @return List: `data` (standard table), `meta` (sidecar list).
standardize_cloudnet_mlh <- function(hh, site_id, exclusions, utc_offset_h, info) {
  source_id <- "cloudnet_eps_mlh"
  d <- flag_exclusions(hh, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, source_id)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL
  data <- dplyr::bind_cols(
    dplyr::tibble(timestamp_start_utc = d$timestamp_start_utc, timestamp_end_utc = d$timestamp_end_utc,
                  pblh_m = d$pblh_m, source_id = source_id, qc_flag = d$qc_flag, uncertainty_m = NA_real_,
                  n_samples = d$n_samples),
    d[, grep("^x_", names(d))]
  )
  meta <- standard_meta(
    site_id, source_id, "turbulence",
    instrument = info$instrument,
    retrieval_algorithm = paste(
      "surface-connected eps >= 1e-4 m2 s-3 from the lowest reliable gate (Manninen et al. 2018, Sect. 3.2.1);",
      "MLH = top of the last turbulent gate; cloud-driven (top-down) mixing not included"
    ),
    height_reference = "agl",
    instrument_altitude_m = info$instrument_altitude_m,
    time_reference = "interval_mean",
    original_resolution = info$original_resolution,
    original_qc = paste(
      "Cloudnet epsilon-lidar: epsilon missing where the signal is too low (no QC/SNR/precipitation field);",
      "precipitation from the Cloudnet classification product"
    ),
    qc_flag_meaning = list(`0` = "all valid profiles measured", `1` = "censored below z_min (more than half of the valid profiles)",
                           `2` = "partly censored (median of the measured profiles)",
                           `NA` = "coverage below min_coverage, or no valid profile"),
    processing_history = paste0(info$processing_history, " Rows in exclusion periods dropped (", ex$n_excluded,
                                "); standardized by standardize_cloudnet_mlh() (spec 07)."),
    licence = "CC BY 4.0 (ACTRIS/Cloudnet)",
    citation = info$citation
  )
  list(data = data, meta = meta)
}
