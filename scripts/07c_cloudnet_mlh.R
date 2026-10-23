# Spec 07 section 3: turbulence-based mixing-layer height from the Cloudnet Doppler-lidar dissipation rate
# (Manninen et al. 2018, surface-connected epsilon >= 1e-4 m2 s-3) for one ICOS station, written as the standard
# source `cloudnet_eps_mlh`.
# Inputs: the epsilon-lidar and classification files downloaded by scripts/07b_cloudnet_download.R into
# <cloudnet_dir> (with their provenance tables). Lowest reliable gate: `z_min_agl_m` of config/measured_sites.csv
# (FI-Hyy: 105 m agl, G 2026-10-23), else the lowest valid gate of each file. Heights above ground = height a.s.l.
# minus the altitude in the Cloudnet file.
# Output:
#   <derived_dir>/cloudnet/<site>/<site>_cloudnet_eps_mlh_profiles.csv.gz  (one row per lidar profile)
#   <standard_dir>/<site>/<site>__cloudnet_eps_mlh.csv.gz + .yml           (half-hours, standard format)
#   <standard_dir>/<site>/<site>__cloudnet_eps_mlh_files.csv                (uuid / PID of every input file)
#   <output_dir>/tables/07c_profile_status_<site>.csv, 07c_censoring_<site>.csv, 07c_altitude_sensitivity_<site>.csv
# Reports: profile status counts (measured, censored weak, lowest gate missing, no signal, precipitation), half-hour
# QC counts and coverage, censoring by hour x season, precipitation flag vs `rain_detected`, and the sensitivity to
# the ground altitude (ICOS site altitude instead of the Cloudnet altitude).
#
# Usage: Rscript scripts/07c_cloudnet_mlh.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

actris_ack <- paste(
  "The author(s) wish to acknowledge the ACTRIS (Aerosol, Clouds, and Trace Gases Research Infrastructure) for its",
  "essential role in providing access to high-quality data, facilities, and technical support that were integral",
  "to the execution of this research."
)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FI-Hyy}.")
sites <- read_sites(resolve_path(cfg$sites_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "cloudnet" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) Cloudnet pairing in {.file measured_sites.csv}.")
meta_site <- sites[sites$site_id == site_id, ]
cn_dir <- resolve_path(cfg$cloudnet_dir, site = site_id)
derived_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "cloudnet", site_id)
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
for (d in c(derived_dir, std_dir, tables_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

opts <- cloudnet_mlh_options(z_min_agl = if (is.na(ms$z_min_agl_m)) NULL else ms$z_min_agl_m)

# --- Input files (provenance tables written by script 07b) ---
read_prov <- function(product) {
  f <- file.path(cn_dir, paste0("cloudnet_files_provenance_", product, ".csv"))
  if (!file.exists(f)) return(NULL)
  p <- readr::read_csv(f, col_types = readr::cols(.default = "c"), progress = FALSE)
  p <- p[p$status %in% c("downloaded", "present"), ]
  p$path <- file.path(cn_dir, p$filename)
  p[file.exists(p$path), ]
}
prov_eps <- read_prov("epsilon-lidar")
prov_cls <- read_prov("classification")
if (is.null(prov_eps) || nrow(prov_eps) == 0) cli::cli_abort("No epsilon-lidar files in {.path {cn_dir}} (run 07b).")
if (is.null(prov_cls)) prov_cls <- prov_eps[0, ]
if (anyDuplicated(prov_eps$measurement_date)) cli::cli_abort("More than one epsilon-lidar file on some days.")
days <- sort(prov_eps$measurement_date)
cls_path <- stats::setNames(prov_cls$path, prov_cls$measurement_date)

cli::cli_h1("Cloudnet epsilon MLH: {ms$provider_site} for {site_id}")
cli::cli_inform(c(
  "{length(days)} epsilon-lidar day{?s} ({days[1]} to {days[length(days)]}); classification on {sum(days %in% names(cls_path))}.",
  "Options: eps_surface = {opts$eps_surface} m2 s-3; z_min = {if (is.null(opts$z_min_agl)) 'lowest valid gate per file' else paste(opts$z_min_agl, 'm agl')}; min_coverage = {opts$min_coverage}; precipitation classes {paste(opts$precip_classes, collapse = ', ')} at or below {opts$precip_top_agl} m agl within +/- {opts$precip_window_s} s."
))

# --- Profiles (and the ground-altitude sensitivity on the same files) ---
alt_icos <- meta_site$altitude_m
prof <- vector("list", length(days))
sens <- vector("list", length(days))
eps_attrs <- vector("list", length(days))
for (k in seq_along(days)) {
  dd <- days[k]
  pe <- prov_eps$path[prov_eps$measurement_date == dd]
  pc <- if (dd %in% names(cls_path)) cls_path[[dd]] else NULL
  p <- cloudnet_mlh_day(pe, pc, opts)
  prof[[k]] <- p
  eps_attrs[[k]] <- c(attr(p, "eps_attrs"), altitude_m = attr(p, "altitude_m"), gate_m = attr(p, "gate_m"))
  # Sensitivity: heights above ground with the ICOS altitude. (a) z_min follows the same gate (z_min - offset);
  # (b) z_min stays at the same height above ground (another gate). Precipitation flags from the base run.
  if (!is.na(alt_icos) && !is.null(opts$z_min_agl)) {
    e <- read_cloudnet_epsilon(pe, ground_altitude_m = alt_icos)
    off <- alt_icos - e$altitude_m
    pa <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, utils::modifyList(opts, list(z_min_agl = opts$z_min_agl - off)))
    pb <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, opts)
    fix <- function(x) {
      rain <- p$status == "precipitation"
      x$status[rain] <- "precipitation"
      x$mlh_m[rain] <- NA_real_
      x$censored[rain] <- FALSE
      dplyr::bind_cols(p[, c("time", "precip", "rain_detected", "cloud_base_agl_m", "n_expected", "dt_s", "file_uuid")],
                       x[, c("status", "mlh_m", "censored", "gap_in_column", "z_min_m")])
    }
    sens[[k]] <- list(a = fix(pa), b = fix(pb))
  }
  if (k %% 100 == 0) cli::cli_inform("{k} / {length(days)} days processed.")
}
prof <- dplyr::bind_rows(prof)
prof$qc_na <- prof$status %in% c("no_signal", "precipitation")
ea <- dplyr::bind_rows(lapply(eps_attrs, function(a) as.data.frame(a[!vapply(a, is.null, logical(1))])))

prof_file <- file.path(derived_dir, paste0(site_id, "_cloudnet_eps_mlh_profiles.csv.gz"))
readr::write_csv(dplyr::mutate(prof[, c("time", "status", "mlh_m", "censored", "gap_in_column", "qc_na", "z_min_m",
                                       "precip", "rain_detected", "cloud_base_agl_m", "n_expected", "dt_s",
                                       "file_uuid")],
                               time = format(.data$time, "%Y-%m-%dT%H:%M:%OS2Z", tz = "UTC")),
                 prof_file)
cli::cli_alert_success("Profiles written: {.path {prof_file}} ({nrow(prof)} rows).")

# --- Half-hours and standard file ---
hh <- cloudnet_mlh_halfhour(prof, opts)
inst <- unique(prov_eps[, c("instrument_name", "instrument_pid")])
z_used <- sort(unique(prof$z_min_m))
files_tab <- dplyr::bind_rows(
  dplyr::mutate(prov_eps, product = "epsilon-lidar"), dplyr::mutate(prov_cls, product = "classification")
)[, c("product", "measurement_date", "filename", "uuid", "pid", "version", "volatile", "checksum_sha256",
      "instrument_name", "instrument_pid", "download_url", "downloaded_utc")]
files_file <- file.path(std_dir, paste0(site_id, "__cloudnet_eps_mlh_files.csv"))
readr::write_csv(files_tab, files_file)
dl_dates <- range(as.Date(substr(files_tab$downloaded_utc, 1, 10)), na.rm = TRUE)
rng <- function(x) paste(sort(unique(stats::na.omit(x))), collapse = ", ")
info <- list(
  instrument = paste0(paste0(inst$instrument_name, " (HALO Photonics StreamLine Doppler lidar; PID ",
                             inst$instrument_pid, ")"), collapse = "; "),
  instrument_altitude_m = paste(rng(ea$altitude_m), collapse = ", "),
  original_resolution = paste0(
    "vertical-stare profiles every ", round(stats::median(prof$dt_s, na.rm = TRUE), 2), " s (epsilon over a ",
    rng(ea$rolling_window_s), " s rolling window, ", rng(ea$ray_accumulation_s), " s rays), about 14-15 min of ",
    "stare per 30 min; ", rng(ea$gate_m), " m gates"
  ),
  citation = paste0(
    "ACTRIS Cloudnet data, site ", ms$provider_site, ": epsilon-lidar (", nrow(prov_eps), " daily files) and ",
    "classification (", nrow(prov_cls), " daily files), instrument ", paste(inst$instrument_name, inst$instrument_pid),
    "; file uuids and PIDs in ", basename(files_file), "; downloaded ", paste(unique(dl_dates), collapse = " to "),
    " from https://cloudnet.fmi.fi; licence CC BY 4.0. ",
    "Epsilon: O'Connor et al. (2010), J. Atmos. Oceanic Technol. 27, 1652-1664, doi:10.1175/2010JTECHA1455.1 ",
    "(OConnor2010). MLH: Manninen et al. (2018) (Manninen2018). ", actris_ack
  ),
  processing_history = paste0(
    "Cloudnet epsilon-lidar files (doppy ", rng(ea$doppy_version), ", CloudnetPy ", rng(ea$cloudnetpy_version),
    ") read with read_cloudnet_epsilon(); heights above ground = height a.s.l. - altitude in the Cloudnet file (",
    rng(ea$altitude_m), " m; ICOS site altitude ", alt_icos, " m, see the sensitivity in the spec 07 report); ",
    "gate top = centre + ", rng(ea$gate_m / 2), " m. Profile MLH (cloudnet_mlh_profiles()): eps_surface = ",
    opts$eps_surface, " m2 s-3, z_min = ", paste(z_used, collapse = ", "), " m agl (gate centre",
    if (!is.null(opts$z_min_agl)) ", fixed for the site in config/measured_sites.csv" else ", per file", "); ",
    "status no_signal (no valid epsilon) and precipitation (Cloudnet classification classes ",
    paste(opts$precip_classes, collapse = ", "), " at or below ", opts$precip_top_agl, " m agl within +/- ",
    opts$precip_window_s, " s) are NA; censored = epsilon below threshold or missing at z_min. ",
    "Half-hours (cloudnet_mlh_halfhour()): expected profiles = median count per half-hour of each daily file, ",
    "coverage = valid / expected (capped at 1), min_coverage = ", opts$min_coverage,
    "; censored if more than half the valid profiles are censored; else median of the measured profiles. ",
    sum(!days %in% names(cls_path)), " day(s) without classification: no precipitation check (x_precip_checked). ",
    "Profile-level output kept in derived/cloudnet (", basename(prof_file), ")."
  )
)
std <- standardize_cloudnet_mlh(hh, site_id, exclusions, meta_site$utc_offset_h, info)
paths <- write_standard_pblh(std$data, std$meta, std_dir)
cli::cli_alert_success("Standard file: {.path {paths[['data']]}} ({nrow(std$data)} half-hours); files table {.path {files_file}}.")

# --- Report ---
season_hour <- function(t) {
  data.frame(season = season_of(t), hour_utc = as.integer(format(t, "%H", tz = "UTC")))
}
cli::cli_h2("Profiles by status")
st_levels <- c("measured", "censored_weak", "lowest_gate_missing", "no_signal", "precipitation")
tab <- as.data.frame(table(status = factor(prof$status, st_levels)))
tab$pct <- round(100 * tab$Freq / sum(tab$Freq), 2)
print_tbl(tab)
cli::cli_inform("Measured profiles that stopped at a missing gate (gap in the column): {sum(prof$gap_in_column)} ({round(100 * mean(prof$gap_in_column[prof$status == 'measured']), 1)} % of measured). Profiles without precipitation check: {sum(is.na(prof$precip))}.")
ps <- cbind(season_hour(prof$time), status = factor(prof$status, st_levels))
st_tab <- as.data.frame(table(ps))
readr::write_csv(st_tab, file.path(tables_dir, paste0("07c_profile_status_", site_id, ".csv")))
st_season <- as.data.frame.matrix(table(ps$season, ps$status))
print(round(100 * st_season / rowSums(st_season), 1))

cli::cli_h2("Precipitation flag (classification) vs rain_detected (profiles with both)")
both <- !is.na(prof$precip) & !is.na(prof$rain_detected)
print(table(precip = prof$precip[both], rain_detected = prof$rain_detected[both]))

cli::cli_h2("Half-hours")
hq <- std$data
cli::cli_inform(c(
  "{nrow(hq)} half-hours with at least one profile; coverage: median {round(stats::median(hq$x_coverage), 2)}, below {opts$min_coverage}: {sum(is.na(hq$qc_flag))}.",
  "Expected profiles per half-hour: median {stats::median(hq$x_n_expected)} (range {paste(range(hq$x_n_expected), collapse = '-')}).",
  "Half-hours without precipitation check (no classification file): {sum(!hq$x_precip_checked)}."
))
print(table(qc_flag = hq$qc_flag, useNA = "ifany"))
cli::cli_inform("MLH [m agl], qc 0 or 2:")
print(round(stats::quantile(hq$pblh_m, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)))
cens <- cbind(season_hour(hq$timestamp_start_utc), cens = hq$qc_flag %in% 1L, valid = !is.na(hq$qc_flag))
cens_tab <- dplyr::summarise(dplyr::group_by(cens, .data$season, .data$hour_utc),
                             n_valid = sum(.data$valid), n_censored = sum(.data$cens),
                             pct_censored = round(100 * sum(.data$cens) / max(1, sum(.data$valid)), 1), .groups = "drop")
readr::write_csv(cens_tab, file.path(tables_dir, paste0("07c_censoring_", site_id, ".csv")))
cli::cli_inform("Share of valid half-hours censored below z_min [%], season x hour (UTC):")
print(with(cens_tab, tapply(pct_censored, list(hour_utc = hour_utc, season = season), sum)))
cli::cli_inform("Median MLH [m agl] (qc 0/2), season x hour (UTC):")
mm <- cbind(season_hour(hq$timestamp_start_utc), pblh = hq$pblh_m)
print(round(with(mm, tapply(pblh, list(hour_utc = hour_utc, season = season), stats::median, na.rm = TRUE))))

if (!all(vapply(sens, is.null, logical(1)))) {
  cli::cli_h2("Sensitivity: ground altitude {alt_icos} m (ICOS) instead of the Cloudnet altitude")
  cmp <- function(v, label) {
    pv <- dplyr::bind_rows(lapply(sens, `[[`, v))
    hv <- cloudnet_mlh_halfhour(pv, utils::modifyList(opts, list(z_min_agl = stats::median(pv$z_min_m, na.rm = TRUE))))
    m <- merge(hh[, c("timestamp_start_utc", "pblh_m", "qc_flag")], hv[, c("timestamp_start_utc", "pblh_m", "qc_flag")],
               by = "timestamp_start_utc", suffixes = c("", "_s"))
    qc_chg <- !identical_na(m$qc_flag, m$qc_flag_s)
    both_v <- !is.na(m$pblh_m) & !is.na(m$pblh_m_s)
    d <- m$pblh_m_s[both_v] - m$pblh_m[both_v]
    data.frame(variant = label, z_min_m = stats::median(pv$z_min_m, na.rm = TRUE), n_halfhours = nrow(m),
               pct_qc_changed = round(100 * mean(qc_chg), 2),
               pct_censoring_changed = round(100 * mean(xor(m$qc_flag %in% 1L, m$qc_flag_s %in% 1L)), 2),
               pct_mlh_changed = round(100 * mean(abs(d) > 1e-6), 2), median_mlh_diff_m = stats::median(d),
               pct_profiles_status_changed = round(100 * mean(pv$status != prof$status), 2))
  }
  identical_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  sens_tab <- rbind(cmp("a", "same gate as z_min (z_min shifted)"), cmp("b", "same z_min above ground (other gate)"))
  print_tbl(sens_tab)
  readr::write_csv(sens_tab, file.path(tables_dir, paste0("07c_altitude_sensitivity_", site_id, ".csv")))
}
