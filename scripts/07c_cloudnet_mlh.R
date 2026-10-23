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
  eps_attrs[[k]] <- c(attr(p, "eps_attrs"), altitude_m = attr(p, "altitude_m"), gate_m = attr(p, "gate_m"))
  # Same epsilon-only profile with a different retrieval, precipitation flags of the base run.
  redo <- function(x) {
    y <- p[, c("time", "precip", "rain_detected", "cloud_base_agl_m", "n_expected", "dt_s", "file_uuid")]
    y$status_eps <- x$status
    y$mlh_eps_m <- x$mlh_m
    y$gap_eps <- x$gap_in_column
    y$z_min_m <- x$z_min_m
    cloudnet_apply_precip(y, p$precip)
  }
  # Sensitivity: bridge a single missing gate when the next gate is turbulent.
  e <- read_cloudnet_epsilon(pe)
  pbr <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, utils::modifyList(opts, list(bridge_single_gap = TRUE)))
  p$status_bridged <- pbr$status
  p$mlh_bridged_m <- pbr$mlh_m
  sens[[k]] <- list(bridge = redo(pbr))
  # Sensitivity: heights above ground with the ICOS altitude. (a) z_min follows the same gate (z_min - offset);
  # (b) z_min stays at the same height above ground (another gate).
  if (!is.na(alt_icos) && !is.null(opts$z_min_agl)) {
    e <- read_cloudnet_epsilon(pe, ground_altitude_m = alt_icos)
    off <- alt_icos - e$altitude_m
    pa <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, utils::modifyList(opts, list(z_min_agl = opts$z_min_agl - off)))
    pb <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, opts)
    sens[[k]]$a <- redo(pa)
    sens[[k]]$b <- redo(pb)
  }
  prof[[k]] <- p
  if (k %% 100 == 0) cli::cli_inform("{k} / {length(days)} days processed.")
}
prof <- dplyr::bind_rows(prof)
prof$qc_na <- prof$status %in% c("no_signal", "precipitation")
ea <- dplyr::bind_rows(lapply(eps_attrs, function(a) as.data.frame(a[!vapply(a, is.null, logical(1))])))

prof_file <- file.path(derived_dir, paste0(site_id, "_cloudnet_eps_mlh_profiles.csv.gz"))
readr::write_csv(dplyr::mutate(prof[, c("time", "status", "mlh_m", "censored", "gap_in_column", "qc_na", "z_min_m",
                                       "status_eps", "mlh_eps_m", "eps_top_ratio", "precip", "precip_category",
                                       "precip_column", "precip_strict", "rain_detected", "cloud_base_agl_m",
                                       "n_expected", "dt_s", "file_uuid")],
                               time = format(.data$time, "%Y-%m-%dT%H:%M:%OS2Z", tz = "UTC")),
                 prof_file)
cli::cli_alert_success("Profiles written: {.path {prof_file}} ({nrow(prof)} rows).")

# --- Half-hours and standard file ---
hh <- cloudnet_mlh_halfhour(prof, opts)
# Sensitivity cases as extra columns of the same half-hours (G 2026-10-23): precipitation rules "column" and
# "strict" (<= precip_top_agl), and single-gap bridging.
add_sens <- function(hh, hv, tag) {
  m <- match(hh$timestamp_start_utc, hv$timestamp_start_utc)
  hh[[paste0("x_pblh_m_", tag)]] <- hv$pblh_m[m]
  hh[[paste0("x_qc_flag_", tag)]] <- hv$qc_flag[m]
  hh
}
hh <- add_sens(hh, cloudnet_mlh_halfhour(cloudnet_apply_precip(prof, prof$precip_column), opts), "precip_column")
hh <- add_sens(hh, cloudnet_mlh_halfhour(cloudnet_apply_precip(prof, prof$precip_strict), opts), "precip_strict")
hh <- add_sens(hh, cloudnet_mlh_halfhour(dplyr::bind_rows(lapply(sens, `[[`, "bridge")), opts), "bridged")
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
  instrument_altitude_m = if (length(unique(ea$altitude_m)) == 1) unique(ea$altitude_m) else rng(ea$altitude_m),
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
    rng(ea$altitude_m), " m). Ground-reference uncertainty +/- ", abs(alt_icos - as.numeric(rng(ea$altitude_m))),
    " m: ICOS site altitude ", alt_icos, " m; the remote-sensing instruments stand on a field about 200 m from the ",
    "main SMEAR II station (Sinclair et al. 2022). Gate top = centre + ", rng(ea$gate_m / 2), " m. ",
    "Profile MLH (cloudnet_mlh_profiles()): eps_surface = ",
    opts$eps_surface, " m2 s-3, z_min = ", paste(z_used, collapse = ", "), " m agl (gate centre",
    if (!is.null(opts$z_min_agl)) ", fixed for the site in config/measured_sites.csv" else ", per file", "); ",
    "status no_signal (no valid epsilon) is NA; precipitation (rule '", opts$precip_rule, "': Cloudnet ",
    "classification rain_detected within +/- ", opts$precip_window_s, " s) is NA; per-profile precipitation ",
    "category kept in the profile output. Sensitivity columns: x_pblh_m_/x_qc_flag_precip_column (rain_detected, ",
    "or classification classes ", paste(opts$precip_classes, collapse = ", "), " from the bottom of the z_min gate ",
    "up to one lidar gate above the gate where the upward search stopped), _precip_strict (classes at or below ",
    opts$precip_top_agl, " m agl), _bridged (single missing gate bridged when the gates below and above are ",
    "turbulent). Censored = epsilon below threshold or missing at z_min. ",
    "Half-hours (cloudnet_mlh_halfhour()): expected profiles = median count per half-hour of each daily file, ",
    "coverage = valid / expected (capped at 1), min_coverage = ", opts$min_coverage,
    "; censored if more than half the valid profiles are censored; else median of the measured profiles, ",
    "qc 3 if more than ", opts$signal_limited_max, " of them ended at a missing gate (signal-limited). ",
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

cli::cli_h2("Precipitation category (profiles with a classification check)")
print(round(100 * prop.table(table(season = season_of(prof$time), category = prof$precip_category), 1), 1))
print(table(category = prof$precip_category, useNA = "ifany"))

cli::cli_h2("Half-hours")
hq <- std$data
cli::cli_inform(c(
  "{nrow(hq)} half-hours with at least one profile; coverage: median {round(stats::median(hq$x_coverage), 2)}, below {opts$min_coverage}: {sum(is.na(hq$qc_flag))}.",
  "Expected profiles per half-hour: median {stats::median(hq$x_n_expected)} (range {paste(range(hq$x_n_expected), collapse = '-')}).",
  "Half-hours without precipitation check (no classification file): {sum(!hq$x_precip_checked)}."
))
print(table(qc_flag = hq$qc_flag, useNA = "ifany"))
cli::cli_inform("MLH [m agl], qc 0, 2 or 3:")
print(round(stats::quantile(hq$pblh_m, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)))

qc_lab <- function(q) factor(ifelse(is.na(q), "NA", q), levels = c("0", "1", "2", "3", "NA"))
for (tz_lab in c("UTC", "local")) {
  off <- if (tz_lab == "UTC") 0 else meta_site$utc_offset_h
  tl <- hq$timestamp_start_utc + off * 3600
  hr <- as.integer(format(tl, "%H", tz = "UTC"))
  cli::cli_inform("qc_flag [% of half-hours] by hour ({tz_lab}{if (off != 0) paste0(', UTC+', off) else ''}) x season (season of the UTC date):")
  for (s in levels(season_of(hq$timestamp_start_utc))) {
    sel <- season_of(hq$timestamp_start_utc) == s
    tb <- table(hour = hr[sel], qc = qc_lab(hq$qc_flag[sel]))
    cat("--", s, "\n")
    print(round(100 * prop.table(tb, 1)))
  }
  qc_long <- as.data.frame(table(season = season_of(hq$timestamp_start_utc), hour = hr, qc = qc_lab(hq$qc_flag)))
  readr::write_csv(qc_long, file.path(tables_dir, paste0("07c_qc_by_hour_season_", tz_lab, "_", site_id, ".csv")))
}
cens <- cbind(season_hour(hq$timestamp_start_utc), cens = hq$qc_flag %in% 1L, valid = !is.na(hq$qc_flag))
cens_tab <- dplyr::summarise(dplyr::group_by(cens, .data$season, .data$hour_utc),
                             n_valid = sum(.data$valid), n_censored = sum(.data$cens),
                             pct_censored = round(100 * sum(.data$cens) / max(1, sum(.data$valid)), 1), .groups = "drop")
readr::write_csv(cens_tab, file.path(tables_dir, paste0("07c_censoring_", site_id, ".csv")))
cli::cli_inform("Median MLH [m agl] (qc 0/2/3), season x hour (UTC):")
mm <- cbind(season_hour(hq$timestamp_start_utc), pblh = hq$pblh_m)
print(round(with(mm, tapply(pblh, list(hour_utc = hour_utc, season = season), stats::median, na.rm = TRUE))))

# --- Precipitation rule sensitivity (half-hours) ---
cli::cli_h2("Precipitation rule sensitivity")
rule_summary <- function(flag, label) {
  pv <- cloudnet_apply_precip(prof, flag)
  hv <- cloudnet_mlh_halfhour(pv, opts)
  s <- season_of(hv$timestamp_start_utc)
  dplyr::bind_rows(lapply(c(levels(s), "all"), function(ss) {
    sel <- if (ss == "all") rep(TRUE, nrow(hv)) else s == ss
    q <- hv$qc_flag[sel]
    data.frame(rule = label, season = ss, pct_profiles_precip = round(100 * mean(pv$status[if (ss == "all") TRUE else season_of(pv$time) == ss] == "precipitation"), 1),
               n_qc0 = sum(q %in% 0L), n_qc1 = sum(q %in% 1L), n_qc2 = sum(q %in% 2L), n_qc3 = sum(q %in% 3L),
               n_na = sum(is.na(q)), median_mlh_m = stats::median(hv$pblh_m[sel], na.rm = TRUE),
               pct_censored = round(100 * sum(q %in% 1L) / max(1, sum(!is.na(q))), 1))
  }))
}
rule_tab <- rbind(rule_summary(prof$precip, "rain_detected (default)"), rule_summary(prof$precip_column, "column"),
                  rule_summary(prof$precip_strict, "strict (<= 300 m)"))
print_tbl(rule_tab)
readr::write_csv(rule_tab, file.path(tables_dir, paste0("07c_precip_rule_sensitivity_", site_id, ".csv")))

# --- Signal-limited profiles ---
cli::cli_h2("Signal-limited profiles (search ended at a missing gate)")
meas <- prof$status == "measured"
sh <- season_hour(prof$time[meas])
lim <- prof$gap_in_column[meas]
cli::cli_inform("{sum(lim)} of {sum(meas)} measured profiles ({round(100 * mean(lim), 1)} %); half-hours with qc 3: {sum(hq$qc_flag %in% 3L)}.")
cli::cli_inform("(a) epsilon in the last turbulent gate before the gap / 1e-4, median, season x hour (UTC):")
print(round(tapply(prof$eps_top_ratio[meas][lim], list(hour_utc = sh$hour_utc[lim], season = sh$season[lim]),
                   stats::median), 1))
cli::cli_inform("    quantiles over all signal-limited profiles:")
print(round(stats::quantile(prof$eps_top_ratio[meas][lim], c(.05, .25, .5, .75, .95)), 2))
cli::cli_inform("(b) median MLH [m agl] of signal-limited / epsilon-limited profiles, season x hour (UTC):")
ml <- prof$mlh_m[meas]
print(round(tapply(ml[lim], list(hour_utc = sh$hour_utc[lim], season = sh$season[lim]), stats::median)))
print(round(tapply(ml[!lim], list(hour_utc = sh$hour_utc[!lim], season = sh$season[!lim]), stats::median)))
cli::cli_inform("    quantiles, signal-limited vs epsilon-limited:")
print(rbind(signal_limited = round(stats::quantile(ml[lim], c(.05, .25, .5, .75, .95))),
            eps_limited = round(stats::quantile(ml[!lim], c(.05, .25, .5, .75, .95)))))
cli::cli_inform("(c) bridging a single missing gate (next gate turbulent):")
b_meas <- prof$status_eps == "measured"
dmb <- prof$mlh_bridged_m[b_meas] - prof$mlh_eps_m[b_meas]
cli::cli_inform("measured profiles whose MLH changes: {sum(dmb > 0)} of {sum(b_meas)} ({round(100 * mean(dmb > 0), 1)} %); increase among changed: median {stats::median(dmb[dmb > 0])} m, p90 {stats::quantile(dmb[dmb > 0], 0.9)} m, max {max(dmb)} m.")
if (!all(vapply(sens, is.null, logical(1)))) {
  pv <- dplyr::bind_rows(lapply(sens, `[[`, "bridge"))
  hv <- cloudnet_mlh_halfhour(pv, opts)
  m <- merge(hh[, c("timestamp_start_utc", "pblh_m", "qc_flag")], hv[, c("timestamp_start_utc", "pblh_m", "qc_flag")],
             by = "timestamp_start_utc", suffixes = c("", "_b"))
  bv <- !is.na(m$pblh_m) & !is.na(m$pblh_m_b)
  d <- m$pblh_m_b[bv] - m$pblh_m[bv]
  cli::cli_inform("half-hours with MLH: {sum(bv)}; changed {sum(abs(d) > 1e-6)} ({round(100 * mean(abs(d) > 1e-6), 1)} %), median change among changed {stats::median(d[abs(d) > 1e-6])} m; qc 3 {sum(m$qc_flag %in% 3L)} -> {sum(m$qc_flag_b %in% 3L)}.")
}

if (!is.null(sens[[1]]$a)) {
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
