# Spec 01: check the input layer for the test site(s).
# Reads ICOS L2 (FLUXES, FLUXNET, METEO) and ERA5 BLH, prints a concise report
# and writes the range/coverage report to <output_dir>/tables/01_check_inputs_<site>.csv and the input summary
# (files, reference period, ERA5 cell; spec 09 A4) to 01_inputs_<site>.csv.
#
# Usage: Rscript scripts/01_check_inputs.R [SITE ...]
#   Without arguments the sites in cfg$test_sites are used.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
site_ids <- if (length(args) > 0) args else unlist(cfg$test_sites)
if (length(site_ids) == 0) {
  cli::cli_abort("No sites: pass station codes as arguments or set {.field test_sites} in {.file config.yml}.")
}

sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

fmt_t <- function(x) format(x, "%Y-%m-%d %H:%M UTC")
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 6)

for (site_id in site_ids) {
  cli::cli_h1("Input check: {site_id}")

  # --- ICOS L2 ---
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
  info <- attr(d, "read_info")

  cli::cli_h2("ICOS L2 files")
  info_print <- dplyr::transmute(
    info, file,
    first_start = fmt_t(first_start_utc), last_start = fmt_t(last_start_utc),
    n_rows_file, n_rows_added, n_rows_outside_ref, missing_cols, optional_cols
  )
  print_tbl(info_print)
  ta_lev <- info$optional_cols[info$file == "METEO"]
  cli::cli_inform("METEO TA profile levels found in the file header: {if (length(ta_lev) && nzchar(ta_lev)) ta_lev else 'none'}.")
  cli::cli_inform(c(
    "Merged:{nrow(d)} half-hours, {fmt_t(min(d$timestamp_start_utc))} to {fmt_t(max(d$timestamp_start_utc))}.",
    "Rows present in: FLUXES {sum(d$src_fluxes)}, FLUXNET {sum(d$src_fluxnet)}, METEO {sum(d$src_meteo)}; in none: {sum(!(d$src_fluxes | d$src_fluxnet | d$src_meteo))}."
  ))

  cli::cli_h2("Range / coverage report (report only, no filtering)")
  ranges <- dplyr::bind_rows(
    icos_default_ranges(),
    dplyr::tibble(variable = "fn_H_F_MDS", lower = NA_real_, upper = NA_real_,
                  abs = FALSE, strict_lower = FALSE) # coverage only, no range in spec
  )
  rr <- icos_range_report(d, ranges)
  print_tbl(rr)

  cli::cli_h2("Site heights")
  print(table(zm_source = d$zm_source, useNA = "ifany"))
  print(table(hc_source = d$hc_source, useNA = "ifany"))
  cli::cli_inform("zm range: {paste(range(d$zm_m, na.rm = TRUE), collapse = ' - ')} m; canopy height range: {paste(round(range(d$canopy_height_m, na.rm = TRUE), 3), collapse = ' - ')} m.")

  cli::cli_h2("Exclusions")
  cli::cli_inform("{sum(d$excl_flag)} of {nrow(d)} records flagged as excluded.")

  cli::cli_h2("METEO PA unit check")
  if ("met_PA_PA" %in% names(d) && any(!is.na(d$met_PA_PA))) {
    pa_kpa <- stats::median(d$met_PA_PA, na.rm = TRUE) / 1000 # value as in the file
    if (pa_kpa > 50 && pa_kpa < 110) {
      cli::cli_alert_success("METEO PA median in file = {round(pa_kpa, 2)}: consistent with kPa (converted to Pa).")
    } else {
      cli::cli_alert_danger("METEO PA median in file = {round(pa_kpa, 2)}: NOT consistent with kPa. Check the unit before using {.field met_PA_PA}.")
    }
  } else {
    cli::cli_alert_info("METEO PA not available.")
  }

  # --- ERA5 ---
  cli::cli_h2("ERA5 BLH")
  era_files <- find_era5_files(site_id, cfg, "single_level", "boundary_layer_height", sites = sites)
  print_tbl(era_files[c("file", "cell_row", "cell_col")])
  era <- read_era5_single_level(era_files, "boundary_layer_height", cells = "center")
  cells <- dplyr::summarise(
    era, .by = c(cell_row, cell_col),
    grid_lat = dplyr::first(grid_lat), grid_lon = dplyr::first(grid_lon),
    grid_dist_km = round(dplyr::first(grid_dist_km), 2),
    first = fmt_t(min(timestamp_utc)), last = fmt_t(max(timestamp_utc)),
    n_hours = dplyr::n(), n_valid = sum(!is.na(blh_m))
  )
  print_tbl(cells)
  print(table(year = lubridate::year(era$timestamp_utc)))

  out <- file.path(tables_dir, paste0("01_check_inputs_", site_id, ".csv"))
  readr::write_csv(dplyr::mutate(rr, site_id = site_id, .before = 1), out)
  cli::cli_alert_success("Range report written to {.path {out}}")

  # Input summary for the network collection (spec 09 A4), long format: one value per item. Reaching this point
  # means the unit checks of the readers passed (they abort otherwise).
  fi <- function(k, what) { x <- info[[what]][info$file == k]; if (length(x)) as.character(x[1]) else NA_character_ }
  cc <- cells[cells$cell_row == 0 & cells$cell_col == 0, ]
  items <- c(
    file_fluxes = fi("FLUXES", "path"), file_fluxnet = fi("FLUXNET", "path"), file_meteo = fi("METEO", "path"),
    reference_start_utc = fmt_t(min(d$timestamp_start_utc)), reference_end_utc = fmt_t(max(d$timestamp_end_utc)),
    n_halfhours = nrow(d), n_excluded = sum(d$excl_flag),
    n_fluxes = sum(d$src_fluxes), n_fluxnet = sum(d$src_fluxnet), n_meteo = sum(d$src_meteo),
    meteo_ta_levels = if (length(ta_lev) && nzchar(ta_lev)) ta_lev else NA_character_,
    zm_min_m = min(d$zm_m, na.rm = TRUE), zm_max_m = max(d$zm_m, na.rm = TRUE),
    hc_min_m = min(d$canopy_height_m, na.rm = TRUE), hc_max_m = max(d$canopy_height_m, na.rm = TRUE),
    units_checked = "ok",
    era5_first_utc = cc$first, era5_last_utc = cc$last, era5_n_hours = cc$n_hours, era5_n_valid = cc$n_valid,
    era5_cell_dist_km = cc$grid_dist_km
  )
  items[c("file_fluxes", "file_fluxnet", "file_meteo")] <- basename(items[c("file_fluxes", "file_fluxnet", "file_meteo")])
  out_i <- file.path(tables_dir, paste0("01_inputs_", site_id, ".csv"))
  readr::write_csv(dplyr::tibble(site_id = site_id, item = names(items), value = unname(items)), out_i)
  cli::cli_alert_success("Input summary written to {.path {out_i}}")
}
