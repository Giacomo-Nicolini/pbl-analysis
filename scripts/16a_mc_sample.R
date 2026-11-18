# Spec 16 (phase 5b): footprint uncertainty from the FFP inputs. Step 16a: strata, stratified sample, perturbation
# draws and chunk list (no FFP).
# Per site (config/footprint_mc_sites.txt, or the sites given): the approved reference year (config/reference_years.csv)
#   of the 08a input table; eligible rows = valid under the etc set for every PBLH case present among pblhmod,
#   pblhmod_clim, era5 (spec-15 common_pe), subset valid, with an FFP result in the 08b chunks of each case;
#   strata zm/L class x day/night (6); n per stratum (170; all rows if fewer); seed per site (mc_site_seed()), stored.
#   Design A: all sampled rows (17 runs each under the etc set plus 4 rsl_off runs of hc and z0, spec 16 amendment 1;
#   built in 16b from mc_oat_runs(); chunks sized on mc_oat_ffp_runs() = 19 distinct FFP inputs). Design B: 50 rows
#   per stratum, 50 joint draws per row (standard deviates stored here), two h options: h_pct30 (30 %) and h_empirical (day: sigma_ln of the
#   daily ln(PBLHmod / measured) at the aerosol sites, median over sites of IQR / 1.349; night: 30 %).
#   The spec-15 values of each sampled row (TA fraction, x_80, h of pblhmod, pblhmod_clim, era5) come from the 08b
#   chunk _hh files (no new runs).
# h_empirical: computed once from the standard files of the primary aerosol sources (config/measured_sites.csv) and
#   cached in <derived_dir>/footprint_mc/h_empirical_sigma.csv (--recompute-h-empirical to redo).
# Output: <derived_dir>/footprint_mc/<site>_draws.rds (sample, weights, base inputs, spec-15 values, B deviates, seed,
#   input spec, h options); <derived_dir>/footprint_mc/chunk_list.csv (site x design x part; rows of the other sites
#   kept); log lines with n per stratum and the chunk estimate.
#
# Usage: Rscript scripts/16a_mc_sample.R [SITE ...] [--n-per-stratum=170] [--n-mc-per-stratum=50] [--n-draws=50]
#          [--sigmav-cv=0.15] [--cpus=2] [--sec-per-run=S] [--year=YYYY] [--period=YYYY-MM] [--recompute-h-empirical]
#   --period (PC smoke test): the 08a input table of that month (<site>_footprint_inputs_<period>.rds) and its 08b
#   chunks, instead of the reference year. Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
known <- c("--n-per-stratum", "--n-mc-per-stratum", "--n-draws", "--sigmav-cv", "--cpus", "--sec-per-run", "--year",
           "--period", "--recompute-h-empirical")
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), known)
if (length(bad)) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
cfg <- get_cfg()
n_per <- as.integer(opt("n-per-stratum", "170"))
n_mc_per <- as.integer(opt("n-mc-per-stratum", "50"))
n_draws <- as.integer(opt("n-draws", "50"))
sigmav_cv <- as.numeric(opt("sigmav-cv", "0.15"))
cpus <- as.integer(opt("cpus", "2"))
sec_per_run <- as.numeric(opt("sec-per-run", as.character(cfg$sec_per_ffp_halfhour)))
period <- opt("period")
spec <- mc_input_spec(sigmav_cv = sigmav_cv)
sites_tab <- read_sites(resolve_path(cfg$sites_file))
site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
if (!length(site_ids)) {
  site_ids <- trimws(readLines(here::here("config", "footprint_mc_sites.txt"), warn = FALSE))
  site_ids <- site_ids[nzchar(site_ids) & !startsWith(site_ids, "#")]
}
if (length(setdiff(site_ids, sites_tab$site_id))) cli::cli_abort("Unknown site{?s} {.val {setdiff(site_ids, sites_tab$site_id)}}.")
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
mc_dir <- file.path(derived, "footprint_mc")
dir.create(mc_dir, recursive = TRUE, showWarnings = FALSE)
cases <- c("pblhmod", "pblhmod_clim", "era5")

# --- h_empirical: sigma_ln of the daily ln(PBLHmod / measured) at the primary aerosol sources ---
he_file <- file.path(mc_dir, "h_empirical_sigma.csv")
if (!file.exists(he_file) || "--recompute-h-empirical" %in% args) {
  ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
  aer <- ms[!is.na(ms$primary_source_id) & ms$method_class %in% "aerosol" & ms$status %in% "used", ]
  std_dir <- resolve_path(cfg$standard_dir)
  he <- dplyr::bind_rows(lapply(seq_len(nrow(aer)), function(k) {
    s <- aer$icos_site[k]
    p <- function(src) file.path(std_dir, s, paste0(s, "__", src, ".csv.gz"))
    need <- c(p(aer$primary_source_id[k]), p("pblhmod_default"), p("era5_blh_center"))
    if (!all(file.exists(need))) {
      cli::cli_alert_warning("h_empirical: {s} skipped (standard file{?s} missing: {basename(need[!file.exists(need)])}).")
      return(NULL)
    }
    x <- mc_daily_log_spread(read_standard_pblh(need[1]), read_standard_pblh(need[2]),
                             list(era5_blh_center = read_standard_pblh(need[3])), sites_tab$lon[sites_tab$site_id == s])
    dplyr::bind_cols(dplyr::tibble(site_id = s, source_id = aer$primary_source_id[k]), x)
  }))
  if (!nrow(he) || all(is.na(he$sd_ln_daily))) cli::cli_abort("h_empirical: no aerosol site with standard files.")
  write_csv_header(he, he_file, c(
    "Spec 16 h_empirical: per primary aerosol source, daily median of ln(PBLHmod / measured) on the comparison set",
    "(apparent solar window of measured_compare_options(), PBLHmod convective); sd_ln_daily = IQR / 1.349 of the daily",
    "values [-]; sigma_ln used = median over sites of sd_ln_daily (applied by day; night: 30 %)."))
}
he <- readr::read_csv(he_file, comment = "#", show_col_types = FALSE, progress = FALSE)
sigma_ln_emp <- stats::median(he$sd_ln_daily, na.rm = TRUE)
cv_emp <- sqrt(exp(sigma_ln_emp^2) - 1) # the lognormal of mc_factor() is parameterised by the CV
h_options <- list(h_pct30 = list(day = 0.30, night = 0.30), h_empirical = list(day = cv_emp, night = 0.30))
cli::cli_inform("h_empirical: sigma_ln {signif(sigma_ln_emp, 3)} (median over {sum(!is.na(he$sd_ln_daily))} aerosol site{?s}: {paste(he$site_id, signif(he$sd_ln_daily, 3), collapse = ', ')}); CV {signif(cv_emp, 3)}.")

ry <- if (is.null(period) && is.null(opt("year"))) read_reference_years(resolve_path(cfg$reference_years_file), site_ids) else NULL
chunks_all <- list()
for (site_id in site_ids) {
  cli::cli_h2("{site_id}")
  fp_dir <- file.path(derived, "footprint", site_id)
  inp_file <- file.path(fp_dir, paste0(site_id, "_footprint_inputs", if (is.null(period)) "" else paste0("_", period), ".rds"))
  if (!file.exists(inp_file)) cli::cli_abort("{site_id}: input table {.path {inp_file}} not found (step 08a).")
  inp <- readRDS(inp_file)
  tab <- inp$tab
  t <- tab$timestamp_start_utc
  year <- if (!is.null(period)) as.integer(substr(period, 1, 4)) else if (!is.null(opt("year"))) as.integer(opt("year")) else ry$year[ry$site_id == site_id]
  sel <- if (!is.null(period)) format(t, "%Y-%m", tz = "UTC") == period else format(t, "%Y", tz = "UTC") == as.character(year)
  cs <- cases[paste0("h_", cases) %in% names(tab)]
  if (!all(c("pblhmod", "era5") %in% cs)) cli::cli_abort("{site_id}: pblhmod and era5 cases needed in the input table.")
  if (length(setdiff(cases, cs))) cli::cli_alert_warning("{site_id}: case{?s} {.val {setdiff(cases, cs)}} not in the input table: row set and source effects without {?it/them}.")
  # spec-15 values from the 08b chunk _hh files of the selected months
  months <- sort(unique(format(t[sel], "%Y-%m", tz = "UTC")))
  cdir <- file.path(fp_dir, "chunks")
  ref <- list()
  for (k in cs) {
    parts <- unlist(lapply(months, function(mo) chunk_parts_on_disk(cdir, site_id, k, mo)))
    files <- file.path(cdir, paste0(site_id, "_", k, "_", parts, "_hh.csv.gz"))
    files <- files[file.exists(files)]
    if (!length(files)) cli::cli_abort("{site_id}: no 08b chunk of case {k} for {paste(months, collapse = ', ')}.")
    x <- dplyr::bind_rows(lapply(files, function(f) readr::read_csv(f, col_select = c("timestamp_start_utc", "h", "x_80", "ffp_ok", "f_target_area"),
                                                                    show_col_types = FALSE, progress = FALSE)))
    m <- match(as.numeric(t), as.numeric(x$timestamp_start_utc))
    ref[[k]] <- list(ta = x$f_target_area[m], x80 = x$x_80[m], h = x$h[m], ok = x$ffp_ok[m] %in% TRUE & is.finite(x$f_target_area[m]))
  }
  sub <- inp$subsets
  vcol <- function(u, k) { v <- sub[[paste(u, "etc", k, sep = "__")]]; if (is.null(v)) rep(FALSE, nrow(tab)) else v %in% TRUE }
  elig <- sel
  for (k in cs) elig <- elig & vcol("valid", k) & ref[[k]]$ok
  qc <- Reduce(`&`, lapply(cs, function(k) vcol("valid_qc", k)))
  qcu <- Reduce(`&`, lapply(cs, function(k) vcol("valid_qc_used", k)))
  st <- mc_stratum(tab$zm_m, tab$L, tab$day)
  pos <- which(elig)
  seed <- mc_site_seed(site_id)
  set.seed(seed)
  smp <- mc_sample(st[pos], n_per)
  i <- pos[smp$i]
  smp$i <- NULL
  smp <- dplyr::bind_cols(
    dplyr::tibble(row_id = seq_along(i), timestamp_start_utc = t[i]), smp,
    dplyr::tibble(ec_height_m = tab$ec_height_m[i], canopy_height_m = tab$canopy_height_m[i], zm_m = tab$zm_m[i],
                  z0_m = tab$z0_m[i], umean = tab$umean[i], L = tab$L[i], sigmav = tab$sigmav[i], ustar = tab$ustar[i],
                  wind_dir = tab$wind_dir[i], h = tab$h_pblhmod[i], day = tab$day[i],
                  ustar_implausible = tab$ustar_implausible[i] %in% TRUE, excluded = tab$excluded[i] %in% TRUE,
                  valid_qc = qc[i], valid_qc_used = qcu[i]))
  for (k in cases) {
    for (q in c("ta", "x80", "h")) smp[[paste0(q, "__", k)]] <- if (k %in% cs) ref[[k]][[q]][i] else NA_real_
  }
  smp$in_b <- mc_subsample(smp, n_mc_per)
  rb <- smp$row_id[smp$in_b]
  zb <- stats::setNames(lapply(rb, function(r) mc_draw_z(n_draws, spec)), as.character(rb))
  ns <- dplyr::count(dplyr::tibble(stratum = factor(smp$stratum, levels = mc_strata_levels()), b = smp$in_b), stratum, b)
  cli::cli_inform("{if (is.null(period)) paste('reference year', year) else paste('period', period)}: {sum(sel)} half-hours, {length(pos)} eligible (etc, valid, cases {paste(cs, collapse = ', ')}); seed {seed}.")
  print(as.data.frame(dplyr::summarise(dplyr::group_by(smp, stratum), n_eligible = dplyr::first(n_stratum), n_A = dplyr::n(),
                                       n_B = sum(in_b), w_stratum = round(dplyr::first(w_stratum), 3), .groups = "drop")), row.names = FALSE)
  draws <- list(site_id = site_id, year = year, period = period, seed = seed, rng_kind = RNGkind(), spec = spec,
                sample = smp, z_b = zb, n_draws = n_draws, h_options = h_options, sigma_ln_empirical = sigma_ln_emp,
                h_empirical_sites = he, m_etc = inp$m_etc, cases = cs, n_eligible = length(pos), inputs_file = basename(inp_file),
                created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  out <- file.path(mc_dir, paste0(site_id, "_draws.rds"))
  write_rds_atomic(draws, out)
  cl <- mc_chunk_list(site_id, nrow(smp), length(rb), runs_a = mc_oat_ffp_runs(mc_oat_runs(spec)),
                      runs_b = n_draws * length(h_options), sec_per_run = sec_per_run, cpus = cpus,
                      target_min = cfg$chunk_target_min)
  cli::cli_inform("{.path {out}}; chunks: {nrow(cl)} ({sum(cl$design == 'A')} A, {sum(cl$design == 'B')} B), {round(sum(cl$cpu_hours), 1)} CPU-hours before removing invalid and duplicated runs ({sec_per_run} s per run).")
  chunks_all[[site_id]] <- cl
}
cl_file <- file.path(mc_dir, "chunk_list.csv")
new <- dplyr::bind_rows(chunks_all)
old <- if (file.exists(cl_file)) readr::read_csv(cl_file, col_types = readr::cols(.default = "c"), progress = FALSE) else NULL
if (!is.null(old)) {
  old <- old[!old$site_id %in% site_ids, ]
  new <- dplyr::bind_rows(dplyr::mutate(old, dplyr::across(c(first_row, last_row, n_rows, n_runs), as.integer),
                                        dplyr::across(c(est_min, cpu_hours), as.numeric)), new)
}
readr::write_csv(new, cl_file)
cli::cli_alert_success("{.path {cl_file}}: {nrow(new)} chunks, {round(sum(new$cpu_hours), 1)} CPU-hours (upper bound), sites {paste(unique(new$site_id), collapse = ', ')}.")
