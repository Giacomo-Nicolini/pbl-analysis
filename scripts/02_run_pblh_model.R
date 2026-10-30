# Spec 02: run PBLHmod v2 for the test site(s) and print a summary.
# Runs: default (= PBLHmod, EC data only: bg_only, neutral auto, constant gamma 0.01 K/m); reference gamma values
# (G 2026-10-30): gamma001 (0.01 K/m, Batchvarova & Gryning 1991 as cited by Kljun et al. 2015; equal to default
# until the EC-only constant is revised) and gamma0025 (0.0025 K/m, legacy ETC code; no literature source found);
# max_B2 (combine max, neutral B2: behaviour before amendment 11); max_auto;
# u* filter sensitivity (spec 05, spec 02 amendment 12): ustar04 / ustar06 (u*/WS limit 0.4 / 0.6) and
# nofilter (u* filter disabled); L consistency (spec 06, spec 02 amendment 13): l_L2 (no check = spec 05),
# l_forcing (L always from the forcing, option D), lcons2 / lcons5 (factor 2 / 5); nofilter_l2 (no u* filter,
# no L check = spec 02). All other cases use the defaults (u* filter on, l_source = "L2_checked").
# Ends with u* filter and L check summaries across cases and a side-by-side table of quantiles per regime.
# Output: <derived_dir>/pblhmod/<site>_pblhmod_<label>.csv.gz
# Regression chain (if the backups exist in <derived_dir>/pblhmod/): nofilter_l2 vs pre_spec05/ (spec-02
# default), l_L2 vs pre_spec06/ (spec-05 default) and l_check_allday (L check also at night = spec 06
# amendment 1) vs pre_spec06a2/, and every case vs pre_spec07/ (outputs before the reference-period rule,
# spec 07 step 0), and every case vs pre_spec11/ (outputs before the gamma_source option, spec 11 B: the existing
# cases must be unchanged); all must be identical (< 1e-9 m, same NA, regime, method).
# Spec 11 B cases (ERA5 gamma from <derived_dir>/gamma/, script 11a): clim (= PBLHmod-clim, standard source
# pblhmod_clim: site x month climatology; G 2026-10-30),
# gamma_era5daily (morning gamma of the day, climatology fallback), gamma_era5clim_zmin1000 (fit layer from 1000 m),
# gamma_era5clim_bgonly (ERA5 climatology in B&G only, N (Eq. B3) from the constant gamma; G 2026-10-30).
#
# Usage: Rscript scripts/02_run_pblh_model.R [SITE ...] [--cases=label1,label2]
#   Without sites the sites in cfg$test_sites are used; --cases runs only the listed cases.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
case_arg <- grep("^--cases=", args, value = TRUE)
args <- setdiff(args, case_arg)
site_ids <- if (length(args) > 0) args else unlist(cfg$test_sites)
if (length(site_ids) == 0) {
  cli::cli_abort("No sites: pass station codes as arguments or set {.field test_sites} in {.file config.yml}.")
}

sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
out_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "pblhmod")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- list(
  default = pblh_options(),
  gamma0025 = pblh_options(gamma = 0.0025, label = "gamma0025"),
  max_B2 = pblh_options(combine = "max", neutral_method = "B2", label = "max_B2"),
  max_auto = pblh_options(combine = "max", neutral_method = "auto", label = "max_auto"),
  ustar04 = pblh_options(ustar_ws_max = 0.4, label = "ustar04"),
  ustar06 = pblh_options(ustar_ws_max = 0.6, label = "ustar06"),
  nofilter = pblh_options(ustar_ws_max = Inf, ustar_max = Inf, label = "nofilter"),
  l_L2 = pblh_options(l_source = "L2", label = "l_L2"),
  l_forcing = pblh_options(l_source = "forcing", label = "l_forcing"),
  lcons2 = pblh_options(l_consistency_factor = 2, label = "lcons2"),
  lcons5 = pblh_options(l_consistency_factor = 5, label = "lcons5"),
  nofilter_l2 = pblh_options(ustar_ws_max = Inf, ustar_max = Inf, l_source = "L2", label = "nofilter_l2"),
  l_check_allday = pblh_options(l_check_daytime_only = FALSE, label = "l_check_allday"),
  clim = pblh_options(gamma_source = "era5_clim", label = "clim"),
  gamma001 = pblh_options(gamma = 0.01, label = "gamma001"),
  gamma_era5daily = pblh_options(gamma_source = "era5_daily", label = "gamma_era5daily"),
  gamma_era5clim_zmin1000 = pblh_options(gamma_source = "era5_clim", label = "gamma_era5clim_zmin1000"),
  gamma_era5clim_bgonly = pblh_options(gamma_source = "era5_clim", gamma_scope = "bg_only",
                                       label = "gamma_era5clim_bgonly")
)
# ERA5 gamma tables per case (file suffix in <derived_dir>/gamma/: "" = fit layer 500-3500 m, "_zmin1000")
gamma_layer <- c(gamma_era5clim_zmin1000 = "_zmin1000")
if (length(case_arg) == 1) {
  keep <- strsplit(sub("^--cases=", "", case_arg), ",")[[1]]
  bad <- setdiff(keep, names(runs))
  if (length(bad) > 0) cli::cli_abort("Unknown case{?s}: {.val {bad}}.")
  runs <- runs[keep]
}
g_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "gamma")
read_gamma_tables <- function(site_id, sfx) {
  p <- file.path(g_dir, paste0(site_id, sfx, c("_daily.csv", "_monthly.csv")))
  if (!all(file.exists(p))) cli::cli_abort("ERA5 gamma tables {.path {p}} not found: run scripts/11a_gamma_era5.R first.")
  list(daily = readr::read_csv(p[1], show_col_types = FALSE, progress = FALSE),
       clim = readr::read_csv(p[2], show_col_types = FALSE, progress = FALSE))
}

# L consistency accounting for one case (rule in l_consistency(), spec 02 amendment 13).
l_report <- function(out, opts) {
  ok <- !out$excl_flag
  chk <- out$L_check %in% c("consistent", "sign", "magnitude")
  rp <- out$flag_L_replaced
  data.frame(
    case = opts$label, l_source = opts$l_source, checked = sum(chk),
    repl_sign = sum(rp & out$L_check %in% "sign"), repl_magnitude = sum(rp & out$L_check %in% "magnitude"),
    replaced = sum(rp), days_with_repl = length(unique(out$day_id[rp])),
    pct_L_from_L2 = round(100 * sum(out$L_source %in% "L2" & ok) / sum(!is.na(out$L) & ok), 2),
    sign_inconsistent = sum(out$flag_sign_inconsistent),
    gt5km = sum(out$pblh_m > 5000, na.rm = TRUE), gt10km = sum(out$pblh_m > 10000, na.rm = TRUE)
  )
}

# Compare a case with a saved reference output (regression chain).
check_reference <- function(new, ref_path, what) {
  if (is.null(new)) return(invisible(NA)) # case not run in this call
  if (!file.exists(ref_path)) {
    cli::cli_alert_info("{what}: reference {.path {ref_path}} not found, check skipped.")
    return(invisible(NA))
  }
  cli::cli_h2(what)
  old <- readr::read_csv(ref_path, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
  same_rows <- nrow(old) == nrow(new) && all(old$timestamp_start_utc == new$timestamp_start_utc)
  cmp <- vapply(c("pblh_m", "h_diag_m", "h_bg_m", "ustar", "L"), function(k) {
    if (!identical(is.na(old[[k]]), is.na(new[[k]]))) return(Inf)
    max(c(0, abs(old[[k]] - new[[k]])), na.rm = TRUE)
  }, numeric(1))
  print_tbl(data.frame(column = names(cmp), max_abs_diff = cmp))
  chr_same <- all(vapply(c("regime", "method"), function(k) identical(old[[k]], new[[k]]), logical(1)))
  ok <- same_rows && chr_same && all(cmp < 1e-9)
  if (ok) {
    cli::cli_alert_success("{what}: identical (max |diff| < 1e-9 m, same NA pattern, same regime/method).")
  } else {
    cli::cli_alert_danger("{what}: differs (rows {same_rows}, regime/method {chr_same}).")
  }
  invisible(ok)
}

# u* filter accounting for one case (rule in ustar_implausible(), spec 02 amendment 12).
ustar_report <- function(d, out, opts) {
  ok <- !d$excl_flag & !is.na(d$USTAR)
  r <- ustar_implausible(d$USTAR, d$WS, opts)
  rej <- out$flag_ustar_implausible
  data.frame(
    case = opts$label, n_ustar = sum(ok),
    ratio_only = sum(ok & r$ratio & !r$cap), cap_only = sum(ok & r$cap & !r$ratio), both = sum(ok & r$ratio & r$cap),
    rejected = sum(rej), pct = round(100 * sum(rej) / sum(ok), 3),
    interpolated = sum(rej & out$flag_ustar_interp), left_na = sum(rej & is.na(out$ustar)),
    gt5km = sum(out$pblh_m > 5000, na.rm = TRUE), gt10km = sum(out$pblh_m > 10000, na.rm = TRUE)
  )
}
probs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 5)
qtab <- function(x, label) {
  q <- stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
  stats::setNames(as.data.frame(t(c(n = sum(!is.na(x)), round(q, 1)))),
                  c("n", paste0("p", probs * 100)))
}
count_tbl <- function(x, name) {
  t <- as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE)
  names(t) <- c(name, "n")
  t$pct <- round(100 * t$n / sum(t$n), 2)
  t
}

for (site_id in site_ids) {
  cli::cli_h1("PBLHmod: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)

  cli::cli_h2("Inputs")
  cli::cli_inform("{nrow(d)} half-hours; {sum(d$excl_flag)} excluded.")
  cli::cli_inform("H_F_MDS_QC distribution (FLUXNET convention: 0 measured, 1 good, 2 medium, 3 poor gap-fill):")
  print_tbl(count_tbl(d$fn_H_F_MDS_QC, "H_F_MDS_QC"))

  outs <- list()
  for (lab in names(runs)) {
    opts <- runs[[lab]]
    cli::cli_h2("Run {.val {opts$label}} (gamma = {opts$gamma} K m-1, stable = {opts$stable_method}, neutral = {opts$neutral_method}, H = {opts$h_source}, combine = {opts$combine})")
    ge <- NULL
    if (opts$gamma_source != "constant") {
      sfx <- if (lab %in% names(gamma_layer)) gamma_layer[[lab]] else ""
      ge <- read_gamma_tables(site_id, sfx)
      cli::cli_inform("ERA5 gamma ({opts$gamma_source}, scope {opts$gamma_scope}, tables {site_id}{sfx}_*.csv).")
    }
    tm <- system.time(out <- run_pblh_model(d, meta, opts, gamma_era5 = ge))
    if (opts$gamma_source != "constant") {
      cli::cli_inform("gamma used in B&G: median {signif(stats::median(out$gamma_used, na.rm = TRUE), 3)} K m-1; source: {paste(names(table(out$gamma_source_used)), table(out$gamma_source_used), sep = ' ', collapse = ', ')}.")
    }
    outs[[opts$label]] <- out
    cli::cli_inform("Run time: {round(tm[['elapsed']], 1)} s.")

    path <- file.path(out_dir, paste0(site_id, "_pblhmod_", opts$label, ".csv.gz"))
    readr::write_csv(out, path)
    cli::cli_alert_success("Written {.path {path}}")

    cli::cli_inform("u* filter (ustar_ws_max = {opts$ustar_ws_max}, ws_min_ratio_test = {opts$ws_min_ratio_test} m s-1, ustar_max = {opts$ustar_max} m s-1):")
    print_tbl(ustar_report(d, out, opts))

    ok <- !out$excl_flag
    cli::cli_inform("Valid pblh_m: {round(100 * mean(!is.na(out$pblh_m)), 2)}% of all rows, {round(100 * mean(!is.na(out$pblh_m[ok])), 2)}% of non-excluded rows.")
    print_tbl(count_tbl(out$regime, "regime"))
    print_tbl(count_tbl(out$method, "method"))

    cli::cli_inform("Quantiles [m]: PBLHmod vs operational PBLH on the same rows (both valid):")
    both <- !is.na(out$pblh_m) & !is.na(d$PBLH)
    rows <- list(all = both)
    for (r in c("stable", "unstable_diag", "convective")) rows[[r]] <- both & out$regime %in% r
    qt <- do.call(rbind, lapply(names(rows), function(r) {
      rbind(
        cbind(subset = r, source = "PBLHmod", qtab(out$pblh_m[rows[[r]]])),
        cbind(subset = r, source = "PBLH_oper", qtab(d$PBLH[rows[[r]]]))
      )
    }))
    print_tbl(qt)
    cli::cli_inform("PBLHmod quantiles on all valid rows:")
    print_tbl(qtab(out$pblh_m))

    flags <- grep("^flag_", names(out), value = TRUE)
    cli::cli_inform("Flag counts (rows):")
    print_tbl(data.frame(flag = flags, n = vapply(flags, function(fl) sum(out[[fl]]), numeric(1))))

    days <- unique(out$day_id[out$is_day %in% TRUE])
    day_share <- function(fl) {
      dflag <- unique(out$day_id[out[[fl]]])
      round(100 * mean(days %in% dflag), 2)
    }
    conv <- out$regime %in% "convective"
    cli::cli_inform(c(
      "Days with daylight: {length(days)}; days with a convective window: {length(unique(out$day_id[conv]))}.",
      "Share of days: flag_no_h0 {day_share('flag_no_h0')}%, flag_bg_stopped_gap {day_share('flag_bg_stopped_gap')}%, flag_no_positive_h {day_share('flag_no_positive_h')}%, flag_polar_day {day_share('flag_polar_day')}%.",
      "flag_sign_inconsistent: {sum(out$flag_sign_inconsistent)} rows = {round(100 * mean(out$flag_sign_inconsistent[conv]), 2)}% of convective rows."
    ))
    print_tbl(count_tbl(out$H_source, "H_source"))
    print_tbl(count_tbl(out$T_source, "T_source"))
    print_tbl(count_tbl(out$p_source, "p_source"))

    if (lab == "default") {
      # Conditions on rejected vs kept rows, to spot rejections in plausible conditions.
      cli::cli_h2("u* filter: conditions on rejected vs kept rows (raw L2 values)")
      ok <- !d$excl_flag & !is.na(d$USTAR)
      r <- ustar_implausible(d$USTAR, d$WS, opts)
      rho <- out$p_Pa / (pblh_constants()$r_d * out$T_K)
      grp <- ifelse(!ok, NA, ifelse(r$ratio & r$cap, "both", ifelse(r$ratio, "ratio", ifelse(r$cap, "cap", "kept"))))
      cond <- dplyr::tibble(grp = grp, ustar = d$USTAR, ws = d$WS, ratio = d$USTAR / d$WS,
                            sw_us = d$W_SIGMA / d$USTAR, ustar_tau = sqrt(abs(d$TAU) / rho))
      cond <- cond[!is.na(cond$grp), ]
      med <- function(x) round(stats::median(x, na.rm = TRUE), 3)
      print_tbl(dplyr::summarise(cond, .by = grp, n = dplyr::n(), ustar_med = med(ustar), ustar_p90 = round(stats::quantile(ustar, 0.9, na.rm = TRUE), 3),
                                 ws_med = med(ws), ws_p90 = round(stats::quantile(ws, 0.9, na.rm = TRUE), 2), ratio_med = med(ratio),
                                 sigw_ustar_med = med(sw_us), ustar_from_tau_med = med(ustar_tau)))
      rr <- cond[cond$grp == "ratio", ]
      cli::cli_inform(c(
        "Ratio-rejected rows by wind speed (WS bins [m s-1]) and u*/WS just above the limit (0.5-0.6):",
        paste(names(table(cut(rr$ws, c(1, 2, 4, 6, 8, Inf), right = FALSE))), table(cut(rr$ws, c(1, 2, 4, 6, 8, Inf), right = FALSE)),
              sep = ": ", collapse = "; "),
        "u*/WS in (0.5, 0.6]: {sum(rr$ratio <= 0.6)} rows, of which WS >= 6 m s-1: {sum(rr$ratio <= 0.6 & rr$ws >= 6)}; median sigma_w/u* on those: {med(rr$sw_us[rr$ratio <= 0.6])}."
      ))

      cli::cli_h2("is_day vs FLUXNET SW_IN_POT > 0 (FLUXNET period)")
      sel <- d$src_fluxnet & !is.na(d$fn_SW_IN_POT) & !is.na(out$is_day)
      pot_day <- d$fn_SW_IN_POT > 0
      mism <- which(sel & (out$is_day != pot_day))
      cli::cli_inform("Agreement: {round(100 * (1 - length(mism) / sum(sel)), 3)}% of {sum(sel)} half-hours; {length(mism)} mismatches.")
      if (length(mism) > 0) {
        # transitions of is_day: index of the first row after a change
        trans <- which(diff(out$is_day) != 0) + 1
        nearest <- trans[findInterval(mism, c(-Inf, (trans[-1] + trans[-length(trans)]) / 2, Inf))]
        dist <- ifelse(mism >= nearest, mism - nearest, nearest - 1 - mism) # 0 = adjacent to the change
        per_trans <- table(nearest)
        cli::cli_inform(c(
          "Mismatches by distance from the nearest is_day change [half-hours, 0 = adjacent]:",
          paste(names(table(dist)), table(dist), sep = ": ", collapse = "; "),
          "Max mismatches at one sunrise/sunset: {max(per_trans)} ({sum(per_trans > 1)} transitions with more than one)."
        ))
      }
    }
  }

  cli::cli_h2("u* filter summary across cases (gt5km / gt10km: pblh_m > 5 / 10 km)")
  print_tbl(do.call(rbind, lapply(names(outs), function(k) ustar_report(d, outs[[k]], runs[[k]]))))

  cli::cli_h2("L consistency summary across cases (checked: L_check consistent / sign / magnitude)")
  print_tbl(do.call(rbind, lapply(names(outs), function(k) l_report(outs[[k]], runs[[k]]))))

  # Default case: where and when L is replaced, and every value > 5 km that remains.
  od <- outs$default
  if (!is.null(od)) {
    rp <- od$flag_L_replaced
    cli::cli_h2("Default: L replacements by regime, day/night, year; days with most replacements")
    cnt <- function(x, name) count_tbl(x, name)
    print_tbl(merge(cnt(od$regime[rp], "regime"), cnt(od$regime[od$L_check %in% c("consistent", "sign", "magnitude")],
                                                      "regime"), by = "regime", suffixes = c("_replaced", "_checked")))
    print_tbl(cnt(ifelse(od$is_day[rp], "day", "night"), "is_day"))
    print_tbl(cnt(format(od$timestamp_start_utc[rp], "%Y"), "year"))
    top <- sort(table(od$day_id[rp]), decreasing = TRUE)[1:10]
    print_tbl(data.frame(day_id = names(top), n_replaced = as.integer(top)))
    big <- od[!is.na(od$pblh_m) & od$pblh_m > 5000, ]
    cli::cli_h2("Default: rows with pblh_m > 5 km ({nrow(big)})")
    if (nrow(big) > 0) {
      print_tbl(dplyr::transmute(big, t_start = format(timestamp_start_utc, "%Y-%m-%d %H:%M"), pblh_m = round(pblh_m),
                                 regime, method, L = signif(L, 3), L_l2 = signif(L_l2, 3), L_forcing = signif(L_forcing, 3),
                                 L_check, wtheta = signif(wtheta, 3), ustar = round(ustar, 3), H_source))
    }
  }

  # Regression chain (cases not run in this call are skipped).
  check_reference(outs$nofilter_l2, file.path(out_dir, "pre_spec05", paste0(site_id, "_pblhmod_default.csv.gz")),
                  "nofilter_l2 vs spec-02 default (pre_spec05)")
  check_reference(outs$l_L2, file.path(out_dir, "pre_spec06", paste0(site_id, "_pblhmod_default.csv.gz")),
                  "l_L2 vs spec-05 default (pre_spec06)")
  check_reference(outs$l_check_allday, file.path(out_dir, "pre_spec06a2", paste0(site_id, "_pblhmod_default.csv.gz")),
                  "l_check_allday vs spec-06 amendment-1 default (pre_spec06a2)")
  # Spec 07 step 0 (reference period = FLUXES range): every case vs its output before the change.
  for (k in names(outs)) {
    ref_k <- file.path(out_dir, "pre_spec07", paste0(site_id, "_pblhmod_", k, ".csv.gz"))
    if (file.exists(ref_k)) check_reference(outs[[k]], ref_k, paste0(k, " vs before the reference period (pre_spec07)"))
  }
  # Spec 11 B: every existing case vs its output before the gamma_source option.
  for (k in names(outs)) {
    ref_k <- file.path(out_dir, "pre_spec11", paste0(site_id, "_pblhmod_", k, ".csv.gz"))
    if (file.exists(ref_k)) check_reference(outs[[k]], ref_k, paste0(k, " vs before gamma_source (pre_spec11)"))
  }
  # gamma001 equals default while the default constant is 0.01 (G 2026-10-30).
  check_reference(outs$gamma001, file.path(out_dir, paste0(site_id, "_pblhmod_default.csv.gz")),
                  "gamma001 vs default (equal while the default gamma is 0.01)")

  cli::cli_h2("Side by side: quantiles per regime [m] (rows valid in all cases and in the operational PBLH)")
  reg <- outs[[1]]$regime # regimes do not depend on combine / neutral_method / gamma
  all_ok <- Reduce(`&`, lapply(outs, function(o) !is.na(o$pblh_m))) & !is.na(d$PBLH)
  series <- c(lapply(outs, `[[`, "pblh_m"), list(PBLH_oper = d$PBLH))
  side <- do.call(rbind, lapply(c("all", "stable", "unstable_diag", "convective"), function(r) {
    sel <- all_ok & (r == "all" | reg %in% r)
    do.call(rbind, lapply(names(series), function(s) {
      q <- stats::quantile(series[[s]][sel], c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
      data.frame(regime = r, case = s, n = sum(sel), p5 = round(q[1]), p50 = round(q[2]), p95 = round(q[3]))
    }))
  }))
  print_tbl(side)
}
