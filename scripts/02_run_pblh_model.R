# Spec 02: run PBLHmod v2 for the test site(s) and print a summary.
# Runs: default (bg_only, neutral auto, gamma 0.01 K/m); gamma0025 (legacy gamma);
# max_B2 (combine max, neutral B2: behaviour before amendment 11); max_auto;
# u* filter sensitivity (spec 05, spec 02 amendment 12): ustar04 / ustar06 (u*/WS limit 0.4 / 0.6) and
# nofilter (filter disabled = previous default). All other cases use the filter defaults.
# Ends with a u* filter summary across cases and a side-by-side table of quantiles per regime.
# Output: <derived_dir>/pblhmod/<site>_pblhmod_<label>.csv.gz
# If <derived_dir>/pblhmod/pre_spec05/<site>_pblhmod_default.csv.gz exists (spec-02 output saved before
# spec 05), the nofilter case is compared with it (must be identical, < 1e-9 m).
#
# Usage: Rscript scripts/02_run_pblh_model.R [SITE ...]
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
out_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "pblhmod")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- list(
  default = pblh_options(),
  gamma0025 = pblh_options(gamma = 0.0025, label = "gamma0025"),
  max_B2 = pblh_options(combine = "max", neutral_method = "B2", label = "max_B2"),
  max_auto = pblh_options(combine = "max", neutral_method = "auto", label = "max_auto"),
  ustar04 = pblh_options(ustar_ws_max = 0.4, label = "ustar04"),
  ustar06 = pblh_options(ustar_ws_max = 0.6, label = "ustar06"),
  nofilter = pblh_options(ustar_ws_max = Inf, ustar_max = Inf, label = "nofilter")
)

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
    tm <- system.time(out <- run_pblh_model(d, meta, opts))
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

  ref <- file.path(out_dir, "pre_spec05", paste0(site_id, "_pblhmod_default.csv.gz"))
  if (file.exists(ref)) {
    cli::cli_h2("nofilter vs spec-02 default saved before spec 05")
    old <- readr::read_csv(ref, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    nf <- outs$nofilter
    same_rows <- nrow(old) == nrow(nf) && all(old$timestamp_start_utc == nf$timestamp_start_utc)
    cmp <- vapply(c("pblh_m", "h_diag_m", "h_bg_m", "ustar", "L"), function(k) {
      same_na <- identical(is.na(old[[k]]), is.na(nf[[k]]))
      if (!same_na) return(Inf)
      max(c(0, abs(old[[k]] - nf[[k]])), na.rm = TRUE)
    }, numeric(1))
    print_tbl(data.frame(column = names(cmp), max_abs_diff = cmp))
    chr_same <- all(vapply(c("regime", "method"), function(k) identical(old[[k]], nf[[k]]), logical(1)))
    if (same_rows && chr_same && all(cmp < 1e-9)) {
      cli::cli_alert_success("nofilter reproduces the spec-02 default (max |diff| < 1e-9 m, same NA pattern, same regime/method).")
    } else {
      cli::cli_alert_danger("nofilter differs from the spec-02 default (rows {same_rows}, regime/method {chr_same}).")
    }
  }

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
