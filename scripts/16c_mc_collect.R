# Spec 16 (phase 5b): footprint uncertainty from the FFP inputs. Step 16c: tables and figures from the 16a draws and
# the 16b runs.
#   A  one-at-a-time: per row and input, the change at -1 / +1 sigma of the TA fraction [pp] and of ln x80 against the
#      row's base run; the half-range as the input's 1-sigma effect; the base run against the spec-15 value of the row
#      (08b chunks; tolerance 0.01 pp, counted and printed). Amendment 1: hc and z0 also with the etc
#      roughness-sublayer test off (validity rsl_off), reported next to the main runs (items <input>|rsl_off,
#      share_rsl_flip) and used for hc and z0 in 16_input_ranking; their share of runs failing the test under etc per
#      site in 16s_validity_flips (close-out 2026-11-19: per site, no network median in the main figure).
#   B  joint draws: per row and h option, quantiles p05 ... p95 and the 1-sigma width (p84 - p16) / 2.
#   C  attribution: squared standardised regression coefficients per stratum (and R^2).
#   Section 4: 1-sigma effects vs the PBLH-source effect (|ERA5 - PBLHmod|, |PBLHmod-clim - PBLHmod|), ratio to the
#      joint width and shares of rows above 1 and 2 widths; ICOS criteria under uncertainty (probability of a pass).
#   Site statistics over all strata use the stratum weights of 16a (share of the eligible rows); network = median over
#   the sites of the site values.
# Output: <output_dir>/tables/16_oat_<site>.csv, 16_mc_<site>.csv, 16_attribution_<site>.csv,
#   16_network_uncertainty.csv, 16_icos_uncertainty.csv (each with a '#' header); figures/16_input_ranking.png,
#   16_input_ranking_x80.png, 16s_mc_width.png, 16s_validity_flips.png. Runs of a period (PC smoke test) are labelled
#   "test, one month" in the figures and headers.
#
# Usage: Rscript scripts/16c_mc_collect.R [SITE ...]  (default: every site with a draws file)
#        Rscript scripts/16c_mc_collect.R --from-tables=DIR [--fp-table=FILE]
#          re-renders 16_input_ranking, 16_input_ranking_x80 and 16s_validity_flips from collected tables in DIR
#          (16_network_uncertainty.csv, 16_oat_<site>.csv; no draws or runs needed; 16s_mc_width is not redrawn).
#          Tables written before w_row was added: row weights from the eligible rows per stratum of the spec-15
#          table 09_network_fp_target_area.csv (set etc, subset valid; default DIR/09_network_fp_target_area.csv).
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), c("--from-tables", "--fp-table"))
if (length(bad)) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
cfg <- get_cfg()
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
mc_dir <- file.path(derived, "footprint_mc")
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
oat_types <- readr::cols(input = "c", fail = "c", validity = "c", stratum = "c", timestamp_start_utc = "c")

# Labels: per site the reference year (or the test period); figures say "reference year per site" and list the years
run_label <- function(lab_site) {
  if (all(grepl("^reference year", lab_site))) "reference year per site" else paste(unique(lab_site), collapse = "; ")
}
sites_text <- function(lab_site) {
  if (all(grepl("^reference year", lab_site))) {
    paste0("Sites (reference year): ", paste(names(lab_site), sub("^reference year ", "", lab_site), collapse = ", "))
  } else {
    paste0("Sites: ", paste0(names(lab_site), " (", lab_site, ")", collapse = ", "))
  }
}

# Figures from the network table and the design-A tables (with w_row); mcs / srcs for 16s_mc_width (NULL: skipped)
mc_figures <- function(tab, oat_all, lab_site, mcs = NULL, srcs = NULL) {
  lab <- run_label(lab_site)
  fs <- mc_flip_shares(oat_all)
  top <- fs[fs$input %in% c("hc", "z0f") & fs$share_rsl > 0, ]
  top <- top[order(-top$share_rsl), ]
  flip_txt <- if (nrow(top)) {
    paste0(" (", paste0(top$site_id, " ", style_of(top$input, "label", "mc_input"), " ", round(100 * top$share_rsl), " %",
                        collapse = ", "), ")")
  } else " (none)"
  cap <- paste(sites_text(lab_site), "- stratum = zm/L class x day/night;",
               "1σ effect = half-range of the one-at-a-time runs at ±1σ; small grey points: sites (median over rows); large point and bar: network median and IQR over sites.",
               "Lines (PBLH source): network median of the stratum of |ERA5 − PBLHmod| and |PBLHmod-clim − PBLHmod| (spec 15).",
               "hc and z0: effects from the runs with the etc roughness-sublayer test (zm ≥ m z0) switched off ('rsl test off'), all other validity",
               "conditions kept. Under the etc set a share of these runs fails that test; the share is per site in 16s_validity_flips",
               paste0(flip_txt, ":"),
               "a failing run leaves the half-hour without a valid footprint, so an effect near zero there is not a zero uncertainty.")
  p <- plot_mc_input_ranking(tab, "ta_pp", paste0("Footprint uncertainty from the FFP inputs: target-area fraction (", lab, ")"), cap)
  save_figure(p, "16_input_ranking", 22, 15)
  p <- plot_mc_input_ranking(tab, "lnx80", paste0("Footprint uncertainty from the FFP inputs: ln x80 (", lab, ")"), cap)
  save_figure(p, "16_input_ranking_x80", 22, 15)
  if (!is.null(mcs) && nrow(mcs)) {
    p <- plot_mc_width(mcs, srcs, "h_pct30", subtitle = paste0("Joint 1σ width vs PBLH-source effect per row (", lab, ")"),
                       caption = paste(sites_text(lab_site), "- values below 0.001 pp drawn at 0.001 pp (log axes)."))
    save_figure(p, "16s_mc_width", 22, 13)
  }
  p <- plot_mc_validity_flips(fs, paste0("Validity flips under one-at-a-time perturbation, per site (", lab, ")"),
                              paste(sites_text(lab_site), "- base rows valid under the etc set; a perturbed run is invalid if any condition of the set",
                                    "fails with the perturbed inputs (z0 enters only here). Share of the runs at −1σ and +1σ together, stratum-weighted",
                                    "over the sampled rows (16a weights); in brackets the share failing the roughness-sublayer test (zm ≥ m z0)."))
  save_figure(p, "16s_validity_flips", 22, 13)
  cli::cli_alert_success("Figures 16_input_ranking, 16_input_ranking_x80, 16s_validity_flips{if (!is.null(mcs)) ', 16s_mc_width'} written ({lab}).")
  print(as.data.frame(dplyr::mutate(fs[fs$share_invalid > 0, ], dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))),
        row.names = FALSE)
  invisible(fs)
}

# --- re-render from collected tables ---
from <- opt("from-tables")
if (!is.null(from)) {
  if (!dir.exists(from)) cli::cli_abort("{.path {from}} not found.")
  tab <- readr::read_csv(file.path(from, "16_network_uncertainty.csv"), comment = "#", show_col_types = FALSE, progress = FALSE)
  fo <- list.files(from, pattern = "^16_oat_.*\\.csv$", full.names = TRUE)
  if (!length(fo)) cli::cli_abort("No 16_oat_<site>.csv in {.path {from}}.")
  lab_site <- character()
  oat_all <- dplyr::bind_rows(lapply(fo, function(f) {
    s <- sub("^16_oat_(.*)\\.csv$", "\\1", basename(f))
    h1 <- readLines(f, n = 1)
    lab_site[s] <<- sub("^[^(]*\\(([^)]*)\\).*$", "\\1", h1)
    dplyr::mutate(readr::read_csv(f, comment = "#", col_types = oat_types, progress = FALSE), site_id = s, .before = 1)
  }))
  if (!"w_row" %in% names(oat_all) || anyNA(oat_all$w_row)) {
    fpt <- opt("fp-table", file.path(from, "09_network_fp_target_area.csv"))
    if (!file.exists(fpt)) cli::cli_abort("Tables without w_row: pass --fp-table= (spec-15 09_network_fp_target_area.csv).")
    ta <- readr::read_csv(fpt, comment = "#", show_col_types = FALSE, progress = FALSE)
    ta <- ta[ta$scope == "site" & ta$set == "etc" & ta$subset == "valid" & ta$pair == "era5_vs_pblhmod" &
               ta$regime != "all" & ta$daynight != "all", ]
    counts <- dplyr::tibble(site_id = ta$site_id, stratum = paste(ta$regime, ta$daynight, sep = "_"), n_stratum = ta$n)
    oat_all <- mc_weights_from_counts(oat_all, counts)
    cli::cli_inform("Row weights from the eligible rows per stratum of {.path {fpt}} (counts match the samples).")
  }
  mc_figures(tab, oat_all, lab_site)
  quit(save = "no")
}

site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
if (!length(site_ids)) site_ids <- sub("_draws\\.rds$", "", list.files(mc_dir, pattern = "_draws\\.rds$"))
if (!length(site_ids)) cli::cli_abort("No draws file in {.path {mc_dir}}: run 16a and 16b first.")

tol_pp <- 0.01
all_sum <- list(); all_oat <- list(); all_mc <- list(); all_src <- list(); all_icos <- list(); labels <- character()
for (site_id in site_ids) {
  draws <- readRDS(file.path(mc_dir, paste0(site_id, "_draws.rds")))
  f <- list.files(file.path(mc_dir, site_id), pattern = "_runs\\.csv\\.gz$", full.names = TRUE)
  if (!length(f)) { cli::cli_alert_warning("{site_id}: no 16b runs, skipped."); next }
  runs <- dplyr::bind_rows(lapply(f, readr::read_csv, show_col_types = FALSE, progress = FALSE,
                                  col_types = readr::cols(input = "c", h_option = "c", fail = "c", run = "c", design = "c",
                                                           validity = "c")))
  x <- mc_join_runs(draws, runs)
  lab <- if (is.null(draws$period)) paste("reference year", draws$year) else paste("test, one month", draws$period)
  labels[site_id] <- lab
  hdr0 <- c(paste0("Spec 16 footprint uncertainty, ", site_id, " (", lab, "); seed ", draws$seed, " (", paste(draws$rng_kind, collapse = ", "), ")."),
            "Inputs perturbed relative to the published L2 values (1 sigma: sigma_v 15 %, u* 10 %, L 20 % via 1/L, umean 5 %, hc 10 %, z0 factor 50 %, WD 7.5 deg, h 30 %); lognormal factors with median 1 and that CV.",
            "Validity: set etc recomputed with the perturbed inputs (ETC z0 multiplier of 08a); invalid runs not computed.",
            "Validity rsl_off (spec 16 amendment 1): the hc and z0 runs of design A again with the etc roughness-sublayer test (zm >= m z0) off, all other conditions kept; rsl_fail = that test fails; valid_run = validity under the run's own set.")
  oat <- mc_oat_table(x)
  write_csv_header(oat, file.path(tdir, paste0("16_oat_", site_id, ".csv")), c(hdr0,
    "Design A: d_ta_pp = TA fraction of the run minus the row's base run [pp]; d_lnx80 = ln(x80 run / x80 base) [-] (x80: crosswind-integrated, ffp_ci_paper(), as 08b); side -1 / +1 = quantiles 0.1587 / 0.8413.",
    "w_row = row weight of 16a (stratum share of the eligible rows / rows sampled in the stratum; sum 1 over the rows of the site).",
    "base_diff_pp: base run minus the spec-15 PBLHmod value of the row (08b chunk) [pp]."))
  b0 <- oat[oat$input == "base", ]
  nbad <- sum(abs(b0$base_diff_pp) > tol_pp, na.rm = TRUE)
  cli::cli_inform("{site_id} ({lab}): {length(unique(oat$row_id))} rows A; base run vs spec-15 value: max |diff| {signif(max(abs(b0$base_diff_pp), na.rm = TRUE), 3)} pp, {nbad} row{?s} above {tol_pp} pp; x80 max |diff| {signif(max(abs(b0$base_diff_x80_m), na.rm = TRUE), 3)} m.")
  if (nbad > 0) cli::cli_alert_danger("{site_id}: base runs differ from the spec-15 values: check the inputs and the FFP code version.")
  eff <- mc_oat_effects(oat)
  src <- mc_source_effects(draws$sample)
  mc <- mc_mc_table(x)
  if (nrow(mc)) {
    write_csv_header(mc, file.path(tdir, paste0("16_mc_", site_id, ".csv")), c(hdr0,
      paste0("Design B: ", draws$n_draws, " joint draws per row (L' = L0 (u*'/u*0)^3 / eps_L); quantiles of the valid draws; width = (p84 - p16) / 2. h options: h_pct30 (30 %), h_empirical (day sigma_ln ",
             signif(draws$sigma_ln_empirical, 3), ", night 30 %).")))
    at <- mc_attribution_table(x, draws$spec$input)
    write_csv_header(at, file.path(tdir, paste0("16_attribution_", site_id, ".csv")), c(hdr0,
      "Design C: share = squared standardised regression coefficient of the response (TA fraction [pp] or ln x80) on the standard deviates of the inputs, within rows (row means removed); the L-u* link is attributed to u*; r2 of the linear model (low r2 = non-linear response: rely on design A)."))
  }
  s <- mc_site_summary(site_id, eff, src, mc)
  all_sum[[site_id]] <- s
  all_oat[[site_id]] <- dplyr::mutate(oat, site_id = site_id, .before = 1)
  all_mc[[site_id]] <- if (nrow(mc)) dplyr::mutate(mc, site_id = site_id, .before = 1)
  all_src[[site_id]] <- dplyr::mutate(src, site_id = site_id, .before = 1)
  all_icos[[site_id]] <- mc_icos_table(x, src, site_id)
  fl <- dplyr::summarise(dplyr::group_by(oat[oat$input != "base" & oat$validity == "etc", ], input),
                         share_invalid = round(mean(!valid_etc), 3), .groups = "drop")
  cli::cli_inform("{site_id}: invalid perturbed runs by input (A, etc set, unweighted): {paste(fl$input, fl$share_invalid, collapse = ', ')}.")
  for (inp in c("hc", "z0f")) { # spec 16 amendment 1: rsl_off next to the main runs
    ro <- s[s$kind == "oat_1sigma" & s$stratum == "all" & s$item %in% c(inp, paste0(inp, "|rsl_off")), ]
    v <- function(it, q) {
      r <- ro[ro$item == it & ro$quantity == q, ]
      if (nrow(r)) paste0(signif(r$median, 3), " (n ", r$n, ")") else "NA"
    }
    cli::cli_inform(paste0(site_id, " ", inp, ": 1-sigma effect, all strata: TA [pp] etc ", v(inp, "ta_pp"), ", rsl_off ",
                           v(paste0(inp, "|rsl_off"), "ta_pp"), "; ln x80 etc ", v(inp, "lnx80"), ", rsl_off ",
                           v(paste0(inp, "|rsl_off"), "lnx80"), "; rsl flips ", round(100 * ro$share_rsl_flip[1]), " %."))
  }
}
if (!length(all_sum)) cli::cli_abort("No site with runs.")
lab_hdr <- if (all(grepl("^reference year", labels))) {
  paste0("reference year per site: ", paste(names(labels), sub("^reference year ", "", labels), collapse = ", "))
} else paste(unique(labels), collapse = "; ")
s <- dplyr::bind_rows(all_sum)
net <- mc_network_summary(s)
tab <- dplyr::bind_rows(s, net)
write_csv_header(tab, file.path(tdir, "16_network_uncertainty.csv"), c(
  paste0("Spec 16 section 4 (", lab_hdr, "): kind oat_1sigma = half-range 1-sigma effect of an input (design A); source_effect = |first minus second| of the spec-15 values (pairs era5_vs_pblhmod, pblhmod_clim_vs_pblhmod);"),
  "joint_width = (p84 - p16) / 2 of the joint draws (design B, per h option); ratio = |source effect| / joint width, share_gt_1s / share_gt_2s = share of rows above 1 and 2 widths.",
  "item <input>|rsl_off = the input's runs with the etc roughness-sublayer test off (spec 16 amendment 1); share_rsl_flip = share of the input's -1 / +1 sigma runs that fail that test (oat_1sigma rows).",
  "quantity ta_pp [pp] or lnx80 [-]; stratum all = stratum-weighted over the rows of the site; scope network = median (p25, p75) over the sites of the site medians, n = sites."))
icos <- dplyr::bind_rows(all_icos)
if (nrow(icos)) write_csv_header(icos, file.path(tdir, "16_icos_uncertainty.csv"), c(
  paste0("Spec 16 section 4 (", lab_hdr, "): ICOS criteria under uncertainty from the joint draws (B rows, stratum-weighted): rebmann_ta90 = share of rows with TA >= 90 %, pass if >= 0.5 (Rebmann et al. 2018);"),
  "icos_lab_ta70 = share of the valid_qc rows with TA >= 70 %, pass if >= 0.7 (ICOS labelling rule, ETC); prob_pass = share of draws that pass; uncertain = 0.1 < prob_pass < 0.9; cases era5, pblhmod_clim: draws of PBLHmod plus the row's spec-15 difference."))

# --- figures ---
mc_figures(tab, dplyr::bind_rows(all_oat), labels, dplyr::bind_rows(all_mc), dplyr::bind_rows(all_src))

# printed summary: input ranking (all strata) and source effect
pr <- tab[tab$stratum == "all" & tab$kind %in% c("oat_1sigma", "source_effect", "joint_width"), ]
pr <- pr[order(pr$scope, pr$site_id, pr$quantity, -pr$median), c("scope", "site_id", "quantity", "kind", "item", "n", "median", "p25", "p75")]
print(as.data.frame(dplyr::mutate(pr, dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
