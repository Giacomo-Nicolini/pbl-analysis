# Footprint uncertainty from the FFP inputs (spec 16, phase 5b; G 2026-11-17).
# Perturbations of the FFP inputs around the published L2 values of a stratified sample of half-hours, one-at-a-time
# (design A) and joint Monte Carlo (design B), with the PBLH-source effect of spec 15 as the reference. FFP model:
# Kljun et al. (2015), `Kljun2015`, run through the same wrapper and settings as step 08b ([ffp_halfhour_run()],
# [ffp_ci_paper()]); validity under perturbation with [footprint_validity()] (set `etc`).
#
# Notes on the FFP code (umean only, z0 = NaN; spec 08): the crosswind-integrated distances x_r depend on zm, h,
# umean and u* only (Kljun et al. 2015 Eq. 6 with umean), so sigma_v, L and the wind direction change x_r only
# through validity; in the 2D field, L enters only through the sign of zm/L in the sigma_y scale constant
# (0.8 unstable, 0.55 stable; the 1e-5 |L / zm| term matters only for |L| >~ 1e4 zm, and |L| > 5000 m is set to
# neutral, `ol_n`), so a sign-preserving perturbation of L acts on the footprint mainly through validity
# (zm/L > -15.5, zm/L <= 1) and the stratum.

# --- inputs and draws -------------------------------------------------------------------------------------------------

#' Input uncertainties of spec 16 (section 2; G 2026-11-17)
#'
#' One standard deviation per input, relative to the published L2 value of the row. Lognormal factors use
#' sigma_ln = sqrt(ln(1 + cv^2)) ([mc_sigma_ln()]): median 1, coefficient of variation `cv`.
#'
#' @param sigmav_cv 1 sigma of sigma_v [-] (0.15; sensitivity 0.10, 0.20).
#' @param h_cv 1 sigma of h [-] (0.30, `h_pct30`).
#' @return Tibble: `input` (`sigmav`, `ustar`, `L`, `umean`, `hc`, `z0f`, `wd`, `h`), `kind` (`lognormal`,
#'   `normal_trunc` = multiplicative normal truncated at +-3 sigma, `additive` [deg]), `sd` [-, or deg for `wd`],
#'   `label`, `note`.
mc_input_spec <- function(sigmav_cv = 0.15, h_cv = 0.30) {
  dplyr::tibble(
    input = c("sigmav", "ustar", "L", "umean", "hc", "z0f", "wd", "h"),
    kind = c("lognormal", "lognormal", "lognormal", "lognormal", "normal_trunc", "lognormal", "additive", "lognormal"),
    sd = c(sigmav_cv, 0.10, 0.20, 0.05, 0.10, 0.50, 7.5, h_cv),
    label = c("σv", "u*", "L", "umean", "hc (d, zm, z0)", "z0 factor", "wind direction", "h"),
    note = c("V_SIGMA, multiplicative", "USTAR, multiplicative", "MO_LENGTH via 1/L; joint draws: L' = L0 (u*'/u*0)^3 / eps_L",
             "WS, multiplicative", "canopy height: d = 2/3 hc, zm = EC height - d, z0 = 0.15 hc",
             "factor on 0.15 in z0 = 0.15 hc (validity only: the FFP runs with umean, z0 = NaN)",
             "WD, additive [deg], wrapped to [0, 360)", "PBLHmod h, multiplicative (h_pct30; h_empirical by day)")
  )
}

#' Log-standard deviation of a lognormal factor with median 1 and a given coefficient of variation
#'
#' For a lognormal factor F = exp(sigma_ln Z), CV(F) = sqrt(exp(sigma_ln^2) - 1); hence sigma_ln = sqrt(ln(1 + cv^2)).
#'
#' @param cv Coefficient of variation [-].
#' @return sigma_ln [-].
mc_sigma_ln <- function(cv) sqrt(log1p(cv^2))

#' Standard normal deviate at a probability, optionally truncated at +-3
#'
#' @param p Probability [0-1].
#' @param trunc Truncate the standard normal at +-3 [logical].
#' @return Quantile [-].
mc_qz <- function(p, trunc = FALSE) {
  if (!trunc) return(stats::qnorm(p))
  stats::qnorm(stats::pnorm(-3) + p * (stats::pnorm(3) - stats::pnorm(-3)))
}

#' Random standard deviates for the inputs (independent; hc truncated at +-3)
#'
#' @param n Number of draws [-].
#' @param spec [mc_input_spec()].
#' @return Matrix n x inputs of standard deviates (column names = `spec$input`).
mc_draw_z <- function(n, spec = mc_input_spec()) {
  z <- vapply(seq_len(nrow(spec)), function(j) {
    if (spec$kind[j] == "normal_trunc") mc_qz(stats::runif(n), trunc = TRUE) else stats::rnorm(n)
  }, numeric(n))
  z <- matrix(z, nrow = n)
  colnames(z) <- spec$input
  z
}

#' Factor (or additive deviation) of an input from its standard deviate
#'
#' @param z Standard deviate(s) [-].
#' @param kind,sd Row of [mc_input_spec()].
#' @return Multiplicative factor [-] (`lognormal`: exp(sigma_ln z); `normal_trunc`: 1 + sd z), or the additive
#'   deviation [deg] (`additive`: sd z).
mc_factor <- function(z, kind, sd) {
  switch(kind,
         lognormal = exp(mc_sigma_ln(sd) * z),
         normal_trunc = 1 + sd * z,
         additive = sd * z,
         cli::cli_abort("Unknown perturbation kind {.val {kind}}."))
}

#' Perturbed FFP inputs
#'
#' Applies the factors to the base values (spec 16 section 2). Heights follow the ETC convention of
#' [footprint_heights()]: d = 2/3 hc, zm = EC height - d, z0 = 0.15 hc times the z0 factor. L: 1/L is multiplied by
#' the L factor (L' = L0 / f_L, sign kept); with `link_L = TRUE` (joint draws) also L' = L0 (u*'/u*0)^3 / f_L
#' (L proportional to u*^3). Wind direction: additive, wrapped to [0, 360).
#'
#' @param base Tibble of base values: `ec_height_m`, `canopy_height_m` [m], `umean` [m s-1], `L` [m], `sigmav`, `ustar`
#'   [m s-1], `wind_dir` [deg], `h` [m].
#' @param f Named list or tibble of factors with the same number of rows (or length 1): `sigmav`, `ustar`, `L`,
#'   `umean`, `hc`, `z0f`, `h` [-] and the additive `wd` [deg]; missing entries = no perturbation.
#' @param link_L Link L to the u* draw [logical].
#' @return Tibble: `zm_m`, `z0_m`, `d_m`, `canopy_height_m`, `ec_height_m` [m], `umean`, `L`, `sigmav`, `ustar`,
#'   `wind_dir`, `h`.
mc_perturb <- function(base, f, link_L = FALSE) {
  g <- function(k, def) if (is.null(f[[k]])) def else f[[k]]
  hc <- base$canopy_height_m * g("hc", 1)
  d <- 2 / 3 * hc
  ustar <- base$ustar * g("ustar", 1)
  L <- base$L / g("L", 1)
  if (link_L) L <- L * (ustar / base$ustar)^3
  dplyr::tibble(
    zm_m = base$ec_height_m - d, z0_m = 0.15 * hc * g("z0f", 1), d_m = d, canopy_height_m = hc,
    ec_height_m = base$ec_height_m, umean = base$umean * g("umean", 1), L = L, sigmav = base$sigmav * g("sigmav", 1),
    ustar = ustar, wind_dir = (base$wind_dir + g("wd", 0)) %% 360, h = base$h * g("h", 1)
  )
}

#' Validity of perturbed inputs under the `etc` set
#'
#' [footprint_validity()] recomputed with the perturbed zm, z0, h, L, u*, umean, sigma_v and wind direction; the ETC
#' z0 multiplier m stays the site's 08a value; the spec-05 u* flag and the exclusion flag of the base row are kept.
#'
#' @param p Output of [mc_perturb()].
#' @param m_etc Site's ETC multiplier [-].
#' @param ustar_implausible,excluded Base-row flags [logical].
#' @return Tibble: `valid_etc` [logical] and `fail` (first failing condition of the set in the order inputs, code
#'   checks, domain, e_rsl, e_el, e_tower, e_stab, e_ustar; `NA` when valid); `valid_norsl` [logical, the `etc` set
#'   without the roughness-sublayer test e_rsl, all other conditions kept: validity of the `rsl_off` runs, spec 16
#'   amendment 1]; `rsl_fail` [logical, e_rsl zm >= m z0 fails].
mc_validity <- function(p, m_etc, ustar_implausible = FALSE, excluded = FALSE) {
  v <- footprint_validity(p$zm_m, p$z0_m, p$ec_height_m, p$h, p$umean, p$L, p$sigmav, p$ustar, p$wind_dir, m_etc,
                          ustar_implausible, excluded)
  conds <- c("inp", "code_h10", "code_ustar01", "code_stab", "code_nonneg", "domain", "e_rsl", "e_el", "e_tower",
             "e_stab", "e_ustar")
  m <- as.matrix(v[conds])
  first <- apply(!m, 1, function(x) if (any(x)) conds[which(x)[1]] else NA_character_)
  keep <- setdiff(conds, "e_rsl")
  norsl <- apply(as.matrix(v[keep]), 1, function(x) all(x %in% TRUE))
  dplyr::tibble(valid_etc = v$valid_etc, fail = first, valid_norsl = norsl, rsl_fail = v$e_rsl %in% FALSE)
}

# --- strata and sample ------------------------------------------------------------------------------------------------

#' Stratum of a half-hour: zm/L class x day/night
#'
#' @param zm Measurement height above d [m]; @param L Obukhov length [m]; @param day Daytime [logical].
#' @return Character `<class>_<day|night>` (classes of [zl_class()]: unstable zm/L < -0.1 <= neutral <= 0.1 <
#'   stable), `NA` where missing.
mc_stratum <- function(zm, L, day) {
  zc <- as.character(zl_class(zm, L))
  dn <- ifelse(day %in% TRUE, "day", ifelse(day %in% FALSE, "night", NA_character_))
  ifelse(is.na(zc) | is.na(dn), NA_character_, paste(zc, dn, sep = "_"))
}

#' The six strata in reporting order
#' @return Character.
mc_strata_levels <- function() paste(rep(c("unstable", "neutral", "stable"), each = 2), c("day", "night"), sep = "_")

#' Stratified random sample of half-hours
#'
#' `n_per` rows per stratum (all rows of a stratum with fewer), drawn with the current RNG state (set the seed
#' before). Weights: the stratum's share of the eligible rows, divided by the rows sampled in it, so that the weights
#' of a site sum to 1 and a weighted statistic represents the eligible rows of the reference year.
#'
#' @param stratum Stratum per eligible row [character].
#' @param n_per Rows per stratum [-].
#' @return Tibble: `i` (row position in `stratum`), `stratum`, `n_stratum` [-, eligible rows], `n_sampled` [-],
#'   `w_stratum` [0-1], `w_row` [0-1], ordered by stratum then position.
mc_sample <- function(stratum, n_per) {
  ok <- !is.na(stratum)
  n_all <- sum(ok)
  lev <- intersect(mc_strata_levels(), unique(stratum[ok]))
  out <- lapply(lev, function(s) {
    pos <- which(stratum %in% s)
    pick <- if (length(pos) <= n_per) pos else sort(pos[sample.int(length(pos), n_per)])
    dplyr::tibble(i = pick, stratum = s, n_stratum = length(pos), n_sampled = length(pick),
                  w_stratum = length(pos) / n_all, w_row = length(pos) / n_all / length(pick))
  })
  dplyr::bind_rows(out)
}

#' Subsample of the design-B rows: the first `n_per` sampled rows of each stratum in a random order
#'
#' @param smp Output of [mc_sample()].
#' @param n_per Rows per stratum [-].
#' @return Logical per row of `smp`.
mc_subsample <- function(smp, n_per) {
  keep <- rep(FALSE, nrow(smp))
  for (s in unique(smp$stratum)) {
    pos <- which(smp$stratum == s)
    keep[pos[sample.int(length(pos), min(n_per, length(pos)))]] <- TRUE
  }
  keep
}

#' Seed of a site (fixed, documented in every output)
#'
#' @param site_id Station code [character].
#' @param base Base seed [integer].
#' @return Integer seed: `base` + sum of the character codes weighted by position.
mc_site_seed <- function(site_id, base = 20261117L) {
  u <- utf8ToInt(site_id)
  as.integer(base + sum(u * seq_along(u)))
}

# --- run tables -------------------------------------------------------------------------------------------------------

#' One-at-a-time run table (design A)
#'
#' Base run plus each input at the 0.1587 and 0.8413 quantiles of its distribution (about -1 and +1 sigma), the others
#' at base: 1 + 2 x 8 = 17 runs per row under the `etc` set. Spec 16 amendment 1 (G 2026-11-18): the runs of the
#' inputs in `rsl_off` (hc, z0) again with the roughness-sublayer test of the `etc` set switched off (validity
#' `rsl_off`, all other conditions kept): 4 more runs, 21 per row. Their FFP inputs equal those of the `etc` runs of
#' the same input (z0: the base run), so they cost an FFP call only where the `etc` run fails e_rsl and hc changes zm.
#'
#' @param spec [mc_input_spec()].
#' @param rsl_off Inputs also run without the e_rsl test [character].
#' @return Tibble: `run` (`base`, `<input>_minus`, `<input>_plus`, `<input>_minus_rsl_off`, `<input>_plus_rsl_off`),
#'   `input` (`NA` for base), `side` (`-1`, `+1`, `0`), `z` [-], `validity` (`etc`, `rsl_off`).
mc_oat_runs <- function(spec = mc_input_spec(), rsl_off = c("hc", "z0f")) {
  p <- c(0.1587, 0.8413)
  r <- dplyr::bind_rows(lapply(seq_len(nrow(spec)), function(j) {
    dplyr::tibble(run = paste0(spec$input[j], c("_minus", "_plus")), input = spec$input[j], side = c(-1L, 1L),
                  z = mc_qz(p, trunc = spec$kind[j] == "normal_trunc"), validity = "etc")
  }))
  ro <- r[r$input %in% rsl_off, ]
  ro$run <- paste0(ro$run, "_rsl_off")
  ro$validity <- "rsl_off"
  dplyr::bind_rows(dplyr::tibble(run = "base", input = NA_character_, side = 0L, z = 0, validity = "etc"), r, ro)
}

#' Upper bound of the distinct FFP input sets of design A per row (chunk sizing)
#'
#' @param oat [mc_oat_runs()].
#' @return Integer: the `etc` runs plus the `rsl_off` runs that change an FFP input (hc; z0 enters validity only).
mc_oat_ffp_runs <- function(oat) {
  as.integer(sum(oat$validity == "etc") + sum(oat$validity == "rsl_off" & oat$input != "z0f"))
}

#' Factors of one run of a row
#'
#' @param spec [mc_input_spec()].
#' @param z Named standard deviates (missing inputs = 0).
#' @param h_sd Log or relative sd of h for this run [-] (`NULL` = `spec`).
#' @return Named list of factors (additive deviation for `wd`).
mc_factors <- function(spec, z, h_sd = NULL) {
  out <- list()
  for (j in seq_len(nrow(spec))) {
    k <- spec$input[j]
    zk <- if (k %in% names(z)) z[[k]] else 0
    sd <- if (k == "h" && !is.null(h_sd)) h_sd else spec$sd[j]
    out[[k]] <- mc_factor(zk, spec$kind[j], sd)
  }
  out
}

#' Perturbed inputs of all runs of a chunk
#'
#' @param smp Sample rows of the chunk (base values: `row_id`, `ec_height_m`, `canopy_height_m`, `umean`, `L`,
#'   `sigmav`, `ustar`, `wind_dir`, `h`, `day`, `ustar_implausible`, `excluded`).
#' @param design `"A"` or `"B"`.
#' @param draws For B: list with `z` (array draws x inputs, per row: a list of matrices named by `row_id`) and
#'   `h_options` (named: sd of h by day / night per option).
#' @param spec [mc_input_spec()].
#' @param m_etc Site's ETC multiplier [-].
#' @return Tibble, one row per run: `row_id`, `design`, `run`, `input`, `side`, `validity` (`etc`, `rsl_off`), `draw`,
#'   `h_option`, the perturbed inputs ([mc_perturb()]), `valid_etc`, `fail`, `valid_norsl`, `rsl_fail` ([mc_validity()]),
#'   `valid_run` (validity under the run's own set: `valid_norsl` for `rsl_off`, else `valid_etc`); for B also the
#'   standard deviates `z_<input>`.
mc_run_inputs <- function(smp, design, draws = NULL, spec = mc_input_spec(), m_etc) {
  out <- list()
  for (r in seq_len(nrow(smp))) {
    b <- smp[r, ]
    if (design == "A") {
      oat <- mc_oat_runs(spec)
      f <- lapply(seq_len(nrow(oat)), function(k) {
        z <- if (is.na(oat$input[k])) list() else stats::setNames(list(oat$z[k]), oat$input[k])
        mc_factors(spec, z)
      })
      fac <- dplyr::bind_rows(lapply(f, dplyr::as_tibble))
      p <- mc_perturb(b[rep(1, nrow(oat)), ], fac, link_L = FALSE)
      x <- dplyr::bind_cols(dplyr::tibble(row_id = b$row_id, design = "A", run = oat$run, input = oat$input,
                                          side = oat$side, validity = oat$validity, draw = NA_integer_,
                                          h_option = NA_character_), p)
    } else {
      zr <- draws$z[[as.character(b$row_id)]]
      x <- dplyr::bind_rows(lapply(names(draws$h_options), function(ho) {
        hsd <- draws$h_options[[ho]][[if (b$day %in% TRUE) "day" else "night"]]
        fac <- dplyr::bind_rows(lapply(seq_len(nrow(zr)), function(k) dplyr::as_tibble(mc_factors(spec, as.list(zr[k, ]), hsd))))
        p <- mc_perturb(b[rep(1, nrow(zr)), ], fac, link_L = TRUE)
        zz <- dplyr::as_tibble(zr)
        names(zz) <- paste0("z_", names(zz))
        dplyr::bind_cols(dplyr::tibble(row_id = b$row_id, design = "B", run = paste0(ho, "_", seq_len(nrow(zr))),
                                       input = NA_character_, side = NA_integer_, validity = "etc",
                                       draw = seq_len(nrow(zr)),
                                       h_option = ho), p, zz)
      }))
    }
    v <- mc_validity(x, m_etc, b$ustar_implausible %in% TRUE, b$excluded %in% TRUE)
    v$valid_run <- ifelse(x$validity == "rsl_off", v$valid_norsl, v$valid_etc)
    out[[r]] <- dplyr::bind_cols(x, v)
  }
  dplyr::bind_rows(out)
}

#' Key of the FFP inputs of a run (identical inputs are computed once)
#'
#' @param x Runs ([mc_run_inputs()]).
#' @return Character key of zm, h, umean, L, sigma_v, u*, wind direction (full precision).
mc_run_key <- function(x) {
  paste(sprintf("%.17g", x$zm_m), sprintf("%.17g", x$h), sprintf("%.17g", x$umean), sprintf("%.17g", x$L),
        sprintf("%.17g", x$sigmav), sprintf("%.17g", x$ustar), sprintf("%.17g", x$wind_dir))
}

# --- h_empirical -------------------------------------------------------------------------------------------------------

#' Spread of ln(PBLHmod / measured) per day at an aerosol site (spec 16, h_empirical)
#'
#' The site's primary aerosol source on its comparison set (apparent solar window of
#' [measured_compare_options()], PBLHmod convective), on the measurement time base ([measured_time_base()]); daily
#' median of ln(PBLHmod / measured); robust sd = IQR / 1.349 of the daily values.
#'
#' @param meas Measured standard table; @param mod PBLHmod standard table (`x_regime`); @param era Named list with
#'   `era5_blh_center` (needed by the time base); @param lon Site longitude [deg E].
#' @param opts [measured_compare_options()].
#' @return One-row tibble: `n_days`, `sd_ln_daily` [-], `median_ln_daily` [-].
mc_daily_log_spread <- function(meas, mod, era, lon, opts = measured_compare_options()) {
  step <- stats::median(as.numeric(meas$timestamp_end_utc) - as.numeric(meas$timestamp_start_utc), na.rm = TRUE)
  tb <- measured_time_base(meas, mod, era, step_s = step)
  ref <- tb$meas[tb$meas$kind %in% "exact", c("timestamp_utc", "pblh_m")]
  pt <- pair_table(tb$mod, ref, "x")
  a <- apparent_solar_hour(pt$timestamp_utc + tb$mid_s, lon)
  pt <- pt[pt$regime %in% "convective" & a >= opts$aerosol_window[1] & a < opts$aerosol_window[2], ]
  if (nrow(pt) == 0) return(dplyr::tibble(n_days = 0L, sd_ln_daily = NA_real_, median_ln_daily = NA_real_))
  dd <- tapply(log(pt$mod / pt$ref), as.Date(pt$timestamp_utc), stats::median)
  dplyr::tibble(n_days = length(dd), sd_ln_daily = unname(stats::IQR(dd) / 1.349), median_ln_daily = stats::median(dd))
}

# --- metrics of the runs -------------------------------------------------------------------------------------------------

#' Half-range 1-sigma effect of an input (design A)
#'
#' @param d_minus,d_plus Change at -1 sigma and +1 sigma [any unit].
#' @return (|d_plus| + |d_minus|) / 2; `NA` if either side is missing (invalid run).
mc_half_range <- function(d_minus, d_plus) (abs(d_plus) + abs(d_minus)) / 2

#' Weighted quantile
#'
#' @param x Values; @param w Weights (>= 0); @param p Probabilities.
#' @return Quantiles (step function of the weighted ECDF); `NA` without values.
mc_wquantile <- function(x, w, p) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(p)))
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(q) x[o][which(cw >= q - 1e-12)[1]], numeric(1))
}

#' Variance attribution by squared standardised regression coefficients (design C)
#'
#' Linear model of the response on the standard deviates of the inputs (independent by construction; the L-u* link of
#' the joint draws is therefore attributed to u*), after removing the row means (within-row variance only). Share of
#' input j = b_j^2 var(z_j) / var(y); the shares add up to about R^2.
#'
#' @param y Response (e.g. TA fraction [pp] or ln x80) per run.
#' @param z Matrix of standard deviates (runs x inputs).
#' @param row Row identifier per run (row means removed).
#' @return Tibble: `input`, `share` [0-1], `r2` [0-1], `n` [-]; `NA` shares when fewer runs than inputs + 2.
mc_attribution <- function(y, z, row) {
  ok <- is.finite(y) & stats::complete.cases(z)
  y <- y[ok]
  z <- z[ok, , drop = FALSE]
  row <- row[ok]
  keep <- apply(z, 2, function(v) stats::var(v) > 0)
  inp <- colnames(z)
  if (length(y) < ncol(z) + 2 || stats::var(y) == 0) {
    return(dplyr::tibble(input = inp, share = NA_real_, r2 = NA_real_, n = length(y)))
  }
  yc <- y - stats::ave(y, row)
  zc <- apply(z[, keep, drop = FALSE], 2, function(v) v - stats::ave(v, row))
  if (stats::var(yc) == 0) return(dplyr::tibble(input = inp, share = 0, r2 = NA_real_, n = length(y)))
  fit <- stats::lm.fit(cbind(1, zc), yc)
  b <- fit$coefficients[-1]
  b[is.na(b)] <- 0
  share <- rep(0, length(inp))
  share[keep] <- b^2 * apply(zc, 2, stats::var) / stats::var(yc)
  r2 <- 1 - sum(fit$residuals^2) / sum((yc - mean(yc))^2)
  dplyr::tibble(input = inp, share = share, r2 = r2, n = length(y))
}

# --- chunks ------------------------------------------------------------------------------------------------------------

#' Chunk list of the Monte Carlo runs
#'
#' One chunk = site x design x part; a part is a block of consecutive sample rows whose estimated duration (all runs,
#' before removing invalid and duplicated runs) stays below `target_min` with `cpus` workers.
#'
#' @param site_id Station code; @param n_rows_a,n_rows_b Rows of designs A and B [-].
#' @param runs_a,runs_b Runs per row [-] (A: 17; B: draws x h options).
#' @param sec_per_run FFP cost [s per run, one core]; @param cpus Workers per chunk [-]; @param target_min Target
#'   chunk duration [min].
#' @return Tibble: `site_id`, `design`, `part`, `first_row`, `last_row` [position in the design's row list],
#'   `n_rows`, `n_runs`, `est_min` [min wall time], `cpu_hours` [h].
mc_chunk_list <- function(site_id, n_rows_a, n_rows_b, runs_a = 17L, runs_b = 100L, sec_per_run = 3, cpus = 2,
                          target_min = 120) {
  one <- function(design, n_rows, runs) {
    if (n_rows == 0) return(NULL)
    per <- max(1L, floor(target_min * 60 * cpus / (runs * sec_per_run)))
    first <- seq(1L, n_rows, by = per)
    last <- pmin(first + per - 1L, n_rows)
    n <- (last - first + 1L) * runs
    dplyr::tibble(site_id = site_id, design = design, part = sprintf("%s%03d", design, seq_along(first)),
                  first_row = as.integer(first), last_row = as.integer(last), n_rows = as.integer(last - first + 1L),
                  n_runs = as.integer(n), est_min = round(n * sec_per_run / cpus / 60, 1),
                  cpu_hours = round(n * sec_per_run / 3600, 2))
  }
  dplyr::bind_rows(one("A", n_rows_a, runs_a), one("B", n_rows_b, runs_b))
}

#' Inputs hash of a Monte Carlo chunk
#'
#' @param runs Runs of the chunk ([mc_run_inputs()]); @param ffp_md5 MD5 of the FFP code files; @param r_levels
#'   Source-area levels; @param grid [ffp_grid()]; @param mask Target-area mask.
#' @return Character (digest).
mc_chunk_hash <- function(runs, ffp_md5, r_levels, grid, mask) {
  digest::digest(list(runs, ffp_md5, r_levels, grid$domain, grid$dx, sum(mask), "mc_v1"))
}

# FFP field metrics of the runs `ii` of `sub` (inputs zm_m, h, umean, L, sigmav, ustar, wind_dir); `env` =
# load_ffp_code() of the worker. Top-level so that workers do not receive a closure with the caller's frame.
mc_ffp_block <- function(ii, sub, env, grid, mask, r_levels) {
  cols <- c("ffp_ok", "f_target_area", "f_domain", "n_nonfinite", "n_retry", "xmax_80")
  dplyr::bind_rows(lapply(ii, function(j) {
    m <- ffp_halfhour_run(env, grid, mask, r_levels, sub$zm_m[j], sub$h[j], sub$umean[j], sub$L[j], sub$sigmav[j],
                          sub$ustar[j], sub$wind_dir[j])$metrics
    for (k in setdiff(cols, names(m))) m[[k]] <- NA
    m[cols]
  }))
}

#' FFP of the runs of a chunk (unique valid inputs only)
#'
#' Runs that are invalid under their own set (`valid_run`: `etc`, or `etc` without e_rsl for the `rsl_off` runs) are
#' not computed (FFP `NA`); runs with identical inputs are computed once
#' (e.g. design B with `h_empirical` at night, equal to `h_pct30`). Distances x_peak, x_50 ... x_90 from
#' [ffp_ci_paper()] as in 08b; field metrics from [ffp_halfhour_run()]. Work spread over `workers` processes in fixed
#' blocks (result independent of `workers`).
#'
#' @param runs Output of [mc_run_inputs()].
#' @param ffp_dir Folder of Kljun's FFP code; @param grid,mask,r_levels As in 08b.
#' @param workers Processes [-]; @param block_runs Runs per work block [-].
#' @return `runs` with `ffp_ok`, `f_target_area`, `f_domain` [-], `n_nonfinite`, `n_retry`, `xmax_80` [m] and
#'   `x_peak`, `x_50`, `x_70`, `x_80`, `x_90` [m]; attribute `n_ffp` (FFP calls made).
mc_run_ffp <- function(runs, ffp_dir, grid, mask, r_levels, workers = 1L, block_runs = 50L) {
  key <- mc_run_key(runs)
  vr <- which(runs$valid_run) # first valid run of each distinct input set (an invalid duplicate never hides it)
  todo <- vr[!duplicated(key[vr])]
  blocks <- split(todo, ceiling(seq_along(todo) / block_runs))
  res <- list()
  if (length(todo) && workers > 1 && length(blocks) > 1) {
    cl <- make_worker_cluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (.Platform$OS.type == "windows") {
      rfiles <- list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)
      parallel::clusterCall(cl, function(fs) { for (f in fs) source(f); NULL }, rfiles)
    }
    parallel::clusterCall(cl, function(d) { assign(".ffp_env", load_ffp_code(d), envir = globalenv()); NULL }, ffp_dir)
    sub <- runs[todo, c("zm_m", "h", "umean", "L", "sigmav", "ustar", "wind_dir")]
    bl <- split(seq_along(todo), ceiling(seq_along(todo) / block_runs))
    res <- parallel::parLapply(cl, bl, function(ii, sub, grid, mask, r_levels) {
      mc_ffp_block(ii, sub, get(".ffp_env", envir = globalenv()), grid, mask, r_levels)
    }, sub = sub, grid = grid, mask = mask, r_levels = r_levels)
  } else if (length(todo)) {
    env <- load_ffp_code(ffp_dir)
    sub <- runs[todo, c("zm_m", "h", "umean", "L", "sigmav", "ustar", "wind_dir")]
    res <- lapply(split(seq_along(todo), ceiling(seq_along(todo) / block_runs)),
                  function(ii) mc_ffp_block(ii, sub, env, grid, mask, r_levels))
  }
  met <- dplyr::bind_rows(res)
  pos <- match(key, key[todo]) # every run gets the result of its first identical valid run
  pos[!runs$valid_run] <- NA
  out <- runs
  for (k in c("ffp_ok", "f_target_area", "f_domain", "n_nonfinite", "n_retry", "xmax_80")) {
    out[[k]] <- if (nrow(met)) met[[k]][pos] else NA
  }
  out$ffp_ok <- out$ffp_ok %in% TRUE
  ci <- ffp_ci_paper(runs$zm_m, runs$h, runs$umean, runs$ustar, r_levels)
  ci[!out$ffp_ok, ] <- NA
  out <- dplyr::bind_cols(out, ci)
  attr(out, "n_ffp") <- length(todo)
  out
}

# --- tables (16c) -------------------------------------------------------------------------------------------------------

#' Weights of the design-B rows
#'
#' The stratum's share of the eligible rows divided by the B rows sampled in it (sum 1 over the B rows of a site);
#' `NA` for rows not in B.
#'
#' @param s Sample of 16a (`stratum`, `w_stratum`, `in_b`).
#' @return Numeric per row of `s`.
mc_b_weights <- function(s) {
  nb <- tapply(s$in_b, s$stratum, sum)
  ifelse(s$in_b, s$w_stratum / as.numeric(nb[s$stratum]), NA_real_)
}

#' Runs of a site joined with the sample
#'
#' @param draws Output of 16a (`<site>_draws.rds`).
#' @param runs Bound `_runs.csv.gz` of the site's chunks (16b).
#' @return `runs` with the sample columns `timestamp_start_utc`, `stratum`, `w_row`, `in_b`, `valid_qc`,
#'   `valid_qc_used`, `day` and the spec-15 values (`ta__<case>`, `x80__<case>`); `ta_pct` = TA fraction [%].
mc_join_runs <- function(draws, runs) {
  s <- draws$sample
  keep <- c("row_id", "timestamp_start_utc", "stratum", "w_stratum", "w_row", "in_b", "valid_qc", "valid_qc_used",
            grep("^(ta|x80|h)__", names(s), value = TRUE))
  x <- dplyr::left_join(runs, s[keep], by = "row_id")
  x$ta_pct <- 100 * x$f_target_area
  x
}

#' Design A table: change of each run against the row's base run
#'
#' @param x Output of [mc_join_runs()] (design A rows).
#' @return Tibble, one row per sampled row x run: `row_id`, `timestamp_start_utc`, `stratum`, `valid_qc`,
#'   `valid_qc_used`, `input` (`base` for the base run), `side` (-1, 0, +1), `validity` (`etc`, `rsl_off`),
#'   `valid_etc`, `fail`, `rsl_fail`, `valid_run`, `ffp_ok`, `ta_pct`
#'   [%], `d_ta_pp` [pp, run minus base], `x80_m` [m], `d_lnx80` [-, ln(x80 run / x80 base)], `base_diff_pp` [pp,
#'   base run minus the spec-15 PBLHmod value of the row], `base_diff_x80_m` [m].
mc_oat_table <- function(x) {
  a <- x[x$design == "A", ]
  b <- a[a$run == "base", ]
  m <- match(a$row_id, b$row_id)
  dplyr::tibble(row_id = a$row_id, timestamp_start_utc = a$timestamp_start_utc, stratum = a$stratum,
                valid_qc = a$valid_qc, valid_qc_used = a$valid_qc_used, input = ifelse(is.na(a$input), "base", a$input),
                side = a$side, validity = a$validity, valid_etc = a$valid_etc, fail = a$fail, rsl_fail = a$rsl_fail,
                valid_run = a$valid_run, ffp_ok = a$ffp_ok, ta_pct = a$ta_pct,
                d_ta_pp = a$ta_pct - b$ta_pct[m], x80_m = a$x_80, d_lnx80 = log(a$x_80 / b$x_80[m]),
                base_diff_pp = b$ta_pct[m] - 100 * a$ta__pblhmod, base_diff_x80_m = b$x_80[m] - a$x80__pblhmod)
}

#' Design A: 1-sigma effect per row and input (half-range)
#'
#' @param oat Output of [mc_oat_table()].
#' @return Tibble: `row_id`, `stratum`, `input`, `validity` (`etc`, `rsl_off`), `eff_ta_pp` [pp], `eff_lnx80` [-]
#'   (`NA` when a side is invalid under the run's set), `rsl_flip` [0-1, share of the input's two runs at -1 / +1
#'   sigma that fail the e_rsl test of the `etc` set: 0, 0.5 or 1].
mc_oat_effects <- function(oat) {
  o <- oat[oat$input != "base", ]
  if (!"validity" %in% names(o)) o$validity <- "etc"
  if (!"rsl_fail" %in% names(o)) o$rsl_fail <- FALSE
  mi <- o[o$side == -1, ]
  pl <- o[o$side == 1, ]
  k <- match(paste(mi$row_id, mi$input, mi$validity), paste(pl$row_id, pl$input, pl$validity))
  dplyr::tibble(row_id = mi$row_id, stratum = mi$stratum, input = mi$input, validity = mi$validity,
                eff_ta_pp = mc_half_range(mi$d_ta_pp, pl$d_ta_pp[k]), eff_lnx80 = mc_half_range(mi$d_lnx80, pl$d_lnx80[k]),
                rsl_flip = (as.numeric(mi$rsl_fail %in% TRUE) + as.numeric(pl$rsl_fail[k] %in% TRUE)) / 2)
}

#' PBLH-source effect of the sampled rows (spec 15, from the 08b chunks)
#'
#' @param s Sample of 16a.
#' @return Tibble: `row_id`, `stratum`, `w_row`, `w_b` ([mc_b_weights()]), and per pair (`era5_vs_pblhmod`, `pblhmod_clim_vs_pblhmod`)
#'   `src_ta_pp__<pair>` [pp, first minus second] and `src_lnx80__<pair>` [-].
mc_source_effects <- function(s) {
  out <- dplyr::tibble(row_id = s$row_id, stratum = s$stratum, w_row = s$w_row, w_b = mc_b_weights(s))
  for (cs in c("era5", "pblhmod_clim")) {
    p <- paste0(cs, "_vs_pblhmod")
    out[[paste0("src_ta_pp__", p)]] <- 100 * (s[[paste0("ta__", cs)]] - s$ta__pblhmod)
    out[[paste0("src_lnx80__", p)]] <- log(s[[paste0("x80__", cs)]] / s$x80__pblhmod)
  }
  out
}

#' Design B: row-level quantiles of the joint draws
#'
#' @param x Output of [mc_join_runs()].
#' @return Tibble per row x h option: `row_id`, `stratum`, `h_option`, `n_draws`, `n_valid` [-], `share_invalid`
#'   [0-1], `ta_p05` ... `ta_p95` [%], `lnx80_p05` ... `lnx80_p95` [-, ln m], `width_ta_pp` [pp, (p84 - p16) / 2],
#'   `width_lnx80` [-], `valid_qc`, `valid_qc_used`.
mc_mc_table <- function(x) {
  b <- x[x$design == "B", ]
  if (!nrow(b)) return(dplyr::tibble())
  p <- c(0.05, 0.16, 0.5, 0.84, 0.95)
  g <- split(seq_len(nrow(b)), paste(b$row_id, b$h_option))
  dplyr::bind_rows(lapply(g, function(ii) {
    d <- b[ii, ]
    ok <- d$valid_etc & d$ffp_ok
    qt <- if (any(ok)) stats::quantile(d$ta_pct[ok], p, names = FALSE) else rep(NA_real_, 5)
    ql <- if (any(ok)) stats::quantile(log(d$x_80[ok]), p, names = FALSE) else rep(NA_real_, 5)
    o <- dplyr::tibble(row_id = d$row_id[1], stratum = d$stratum[1], h_option = d$h_option[1], n_draws = nrow(d),
                       n_valid = sum(ok), share_invalid = mean(!d$valid_etc))
    for (j in seq_along(p)) o[[sprintf("ta_p%02d", round(100 * p[j]))]] <- qt[j]
    for (j in seq_along(p)) o[[sprintf("lnx80_p%02d", round(100 * p[j]))]] <- ql[j]
    o$width_ta_pp <- (qt[4] - qt[2]) / 2
    o$width_lnx80 <- (ql[4] - ql[2]) / 2
    o$valid_qc <- d$valid_qc[1]
    o$valid_qc_used <- d$valid_qc_used[1]
    o
  }))
}

#' Design C: attribution per stratum and h option
#'
#' @param x Output of [mc_join_runs()].
#' @param inputs Inputs [character] (`z_<input>` columns of the B runs).
#' @return Tibble: `stratum` (six strata and `all`), `h_option`, `response` (`ta_pp`, `lnx80`), `input`, `share`,
#'   `r2`, `n` (see [mc_attribution()]).
mc_attribution_table <- function(x, inputs = mc_input_spec()$input) {
  b <- x[x$design == "B" & x$valid_etc & x$ffp_ok, ]
  if (!nrow(b)) return(dplyr::tibble())
  z <- as.matrix(b[paste0("z_", inputs)])
  colnames(z) <- inputs
  out <- list()
  for (ho in unique(b$h_option)) for (st in c(mc_strata_levels(), "all")) {
    i <- b$h_option == ho & (st == "all" | b$stratum == st)
    if (!any(i)) next
    for (resp in c("ta_pp", "lnx80")) {
      y <- if (resp == "ta_pp") b$ta_pct[i] else log(b$x_80[i])
      out[[length(out) + 1]] <- dplyr::bind_cols(dplyr::tibble(stratum = st, h_option = ho, response = resp),
                                                 mc_attribution(y, z[i, , drop = FALSE], b$row_id[i]))
    }
  }
  dplyr::bind_rows(out)
}

#' Summary statistics of a per-row quantity by stratum (and all strata, row-weighted)
#'
#' @param v Values; @param stratum Stratum per value; @param w Row weights (for `all`).
#' @return Tibble: `stratum`, `n`, `median`, `p25`, `p75`.
mc_stratum_stats <- function(v, stratum, w) {
  one <- function(st) {
    i <- !is.na(v) & (st == "all" | stratum == st)
    q <- if (st == "all") mc_wquantile(v[i], w[i], c(0.5, 0.25, 0.75)) else
      if (any(i)) stats::quantile(v[i], c(0.5, 0.25, 0.75), names = FALSE) else rep(NA_real_, 3)
    dplyr::tibble(stratum = st, n = sum(i), median = q[1], p25 = q[2], p75 = q[3])
  }
  dplyr::bind_rows(lapply(c(mc_strata_levels(), "all"), one))
}

#' Section 4 summaries of one site
#'
#' One-sigma effects of the inputs (A), the PBLH-source effect (|first minus second| of the spec-15 values), the joint
#' 1-sigma width (B, per h option), and the ratio |source effect| / joint width with the shares of rows above 1 and 2
#' joint widths. By stratum, and `all` (row-weighted, stratum weights of 16a).
#'
#' @param site_id Station code; @param eff Output of [mc_oat_effects()]; @param src Output of
#'   [mc_source_effects()] (with `w_b`); @param mc Output of [mc_mc_table()].
#' @return Tibble: `kind` (`oat_1sigma`, `source_effect`, `joint_width`, `ratio`), `scope` (`site`), `site_id`,
#'   `stratum`, `quantity` (`ta_pp`, `lnx80`), `item` (input, `input|rsl_off` for the runs without the e_rsl test,
#'   pair, h option, or `pair|h option`), `n`, `median`, `p25`, `p75`, `share_gt_1s`, `share_gt_2s` [0-1, ratio rows
#'   only], `share_rsl_flip` [0-1, oat_1sigma rows only: share of the input's -1 / +1 sigma runs that fail e_rsl;
#'   `all` row-weighted].
mc_site_summary <- function(site_id, eff, src, mc) {
  w_of <- function(ids) src$w_row[match(ids, src$row_id)]
  out <- list()
  add <- function(kind, quantity, item, v, stratum, w, extra = NULL) {
    s <- mc_stratum_stats(v, stratum, w)
    s <- dplyr::bind_cols(dplyr::tibble(kind = kind, scope = "site", site_id = site_id, quantity = quantity, item = item), s)
    if (!is.null(extra)) s <- dplyr::bind_cols(s, extra)
    out[[length(out) + 1]] <<- s
  }
  if (!"validity" %in% names(eff)) eff$validity <- "etc"
  if (!"rsl_flip" %in% names(eff)) eff$rsl_flip <- 0
  lev <- c(mc_strata_levels(), "all")
  for (vs in unique(eff$validity)) for (inp in unique(eff$input[eff$validity == vs])) {
    e <- eff[eff$input == inp & eff$validity == vs, ]
    w <- w_of(e$row_id)
    item <- if (vs == "etc") inp else paste0(inp, "|", vs)
    fl <- vapply(lev, function(st) {
      i <- !is.na(e$rsl_flip) & (st == "all" | e$stratum == st)
      if (!any(i)) return(NA_real_)
      if (st == "all") sum(w[i] * e$rsl_flip[i]) / sum(w[i]) else mean(e$rsl_flip[i])
    }, numeric(1))
    ex <- dplyr::tibble(share_rsl_flip = unname(fl))
    add("oat_1sigma", "ta_pp", item, e$eff_ta_pp, e$stratum, w, extra = ex)
    add("oat_1sigma", "lnx80", item, e$eff_lnx80, e$stratum, w, extra = ex)
  }
  pairs <- c("era5_vs_pblhmod", "pblhmod_clim_vs_pblhmod")
  for (p in pairs) {
    add("source_effect", "ta_pp", p, abs(src[[paste0("src_ta_pp__", p)]]), src$stratum, src$w_row)
    add("source_effect", "lnx80", p, abs(src[[paste0("src_lnx80__", p)]]), src$stratum, src$w_row)
  }
  if (nrow(mc)) {
    wst <- src$w_b[match(mc$row_id, src$row_id)] # weights of the B rows ([mc_b_weights()])
    for (ho in unique(mc$h_option)) {
      m <- mc[mc$h_option == ho, ]
      w <- wst[mc$h_option == ho]
      add("joint_width", "ta_pp", ho, m$width_ta_pp, m$stratum, w)
      add("joint_width", "lnx80", ho, m$width_lnx80, m$stratum, w)
      for (p in pairs) for (q in c("ta_pp", "lnx80")) {
        sv <- abs(src[[paste0(if (q == "ta_pp") "src_ta_pp__" else "src_lnx80__", p)]])[match(m$row_id, src$row_id)]
        wd <- if (q == "ta_pp") m$width_ta_pp else m$width_lnx80
        r <- sv / wd
        r[!is.finite(r)] <- NA
        sh <- function(st, k) {
          i <- !is.na(sv) & !is.na(wd) & (st == "all" | m$stratum == st)
          if (!any(i)) return(NA_real_)
          if (st == "all") sum(w[i] * (sv[i] > k * wd[i])) / sum(w[i]) else mean(sv[i] > k * wd[i])
        }
        lev <- c(mc_strata_levels(), "all")
        add("ratio", q, paste0(p, "|", ho), r, m$stratum, w,
            extra = dplyr::tibble(share_gt_1s = vapply(lev, sh, numeric(1), k = 1),
                                  share_gt_2s = vapply(lev, sh, numeric(1), k = 2)))
      }
    }
  }
  x <- dplyr::bind_rows(out)
  for (k in c("share_gt_1s", "share_gt_2s", "share_rsl_flip")) if (!k %in% names(x)) x[[k]] <- NA_real_
  x
}

#' Network rows of the section 4 summaries
#'
#' Per kind x stratum x quantity x item: median (and p25, p75) over the sites of the site medians; with `all`, the
#' site values are already stratum-weighted, so the network `all` is the median over sites of the stratum-weighted
#' site medians.
#'
#' @param s Bound outputs of [mc_site_summary()].
#' @return Tibble with the columns of [mc_site_summary()], `scope = "network"`, `site_id = "network"`, `n` = sites.
mc_network_summary <- function(s) {
  g <- split(seq_len(nrow(s)), paste(s$kind, s$stratum, s$quantity, s$item, sep = "\r"))
  dplyr::bind_rows(lapply(g, function(ii) {
    d <- s[ii, ]
    v <- d$median[!is.na(d$median)]
    q <- if (length(v)) stats::quantile(v, c(0.5, 0.25, 0.75), names = FALSE) else rep(NA_real_, 3)
    dplyr::tibble(kind = d$kind[1], scope = "network", site_id = "network", quantity = d$quantity[1], item = d$item[1],
                  stratum = d$stratum[1], n = length(v), median = q[1], p25 = q[2], p75 = q[3],
                  share_gt_1s = stats::median(d$share_gt_1s, na.rm = TRUE), share_gt_2s = stats::median(d$share_gt_2s, na.rm = TRUE),
                  share_rsl_flip = if (all(is.na(d$share_rsl_flip))) NA_real_ else stats::median(d$share_rsl_flip, na.rm = TRUE))
  }))
}

#' ICOS criteria under uncertainty (design B)
#'
#' Per draw: the stratum-weighted share of the B rows with TA >= 90 % (Rebmann et al. 2018: pass if >= 50 %) and, among
#' the `valid_qc` rows, the share with TA >= 70 % (ICOS labelling rule, ETC: pass if >= 70 %); the probability of a pass
#' is the share of draws that pass. Case `pblhmod` from the draws; `era5` and `pblhmod_clim` approximated by adding the
#' row's spec-15 difference (case minus PBLHmod) to every draw (flagged in `method`). Runs invalid under perturbation
#' leave the row out of that draw.
#'
#' @param x Output of [mc_join_runs()]; @param src Output of [mc_source_effects()]; @param site_id Station code.
#' @return Tibble: `site_id`, `h_option`, `case`, `criterion` (`rebmann_ta90`, `icos_lab_ta70`), `threshold_share`,
#'   `share_base` [0-1, spec-15 values of the B rows], `pass_base`, `share_p05`, `share_p50`, `share_p95` [0-1],
#'   `prob_pass` [0-1], `uncertain` [logical, 0.1 < prob_pass < 0.9], `n_rows`, `n_draws`, `method`.
mc_icos_table <- function(x, src, site_id) {
  b <- x[x$design == "B", ]
  if (!nrow(b)) return(dplyr::tibble())
  smp_b <- unique(b[c("row_id", "stratum", "valid_qc")])
  smp_b$w <- src$w_b[match(smp_b$row_id, src$row_id)]
  base <- unique(b[c("row_id", "valid_qc", "ta__pblhmod", "ta__pblhmod_clim", "ta__era5")])
  out <- list()
  for (ho in unique(b$h_option)) for (cs in c("pblhmod", "pblhmod_clim", "era5")) {
    d <- b[b$h_option == ho, ]
    shift <- if (cs == "pblhmod") 0 else 100 * (d[[paste0("ta__", cs)]] - d$ta__pblhmod)
    if (all(is.na(shift))) next
    ta <- ifelse(d$valid_etc & d$ffp_ok, d$ta_pct + shift, NA_real_)
    w <- smp_b$w[match(d$row_id, smp_b$row_id)]
    base_ta <- 100 * base[[paste0("ta__", cs)]]
    base_w <- smp_b$w[match(base$row_id, smp_b$row_id)]
    base_qc <- base$valid_qc
    for (crit in c("rebmann_ta90", "icos_lab_ta70")) {
      thr <- if (crit == "rebmann_ta90") c(90, 0.5) else c(70, 0.7)
      sel <- if (crit == "rebmann_ta90") rep(TRUE, nrow(d)) else d$valid_qc %in% TRUE
      per_draw <- vapply(split(seq_len(nrow(d))[sel], d$draw[sel]), function(ii) {
        ok <- !is.na(ta[ii])
        if (!any(ok)) return(NA_real_)
        sum(w[ii][ok] * (ta[ii][ok] >= thr[1])) / sum(w[ii][ok])
      }, numeric(1))
      bsel <- if (crit == "rebmann_ta90") !is.na(base_ta) else !is.na(base_ta) & base_qc %in% TRUE
      sb <- if (any(bsel)) sum(base_w[bsel] * (base_ta[bsel] >= thr[1])) / sum(base_w[bsel]) else NA_real_
      pd <- per_draw[!is.na(per_draw)]
      pp <- if (length(pd)) mean(pd >= thr[2]) else NA_real_
      q <- if (length(pd)) stats::quantile(pd, c(0.05, 0.5, 0.95), names = FALSE) else rep(NA_real_, 3)
      out[[length(out) + 1]] <- dplyr::tibble(
        site_id = site_id, h_option = ho, case = cs, criterion = crit, threshold_share = thr[2],
        share_base = sb, pass_base = sb >= thr[2], share_p05 = q[1], share_p50 = q[2], share_p95 = q[3], prob_pass = pp,
        uncertain = !is.na(pp) & pp > 0.1 & pp < 0.9, n_rows = length(unique(d$row_id[sel])), n_draws = length(pd),
        method = if (cs == "pblhmod") "draws" else "draws of PBLHmod + row difference case - PBLHmod (spec 15)")
    }
  }
  dplyr::bind_rows(out)
}
