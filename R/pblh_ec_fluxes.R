# Spec 18: where the PBLH source matters for EC fluxes (surface layer, u* filter, storage, sigma_v, SCF).
# Analysis only: no PBLHmod or FFP run. Rows = every half-hour of the station's L2 record from its labelling date
# (first valid pblhmod_default value, as network_record_start()), exclusions removed. h sources: standard files
# era5_blh_center (instantaneous, interpolated to the half-hour mid as the footprint cases, era5_at_halfhour_mid()),
# pblhmod_default, pblhmod_clim, and the site's primary measured source (exact values; censored ones as their own
# class, never imputed).
#
# Circularity (spec 18, kept in every result): PBLHmod's night-time h is computed from u* and L (spec 02, B1), so with
# PBLHmod "zm/h > 0.1" is close to "low u*". ERA5 is the main source of sections 1-3.
#
# Panofsky et al. (1977, Boundary-Layer Meteorol. 11, 355-361): sigma_v / u* = (12 - 0.5 h / L)^(1/3) in the
# unstable surface layer. Not in docs/references/*.bib (checked 2026-11-27): to be added by G after checking.

pec_sources <- function() c("era5_blh_center", "pblhmod_default", "pblhmod_clim")

#' Missing requirements of a spec-18 part
#'
#' @param need Named list of conditions (one per variable the part needs; `TRUE` = available).
#' @return Names of the unmet conditions [character]; empty = the part can run.
pec_missing <- function(need) names(need)[!vapply(need, isTRUE, logical(1))]

#' Any finite value in a column of the rows
#'
#' @param r Rows; @param k Column name.
#' @return Logical: the column exists and has a finite value.
pec_has <- function(r, k) k %in% names(r) && any(is.finite(r[[k]]))

#' Night id of each half-hour (sunset to sunrise across midnight)
#'
#' Night rows grouped by the local solar day of (half-hour mid - 12 h): a night from sunset to sunrise falls in one
#' group across midnight; during polar night every solar day (noon to noon) is a group; during polar day there are no
#' night rows.
#'
#' @param t_start Half-hour starts [POSIXct UTC]; @param night Night row [logical, `SW_IN_POT == 0`].
#' @param lon Station longitude [deg E] (local solar time = UTC + lon / 15 h).
#' @return Integer night id (days since 1970-01-01 of the evening's solar date), `NA` on day rows.
night_groups <- function(t_start, night, lon) {
  solar <- as.numeric(t_start) + 900 + lon / 15 * 3600 - 12 * 3600
  id <- as.integer(floor(solar / 86400))
  id[!(night %in% TRUE)] <- NA_integer_
  id
}

#' u* bin of spec 18 (u* / u*_thr)
#'
#' @param ratio u* / u*_thr [-].
#' @return Character: `1-1.5`, `1.5-2`, `>2` (left-closed), `NA` below 1 or missing.
ustar_ratio_bin <- function(ratio) {
  as.character(cut(ratio, c(1, 1.5, 2, Inf), labels = c("1-1.5", "1.5-2", ">2"), right = FALSE))
}

#' Panofsky et al. (1977) sigma_v / u* in the unstable surface layer
#'
#' @param h Boundary-layer height [m]; @param L Obukhov length [m, negative unstable].
#' @return Predicted sigma_v / u* [-] = (12 - 0.5 h / L)^(1/3) (Panofsky et al. 1977); `NA` where L >= 0 or h <= 0.
panofsky_sv_ustar <- function(h, L) {
  ok <- is.finite(h) & is.finite(L) & L < 0 & h > 0
  out <- rep(NA_real_, length(h))
  out[ok] <- (12 - 0.5 * h[ok] / L[ok])^(1 / 3)
  out
}

#' Agreement of the two u*-filter definitions (QC flag vs threshold)
#'
#' @param below_thr USTAR below the VUT threshold [logical]; @param qc NEE_VUT_REF_QC [-, 0 measured].
#' @return One-row tibble: `n` (rows with both), `n_both`, `n_thr_only`, `n_qc_only`, `n_neither`, `agreement` [-].
ustar_filter_agreement <- function(below_thr, qc) {
  ok <- !is.na(below_thr) & !is.na(qc)
  a <- below_thr[ok]
  b <- qc[ok] > 0
  dplyr::tibble(n = sum(ok), n_both = sum(a & b), n_thr_only = sum(a & !b), n_qc_only = sum(!a & b),
                n_neither = sum(!a & !b), agreement = if (sum(ok)) mean(a == b) else NA_real_)
}

#' Storage profile and SC check of a site
#'
#' ICOS L2 NEE = FC + SC; `SC_FLAG` 0 = profile (checked 2026-11-27 at BE-Bra, FI-Hyy, CH-Dav: rows with flag 1 have
#' SC equal to `SC_1P`, the single-point estimate, and FR-Gri, without profile levels, has flag 1 only).
#'
#' @param header FLUXES column names; @param d Rows with `FC`, `NEE`, `SC`, `SC_FLAG` (any may be absent).
#' @param tol Tolerance of SC = NEE - FC [umol m-2 s-1] (values are rounded to 2 decimals in the files).
#' @return One-row tibble: `n_profile_levels` (`CO2_DRY_<n>` columns), `has_sc`, `n_sc_profile` (SC with flag 0),
#'   `n_sc_flag1`, `n_sc_flag_other`, `sc_check_n`, `sc_check_share` (share with |NEE - FC - SC| <= tol),
#'   `storage` [logical: profile levels and SC flag-0 rows].
storage_profile_info <- function(header, d, tol = 0.011) {
  lev <- length(grep("^CO2_DRY_[0-9]+$", header))
  has_sc <- all(c("SC", "FC", "NEE") %in% names(d))
  fl <- if ("SC_FLAG" %in% names(d)) d$SC_FLAG else rep(NA_real_, nrow(d))
  sc_ok <- if (has_sc) is.finite(d$SC) else rep(FALSE, nrow(d))
  ck <- if (has_sc) is.finite(d$SC) & is.finite(d$FC) & is.finite(d$NEE) else rep(FALSE, nrow(d))
  dplyr::tibble(n_profile_levels = lev, has_sc = has_sc, n_sc_profile = sum(sc_ok & fl %in% 0),
                n_sc_flag1 = sum(sc_ok & fl %in% 1), n_sc_flag_other = sum(sc_ok & !(fl %in% c(0, 1))),
                sc_check_n = sum(ck),
                sc_check_share = if (sum(ck)) mean(abs(d$NEE[ck] - d$FC[ck] - d$SC[ck]) <= tol) else NA_real_,
                storage = lev > 0 && sum(sc_ok & fl %in% 0) > 0)
}

#' Measured standard values on the half-hour grid
#'
#' @param meas Measured standard table; @param t_start Half-hour starts [POSIXct UTC].
#' @return Tibble per half-hour: `h_meas` [m, exact values only], `meas_kind` (`exact`, `below`, `lower_bound`, `NA`).
#'   A measured interval covering the half-hour (start <= t < end, e.g. hourly SIRTA rows) gives its value.
measured_on_halfhours <- function(meas, t_start) {
  mv <- measured_values(meas)
  st <- as.numeric(meas$timestamp_start_utc)
  en <- as.numeric(meas$timestamp_end_utc)
  o <- order(st)
  st <- st[o]; en <- en[o]; mv <- mv[o, ]
  t <- as.numeric(t_start)
  i <- findInterval(t, st)
  i[i == 0] <- NA_integer_
  ok <- !is.na(i) & t < en[pmax(i, 1)]
  i[!ok] <- NA_integer_
  kind <- mv$kind[i]
  dplyr::tibble(h_meas = ifelse(kind %in% "exact", mv$pblh_m[i], NA_real_), meas_kind = kind)
}

#' Rows of spec 18 for one site
#'
#' Reads the FLUXES and FLUXNET HH files (only the columns spec 18 needs that are present), the h standard files and
#' the ONEFlux VUT u* thresholds; removes excluded periods (counted) and the rows before the labelling date.
#'
#' @param site_id Station; @param cfg Configuration ([get_cfg()]); @param sites,sites_tv,exclusions Metadata
#'   ([read_sites()], [read_sites_timevarying()], [read_exclusions()]); @param ref_years Reference years
#'   ([read_reference_years()], may be `NULL`); @param ms Measured pairings ([read_measured_sites()], may be `NULL`).
#' @return List: `rows` (tibble per half-hour: `timestamp_start_utc`, `year_local`, `ref_year` [logical], `season`,
#'   `night` [logical, SW_IN_POT = 0], `night_id`, `zm_m` [m, EC height - 2/3 hc], `ec_height_m`, `zl` (zm/L class),
#'   `L` [m], `ustar` [m s-1, L2], `ustar_fn` [m s-1, FLUXNET], `ustar_thr` [m s-1], `below_thr` [logical],
#'   `ustar_ratio` [-], `nee_qc`, `nee_vut`, `reco_nt`, `fc`, `nee_l2`, `nee_data_flag`, `sc`, `sc_flag`
#'   [umol m-2 s-1 / -], `ta_f_k` [K], `v_sigma`, `ws` [m s-1], `*_scf` [-], `h_<source>` [m], `h_meas` [m],
#'   `meas_kind`), `vars` (one-row availability table), `meas_src` (measured source id or `NA`), `meas_type`
#'   (`aerosol`, `turbulence` or `NA`).
pblh_ec_site_rows <- function(site_id, cfg, sites, sites_tv, exclusions, ref_years = NULL, ms = NULL) {
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) cli::cli_abort("Site {.val {site_id}} not in the site metadata.")
  off <- meta$utc_offset_h
  files <- find_icos_files(site_id, cfg)
  if (is.null(files$fluxes) || is.null(files$fluxnet)) cli::cli_abort("{site_id}: FLUXES or FLUXNET HH file missing.")
  tc <- c("TIMESTAMP_START", "TIMESTAMP_END")
  want_fx <- c("USTAR", "MO_LENGTH", "WS", "V_SIGMA", "W_SIGMA", "FC", "NEE", "SC", "SC_FLAG", "NEE_DATA_FLAG",
               "FC_SSITC_TEST", "H_SCF_STAT", "LE_SCF_STAT", "FC_SCF_STAT", "TAU_SCF_STAT")
  want_fn <- c("NEE_VUT_REF", "NEE_VUT_REF_QC", "RECO_NT_VUT_REF", "TA_F", "SW_IN_POT", "NIGHT", "USTAR")
  hfx <- icos_zip_header(files$fluxes)
  hfn <- icos_zip_header(files$fluxnet)
  fx <- read_icos_csv(files$fluxes, c(tc, intersect(want_fx, hfx)), off, label = "FLUXES")
  fn <- read_icos_csv(files$fluxnet, c(tc, intersect(want_fn, hfn)), off, label = "FLUXNET")
  col <- function(d, k) if (k %in% names(d)) d[[k]] else rep(NA_real_, nrow(d))
  m <- match(as.numeric(fx$timestamp_start_utc), as.numeric(fn$timestamp_start_utc))
  t <- fx$timestamp_start_utc
  hs <- get_site_heights(site_id, t, sites, sites_tv)
  r <- dplyr::tibble(
    timestamp_start_utc = t,
    year_local = as.integer(substr(fx$TIMESTAMP_START, 1, 4)),
    ec_height_m = hs$zm_m, canopy_height_m = hs$canopy_height_m,
    zm_m = hs$zm_m - 2 / 3 * hs$canopy_height_m,
    L = col(fx, "MO_LENGTH"), ustar = col(fx, "USTAR"), ws = col(fx, "WS"), v_sigma = col(fx, "V_SIGMA"),
    fc = col(fx, "FC"), nee_l2 = col(fx, "NEE"), sc = col(fx, "SC"), sc_flag = col(fx, "SC_FLAG"),
    nee_data_flag = col(fx, "NEE_DATA_FLAG"), fc_ssitc = col(fx, "FC_SSITC_TEST"),
    h_scf = col(fx, "H_SCF_STAT"), le_scf = col(fx, "LE_SCF_STAT"), fc_scf = col(fx, "FC_SCF_STAT"),
    tau_scf = col(fx, "TAU_SCF_STAT"),
    nee_vut = col(fn, "NEE_VUT_REF")[m], nee_qc = col(fn, "NEE_VUT_REF_QC")[m], reco_nt = col(fn, "RECO_NT_VUT_REF")[m],
    ta_f_k = col(fn, "TA_F_K")[m], sw_in_pot = col(fn, "SW_IN_POT")[m], ustar_fn = col(fn, "USTAR")[m]
  )
  # year_local of added grid rows (no TIMESTAMP_START): from the UTC time + offset
  na_y <- is.na(r$year_local)
  r$year_local[na_y] <- as.integer(format(t[na_y] + off * 3600, "%Y", tz = "UTC"))
  vars <- dplyr::tibble(site_id = site_id)
  for (k in want_fx) vars[[paste0("n_", k)]] <- if (k %in% hfx) sum(!is.na(fx[[k]])) else NA_integer_
  for (k in want_fn) vars[[paste0("n_fn_", k)]] <- if (k %in% hfn) sum(!is.na(col(fn, if (k == "TA_F") "TA_F_K" else k))) else NA_integer_
  vars <- dplyr::bind_cols(vars, storage_profile_info(hfx, dplyr::tibble(FC = r$fc, NEE = r$nee_l2, SC = r$sc,
                                                                          SC_FLAG = r$sc_flag)))
  # exclusions
  ex <- flag_exclusions(r, site_id, exclusions, off)
  vars$n_excluded <- sum(ex$excl_flag)
  r <- r[!ex$excl_flag, ]
  # h sources
  sd <- file.path(resolve_path(cfg$standard_dir), site_id)
  rs <- function(src) {
    p <- standard_paths(sd, site_id, src)[["data"]]
    if (!file.exists(p)) cli::cli_abort("{site_id}: standard file {.path {p}} missing.")
    read_standard_pblh(p)
  }
  pm <- rs("pblhmod_default")
  r$h_pblhmod_default <- pm$pblh_m[match(as.numeric(r$timestamp_start_utc), as.numeric(pm$timestamp_start_utc))]
  cl <- rs("pblhmod_clim")
  r$h_pblhmod_clim <- cl$pblh_m[match(as.numeric(r$timestamp_start_utc), as.numeric(cl$timestamp_start_utc))]
  e5 <- rs("era5_blh_center")
  r$h_era5_blh_center <- era5_at_halfhour_mid(r$timestamp_start_utc, e5$timestamp_start_utc, e5$pblh_m)
  mref <- if (!is.null(ms)) measured_reference(ms, site_id) else NULL
  meas_src <- if (is.null(mref)) NA_character_ else mref$src
  meas_type <- if (is.null(mref)) NA_character_ else mref$type
  if (!is.na(meas_src)) {
    mh <- measured_on_halfhours(rs(meas_src), r$timestamp_start_utc)
    r$h_meas <- mh$h_meas
    r$meas_kind <- mh$meas_kind
  } else {
    r$h_meas <- NA_real_
    r$meas_kind <- NA_character_
  }
  # labelling date: first valid pblhmod_default value
  first <- which(is.finite(r$h_pblhmod_default))[1]
  if (is.na(first)) cli::cli_abort("{site_id}: no valid pblhmod_default value.")
  vars$n_before_labelling <- first - 1L
  r <- r[first:nrow(r), ]
  # u* thresholds (ONEFlux VUT, BIF; year in local standard time)
  thr_src <- "VUT"
  thr <- tryCatch({
    th <- oneflux_ustar_ref(read_oneflux_bif(find_oneflux_bif(site_id, cfg, hh_path = files$fluxnet), "GRP_FLUX_REF",
                                             hh_path = files$fluxnet), "VUT")
    th$ustar_threshold[match(r$year_local, th$year)]
  }, error = function(e) {
    thr_src <<- paste0("none (", conditionMessage(e), ")")
    rep(NA_real_, nrow(r))
  })
  r$ustar_thr <- thr
  r$below_thr <- ifelse(is.finite(r$ustar_fn) & is.finite(thr), r$ustar_fn < thr, NA)
  r$ustar_ratio <- r$ustar_fn / thr
  # strata
  r$night <- r$sw_in_pot %in% 0
  r$night[is.na(r$sw_in_pot)] <- NA
  r$night_id <- night_groups(r$timestamp_start_utc, r$night, meta$lon)
  r$season <- as.character(season_of(r$timestamp_start_utc))
  r$zl <- zl_class(r$zm_m, r$L)
  ry <- if (!is.null(ref_years)) ref_years$year[ref_years$site_id == site_id] else integer()
  r$ref_year <- if (length(ry) == 1) r$year_local == ry else NA
  ag <- ustar_filter_agreement(r$below_thr[r$night %in% TRUE], r$nee_qc[r$night %in% TRUE])
  names(ag) <- paste0("ustar_agree_night_", names(ag))
  vars <- dplyr::bind_cols(vars, ag, dplyr::tibble(
    ustar_thr_source = thr_src, n_rows = nrow(r), first_utc = min(r$timestamp_start_utc),
    last_utc = max(r$timestamp_start_utc), labelling_date = as.Date(r$timestamp_start_utc[1] + off * 3600),
    reference_year = if (length(ry) == 1) ry else NA_integer_, zm_median_m = med_fin(r$zm_m),
    ec_height_median_m = med_fin(r$ec_height_m), canopy_height_median_m = med_fin(r$canopy_height_m),
    igbp = if ("igbp" %in% names(meta)) meta$igbp else NA_character_, measured_source = meas_src, measured_type = meas_type,
    n_h_era5 = sum(is.finite(r$h_era5_blh_center)), n_h_pblhmod = sum(is.finite(r$h_pblhmod_default)),
    n_h_clim = sum(is.finite(r$h_pblhmod_clim)), n_h_meas_exact = sum(is.finite(r$h_meas))))
  list(rows = r, vars = vars, meas_src = meas_src, meas_type = meas_type)
}

# --- Section 1 ------------------------------------------------------------------------------------------------------

#' Strata of spec 18 (day/night x zm/L class x season, with "all" margins)
#'
#' @param r Rows ([pblh_ec_site_rows()]).
#' @param seasons Include the season dimension [logical].
#' @return List of lists: `daynight`, `zl`, `season`, `i` (logical row index).
pec_strata <- function(r, seasons = TRUE) {
  dn <- ifelse(r$night %in% TRUE, "night", ifelse(r$night %in% FALSE, "day", NA))
  out <- list()
  for (a in c("all", "day", "night")) for (z in c("all", "unstable", "neutral", "stable"))
    for (s in if (seasons) c("all", "DJF", "MAM", "JJA", "SON") else "all") {
      i <- (a == "all" | dn %in% a) & (z == "all" | r$zl %in% z) & (s == "all" | r$season %in% s)
      out[[length(out) + 1]] <- list(daynight = a, zl = z, season = s, i = i)
    }
  out
}

#' Section 1: share of half-hours above the surface layer (zm/h above a threshold), by source
#'
#' @param r Rows ([pblh_ec_site_rows()]); @param site_id Station; @param thresholds zm/h thresholds [-] (0.1 main).
#' @param meas_src Measured source id or `NA` (own rows only; censored values counted, not used). Pass only a
#'   turbulence source (Cloudnet, valid at night); aerosol MLH at night sits on a floor, not a stable-layer depth.
#' @return Tibble: `site_id`, `source`, `row_basis` (`common` = rows with h of all three model sources, `own`),
#'   `daynight`, `zl`, `season`, `n`, `share_gt_0.05`, `share_gt_0.1`, `share_gt_0.2` [-], `zm_h_median`,
#'   `zm_h_p90` [-], `n_censored_below`, `n_lower_bound` (measured: censored rows in the stratum, not in the shares).
sl_share_table <- function(r, site_id, thresholds = c(0.05, 0.1, 0.2), meas_src = NA_character_) {
  pos <- function(x) is.finite(x) & x > 0
  srcs <- pec_sources()
  common <- Reduce(`&`, lapply(srcs, function(s) pos(r[[paste0("h_", s)]])))
  st <- pec_strata(r)
  hs <- lapply(stats::setNames(srcs, srcs), function(s) r[[paste0("h_", s)]])
  if (!is.na(meas_src)) hs[[meas_src]] <- r$h_meas
  rows <- list()
  for (s in names(hs)) {
    zmh <- r$zm_m / hs[[s]]
    okh <- pos(hs[[s]]) & is.finite(r$zm_m)
    for (b in if (s %in% srcs) c("common", "own") else "own") {
      base <- okh & (if (b == "common") common else TRUE)
      for (g in st) {
        i <- base & g$i
        z <- zmh[i]
        x <- dplyr::tibble(site_id = site_id, source = s, row_basis = b, daynight = g$daynight, zl = g$zl,
                           season = g$season, n = length(z))
        for (th in thresholds) x[[paste0("share_gt_", th)]] <- if (length(z)) mean(z > th) else NA_real_
        x$zm_h_median <- q_fin(z, 0.5)
        x$zm_h_p90 <- q_fin(z, 0.9)
        x$n_censored_below <- if (identical(s, meas_src)) sum(g$i & r$meas_kind %in% "below") else NA_integer_
        x$n_lower_bound <- if (identical(s, meas_src)) sum(g$i & r$meas_kind %in% "lower_bound") else NA_integer_
        rows[[length(rows) + 1]] <- x
      }
    }
  }
  dplyr::bind_rows(rows)
}

# --- Section 2 ------------------------------------------------------------------------------------------------------

#' Section 2: u*-filtered x above the surface layer (2 x 2) on night rows with all three model h
#'
#' @param r Rows; @param site_id Station; @param threshold zm/h threshold [-].
#' @return Tibble per source x filter definition (`threshold` = USTAR (FLUXNET) below the VUT threshold of the year;
#'   `qc` = NEE_VUT_REF_QC > 0): `n`, `n_filt_above`, `n_filt_below`, `n_pass_above`, `n_pass_below`,
#'   `share_above` [-], `share_filtered`, `share_above_removed` (of the above-surface-layer rows, share u*-filtered),
#'   `share_pass_above` (of the u*-passed rows, share above the surface layer).
ustar_cross_table <- function(r, site_id, threshold = 0.1) {
  pos <- function(x) is.finite(x) & x > 0
  srcs <- pec_sources()
  common <- Reduce(`&`, lapply(srcs, function(s) pos(r[[paste0("h_", s)]])))
  base <- common & r$night %in% TRUE & is.finite(r$zm_m)
  defs <- list(threshold = r$below_thr, qc = ifelse(is.na(r$nee_qc), NA, r$nee_qc > 0))
  rows <- list()
  for (s in srcs) for (dname in names(defs)) {
    f <- defs[[dname]]
    i <- base & !is.na(f)
    above <- (r$zm_m / r[[paste0("h_", s)]])[i] > threshold
    fi <- f[i]
    n <- function(a, b) sum(a & b)
    x <- dplyr::tibble(site_id = site_id, source = s, filter_def = dname, n = sum(i),
                       n_filt_above = n(fi, above), n_filt_below = n(fi, !above),
                       n_pass_above = n(!fi, above), n_pass_below = n(!fi, !above))
    x$share_above <- if (x$n) (x$n_filt_above + x$n_pass_above) / x$n else NA_real_
    x$share_filtered <- if (x$n) (x$n_filt_above + x$n_filt_below) / x$n else NA_real_
    x$share_above_removed <- with(x, if (n_filt_above + n_pass_above) n_filt_above / (n_filt_above + n_pass_above) else NA_real_)
    x$share_pass_above <- with(x, if (n_pass_above + n_pass_below) n_pass_above / (n_pass_above + n_pass_below) else NA_real_)
    rows[[length(rows) + 1]] <- x
  }
  dplyr::bind_rows(rows)
}

# --- Section 3 ------------------------------------------------------------------------------------------------------

#' Night-block bootstrap of a statistic
#'
#' Resamples whole nights with replacement (all rows of a drawn night, a night drawn k times enters k times).
#'
#' @param night_id Night id per row; @param stat Function of a row index vector returning a numeric vector.
#' @param n_boot Replicates; @param seed Seed; @param probs Interval quantiles.
#' @return List: `est` (stat on all rows), `lo`, `hi` (quantiles of the replicates, `NA` where fewer than half the
#'   replicates are finite), `reps` (matrix replicates x statistics).
night_block_boot <- function(night_id, stat, n_boot = 1000L, seed = 20261127L, probs = c(0.05, 0.95)) {
  idx <- split(seq_along(night_id), night_id)
  est <- stat(seq_along(night_id))
  if (length(idx) < 2) {
    return(list(est = est, lo = rep(NA_real_, length(est)), hi = rep(NA_real_, length(est)), reps = NULL))
  }
  set.seed(seed)
  reps <- t(vapply(seq_len(n_boot), function(b) {
    pick <- sample.int(length(idx), length(idx), replace = TRUE)
    stat(unlist(idx[pick], use.names = FALSE))
  }, numeric(length(est))))
  if (length(est) == 1) reps <- matrix(reps, ncol = 1)
  qq <- function(v, p) if (sum(is.finite(v)) >= n_boot / 2) unname(stats::quantile(v, p, na.rm = TRUE)) else NA_real_
  list(est = est, lo = apply(reps, 2, qq, probs[1]), hi = apply(reps, 2, qq, probs[2]), reps = reps)
}

#' Section 3: accepted night-time NEE (or FC) residual above vs within the surface layer
#'
#' Rows: night, u* passed (USTAR >= VUT threshold), measured (`NEE_VUT_REF_QC == 0`), with `RECO_NT_VUT_REF` and
#' `TA_F`; for `flux = "FC"` also FC and a profile SC (`SC_FLAG == 0`). Residual r = flux - RECO_NT [umol m-2 s-1];
#' relative r / RECO_NT on rows with RECO_NT > `reco_min`. Cells: u* bin ([ustar_ratio_bin()]) x TA_F tercile (of the
#' flux's base rows); in each cell the difference of the medians, zm/h > 0.1 minus zm/h <= 0.1, valid with at least
#' `min_n` rows on each side. Pooled values: mean of the valid cell differences weighted by the smaller side's n, per
#' u* bin (over the terciles) and overall. RECO_NT is fitted on the same night data: a bias is diluted, not created.
#'
#' @param r Rows; @param site_id Station; @param sources Sources (ERA5 first); @param fluxes `NEE`, `FC`.
#' @param threshold zm/h threshold; @param min_n Minimum rows per side; @param reco_min RECO_NT floor of the relative
#'   form [umol m-2 s-1]; @param n_boot,seed Bootstrap ([night_block_boot()], 90 % interval).
#' @return Tibble: `site_id`, `source`, `flux`, `measure` (`abs` [umol m-2 s-1], `rel` [-]), `stratum` (`cell`,
#'   `ustar_bin`, `all`), `ustar_bin`, `ta_tercile`, `n_above`, `n_below`, `n_nights`, `median_above`, `median_below`,
#'   `diff`, `lo90`, `hi90`, `valid` [logical].
nee_residual_table <- function(r, site_id, sources = c("era5_blh_center", "pblhmod_default"), fluxes = c("NEE", "FC"),
                               threshold = 0.1, min_n = 30L, reco_min = 0.5, n_boot = 1000L, seed = 20261127L) {
  out <- list()
  base0 <- r$night %in% TRUE & r$below_thr %in% FALSE & r$nee_qc %in% 0 & is.finite(r$reco_nt) & is.finite(r$ta_f_k) &
    is.finite(r$zm_m) & !is.na(r$night_id)
  for (fl in fluxes) {
    y <- if (fl == "NEE") r$nee_vut else r$fc
    base_f <- base0 & is.finite(y) & (fl == "NEE" | r$sc_flag %in% 0)
    if (sum(base_f) == 0) next
    tq <- stats::quantile(r$ta_f_k[base_f], c(1 / 3, 2 / 3), names = FALSE)
    terc <- as.character(cut(r$ta_f_k, c(-Inf, tq, Inf), labels = c("T1", "T2", "T3")))
    ub <- ustar_ratio_bin(r$ustar_ratio)
    for (s in sources) {
      h <- r[[paste0("h_", s)]]
      i_all <- which(base_f & is.finite(h) & h > 0 & !is.na(ub))
      if (!length(i_all)) next
      d <- data.frame(res = y[i_all] - r$reco_nt[i_all], rel = ifelse(r$reco_nt[i_all] > reco_min,
                                                                        (y[i_all] - r$reco_nt[i_all]) / r$reco_nt[i_all], NA_real_),
                      above = (r$zm_m[i_all] / h[i_all]) > threshold, ub = ub[i_all], tt = terc[i_all],
                      night = r$night_id[i_all])
      cells <- expand.grid(ub = c("1-1.5", "1.5-2", ">2"), tt = c("T1", "T2", "T3"), stringsAsFactors = FALSE)
      cell_of <- match(paste(d$ub, d$tt), paste(cells$ub, cells$tt))
      for (meas in c("abs", "rel")) {
        v <- if (meas == "abs") d$res else d$rel
        stat_fun <- function(ix) {
          ix <- ix[is.finite(v[ix])]
          ca <- cell_of[ix]; ab <- d$above[ix]; vv <- v[ix]
          dif <- na <- nb <- rep(NA_real_, nrow(cells))
          for (k in seq_len(nrow(cells))) {
            a <- vv[ca == k & ab]; b <- vv[ca == k & !ab]
            na[k] <- length(a); nb[k] <- length(b)
            if (length(a) >= min_n && length(b) >= min_n) dif[k] <- stats::median(a) - stats::median(b)
          }
          w <- pmin(na, nb)
          pool <- function(sel) { ok <- sel & is.finite(dif); if (any(ok)) sum(dif[ok] * w[ok]) / sum(w[ok]) else NA_real_ }
          c(dif, vapply(c("1-1.5", "1.5-2", ">2"), function(u) pool(cells$ub == u), numeric(1)), pool(rep(TRUE, nrow(cells))))
        }
        bt <- night_block_boot(d$night, stat_fun, n_boot, seed)
        lab <- dplyr::bind_rows(
          dplyr::tibble(stratum = "cell", ustar_bin = cells$ub, ta_tercile = cells$tt),
          dplyr::tibble(stratum = "ustar_bin", ustar_bin = c("1-1.5", "1.5-2", ">2"), ta_tercile = "all"),
          dplyr::tibble(stratum = "all", ustar_bin = "all", ta_tercile = "all"))
        sel_rows <- function(j) {
          ok <- is.finite(v)
          if (lab$stratum[j] == "cell") ok <- ok & d$ub == lab$ustar_bin[j] & d$tt == lab$ta_tercile[j]
          if (lab$stratum[j] == "ustar_bin") ok <- ok & d$ub == lab$ustar_bin[j]
          ok
        }
        x <- dplyr::bind_rows(lapply(seq_len(nrow(lab)), function(j) {
          ok <- sel_rows(j)
          dplyr::tibble(n_above = sum(ok & d$above), n_below = sum(ok & !d$above), n_nights = length(unique(d$night[ok])),
                        median_above = q_fin(v[ok & d$above], 0.5), median_below = q_fin(v[ok & !d$above], 0.5))
        }))
        out[[length(out) + 1]] <- dplyr::bind_cols(
          dplyr::tibble(site_id = site_id, source = s, flux = fl, measure = meas), lab, x,
          dplyr::tibble(diff = unname(bt$est), lo90 = unname(bt$lo), hi90 = unname(bt$hi), valid = is.finite(unname(bt$est))))
      }
    }
  }
  dplyr::bind_rows(out)
}

# --- Section 4 ------------------------------------------------------------------------------------------------------

#' Section 4 per night: storage accumulation and morning release against h
#'
#' Profile SC only (`SC_FLAG == 0`). Accumulation = sum of positive SC x 1800 s over the night's rows [umol m-2],
#' with at least `min_cov` of the night's rows valid; night h = median over the night's rows. Release = sum of
#' negative SC x 1800 s over the 8 half-hours from the first day row after the night [umol m-2] (at least 6 valid);
#' morning growth of h = h at the last of these rows minus h at the first [m] (PBLHmod: its B&G growth; ERA5: the
#' interpolated hourly values).
#'
#' @param r Rows; @param site_id Station; @param min_cov Minimum share of valid SC rows in a night [-].
#' @return Tibble per night: `site_id`, `night_id`, `date` [Date, evening], `ref_year`, `n_rows`, `cov`, `accum`,
#'   `h_night_<source>` [m], `release`, `dh_<source>` [m], `zm_m`.
storage_nights <- function(r, site_id, min_cov = 0.8) {
  srcs <- c("era5_blh_center", "pblhmod_default")
  scv <- ifelse(r$sc_flag %in% 0, r$sc, NA_real_)
  ids <- sort(unique(stats::na.omit(r$night_id)))
  if (!length(ids)) return(dplyr::tibble())
  rows <- split(seq_len(nrow(r)), r$night_id)
  out <- lapply(names(rows), function(k) {
    i <- rows[[k]]
    v <- scv[i]
    cov <- mean(is.finite(v))
    x <- dplyr::tibble(site_id = site_id, night_id = as.integer(k), date = as.Date(as.integer(k), origin = "1970-01-01"),
                       ref_year = any(r$ref_year[i] %in% TRUE), n_rows = length(i), cov = cov,
                       accum = if (cov >= min_cov) sum(pmax(v, 0), na.rm = TRUE) * 1800 else NA_real_,
                       zm_m = med_fin(r$zm_m[i]))
    for (s in srcs) x[[paste0("h_night_", s)]] <- med_fin(r[[paste0("h_", s)]][i])
    j0 <- max(i) + 1L
    jj <- j0:(j0 + 7L)
    okm <- all(jj <= nrow(r)) && all(r$night[jj] %in% FALSE) &&
      isTRUE(as.numeric(r$timestamp_start_utc[max(jj)]) - as.numeric(r$timestamp_start_utc[j0]) == 7 * 1800)
    vm <- if (okm) scv[jj] else rep(NA_real_, 8)
    x$release <- if (sum(is.finite(vm)) >= 6) sum(pmin(vm, 0), na.rm = TRUE) * 1800 else NA_real_
    for (s in srcs) {
      h <- r[[paste0("h_", s)]]
      x[[paste0("dh_", s)]] <- if (okm) h[max(jj)] - h[j0] else NA_real_
    }
    x
  })
  dplyr::bind_rows(out)
}

#' Section 4 site summary
#'
#' @param r Rows; @param nights Output of [storage_nights()]; @param site_id Station; @param threshold zm/h threshold.
#' @param h_breaks Bin edges of the night h [m].
#' @return Tibble long: `site_id`, `part` (`accumulation`, `release`, `accumulation_bin`, `sc_share`), `source`,
#'   `stratum` (bin or zm/h class), `n`, `statistic` (`spearman_rho`, `median`), `value` [-, or umol m-2 for the bins].
#'   sc_share: median SC / NEE (L2 NEE, NEE_DATA_FLAG 0, profile SC) on night rows with |NEE| > 0.5 umol m-2 s-1.
storage_h_table <- function(r, nights, site_id, threshold = 0.1, h_breaks = c(0, 50, 100, 200, 400, 800, 1600, Inf)) {
  srcs <- c("era5_blh_center", "pblhmod_default")
  out <- list()
  rho <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) >= 10) unname(suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))) else NA_real_
  }
  add <- function(...) out[[length(out) + 1]] <<- dplyr::tibble(site_id = site_id, ...)
  for (s in srcs) {
    if (nrow(nights)) {
      hn <- nights[[paste0("h_night_", s)]]
      add(part = "accumulation", source = s, stratum = "all", n = sum(is.finite(hn) & is.finite(nights$accum)),
          statistic = "spearman_rho", value = rho(nights$accum, hn))
      add(part = "release", source = s, stratum = "all",
          n = sum(is.finite(nights$release) & is.finite(nights[[paste0("dh_", s)]])), statistic = "spearman_rho",
          value = rho(nights$release, nights[[paste0("dh_", s)]]))
      b <- cut(hn, h_breaks, right = FALSE, dig.lab = 5)
      for (lv in levels(b)) {
        ok <- b %in% lv & is.finite(nights$accum)
        add(part = "accumulation_bin", source = s, stratum = lv, n = sum(ok), statistic = "median",
            value = q_fin(nights$accum[ok], 0.5))
      }
    }
    h <- r[[paste0("h_", s)]]
    ok <- r$night %in% TRUE & r$sc_flag %in% 0 & is.finite(r$sc) & is.finite(r$nee_l2) & r$nee_data_flag %in% 0 &
      abs(r$nee_l2) > 0.5 & is.finite(h) & h > 0 & is.finite(r$zm_m)
    cls <- ifelse((r$zm_m / h) > threshold, "zm_h_gt_0.1", "zm_h_le_0.1")
    for (k in c("zm_h_le_0.1", "zm_h_gt_0.1")) {
      i <- ok & cls %in% k
      add(part = "sc_share", source = s, stratum = k, n = sum(i), statistic = "median",
          value = q_fin(r$sc[i] / r$nee_l2[i], 0.5))
    }
  }
  dplyr::bind_rows(out)
}

# --- Section 5 ------------------------------------------------------------------------------------------------------

#' Section 5: Panofsky test of sigma_v / u* with each h source
#'
#' Rows: day (SW_IN_POT > 0), -zm/L in [0.05, 2], valid V_SIGMA and USTAR (L2), zm/h < 0.1 for the source.
#' ln ratio = ln(observed / predicted) with the prediction of [panofsky_sv_ustar()]. Ranking within the site on common
#' rows (all ranked sources valid): median |ln ratio - the source's median ln ratio| (the site-specific offset removed),
#' rank 1 = smallest; sets `models` (three model sources) and `with_measured` (plus the measured source, exact values).
#'
#' @param r Rows; @param site_id Station; @param meas_src Measured source or `NA`; @param threshold zm/h threshold.
#' @return List: `stats` (tibble `site_id`, `source`, `row_basis` (`own`, `common_models`, `common_with_measured`), `n`,
#'   `median_log_ratio`, `median_abs_log_ratio`, `median_abs_dev` [-], `slope` [-, OLS of ln(sigma_v/u*) on
#'   ln(12 - 0.5 h/L); 1/3 expected], `rank`) and `bins` (tibble `site_id`, `source`, `bin` (-h/L), `x_mid` [-],
#'   `n`, `sv_ustar_median` [-]).
panofsky_table <- function(r, site_id, meas_src = NA_character_, threshold = 0.1,
                           breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) {
  srcs <- pec_sources()
  hs <- lapply(stats::setNames(srcs, srcs), function(s) r[[paste0("h_", s)]])
  if (!is.na(meas_src)) hs[[meas_src]] <- r$h_meas
  zl <- r$zm_m / r$L
  base <- r$sw_in_pot > 0 & is.finite(zl) & -zl >= 0.05 & -zl <= 2 & is.finite(r$v_sigma) & is.finite(r$ustar) &
    r$ustar > 0 & r$v_sigma > 0
  base[is.na(base)] <- FALSE
  obs <- r$v_sigma / r$ustar
  ok_s <- lapply(hs, function(h) base & is.finite(h) & h > 0 & (r$zm_m / h) < threshold)
  lr_s <- lapply(hs, function(h) log(obs / panofsky_sv_ustar(h, r$L)))
  one <- function(s, i, basis) {
    lr <- lr_s[[s]][i]
    xx <- log(12 - 0.5 * hs[[s]][i] / r$L[i])
    yy <- log(obs[i])
    slope <- if (sum(i) >= 10) unname(stats::coef(stats::lm(yy ~ xx))[2]) else NA_real_
    dplyr::tibble(site_id = site_id, source = s, row_basis = basis, n = sum(i), median_log_ratio = q_fin(lr, 0.5),
                  median_abs_log_ratio = q_fin(abs(lr), 0.5),
                  median_abs_dev = if (length(lr)) stats::median(abs(lr - stats::median(lr))) else NA_real_,
                  slope = slope)
  }
  st <- list()
  for (s in names(hs)) st[[length(st) + 1]] <- dplyr::mutate(one(s, ok_s[[s]], "own"), rank = NA_integer_)
  sets <- list(common_models = srcs)
  if (!is.na(meas_src)) sets$common_with_measured <- c(srcs, meas_src)
  for (b in names(sets)) {
    i <- Reduce(`&`, ok_s[sets[[b]]])
    x <- dplyr::bind_rows(lapply(sets[[b]], function(s) one(s, i, b)))
    x$rank <- if (sum(i) >= 10) as.integer(rank(x$median_abs_dev, ties.method = "min")) else NA_integer_
    st[[length(st) + 1]] <- x
  }
  bins <- dplyr::bind_rows(lapply(names(hs), function(s) {
    i <- ok_s[[s]]
    x <- -hs[[s]][i] / r$L[i]
    b <- cut(x, breaks, right = FALSE)
    dplyr::bind_rows(lapply(levels(b), function(lv) {
      k <- b %in% lv
      dplyr::tibble(site_id = site_id, source = s, bin = lv, x_mid = q_fin(x[k], 0.5), n = sum(k),
                    sv_ustar_median = q_fin(obs[i][k], 0.5))
    }))
  }))
  list(stats = dplyr::bind_rows(st), bins = bins)
}

# --- Section 6 ------------------------------------------------------------------------------------------------------

#' Section 6: spectral correction factors against zm/h and -h/L in unstable conditions (descriptive)
#'
#' The low-frequency part of these EddyPro factors is analytic and contains no h: a relation reflects their shared
#' dependence on wind speed, zm and L, not a measured loss.
#'
#' @param r Rows; @param site_id Station.
#' @return Tibble: `site_id`, `source`, `x_var` (`zm_h`, `neg_h_L`), `bin`, `scf_var` (`H_SCF_STAT`, `LE_SCF_STAT`,
#'   `FC_SCF_STAT`, `TAU_SCF_STAT`), `n`, `median` [-], `x_mid` [-].
scf_h_table <- function(r, site_id, sources = c("era5_blh_center", "pblhmod_default"),
                        zmh_breaks = c(0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, Inf),
                        hl_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, Inf)) {
  scf <- c(H_SCF_STAT = "h_scf", LE_SCF_STAT = "le_scf", FC_SCF_STAT = "fc_scf", TAU_SCF_STAT = "tau_scf")
  base <- r$zl %in% "unstable" & is.finite(r$L) & r$L < 0
  out <- list()
  for (s in sources) {
    h <- r[[paste0("h_", s)]]
    ok <- base & is.finite(h) & h > 0
    xs <- list(zm_h = list(x = r$zm_m / h, br = zmh_breaks), neg_h_L = list(x = -h / r$L, br = hl_breaks))
    for (xv in names(xs)) {
      b <- cut(xs[[xv]]$x, xs[[xv]]$br, right = FALSE, dig.lab = 4)
      for (k in names(scf)) {
        v <- r[[scf[[k]]]]
        if (all(is.na(v))) next
        for (lv in levels(b)) {
          i <- ok & b %in% lv & is.finite(v)
          out[[length(out) + 1]] <- dplyr::tibble(site_id = site_id, source = s, x_var = xv, bin = lv, scf_var = k,
                                                  n = sum(i), median = q_fin(v[i], 0.5),
                                                  x_mid = q_fin(xs[[xv]]$x[i], 0.5))
        }
      }
    }
  }
  dplyr::bind_rows(out)
}

#' Headers of the spec-18 tables
#'
#' @return Named list of character vectors (one per table).
pec_table_headers <- function() {
  circ <- "Circularity: PBLHmod's night-time h is computed from u* and L (spec 02, B1), so with PBLHmod zm/h > 0.1 is close to low u*; ERA5 is the main source."
  rows <- "Rows: every half-hour of the L2 record from the labelling date (first valid pblhmod_default value), exclusions removed; zm = EC height - 2/3 hc (config/sites_timevarying.csv); night = SW_IN_POT = 0; zl = zm/L class (unstable < -0.1 <= neutral <= 0.1 < stable); season from the UTC month."
  hsrc <- "h: era5_blh_center (instantaneous, interpolated to the half-hour mid), pblhmod_default, pblhmod_clim, the measured primary source (exact values only; censored values counted apart, never imputed)."
  list(
    variables = c("Spec 18: variable availability per site (n_<VAR> = non-missing values in the FLUXES file, n_fn_<VAR> in the FLUXNET HH file, NA = column absent).",
                  "storage = CO2 profile levels (CO2_DRY_<n>) and SC with SC_FLAG 0 (profile; flag 1 = single point, SC = SC_1P); sc_check_share = share of rows with |NEE - FC - SC| <= 0.011 umol m-2 s-1.",
                  "ustar_agree_night_*: night rows, u*-filtered by threshold (FLUXNET USTAR < ONEFlux VUT threshold of the local year, BIF) vs NEE_VUT_REF_QC > 0.", rows),
    sl_share = c("Spec 18 section 1: share of half-hours with zm/h above 0.05, 0.1 (main), 0.2; median and p90 of zm/h [-].",
                 "row_basis common = rows with h of all three model sources; own = each source on its rows.", rows, hsrc, circ),
    ustar_cross = c("Spec 18 section 2: night rows with all three model h; u*-filtered x above the surface layer (zm/h > 0.1).",
                    "filter_def threshold = FLUXNET USTAR < VUT threshold of the local year; qc = NEE_VUT_REF_QC > 0 (gap-filled, includes other removals).",
                    "share_above_removed = filtered / all above; share_pass_above = above / all u*-passed.", rows, circ),
    nee_residual = c("Spec 18 section 3: night, u*-passed (USTAR >= VUT threshold), measured (NEE_VUT_REF_QC = 0) rows. r = flux - RECO_NT_VUT_REF [umol m-2 s-1] (measure abs); rel = r / RECO_NT on rows with RECO_NT > 0.5.",
                     "flux FC: FC (L2) at rows with a profile SC (SC_FLAG 0). diff = median(r | zm/h > 0.1) - median(r | zm/h <= 0.1) in each u* bin (u*/u*_thr) x TA_F tercile cell, valid with >= 30 rows each side; ustar_bin / all = mean of the valid cell differences weighted by min(n_above, n_below).",
                     "lo90, hi90 = night-block bootstrap (whole nights resampled, 1000 replicates, seed 20261127). RECO_NT is fitted on the same night data, so a bias would be diluted, not created.", rows, circ),
    storage_h = c("Spec 18 section 4 (sites with a CO2 profile; SC_FLAG 0 only): accumulation = sum of positive SC x 1800 s per night [umol m-2] (>= 80 % of the night's rows valid) vs the night median h; release = sum of negative SC x 1800 s over the 4 h from the first day row [umol m-2] (>= 6 of 8 valid) vs the growth of h over the same rows [m].",
                  "spearman_rho per site and source; accumulation_bin = median accumulation per night-h bin [m]; sc_share = median SC / NEE (L2 NEE, NEE_DATA_FLAG 0) on night rows with |NEE| > 0.5 umol m-2 s-1, by zm/h class.", rows, circ),
    storage_nights = c("Spec 18 section 4: one row per night (sunset to sunrise; night_id = evening's local solar date). accum, release [umol m-2]; h_night_<source> [m, median]; dh_<source> [m, growth over the release hours].", rows),
    panofsky = c("Spec 18 section 5: Panofsky et al. (1977): sigma_v / u* = (12 - 0.5 h / L)^(1/3). Rows: day (SW_IN_POT > 0), 0.05 <= -zm/L <= 2, V_SIGMA and USTAR (L2) valid, zm/h < 0.1 for the source.",
                 "median_log_ratio = median ln(observed / predicted); median_abs_dev = median |ln ratio - its median| (site offset removed); slope = OLS slope of ln(sigma_v/u*) on ln(12 - 0.5 h/L) (1/3 expected).",
                 "rank (1 = best, by median_abs_dev) on common rows: common_models = three model sources; common_with_measured = plus the measured source (exact values).", rows, circ),
    panofsky_bins = c("Spec 18 section 5: binned median sigma_v / u* against -h/L per source (rows as 18_panofsky, own rows); x_mid = median -h/L in the bin."),
    scf_h = c("Spec 18 section 6 (descriptive): total spectral correction factors (EddyPro *_SCF_STAT) against zm/h and -h/L, unstable rows (zm/L < -0.1); binned medians.",
              "The low-frequency part of these factors is analytic and contains no h: any relation reflects their shared dependence on wind speed, zm and L, not a measured loss.", rows)
  )
}
