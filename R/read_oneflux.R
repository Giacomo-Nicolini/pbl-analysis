# ONEFlux / FLUXNET BIF files (ICOS ETC; spec 08, G 2026-10-26).
# `ICOSETC_<site>_FLUXNET_BIF_<years>_v<version>.csv` replaces the AUXNEE file: long format with SITE_ID, GROUP_ID,
# VARIABLE_GROUP, VARIABLE, DATAVALUE; one GROUP_ID per record of a variable group.

#' Find the ONEFlux BIF file of a station
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `icos_l2_dir`.
#' @return Path [character]; aborts if none or more than one.
find_oneflux_bif <- function(site_id, cfg) {
  dir <- file.path(resolve_path(cfg$icos_l2_dir), paste0(site_id, "_output_L2I_2026"))
  f <- list.files(dir, pattern = paste0("^ICOSETC_", gsub("-", "\\\\-", site_id), "_FLUXNET_BIF_.*\\.csv$"),
                  full.names = TRUE)
  if (length(f) != 1) cli::cli_abort("{site_id}: {length(f)} ONEFlux BIF file{?s} in {.path {dir}} (expected 1).")
  f
}

# Product years and version of a FLUXNET file name, e.g.
# "ICOSETC_FI-Hyy_FLUXNET_FLUXMET_HH_2018-2026_v1.3_r1.zip" -> "2018-2026_v1.3_r1".
fluxnet_product_tag <- function(x) {
  m <- regmatches(basename(x), regexpr("[0-9]{4}-[0-9]{4}_v[0-9.]+_r[0-9]+", basename(x)))
  if (length(m) == 0) NA_character_ else m
}

#' Read a ONEFlux BIF file
#'
#' Reads the long-format BIF file and returns one variable group as a wide table
#' (one row per GROUP_ID), or all groups as a named list. If `hh_path` is given,
#' checks that the product years and version in `GRP_ONEFLUX/PRODUCT_NAME` match
#' those of the half-hourly FLUXNET file and aborts otherwise (the thresholds must
#' belong to the product actually used; G 2026-10-26).
#'
#' @param path Path of the BIF csv [character].
#' @param group `VARIABLE_GROUP` to return (e.g. `"GRP_FLUX_REF"`), or `NULL` for all.
#' @param hh_path Path of the FLUXNET FLUXMET HH file to check against, or `NULL`.
#' @return Tibble (`SITE_ID`, `GROUP_ID`, one character column per `VARIABLE`), or
#'   a named list of such tibbles. Attribute `product_name` [character].
read_oneflux_bif <- function(path, group = NULL, hh_path = NULL) {
  b <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
  need <- c("SITE_ID", "GROUP_ID", "VARIABLE_GROUP", "VARIABLE", "DATAVALUE")
  miss <- setdiff(need, names(b))
  if (length(miss) > 0) cli::cli_abort("{.path {basename(path)}}: missing column{?s} {.field {miss}}.")
  pn <- b$DATAVALUE[b$VARIABLE_GROUP == "GRP_ONEFLUX" & b$VARIABLE == "PRODUCT_NAME"]
  pn <- if (length(pn)) pn[1] else NA_character_
  if (!is.null(hh_path)) {
    a <- fluxnet_product_tag(pn)
    h <- fluxnet_product_tag(hh_path)
    if (is.na(a) || is.na(h) || a != h) {
      cli::cli_abort(c("{.path {basename(path)}}: ONEFlux product {.val {pn}} does not match the HH file {.val {basename(hh_path)}}.",
                       "i" = "Product years/version {.val {a}} vs {.val {h}}: thresholds of another processing run."))
    }
  }
  wide <- function(g) {
    x <- b[b$VARIABLE_GROUP == g, ]
    if (nrow(x) == 0) cli::cli_abort("{.path {basename(path)}}: no group {.val {g}}.")
    if (anyDuplicated(x[c("GROUP_ID", "VARIABLE")])) cli::cli_abort("{g}: duplicated VARIABLE within a GROUP_ID.")
    w <- tidy_pivot_wider(x[c("SITE_ID", "GROUP_ID", "VARIABLE", "DATAVALUE")])
    w[order(as.numeric(w$GROUP_ID)), ]
  }
  out <- if (is.null(group)) {
    stats::setNames(lapply(sort(unique(b$VARIABLE_GROUP)), wide), sort(unique(b$VARIABLE_GROUP)))
  } else {
    wide(group)
  }
  attr(out, "product_name") <- pn
  out
}

# Long (SITE_ID, GROUP_ID, VARIABLE, DATAVALUE) -> wide, one row per GROUP_ID (base R, no tidyr).
tidy_pivot_wider <- function(x) {
  ids <- unique(x$GROUP_ID)
  vars <- unique(x$VARIABLE)
  out <- dplyr::tibble(SITE_ID = x$SITE_ID[match(ids, x$GROUP_ID)], GROUP_ID = ids)
  for (v in vars) {
    s <- x[x$VARIABLE == v, ]
    out[[v]] <- s$DATAVALUE[match(ids, s$GROUP_ID)]
  }
  out
}

#' Reference u* thresholds of the ONEFlux NEE product
#'
#' From `GRP_FLUX_REF` (output of [read_oneflux_bif()]): the threshold used for the
#' reference NEE (`REF_VARIABLE`, `REF_VERSION`, `REF_TIME_RESOLUTION`). VUT: one
#' threshold per `REF_YEAR`; CUT: one value for all years (`year = NA`).
#'
#' @param flux_ref Wide `GRP_FLUX_REF` table.
#' @param version `"VUT"` or `"CUT"`.
#' @param variable Reference variable [character], default `"NEE"`.
#' @param time_resolution `"HH"` (half-hourly product) or other.
#' @return Tibble: `year` [integer; `NA` for CUT], `ustar_threshold` [m s-1],
#'   `ustar_percentile` [-].
oneflux_ustar_ref <- function(flux_ref, version = c("VUT", "CUT"), variable = "NEE", time_resolution = "HH") {
  version <- match.arg(version)
  col <- function(k) if (k %in% names(flux_ref)) flux_ref[[k]] else rep(NA_character_, nrow(flux_ref))
  sel <- col("REF_VARIABLE") == variable & col("REF_VERSION") == version & col("REF_TIME_RESOLUTION") == time_resolution
  sel[is.na(sel)] <- FALSE
  x <- flux_ref[sel, ]
  if (nrow(x) == 0) cli::cli_abort("No {variable} {version} {time_resolution} reference threshold in GRP_FLUX_REF.")
  out <- dplyr::tibble(
    year = suppressWarnings(as.integer(if ("REF_YEAR" %in% names(x)) x$REF_YEAR else NA)),
    ustar_threshold = as.numeric(x$REF_USTAR_THRESHOLD),
    ustar_percentile = suppressWarnings(as.numeric(if ("REF_USTAR_PERCENTILE" %in% names(x)) x$REF_USTAR_PERCENTILE else NA))
  )
  if (version == "VUT" && (anyNA(out$year) || anyDuplicated(out$year))) {
    cli::cli_abort("VUT thresholds: missing or duplicated REF_YEAR.")
  }
  if (version == "CUT" && nrow(out) != 1) cli::cli_abort("CUT: expected one threshold, found {nrow(out)}.")
  out[order(out$year), ]
}

#' Check the ONEFlux u* thresholds against the FLUXNET NEE QC
#'
#' On the half-hours with QC-good FLUXES NEE (`NEE_DATA_FLAG == 0`): share of the
#' u*-filtered rows ([oneflux_ustar_filtered()]) that are gap-filled in FLUXNET
#' (`NEE_VUT_REF_QC > 0`; expected about 100 %), and share of the gap-filled rows
#' (the exclusions of the spec-08 subset `valid_qc_used`) explained by u*.
#'
#' @param d Table with `USTAR` [m s-1], `NEE`, `NEE_DATA_FLAG`, `fn_NEE_VUT_REF_QC`,
#'   `fn_NIGHT`, `TIMESTAMP_START` (local standard time).
#' @param thr Output of [oneflux_ustar_ref()] (VUT).
#' @return Tibble per year and `period` (`day`, `night`, `all`): `n_good`,
#'   `n_ustar_filtered`, `n_gapfilled`, `pct_ustar_filtered_gapfilled`,
#'   `pct_gapfilled_explained_by_ustar` [%].
oneflux_ustar_check <- function(d, thr) {
  good <- !is.na(d$NEE) & d$NEE_DATA_FLAG %in% 0 & !is.na(d$fn_NEE_VUT_REF_QC)
  x <- d[good, ]
  yr <- as.integer(substr(x$TIMESTAMP_START, 1, 4))
  uf <- oneflux_ustar_filtered(x$USTAR, yr, thr) %in% TRUE
  gf <- x$fn_NEE_VUT_REF_QC > 0
  per <- ifelse(x$fn_NIGHT %in% 1, "night", "day")
  f <- function(i) dplyr::tibble(n_good = length(i), n_ustar_filtered = sum(uf[i]), n_gapfilled = sum(gf[i]),
                                 pct_ustar_filtered_gapfilled = round(100 * sum(uf[i] & gf[i]) / max(1, sum(uf[i])), 2),
                                 pct_gapfilled_explained_by_ustar = round(100 * sum(uf[i] & gf[i]) / max(1, sum(gf[i])), 2))
  rows <- list()
  for (y in sort(unique(yr))) {
    for (p in c("day", "night", "all")) {
      i <- which(yr == y & (p == "all" | per == p))
      rows[[length(rows) + 1]] <- dplyr::bind_cols(dplyr::tibble(year = y, period = p), f(i))
    }
  }
  for (p in c("day", "night", "all")) {
    i <- which(p == "all" | per == p)
    rows[[length(rows) + 1]] <- dplyr::bind_cols(dplyr::tibble(year = NA_integer_, period = p), f(i))
  }
  dplyr::bind_rows(rows)
}

#' Rows removed by the ONEFlux u* filter
#'
#' ONEFlux removes NEE where u* is below the threshold of the year, by day and by
#' night (G 2026-10-26).
#'
#' @param ustar Friction velocity [m s-1].
#' @param year Year of the half-hour (`TIMESTAMP_START`, local standard time) [integer].
#' @param thr Output of [oneflux_ustar_ref()].
#' @return Logical: `TRUE` = below the threshold (u*-filtered); `NA` if u* or the
#'   threshold of the year is missing.
oneflux_ustar_filtered <- function(ustar, year, thr) {
  t <- if (all(is.na(thr$year))) rep(thr$ustar_threshold[1], length(ustar)) else thr$ustar_threshold[match(year, thr$year)]
  ustar < t
}
