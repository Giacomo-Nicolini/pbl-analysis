# Footprint chunks (spec 08 section 3, spec 09 A5): one chunk = one site x PBLH case x month (or half-month when the
# estimated duration exceeds the target). Tasks are built from the 08a input table; a chunk runs the FFP per half-hour
# (Kljun's code, fixed grid), writes the per-half-hour metrics (1 m fields), the summed fields as 5 m block sums
# (storage c5) and the climatologies from them. Within a chunk the half-hours are processed in fixed blocks (one day)
# by the workers; block sums are added in block order, so the result does not depend on the number of workers.

#' Chunk part labels and time ranges
#'
#' @param month Month `"YYYY-MM"` [character].
#' @param split `FALSE`: the whole month (`"YYYY-MM"`); `TRUE`: halves `"YYYY-MMa"` (days 1-15) and `"YYYY-MMb"`
#'   (day 16 to the end).
#' @return Tibble: `part`, `start`, `end` [POSIXct UTC, half-hour starts in `[start, end)`].
chunk_parts <- function(month, split = FALSE) {
  m0 <- as.POSIXct(paste0(month, "-01"), tz = "UTC")
  m1 <- as.POSIXct(seq(as.Date(m0), by = "month", length.out = 2)[2], tz = "UTC")
  if (!split) return(dplyr::tibble(part = month, start = m0, end = m1))
  mid <- m0 + 15 * 86400
  dplyr::tibble(part = paste0(month, c("a", "b")), start = c(m0, mid), end = c(mid, m1))
}

#' Rows for which the FFP is computed
#'
#' @param inp Footprint input object of `scripts/08a_footprint_flags.R` (`tab`, `subsets`).
#' @param case PBLH case [character].
#' @param sets Validity sets.
#' @return Logical per row of `inp$tab`: valid in at least one requested set.
footprint_run_rows <- function(inp, case, sets) {
  Reduce(`|`, lapply(sets, function(s) {
    v <- inp$subsets[[paste("valid", s, case, sep = "__")]]
    if (is.null(v)) rep(FALSE, nrow(inp$tab)) else v
  }))
}

#' Build one footprint chunk task
#'
#' Rows of `case` valid in at least one set within `[start, end)`, and the climatology groups
#' `set|subset|day or night|rowset`: `own` (rows valid for the case), `common_pe` (also valid for `pblhmod` and
#' `era5`; every case except `measured` and `bound`), `common_pme` (also valid for `measured`; cases `pblhmod`,
#' `pblhmod_clim`, `era5`, `measured`). Common groups are restricted to the case's own valid rows (no change for
#' `pblhmod` and `era5`, where they are implied).
#'
#' @param inp Footprint input object (08a).
#' @param case PBLH case.
#' @param part Part label (see [chunk_parts()]).
#' @param start,end Part range [POSIXct UTC].
#' @param sets Validity sets.
#' @param limit Testing only: at most this many half-hours [integer] (`NA`: all).
#' @return List (`site_id`, `case`, `part`, `rows`, `grp`), or `NULL` without rows.
footprint_task <- function(inp, case, part, start, end, sets = c("kljun", "kljun_norsl", "etc", "none"), limit = NA) {
  tab <- inp$tab
  sub_cols <- inp$subsets
  run <- footprint_run_rows(inp, case, sets)
  idx <- which(run & tab$timestamp_start_utc >= start & tab$timestamp_start_utc < end)
  if (!is.na(limit)) idx <- utils::head(idx, limit)
  if (length(idx) == 0) return(NULL)
  rows <- tab[idx, c("timestamp_start_utc", "zm_m", "umean", "L", "sigmav", "ustar", "wind_dir", "day", paste0("h_", case))]
  names(rows)[ncol(rows)] <- "h"
  get <- function(u, s, k) {
    v <- sub_cols[[paste(u, s, k, sep = "__")]]
    if (is.null(v)) rep(FALSE, nrow(tab)) else v
  }
  has_measured <- "h_measured" %in% names(tab)
  grp <- list()
  for (s in sets) for (u in c("valid", "valid_qc", "valid_qc_used", "valid_qc_h")) {
    own <- get(u, s, case)[idx]
    for (dn in c("day", "night")) {
      dmask <- if (dn == "day") rows$day %in% TRUE else rows$day %in% FALSE
      grp[[paste(s, u, dn, "own", sep = "|")]] <- own & dmask
      pe <- (get(u, s, "pblhmod") & get(u, s, "era5"))[idx]
      if (!case %in% c("measured", "bound")) grp[[paste(s, u, dn, "common_pe", sep = "|")]] <- pe & own & dmask
      if (has_measured && case %in% c("pblhmod", "pblhmod_clim", "era5", "measured")) {
        grp[[paste(s, u, dn, "common_pme", sep = "|")]] <- pe & get(u, s, "measured")[idx] & own & dmask
      }
    }
  }
  grp <- grp[vapply(grp, any, logical(1))]
  list(site_id = inp$site_id, case = case, part = part, rows = rows, grp = grp)
}

#' Chunk list with duration estimates
#'
#' One chunk per case x month of the year; a month whose estimated duration exceeds `target_min` is split into two
#' halves ([chunk_parts()]).
#'
#' @param inp Footprint input object (08a, full period).
#' @param year Reference year [integer].
#' @param cases PBLH cases.
#' @param sets Validity sets.
#' @param sec_per_hh FFP cost [s per half-hour, one core].
#' @param cpus Workers per chunk [-].
#' @param target_min Target chunk duration [min].
#' @return Tibble: `site_id`, `case`, `part`, `start`, `end` [character, UTC], `n_halfhours`, `est_min` [min wall
#'   time with `cpus` workers], `cpu_hours` [h].
make_chunk_list <- function(inp, year, cases, sets = c("kljun", "kljun_norsl", "etc", "none"), sec_per_hh = 2,
                            cpus = 6, target_min = 120) {
  t <- inp$tab$timestamp_start_utc
  t_mo <- format(t, "%Y-%m")
  months <- sprintf("%d-%02d", as.integer(year), 1:12)
  out <- list()
  for (k in cases) {
    run <- footprint_run_rows(inp, k, sets)
    for (mo in months) {
      n_mo <- sum(run & t_mo == mo)
      if (n_mo == 0) next
      est <- n_mo * sec_per_hh / cpus / 60
      pp <- chunk_parts(mo, split = est > target_min)
      for (j in seq_len(nrow(pp))) {
        n <- sum(run & t >= pp$start[j] & t < pp$end[j])
        if (n == 0) next
        out[[length(out) + 1]] <- dplyr::tibble(site_id = inp$site_id, case = k, part = pp$part[j],
                                                start = format(pp$start[j], "%Y-%m-%d %H:%M"),
                                                end = format(pp$end[j], "%Y-%m-%d %H:%M"), n_halfhours = n,
                                                est_min = round(n * sec_per_hh / cpus / 60, 1),
                                                cpu_hours = round(n * sec_per_hh / 3600, 2))
      }
    }
  }
  if (length(out) == 0) {
    return(dplyr::tibble(site_id = character(), case = character(), part = character(), start = character(),
                         end = character(), n_halfhours = integer(), est_min = numeric(), cpu_hours = numeric()))
  }
  dplyr::bind_rows(out)
}

# FFP and metrics for the rows `idx` of a task: per-half-hour metrics (1 m field, smoothed) and the 5 m block sums of
# the unsmoothed fields per group. `env` = load_ffp_code() of the worker.
ffp_chunk_block <- function(idx, tk, env, grid, mask, r_levels, block) {
  rows <- tk$rows
  bm <- vector("list", length(idx))
  sums <- list()
  cnt <- list()
  for (j in seq_along(idx)) {
    i <- idx[j]
    fld <- function() ffp_field_halfhour(env, grid, rows$zm_m[i], rows$h[i], rows$umean[i], rows$L[i], rows$sigmav[i],
                                         rows$ustar[i], rows$wind_dir[i])
    f <- fld()
    if (is.null(f)) {
      bm[[j]] <- dplyr::tibble(ffp_ok = FALSE)
      next
    }
    # Guard (spec 08, 2026-10-27): a field with non-finite cells is recomputed once; cells still non-finite are set
    # to 0 and counted, so they never enter the sums.
    n_retry <- 0L
    if (any(!is.finite(f))) { f <- fld(); n_retry <- 1L }
    n_nonfinite <- sum(!is.finite(f))
    f[!is.finite(f)] <- 0
    met <- function() ffp_field_metrics(ffp_smooth(f), grid, mask, r_levels, keep_region = NULL)$metrics
    m <- met()
    if (!all(is.finite(c(m$f_domain, m$f_target_area)))) { m <- met(); n_retry <- n_retry + 1L }
    bm[[j]] <- dplyr::bind_cols(dplyr::tibble(ffp_ok = TRUE, n_nonfinite = n_nonfinite, n_retry = n_retry), m)
    f5 <- ffp_block_sums(f, block)
    for (g in names(tk$grp)) {
      if (!tk$grp[[g]][i]) next
      if (is.null(sums[[g]])) { sums[[g]] <- f5; cnt[[g]] <- 1L } else { sums[[g]] <- sums[[g]] + f5; cnt[[g]] <- cnt[[g]] + 1L }
    }
  }
  list(bm = bm, sums = sums, cnt = cnt)
}

#' Run one footprint chunk
#'
#' FFP per half-hour, in blocks of `block_rows` half-hours spread over `workers` processes (fork on Linux, PSOCK on
#' Windows); block sums added in block order (result independent of `workers`). Writes
#' `<chunk_dir>/<site>_<case>_<part>_hh.csv.gz` (per-half-hour distances and field metrics at 1 m), `_sums.rds`
#' (5 m block sums per group, float32), `_clim.rds` (day, night and all climatologies per group from the block sums,
#' as the annual ones; 80 % regions packed, [ffp_region_pack()]) and `_meta.rds` (hash, written last: a chunk
#' without it is recomputed).
#'
#' @param tk Task ([footprint_task()]) with `hash`.
#' @param ffp_dir Folder of Kljun's FFP code.
#' @param grid [ffp_grid()]; @param mask Target-area mask; @param r_levels Source-area levels.
#' @param chunk_dir Output folder.
#' @param workers Processes [-].
#' @param block Storage block size [cells] (5 = c5).
#' @param block_rows Half-hours per work block [-].
#' @return Message [character].
run_footprint_chunk <- function(tk, ffp_dir, grid, mask, r_levels, chunk_dir, workers = 1L, block = 5L,
                                block_rows = 48L) {
  t_start <- Sys.time()
  rows <- tk$rows
  n <- nrow(rows)
  a <- ffp_ci_paper(rows$zm_m, rows$h, rows$umean, rows$ustar, r_levels)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / block_rows))
  sums <- list()
  cnt <- list()
  bm <- list()
  add <- function(res) {
    for (g in names(res$sums)) {
      if (is.null(sums[[g]])) { sums[[g]] <<- res$sums[[g]]; cnt[[g]] <<- res$cnt[[g]] } else {
        sums[[g]] <<- sums[[g]] + res$sums[[g]]
        cnt[[g]] <<- cnt[[g]] + res$cnt[[g]]
      }
    }
    bm <<- c(bm, res$bm)
  }
  if (workers > 1 && length(blocks) > 1) {
    cl <- if (.Platform$OS.type == "windows") parallel::makePSOCKcluster(workers) else parallel::makeForkCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (.Platform$OS.type == "windows") {
      rfiles <- list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)
      parallel::clusterCall(cl, function(fs) { for (f in fs) source(f); NULL }, rfiles)
    }
    parallel::clusterCall(cl, function(d) { assign(".ffp_env", load_ffp_code(d), envir = globalenv()); NULL }, ffp_dir)
    # rounds of `workers` blocks, reduced in block order (bounded memory on the main process)
    for (r0 in seq(1, length(blocks), by = workers)) {
      ii <- r0:min(r0 + workers - 1, length(blocks))
      res <- parallel::parLapply(cl, blocks[ii], function(idx, tk, grid, mask, r_levels, block) {
        ffp_chunk_block(idx, tk, get(".ffp_env", envir = globalenv()), grid, mask, r_levels, block)
      }, tk = tk, grid = grid, mask = mask, r_levels = r_levels, block = block)
      for (x in res) add(x)
    }
  } else {
    env <- load_ffp_code(ffp_dir)
    for (idx in blocks) add(ffp_chunk_block(idx, tk, env, grid, mask, r_levels, block))
  }
  hh <- dplyr::bind_cols(dplyr::tibble(timestamp_start_utc = rows$timestamp_start_utc, case = tk$case, h = rows$h), a,
                         dplyr::bind_rows(bm))
  dim_full <- c(length(grid$y), length(grid$x))
  part <- list(n = unlist(cnt), sums = lapply(sums, ffp_block_upsample, dim_full = dim_full, k = block))
  clim <- ffp_clims_from_parts(list(part), grid, mask, r_levels)
  b <- file.path(chunk_dir, paste0(tk$site_id, "_", tk$case, "_", tk$part))
  readr::write_csv(hh, paste0(b, "_hh.csv.gz"))
  n_cells <- prod(dim_full)
  saveRDS(lapply(clim, function(c0) { c0$region80 <- ffp_region_pack(c0$region80, n_cells); c0 }), paste0(b, "_clim.rds"))
  ffp_sums_write(sums, cnt, paste0(b, "_sums.rds"), block = block, dim_full = dim_full)
  secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  saveRDS(list(hash = tk$hash, n = n, seconds = secs, workers = workers, block = block), paste0(b, "_meta.rds"))
  sprintf("%s %s %s %s: %d half-hours, %d climatologies, %.0f s, %d worker(s)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          tk$site_id, tk$case, tk$part, n, length(clim), secs, as.integer(workers))
}

#' Pack / unpack a source region
#'
#' Regions are sets of grid-cell indices (e.g. the 80 % source area). In chunk climatologies they are stored as a
#' packed bitmask over the grid (lossless; spec 09 A5): a month's region as integer indices takes up to about 100 MB
#' per chunk (spec-08 test month), the bitmask a few tens of kB after compression.
#'
#' @param idx Cell indices [integer]; @param n_cells Number of grid cells.
#' @param x Packed region, or an integer index vector (files of spec 08), returned unchanged.
#' @return [ffp_region_pack()]: list (`n_cells`, `bits` [raw]); [ffp_region_unpack()]: sorted indices [integer].
ffp_region_pack <- function(idx, n_cells) {
  b <- logical(n_cells + (-n_cells) %% 8)
  b[idx] <- TRUE
  list(n_cells = n_cells, bits = packBits(b))
}

#' @rdname ffp_region_pack
ffp_region_unpack <- function(x) {
  if (!is.list(x)) return(x)
  which(as.logical(rawToBits(x$bits))[seq_len(x$n_cells)])
}

#' Climatologies of every group from summed fields
#'
#' Day, night and all (= day + night) per `set|subset|rowset`, from the sums of all parts ([ffp_clim_from_sums()]).
#'
#' @param parts List of [ffp_sums_read()]-like lists (`n`, `sums` on the full grid).
#' @inheritParams ffp_field_metrics
#' @return List of climatologies (`set`, `subset`, `daynight`, `rowset`, `n`, `metrics`, `region80`).
ffp_clims_from_parts <- function(parts, grid, mask, r = c(0.5, 0.7, 0.8, 0.9)) {
  gk <- unique(unlist(lapply(parts, function(p) names(p$sums))))
  out <- list()
  for (key in unique(sub("\\|(day|night)\\|", "|", gk))) {
    p <- strsplit(key, "\\|")[[1]]
    for (dn in c("day", "night", "all")) {
      k <- paste(p[1], p[2], if (dn == "all") c("day", "night") else dn, p[3], sep = "|")
      res <- ffp_clim_from_sums(parts, k, grid, mask, r)
      if (is.null(res)) next
      out[[length(out) + 1]] <- list(set = p[1], subset = p[2], daynight = dn, rowset = p[3], n = res$n,
                                     metrics = res$metrics, region80 = res$region)
    }
  }
  out
}

#' Chunk parts on disk for a site, case and month
#'
#' @param chunk_dir Chunk folder; @param site_id Station; @param case PBLH case; @param month `"YYYY-MM"`.
#' @return Part labels with a complete chunk (`_meta.rds` present): the month, its two halves, or none.
chunk_parts_on_disk <- function(chunk_dir, site_id, case, month) {
  has <- function(p) file.exists(file.path(chunk_dir, paste0(site_id, "_", case, "_", p, "_meta.rds")))
  if (has(month)) return(month)
  halves <- paste0(month, c("a", "b"))
  if (all(vapply(halves, has, logical(1)))) return(halves)
  character()
}
