# Spec 09 A5: c5 block sums, chunk parts and list, chunk tasks, footprint cases at extra / other sites.

utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("ffp_block_sums: mass conservation, tower cell at a block centre; upsample keeps the sum", {
  set.seed(1)
  n <- 2001
  m <- matrix(stats::runif(n * n), n, n)
  b <- ffp_block_sums(m, 5)
  expect_equal(dim(b), c(401, 401))
  expect_equal(sum(b), sum(m))
  # the block containing the tower cell (1001, 1001) is centred on it: cells 999-1003 in both directions
  one <- matrix(0, n, n)
  one[999:1003, 999:1003] <- 1
  b1 <- ffp_block_sums(one, 5)
  expect_equal(b1[201, 201], 25)
  expect_equal(sum(b1), 25)
  # edge blocks hold 3 x 5 real cells (2 padded rows): the sum is kept on upsampling
  up <- ffp_block_upsample(b, dim(m), 5)
  expect_equal(dim(up), dim(m))
  expect_equal(sum(up), sum(m))
  expect_equal(sum(up[1:3, 1:3]), sum(m[1:3, 1:3]))         # corner block: 3 x 3 real cells
  expect_equal(up[1001, 1001], mean(m[999:1003, 999:1003])) # piecewise constant
  expect_error(ffp_block_sums(matrix(0, 100, 100), 5), "blocks")
})

test_that("c5 sums: write / read round trip, old full-grid files still readable, months add up", {
  g <- ffp_grid(c(-50, 50, -50, 50), 1)
  set.seed(2)
  f <- lapply(1:4, function(i) matrix(stats::runif(101 * 101), 101, 101))
  p <- withr::local_tempfile(fileext = ".rds")
  ffp_sums_write(list(`a|valid|day|own` = f[[1]] + f[[2]]), list(`a|valid|day|own` = 2L), p, block = 5)
  r <- ffp_sums_read(p)
  expect_equal(sum(r$sums[[1]]), sum(f[[1]] + f[[2]]), tolerance = 1e-6)
  expect_equal(r$sums[[1]], ffp_block_upsample(ffp_block_sums(f[[1]] + f[[2]], 5), c(101, 101), 5), tolerance = 1e-6)
  # block sums computed in the chunk (dim_full given) = coarsened by the writer
  p2 <- withr::local_tempfile(fileext = ".rds")
  ffp_sums_write(list(`a|valid|day|own` = ffp_block_sums(f[[1]] + f[[2]], 5)), list(`a|valid|day|own` = 2L), p2,
                 block = 5, dim_full = c(101, 101))
  expect_equal(ffp_sums_read(p2)$sums, r$sums)
  # full-grid file (spec 08 chunks)
  p3 <- withr::local_tempfile(fileext = ".rds")
  ffp_sums_write(list(`a|valid|day|own` = f[[3]]), list(`a|valid|day|own` = 1L), p3)
  expect_equal(ffp_sums_read(p3)$sums[[1]], f[[3]], tolerance = 1e-6)
  # annual = sum of months = sum of half-months (block sums are linear)
  whole <- ffp_block_sums(f[[1]] + f[[2]] + f[[3]] + f[[4]], 5)
  halves <- ffp_block_sums(f[[1]] + f[[2]], 5) + ffp_block_sums(f[[3]] + f[[4]], 5)
  expect_equal(halves, whole)
})

test_that("chunk_parts: whole month or halves", {
  w <- chunk_parts("2024-02")
  expect_equal(w$part, "2024-02")
  expect_equal(w$end, utc("2024-03-01"))
  h <- chunk_parts("2024-02", split = TRUE)
  expect_equal(h$part, c("2024-02a", "2024-02b"))
  expect_equal(h$start, utc(c("2024-02-01", "2024-02-16")))
  expect_equal(h$end, utc(c("2024-02-16", "2024-03-01")))
})

# Synthetic 08a input object: 2 months, cases pblhmod, era5, pblhmod_clim, measured
fake_inp <- function() {
  t <- utc("2024-01-01") + 1800 * (0:(48 * 60 - 1)) # 2024-01-01 to 2024-02-29
  n <- length(t)
  tab <- dplyr::tibble(timestamp_start_utc = t, zm_m = 20, umean = 3, L = -100, sigmav = 0.5, ustar = 0.4,
                       wind_dir = 180, day = (as.integer(format(t, "%H")) %in% 8:16),
                       h_pblhmod = 500, h_era5 = 600, h_pblhmod_clim = 550, h_measured = 700)
  v <- function(p) rep_len(p, n)
  sub <- list()
  for (s in c("etc", "none")) for (u in c("valid", "valid_qc", "valid_qc_used", "valid_qc_h")) {
    sub[[paste(u, s, "pblhmod", sep = "__")]] <- v(c(TRUE, TRUE, TRUE, FALSE))
    sub[[paste(u, s, "era5", sep = "__")]] <- v(c(TRUE, TRUE, FALSE, TRUE))
    sub[[paste(u, s, "pblhmod_clim", sep = "__")]] <- v(c(TRUE, FALSE, TRUE, TRUE))
    sub[[paste(u, s, "measured", sep = "__")]] <- v(c(TRUE, FALSE, FALSE, FALSE))
  }
  list(site_id = "XX-Tst", tab = tab, subsets = dplyr::as_tibble(sub))
}

test_that("footprint_task: rows of the part, own and common groups restricted to the case's valid rows", {
  inp <- fake_inp()
  pp <- chunk_parts("2024-01")
  tk <- footprint_task(inp, "pblhmod_clim", pp$part, pp$start, pp$end, sets = c("etc", "none"))
  expect_equal(nrow(tk$rows), 31 * 48 * 3 / 4)
  expect_true(all(format(tk$rows$timestamp_start_utc, "%Y-%m") == "2024-01"))
  expect_true(all(tk$rows$h == 550))
  own <- tk$grp[["none|valid|day|own"]] | tk$grp[["none|valid|night|own"]]
  pe <- tk$grp[["none|valid|day|common_pe"]] | tk$grp[["none|valid|night|common_pe"]]
  pme <- tk$grp[["none|valid|day|common_pme"]] | tk$grp[["none|valid|night|common_pme"]]
  expect_true(all(own))
  # common_pe = pblhmod & era5 & clim valid: pattern position 1 only; common_pme: also measured
  expect_equal(sum(pe), 31 * 48 / 4)
  expect_equal(pme, pe)
  # measured: no common_pe group; era5_bilinear (not in the table): no task
  tm <- footprint_task(inp, "measured", pp$part, pp$start, pp$end, sets = "none")
  expect_false(any(grepl("common_pe", names(tm$grp))))
  expect_true(any(grepl("common_pme", names(tm$grp))))
  expect_null(footprint_task(inp, "pblhmod", "2025-01", utc("2025-01-01"), utc("2025-02-01")))
  expect_equal(nrow(footprint_task(inp, "pblhmod", pp$part, pp$start, pp$end, sets = "none", limit = 5)$rows), 5)
})

test_that("make_chunk_list: estimates, months split above the target, empty months skipped", {
  inp <- fake_inp()
  cl <- make_chunk_list(inp, 2024, c("pblhmod", "era5"), sets = "none", sec_per_hh = 2, cpus = 6, target_min = 120)
  expect_equal(nrow(cl), 4) # 2 cases x 2 months (no data from March)
  expect_equal(cl$n_halfhours[cl$case == "pblhmod" & cl$part == "2024-01"], 31 * 48 * 3 / 4)
  expect_equal(cl$est_min[1], round(cl$n_halfhours[1] * 2 / 6 / 60, 1))
  expect_equal(cl$cpu_hours[1], round(cl$n_halfhours[1] * 2 / 3600, 2))
  cl2 <- make_chunk_list(inp, 2024, "pblhmod", sets = "none", sec_per_hh = 60, cpus = 6, target_min = 120)
  expect_equal(cl2$part, c("2024-01a", "2024-01b", "2024-02a", "2024-02b"))
  expect_equal(sum(cl2$n_halfhours), sum(cl$n_halfhours[cl$case == "pblhmod"]))
  expect_equal(nrow(make_chunk_list(inp, 2023, "pblhmod", sets = "none")), 0)
})

test_that("footprint_h_cases: main cases everywhere, supplement at extra sites, aerosol measured on its comparable set", {
  t <- utc("2022-09-01 00:00:00") + 1800 * (0:3)
  std <- function(v, q = 0L, tt = t, dt = 1800) dplyr::tibble(timestamp_start_utc = tt, timestamp_end_utc = tt + dt,
                                                             pblh_m = v, qc_flag = q)
  era <- dplyr::tibble(timestamp_start_utc = utc("2022-09-01 00:00:00") + 3600 * (0:2), pblh_m = c(100, 200, 300))
  ref <- list(pblhmod_gamma001 = std(c(1, 2, 3, 4)), pblhmod_gamma0025 = std(c(5, 6, 7, 8)))
  x <- footprint_h_cases(t, std(c(400, 500, NA, 600)), era, era, pblhmod_clim = std(c(410, 510, 610, NA)),
                         reference = ref, extra = FALSE)
  expect_equal(names(x), c("pblhmod", "pblhmod_clim", "era5"))
  x2 <- footprint_h_cases(t, std(c(400, 500, NA, 600)), era, era, pblhmod_clim = std(c(410, 510, 610, NA)),
                          reference = ref, extra = TRUE, oper = std(c(9, 9, 9, 9)))
  expect_equal(names(x2), c("pblhmod", "pblhmod_clim", "era5", "era5_bilinear", "pblhmod_x0.8", "pblhmod_x1.2",
                            "pblhmod_gamma001", "pblhmod_gamma0025", "oper"))
  expect_equal(x2$pblhmod_gamma0025, c(5, 6, 7, 8))
  # aerosol source, hourly interval means: both half-hours of an hour, only on the comparable set, no bound
  hourly <- std(c(1000, 1200), q = c(0L, 1L), tt = utc("2022-09-01 00:00:00") + 3600 * (0:1), dt = 3600)
  xa <- footprint_h_cases(t, std(c(400, 500, NA, 600)), era, era, measured = hourly, measured_type = "aerosol",
                          comparable = c(TRUE, FALSE, TRUE, TRUE), extra = FALSE)
  expect_equal(xa$measured, c(1000, NA, 1200, 1200))
  expect_false("bound" %in% names(xa))
  expect_error(footprint_h_cases(t, std(1:4), era, era, measured = hourly, measured_type = "aerosol"), "comparable")
  # turbulence source: exact (qc 0/2) and bound (qc 3), matched by interval
  xt <- footprint_h_cases(t, std(1:4), era, era, measured = std(c(500, NA, 800, 900), c(0L, 1L, 3L, 2L)), extra = FALSE)
  expect_equal(xt$measured, c(500, NA, NA, 900))
  expect_equal(xt$bound, c(NA, NA, 800, NA))
})

test_that("measured_primary_source and footprint_extra_sites read the config files", {
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  expect_equal(measured_primary_source(ms, "FR-Gri"), list(source_id = "sirta_mlh_chm15k", type = "aerosol"))
  expect_equal(measured_primary_source(ms, "FI-Hyy"), list(source_id = "cloudnet_eps_mlh", type = "turbulence"))
  expect_equal(measured_primary_source(ms, "CH-Dav")$type, "aerosol")
  expect_null(measured_primary_source(ms, "SE-Nor"))
  expect_null(measured_primary_source(ms, "BE-Bra"))
  ex <- footprint_extra_sites(here::here("config", "footprint_extra_sites.txt"))
  expect_equal(ex, c("BE-Bra", "FI-Hyy", "FR-Gri", "CH-Dav", "SE-Svb", "IT-Cp2"))
  expect_true(all(ex %in% read_sites(here::here("config", "sites.csv"))$site_id))
})

test_that("run_footprint_chunk: block size and workers do not change the result (skipped without the Kljun code)", {
  ffp_dir <- here::here("legacy")
  skip_if_not(file.exists(file.path(ffp_dir, "calc_footprint_FFP_climatology.R")), "Kljun FFP code not available")
  inp <- fake_inp()
  inp$tab$L <- rep_len(c(-100, -30, 200, 500), nrow(inp$tab))
  inp$tab$wind_dir <- rep_len(c(180, 270, 90), nrow(inp$tab))
  pp <- chunk_parts("2024-01")
  tk <- footprint_task(inp, "pblhmod", pp$part, pp$start, pp$end, sets = "none", limit = 6)
  tk$hash <- "test"
  g <- ffp_grid(c(-200, 200, -200, 200), 1) # 401 + 4 = 405 = 81 x 5
  mask <- matrix(FALSE, 401, 401)
  mask[150:250, 150:250] <- TRUE
  attr(mask, "area_m2") <- sum(mask)
  run <- function(dir, ...) {
    dir.create(dir)
    suppressWarnings(run_footprint_chunk(tk, ffp_dir, g, mask, c(0.5, 0.7, 0.8, 0.9), dir, ...))
    list(hh = readr::read_csv(file.path(dir, "XX-Tst_pblhmod_2024-01_hh.csv.gz"), show_col_types = FALSE),
         sums = ffp_sums_read(file.path(dir, "XX-Tst_pblhmod_2024-01_sums.rds")),
         clim = readRDS(file.path(dir, "XX-Tst_pblhmod_2024-01_clim.rds")))
  }
  tmp <- withr::local_tempdir()
  a <- run(file.path(tmp, "a"), block_rows = 48L)
  b <- run(file.path(tmp, "b"), block_rows = 2L)
  expect_equal(a$hh, b$hh)
  expect_equal(a$sums, b$sums, tolerance = 1e-6)
  expect_equal(names(a$sums$n), names(tk$grp)[vapply(tk$grp, any, logical(1))])
  expect_true(all(vapply(a$clim, function(c0) c0$n, 1) > 0))
  if (nzchar(Sys.getenv("PBLH_TEST_PARALLEL"))) { # two worker processes: only where several processes are allowed
    w <- run(file.path(tmp, "w"), block_rows = 2L, workers = 2L)
    expect_equal(w$sums, b$sums)
    expect_equal(w$hh, b$hh)
  }
})

test_that("ffp_region_pack / unpack are lossless; index vectors (spec 08 files) pass unchanged", {
  idx <- c(1L, 7L, 8L, 9L, 4004001L)
  p <- ffp_region_pack(idx, 2001L * 2001L)
  expect_equal(ffp_region_unpack(p), idx)
  expect_equal(ffp_region_unpack(ffp_region_pack(integer(), 10L)), integer())
  expect_equal(ffp_region_unpack(idx), idx)
})
