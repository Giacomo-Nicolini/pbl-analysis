# Campbell Scientific SkyVUE PRO (CS135) ceilometer with the built-in mixing-layer height (spec 14 Part A; DE-Tha, TU
# Dresden). Secondary source: shown and tabulated, never pooled with the E-PROFILE sites. Measured data are not
# released; nothing from the files goes into tracked files except variable names and formats.
#
# Manual: Campbell Scientific, SkyVUE PRO (CS135) LIDAR Ceilometer, Product manual, Revision 04/2026 (cited as
# `CampbellSkyVUEPRO2026` below; section and page numbers of that revision):
# - Sect. 6.3.6 "MESSAGE 006 (profile, sky condition, mixing layer heights)", pp. 87-91: seven lines per message.
#   Line 1 SOH "CS" ID(1) OS(3) N(3) STX; line 2 S(1) WA(1) tr(3) h1-h4(5 each) flags(12); line 3 sky condition;
#   line 4 scale(5) res(2) n(4) energy(3) lt(3) ti(2) bl(4) pulse(4) rate(2) sum(3); line 5 mh1 q1 mh2 q2 mh3 q3 (5
#   characters each; "/////" = no layer; all "/////" if the MLH option is not active); line 6 the profile, n x 5
#   hex characters; line 7 ETX CRC-16(4) EOT.
# - Sect. 2.1.3 "Backscatter profile reporting", p. 6: 5 hex characters per gate = signed 20-bit two's complement
#   integer; values > 2^19 - 1 are negative (subtract 2^20); times 1e-8 gives the attenuated backscatter coefficient
#   [m-1 sr-1]; "the values are scaled by the Attenuated_SCALE parameter" (line 4 `scale`, %, default 100); the last
#   two gates are always zero; the profile is not corrected for tilt.
# - Sect. 2.1.4 "Mixing layer height", pp. 6-7: gradient method after the KNMI operational algorithm (de Haij et al.
#   2007, KNMI WR 2007-01); quality factor 1-3, "with 3 indicating the highest confidence".
# - Sect. 6.2 "Checksums used in SkyVue PRO messages", p. 68: CRC-16 (CCITT polynomial 0x1021, initial 0xFFFF,
#   output XOR 0xFFFF) over all characters after SOH up to and including ETX; C code in the manual.
# Site facts (M. Hehn, TU Dresden, 2026-11-11): 113 m from the flux tower at bearing 275 deg, base 0.5 m above
# ground, not tilted; profiles from 2025-02-06, MLH from 2025-06-25; time UTC, timestamp = end of the interval.
# Daily files E1227_<ddmmyy>_Midnight.{csv,txt,nc}; the .nc files are empty (ignored).

skyvue_missing <- "/////"

#' Parse SkyVUE CSV numbers (decimal comma, `/////` missing)
#'
#' @param x Character values of the SkyVUE CSV export (e.g. `"2802,00000425904"`, `"/////"`) [character].
#' @return Numeric; `NA` for `/////`, empty strings and `NA`.
skyvue_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c(skyvue_missing, "")] <- NA_character_
  suppressWarnings(as.numeric(sub(",", ".", x, fixed = TRUE)))
}

#' CRC-16 of SkyVUE messages, vectorised over records
#'
#' `CampbellSkyVUEPRO2026`, Sect. 6.2: CRC-16 with polynomial 0x1021, initial value 0xFFFF and the output XOR'd with
#' 0xFFFF, over the characters after SOH up to and including ETX. Table-driven form of the manual's bitwise C code
#' (identical result).
#'
#' @param bytes Integer matrix, one record per row, the bytes after SOH up to and including ETX; records shorter than
#'   the matrix are padded at the end and given in `len` [-].
#' @param len Number of valid bytes per row [-]; default all columns.
#' @return Integer CRC per record (0-65535).
skyvue_crc16 <- function(bytes, len = rep(ncol(bytes), nrow(bytes))) {
  tab <- vapply(0:255, function(i) {
    crc <- bitwShiftL(i, 8L)
    for (j in 1:8) {
      crc <- if (bitwAnd(crc, 0x8000L) != 0L) bitwXor(bitwAnd(bitwShiftL(crc, 1L), 0xFFFFL), 0x1021L) else
        bitwAnd(bitwShiftL(crc, 1L), 0xFFFFL)
    }
    crc
  }, integer(1))
  crc <- rep(0xFFFFL, nrow(bytes))
  for (j in seq_len(ncol(bytes))) {
    use <- len >= j
    idx <- bitwXor(bitwShiftR(crc[use], 8L), bytes[use, j])
    crc[use] <- bitwXor(bitwAnd(bitwShiftL(crc[use], 8L), 0xFFFFL), tab[idx + 1L])
  }
  bitwXor(crc, 0xFFFFL)
}

#' Decode SkyVUE hex profiles
#'
#' `CampbellSkyVUEPRO2026`, Sect. 2.1.3: each gate is 5 hex characters, a signed 20-bit two's complement integer
#' (values above 2^19 - 1 minus 2^20); times 1e-8 = attenuated backscatter [m-1 sr-1]. The scale parameter (line 4,
#' %) is undone as `value * 100 / scale`: at the default 100 % (all DE-Tha records checked so far) the factor is 1;
#' the direction for other values is our reading of "scaled by the Attenuated_SCALE parameter" (flagged, spec 14).
#'
#' @param hex Profile lines (line 6), one per record [character].
#' @param scale Scale parameter of each record [%].
#' @param gates Gates to decode (1-based) [-]; default all.
#' @return Numeric matrix, records x gates, attenuated backscatter [m-1 sr-1].
skyvue_decode_profile <- function(hex, scale = rep(100, length(hex)), gates = NULL) {
  n <- nchar(hex)
  if (any(n %% 5 != 0)) cli::cli_abort("SkyVUE profile: line length not a multiple of 5 hex characters.")
  if (is.null(gates)) gates <- seq_len(max(n) / 5)
  out <- matrix(NA_real_, length(hex), length(gates))
  for (k in seq_along(gates)) {
    g <- gates[k]
    s <- substr(hex, 5 * g - 4, 5 * g)
    v <- strtoi(s, 16L)
    v <- ifelse(!is.na(v) & v > 524287L, v - 1048576L, v) # two's complement, 20 bits
    out[, k] <- v * 1e-8 * 100 / scale
  }
  out
}

# Split a field string by single spaces, keeping fixed field counts.
skyvue_fields <- function(x, n, what) {
  f <- strsplit(trimws(x), " +")
  bad <- lengths(f) != n
  if (any(bad)) cli::cli_abort("SkyVUE message: {sum(bad)} {what} line{?s} without {n} fields (first record {which(bad)[1]}).")
  matrix(unlist(f), ncol = n, byrow = TRUE)
}

#' Read a SkyVUE PRO message-006 text file (one day)
#'
#' Splits the file into messages at SOH/EOT, checks the message number (006), the seven-line layout
#' (`CampbellSkyVUEPRO2026`, Sect. 6.3.6) and the CRC-16 of every record (Sect. 6.2). The file carries no time stamp:
#' records are paired with the CSV export by order ([read_skyvue_day()]).
#'
#' @param path Path of `E1227_<ddmmyy>_Midnight.txt` [character].
#' @param keep_profile Keep the hex profile lines (`profile_hex`, about 10 kB per record) [logical].
#' @return Tibble, one row per record in file order: `detection_status` (0-6, `/` as `NA`), `warning_alarm` (`0`, `W`,
#'   `A`), `window_tr_pct` [%], `h1_m`..`h4_m` [m, line 2 heights], `status_flags` (12 hex characters), `scale_pct`
#'   [%], `resolution_m` [m], `n_gates` [-], `laser_energy_pct` [%], `laser_temp_c` [degC], `tilt_deg` [deg],
#'   `background_mv` [mV], `pulses_k` [x 1000], `sample_rate_mhz` [MHz], `bs_sum` [-, line 4 sum], `mlh1_m`..`mlh3_m`
#'   [m above the instrument; `NA` = `/////`], `mlh1_q`..`mlh3_q` [quality 1-3, 3 = highest confidence; 0 where no
#'   layer], `crc_ok` [logical], optionally `profile_hex`. Attribute `skyvue_txt_meta`: `file`, `n_records`,
#'   `n_crc_fail`, `message_id`, `record_bytes` (distinct lengths).
read_skyvue_txt <- function(path, keep_profile = FALSE) {
  if (!file.exists(path)) cli::cli_abort("SkyVUE file {.path {path}} not found.")
  b <- readBin(path, "raw", n = file.size(path))
  soh <- which(b == as.raw(1L))
  etx <- which(b == as.raw(3L))
  eot <- which(b == as.raw(4L))
  if (length(soh) == 0) cli::cli_abort("{.path {path}}: no SOH character (not a SkyVUE message file).")
  if (length(etx) != length(soh) || length(eot) != length(soh) || any(etx <= soh) || any(eot != etx + 5L)) {
    cli::cli_abort("{.path {path}}: SOH / ETX / EOT do not form {length(soh)} complete messages (ETX {length(etx)}, EOT {length(eot)}).")
  }
  if (length(soh) > 1 && any(soh[-1] < eot[-length(eot)])) cli::cli_abort("{.path {path}}: overlapping messages.")
  n_rec <- length(soh)
  len <- etx - soh # bytes after SOH up to and including ETX
  ml <- max(len)
  mat <- matrix(0L, n_rec, ml)
  for (j in seq_len(ml)) {
    use <- len >= j
    mat[use, j] <- as.integer(b[soh[use] + j])
  }
  crc_calc <- skyvue_crc16(mat, len)
  crc_txt <- vapply(seq_len(n_rec), function(i) rawToChar(b[(etx[i] + 1L):(etx[i] + 4L)]), "")
  crc_ok <- !is.na(strtoi(crc_txt, 16L)) & strtoi(crc_txt, 16L) == crc_calc
  # message text without SOH ... ETX, split into lines
  txt <- vapply(seq_len(n_rec), function(i) rawToChar(b[(soh[i] + 1L):(etx[i] - 1L)]), "")
  rm(b, mat)
  lines <- strsplit(txt, "\r\n", fixed = TRUE)
  nl <- lengths(lines)
  # Lines: 1 header (ends with STX), 2 status, 3 sky condition, 4 profile header, 5 MLH, 6 profile
  if (any(nl != 6)) cli::cli_abort("{.path {path}}: {sum(nl != 6)} message{?s} without the 6 text lines of message 006 (first record {which(nl != 6)[1]}).")
  L <- matrix(unlist(lines), ncol = 6, byrow = TRUE)
  head1 <- sub("\002$", "", L[, 1])
  if (any(!grepl("^CS.{4}006$", head1))) {
    cli::cli_abort("{.path {path}}: not message 006 ({.val {unique(head1[!grepl('^CS.{4}006$', head1)])[1]}}).")
  }
  s2 <- skyvue_fields(L[, 2], 7, "status")
  s4 <- skyvue_fields(L[, 4], 10, "profile header")
  s5 <- skyvue_fields(L[, 5], 6, "mixing layer")
  num <- function(x) skyvue_num(x)
  ds <- substr(s2[, 1], 1, 1)
  out <- dplyr::tibble(
    detection_status = suppressWarnings(as.integer(ifelse(ds == "/", NA, ds))),
    warning_alarm = substr(s2[, 1], 2, 2), window_tr_pct = num(s2[, 2]),
    h1_m = num(s2[, 3]), h2_m = num(s2[, 4]), h3_m = num(s2[, 5]), h4_m = num(s2[, 6]), status_flags = s2[, 7],
    scale_pct = num(s4[, 1]), resolution_m = num(s4[, 2]), n_gates = as.integer(num(s4[, 3])),
    laser_energy_pct = num(s4[, 4]), laser_temp_c = num(s4[, 5]), tilt_deg = num(s4[, 6]),
    background_mv = num(s4[, 7]), pulses_k = num(s4[, 8]), sample_rate_mhz = num(s4[, 9]), bs_sum = num(s4[, 10]),
    mlh1_m = num(s5[, 1]), mlh1_q = as.integer(num(s5[, 2])), mlh2_m = num(s5[, 3]), mlh2_q = as.integer(num(s5[, 4])),
    mlh3_m = num(s5[, 5]), mlh3_q = as.integer(num(s5[, 6])),
    crc_ok = crc_ok
  )
  bad_len <- nchar(L[, 6]) != 5L * out$n_gates
  if (any(bad_len)) cli::cli_abort("{.path {path}}: {sum(bad_len)} profile line{?s} not n_gates x 5 hex characters.")
  if (keep_profile) out$profile_hex <- L[, 6] else attr(out, "profile_head_hex") <- substr(L[, 6], 1, 100)
  attr(out, "skyvue_txt_meta") <- list(file = basename(path), n_records = n_rec, n_crc_fail = sum(!crc_ok),
                                       message_id = unique(head1), record_bytes = sort(unique(len)))
  out
}

#' Read a SkyVUE PRO CSV export (one day)
#'
#' The CSV uses decimal commas and `/////` for missing values; the backscatter-profile field mixes decimal commas
#' with comma separators (`3,65E-06,9,7E-06,...`) and is not parsed here (spec 14), except for the first gates by
#' [skyvue_csv_profile_head()] for the decode check.
#'
#' @param path Path of `E1227_<ddmmyy>_Midnight.csv` [character].
#' @return Tibble in file order: `time_end_utc` [POSIXct UTC, end of the record interval], `cloud_base1_m` ..
#'   `cloud_base4_m`, `vertical_visibility_m`, `highest_signal_m` [m], `ceil_status`, `mlh1_m`..`mlh3_m` [m],
#'   `mlh1_q`..`mlh3_q` [1-3], `ceil_resolution_m` [m], `profile_csv` (the raw profile field, first 600 characters).
read_skyvue_csv <- function(path) {
  if (!file.exists(path)) cli::cli_abort("SkyVUE file {.path {path}} not found.")
  d <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), na = character(),
                       progress = FALSE, name_repair = "minimal")
  col <- function(pattern) {
    k <- grep(pattern, names(d))
    if (length(k) != 1) cli::cli_abort("{.path {path}}: {length(k)} column{?s} match {.val {pattern}}.")
    d[[k]]
  }
  t <- as.POSIXct(col("^UTC Timestamp"), format = "%d.%m.%Y %H:%M:%S", tz = "UTC")
  if (anyNA(t)) cli::cli_abort("{.path {path}}: {sum(is.na(t))} time stamp{?s} not in dd.mm.yyyy HH:MM:SS.")
  dplyr::tibble(
    time_end_utc = t,
    cloud_base1_m = skyvue_num(col("^Cloud Base 1\\(")), cloud_base2_m = skyvue_num(col("^Cloud Base 2\\(")),
    cloud_base3_m = skyvue_num(col("^Cloud Base 3\\(")), cloud_base4_m = skyvue_num(col("^Cloud Base 4\\(")),
    vertical_visibility_m = skyvue_num(col("^Vertical Visibility\\(")),
    highest_signal_m = skyvue_num(col("^Highest Signal\\(")), ceil_status = col("^Ceil Status\\("),
    mlh1_m = skyvue_num(col("^Mixing Layer 1\\(")), mlh2_m = skyvue_num(col("^Mixing Layer 2\\(")),
    mlh3_m = skyvue_num(col("^Mixing Layer 3\\(")),
    mlh1_q = as.integer(skyvue_num(col("^Confidence 1\\("))), mlh2_q = as.integer(skyvue_num(col("^Confidence 2\\("))),
    mlh3_q = as.integer(skyvue_num(col("^Confidence 3\\("))),
    ceil_resolution_m = skyvue_num(col("^Ceil Resolution\\(")),
    profile_csv = substr(col("^Backscatter Profile\\("), 1, 600)
  )
}

#' First gates of the SkyVUE CSV profile field, where the parse is unambiguous
#'
#' Values are written as `m[,f]E±ee` (decimal comma) or `0`, separated by commas, so `0` followed by `,6E-08` cannot
#' be told from `0,6E-08` (checked on the DE-Tha sample of 2026-07-15). A non-zero value in E notation never has a
#' leading mantissa digit 0, so the first `n` values are taken while every token matches `-?[1-9]...E±ee`; from the
#' first zero (or ambiguous) token on, the record is `NA`.
#'
#' @param x Raw profile field (CSV) [character].
#' @param n Number of gates [-].
#' @return Numeric matrix, records x `n` [CSV unit, equal to the decoded value at scale 100].
skyvue_csv_profile_head <- function(x, n = 20L) {
  out <- matrix(NA_real_, length(x), n)
  re <- "^(-?[1-9][0-9]*(,[0-9]+)?E[+-][0-9]+)(,|$)"
  rest <- x
  ok <- rep(TRUE, length(x))
  for (k in seq_len(n)) {
    m <- regmatches(rest, regexec(re, rest))
    hit <- ok & lengths(m) > 0
    v <- rep(NA_character_, length(x))
    v[hit] <- vapply(m[hit], `[`, "", 2)
    out[hit, k] <- as.numeric(sub(",", ".", v[hit], fixed = TRUE))
    rest[hit] <- substr(rest[hit], nchar(vapply(m[hit], `[`, "", 1)) + 1L, nchar(rest[hit]))
    ok <- hit
  }
  out
}

#' Read one SkyVUE day (CSV + TXT) and check the pairing
#'
#' Spec 14 A1: time stamps and the three mixing layers with their confidence from the CSV; the profile header and
#' (optionally) the profile from the TXT, paired by record order. Checks, each aborting with the first failing record:
#' equal record counts; for every record with a valid CRC, the CSV mixing layers and confidences equal the TXT line 5
#' (`/////` = `NA`, confidence 0 when no layer), and the CSV cloud base 1 (rounded to 1 m) equals TXT `h1` where the
#' detection status is 1-4; the decoded profile (first `n_check` gates) equals the CSV values where the CSV parse is
#' unambiguous (relative difference <= 1 %, CSV values carry 3-4 significant digits). Records with a CRC mismatch
#' are kept with `crc_ok = FALSE` and their mixing layers set to `NA` (counted in the attribute).
#'
#' @param csv,txt Paths of the day's `.csv` and `.txt` files [character].
#' @param keep_profile Keep the hex profile lines [logical].
#' @param n_check Profile gates compared with the CSV [-].
#' @return Tibble of [read_skyvue_txt()] with `time_end_utc` [POSIXct UTC] first. Attribute `skyvue_meta`: `csv`,
#'   `txt`, `n_records`, `n_crc_fail`, `n_profile_records_checked`, `max_rel_diff_profile` [-], `scale_values` [%],
#'   `interval_s_median` [s].
read_skyvue_day <- function(csv, txt, keep_profile = FALSE, n_check = 20L) {
  a <- read_skyvue_csv(csv)
  b <- read_skyvue_txt(txt, keep_profile = keep_profile)
  if (nrow(a) != nrow(b)) cli::cli_abort("SkyVUE pairing: {basename(csv)} has {nrow(a)} records, {basename(txt)} {nrow(b)}.")
  first_bad <- function(ok, what) {
    if (!all(ok)) {
      i <- which(!ok)[1]
      cli::cli_abort("SkyVUE pairing ({basename(csv)} / {basename(txt)}): {what} differ on {sum(!ok)} record{?s}; first record {i} ({format(a$time_end_utc[i], '%Y-%m-%d %H:%M:%S')} UTC).")
    }
  }
  same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) < 0.5)
  ok_crc <- b$crc_ok
  for (k in 1:3) {
    first_bad(!ok_crc | same(a[[paste0("mlh", k, "_m")]], b[[paste0("mlh", k, "_m")]]), paste("mixing layer", k, "heights"))
    qa <- dplyr::coalesce(a[[paste0("mlh", k, "_q")]], 0L)
    qb <- dplyr::coalesce(b[[paste0("mlh", k, "_q")]], 0L)
    first_bad(!ok_crc | qa == qb, paste("mixing layer", k, "confidences"))
  }
  cl <- b$detection_status %in% 1:4
  first_bad(!ok_crc | !cl | same(round(a$cloud_base1_m), b$h1_m), "cloud base 1 (detection status 1-4)")
  hex <- if (keep_profile) b$profile_hex else attr(b, "profile_head_hex")
  dec <- skyvue_decode_profile(substr(hex, 1, 5 * n_check), b$scale_pct, gates = seq_len(n_check))
  csvp <- skyvue_csv_profile_head(a$profile_csv, n_check)
  cmp <- !is.na(csvp) & ok_crc
  rel <- abs(dec - csvp) / pmax(abs(csvp), 1e-12)
  rel[!cmp] <- NA_real_
  bad_rec <- apply(rel > 0.01, 1, function(r) any(r, na.rm = TRUE))
  first_bad(!bad_rec, paste0("decoded profile and CSV profile (first ", n_check, " gates)"))
  for (k in 1:3) {
    b[[paste0("mlh", k, "_m")]][!ok_crc] <- NA_real_
    b[[paste0("mlh", k, "_q")]][!ok_crc] <- NA_integer_
  }
  out <- dplyr::bind_cols(dplyr::tibble(time_end_utc = a$time_end_utc), b)
  attr(out, "skyvue_meta") <- list(
    csv = basename(csv), txt = basename(txt), n_records = nrow(out), n_crc_fail = sum(!ok_crc),
    n_profile_records_checked = sum(rowSums(cmp) > 0),
    max_rel_diff_profile = if (any(!is.na(rel))) max(rel, na.rm = TRUE) else NA_real_,
    scale_values = sort(unique(b$scale_pct)),
    interval_s_median = if (nrow(out) > 1) stats::median(diff(as.numeric(out$time_end_utc))) else NA_real_
  )
  out
}

#' Find the SkyVUE day file pairs under a folder
#'
#' @param dir Folder searched recursively [character].
#' @return Tibble: `day` [Date, from `<ddmmyy>`], `csv`, `txt` [paths]; days with only one of the two files are
#'   reported and dropped; the empty `.nc` files are ignored.
find_skyvue_days <- function(dir) {
  f <- list.files(dir, pattern = "_Midnight\\.(csv|txt)$", recursive = TRUE, full.names = TRUE)
  base <- sub("\\.(csv|txt)$", "", f)
  ext <- sub("^.*\\.", "", f)
  k <- unique(base)
  csv <- f[ext == "csv"][match(k, base[ext == "csv"])]
  txt <- f[ext == "txt"][match(k, base[ext == "txt"])]
  d <- dplyr::tibble(key = k, csv = csv, txt = txt)
  half <- is.na(d$csv) | is.na(d$txt)
  if (any(half)) cli::cli_alert_warning("SkyVUE: {sum(half)} day{?s} with only one of .csv / .txt dropped: {.file {basename(d$key[half])}}.")
  d <- d[!half, ]
  dd <- regmatches(basename(d$key), regexpr("_[0-9]{6}_", basename(d$key)))
  d$day <- as.Date(gsub("_", "", dd), format = "%d%m%y")
  d <- d[order(d$day), c("day", "csv", "txt")]
  if (anyDuplicated(d$day)) cli::cli_abort("SkyVUE: more than one file pair for {.val {format(d$day[duplicated(d$day)])}}.")
  d
}

#' Aggregate SkyVUE records to half-hours (all three layers, no layer chosen)
#'
#' Spec 14 A2: half-hours [t, t + 30 min) UTC; a record belongs to the half-hour containing its time stamp minus 1 s
#' (the stamp is the end of the ~8.6 s record interval, so a record ending exactly on the half-hour belongs to the
#' previous one). Records with a CRC mismatch are not used for the layers. `pblh_m` stays `NA` until G decides a
#' layer rule.
#'
#' @param r Output of [read_skyvue_day()] (one or more days, bound).
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `n_samples` [-, records with a valid
#'   CRC], and per layer k = 1-3 `x_layer<k>_m` [m above the instrument, median of the records with a layer],
#'   `x_layer<k>_conf` [median confidence 1-3 of those records], `x_layer<k>_n` [-, records with layer k],
#'   `x_layer<k>_share_conf3` [0-1, share of those records with confidence 3]; `x_n_records` [-, all records],
#'   `x_cloud_share` [0-1, records with detection status 1-4], `x_window_tr_pct` [%, mean window transmission].
skyvue_halfhour <- function(r) {
  t <- as.numeric(r$time_end_utc) - 1
  hh <- floor(t / 1800) * 1800
  grid <- seq(min(hh), max(hh), by = 1800)
  k <- match(hh, grid)
  n <- length(grid)
  ok <- r$crc_ok %in% TRUE
  agg <- function(x, keep, f) {
    out <- rep(NA_real_, n)
    if (any(keep)) {
      a <- tapply(x[keep], k[keep], f)
      out[as.integer(names(a))] <- a
    }
    out
  }
  st <- as.POSIXct(grid, origin = "1970-01-01", tz = "UTC")
  out <- dplyr::tibble(timestamp_start_utc = st, timestamp_end_utc = st + 1800, n_samples = tabulate(k[ok], nbins = n))
  for (j in 1:3) {
    h <- r[[paste0("mlh", j, "_m")]]
    q <- r[[paste0("mlh", j, "_q")]]
    use <- ok & !is.na(h)
    out[[paste0("x_layer", j, "_m")]] <- agg(h, use, stats::median)
    out[[paste0("x_layer", j, "_conf")]] <- agg(q, use & !is.na(q), stats::median)
    out[[paste0("x_layer", j, "_n")]] <- tabulate(k[use], nbins = n)
    out[[paste0("x_layer", j, "_share_conf3")]] <- agg(q %in% 3L, use, mean)
  }
  out$x_n_records <- tabulate(k, nbins = n)
  out$x_cloud_share <- agg(r$detection_status %in% 1:4, ok, mean)
  out$x_window_tr_pct <- agg(r$window_tr_pct, ok & !is.na(r$window_tr_pct), mean)
  out
}

#' Standardize the SkyVUE PRO built-in MLH (`tud_skyvue_mlh`, secondary)
#'
#' Half-hourly interval values from [skyvue_halfhour()], `method_class = aerosol`, all three layers kept as `x_`
#' columns, `pblh_m` = `NA` and `qc_flag` = `NA` on every row until G decides a layer rule (spec 14 A2). Exclusion
#' periods of the ICOS site and of the source are applied and counted. Internal use only.
#'
#' @param h Output of [skyvue_halfhour()].
#' @param site_id ICOS station [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h ICOS site UTC offset [h] (for the exclusion periods).
#' @param measured_exclusions Output of [read_measured_exclusions()]; `NULL` = none (reported).
#' @param days_meta Tibble of per-day reading results (`day`, `n_records`, `n_crc_fail`, `max_rel_diff_profile`)
#'   [-], summarised in `processing_history`.
#' @return List: `data` (standard table), `meta` (sidecar), `n_excluded`, `n_excluded_measured` [-].
standardize_skyvue_mlh <- function(h, site_id, exclusions, utc_offset_h, measured_exclusions = NULL, days_meta = NULL) {
  src <- "tud_skyvue_mlh"
  d <- dplyr::tibble(timestamp_start_utc = h$timestamp_start_utc, timestamp_end_utc = h$timestamp_end_utc,
                     pblh_m = NA_real_, source_id = src, qc_flag = NA_integer_, uncertainty_m = NA_real_,
                     n_samples = h$n_samples)
  d <- dplyr::bind_cols(d, h[grep("^x_", names(h), value = TRUE)])
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, src)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL
  mx <- apply_measured_exclusions(d, site_id, src, measured_exclusions)
  d <- mx$data
  dm <- if (is.null(days_meta)) "" else paste0(
    nrow(days_meta), " day file pairs ", format(min(days_meta$day)), " to ", format(max(days_meta$day)), ", ",
    sum(days_meta$n_records), " records, ", sum(days_meta$n_crc_fail), " with a CRC mismatch (layers not used), ",
    "profile decode vs CSV: max relative difference ", signif(max(days_meta$max_rel_diff_profile, na.rm = TRUE), 2),
    "; ")
  meta <- standard_meta(
    site_id, src, "aerosol",
    instrument = paste("Campbell Scientific SkyVUE PRO (CS135) lidar ceilometer, TU Dresden, 113 m from the DE-Tha",
                       "flux tower at bearing 275 deg, base 0.5 m above ground, not tilted (site, 2026-11-11);",
                       "profiles from 2025-02-06, MLH option from 2025-06-25"),
    retrieval_algorithm = paste("Built-in mixing-layer height of the SkyVUE PRO (gradient method after the KNMI",
                                "operational algorithm, de Haij et al. 2007; CampbellSkyVUEPRO2026, Sect. 2.1.4), up to",
                                "three layers with a quality factor 1-3 (3 = highest confidence); computed by the",
                                "instrument, not retrieved in this project. No layer chosen: pblh_m is NA until G",
                                "decides a layer rule (spec 14)"),
    height_reference = "above the instrument (0.5 m above ground; heights used as is)",
    instrument_altitude_m = "unknown",
    time_reference = "interval_mean",
    original_resolution = "about 8.6 s records (message 006, 2048 gates of 5 m); time stamp = end of the record interval, UTC",
    original_qc = "instrument quality factor 1-3 per layer; CRC-16 per record (records with a mismatch not used)",
    qc_flag_meaning = list(`NA` = "no layer rule yet: pblh_m is NA on every row (spec 14)"),
    processing_history = paste0(
      dm, "read by read_skyvue_day() (CSV time stamps and layers, TXT message 006 paired by record order and ",
      "checked); half-hourly median of each layer by skyvue_halfhour(); rows in exclusion periods dropped (",
      ex$n_excluded, "); ", mx$history, "; standardized by standardize_skyvue_mlh() (spec 14)."),
    contact = "TU Dresden (DE-Tha)",
    licence = "internal use only: measured data are not released (spec 14; provider TU Dresden)",
    citation = "CampbellSkyVUEPRO2026 (manual, Sect. 2.1.4, 6.3.6); de Haij et al. 2007 (KNMI WR 2007-01)",
    data_owner = "TU Dresden",
    processor = "instrument firmware (built-in MLH option); aggregated in this project",
    known_limitations = "secondary source: manufacturer algorithm, never pooled with the E-PROFILE sites"
  )
  list(data = d, meta = meta, n_excluded = ex$n_excluded, n_excluded_measured = mx$n_excluded)
}
