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

#' Read a SkyVUE PRO message file (message 006 or 004; one day)
#'
#' Splits the file into messages at SOH/EOT, checks the message number, the line layout (`CampbellSkyVUEPRO2026`,
#' Sect. 6.3.6) and the CRC-16 of every record (Sect. 6.2). The file carries no time stamp: records are paired with
#' the CSV export by order ([read_skyvue_day()]).
#'
#' Message 006: seven lines (header, status, sky condition, profile header, mixing layers, profile, ETX line).
#' Message 004 (DE-Tha 2025-02-07 to 2025-03-11, before the MLH option; close-out 2026-11-19): six lines, the same
#' without the mixing-layer line; read with the 006 field layouts of the other lines (the field counts are checked,
#' so a different layout aborts). Its mixing layers are `NA` ("no MLH option"), not an error. The message type is
#' read per record.
#'
#' Incomplete messages (fragments) are dropped, counted in `n_truncated_dropped` and reported with a warning (byte
#' position and the number of complete records before them); the rest of the day is kept. A fragment is a SOH whose
#' message has neither ETX nor EOT before the next SOH (mid-file, probably a logger restart; fix-14a June 2026:
#' DE-Tha 2026-06-21, -22, -28) or before the end of the file (truncated last message, close-out 2026-11-19; there a
#' dangling ETX is allowed). Any other layout (an ETX without its EOT mid-file, two ETX or EOT in one message, EOT not
#' 5 bytes after ETX) aborts. An abort after the fragment scan carries the partial meta (`n_truncated_dropped`) in
#' the condition field `skyvue_partial`.
#'
#' @param path Path of `E1227_<ddmmyy>_Midnight.txt` [character].
#' @param keep_profile Keep the hex profile lines (`profile_hex`, about 10 kB per record) [logical].
#' @return Tibble, one row per record in file order: `message_id` (`"006"`, `"004"`), `detection_status` (0-6, `/` as
#'   `NA`), `warning_alarm` (`0`, `W`,
#'   `A`), `window_tr_pct` [%], `h1_m`..`h4_m` [m, line 2 heights], `status_flags` (12 hex characters), `scale_pct`
#'   [%], `resolution_m` [m], `n_gates` [-], `laser_energy_pct` [%], `laser_temp_c` [degC], `tilt_deg` [deg],
#'   `background_mv` [mV], `pulses_k` [x 1000], `sample_rate_mhz` [MHz], `bs_sum` [-, line 4 sum], `mlh1_m`..`mlh3_m`
#'   [m above the instrument; `NA` = `/////` and in message 004], `mlh1_q`..`mlh3_q` [quality 1-3, 3 = highest
#'   confidence; 0 where no layer, `NA` in message 004], `crc_ok` [logical], optionally `profile_hex`. Attribute
#'   `skyvue_txt_meta`: `file`, `n_records`, `n_crc_fail`, `message_id` (distinct), `n_004`, `n_truncated_dropped`
#'   (fragments dropped, last or mid-file), `fragment_after` (complete records before each fragment) [-],
#'   `fragment_byte` (byte position of its SOH) [-], `record_bytes` (distinct lengths).
read_skyvue_txt <- function(path, keep_profile = FALSE) {
  if (!file.exists(path)) cli::cli_abort("SkyVUE file {.path {path}} not found.")
  b <- readBin(path, "raw", n = file.size(path))
  soh <- which(b == as.raw(1L))
  etx <- which(b == as.raw(3L))
  eot <- which(b == as.raw(4L))
  if (length(soh) == 0) cli::cli_abort("{.path {path}}: no SOH character (not a SkyVUE message file).")
  # every ETX / EOT belongs to the message of the last SOH before it (its window up to the next SOH)
  ie <- findInterval(etx, soh)
  io <- findInterval(eot, soh)
  if (any(ie == 0L) || any(io == 0L)) cli::cli_abort("{.path {path}}: ETX / EOT before the first SOH.")
  ne <- tabulate(ie, length(soh))
  no <- tabulate(io, length(soh))
  is_last <- seq_along(soh) == length(soh)
  frag <- no == 0L & (ne == 0L | (is_last & ne == 1L))
  complete <- ne == 1L & no == 1L
  if (any(!frag & !complete)) {
    k <- which(!frag & !complete)[1]
    cli::cli_abort("{.path {path}}: SOH / ETX / EOT do not form complete messages: {sum(!frag & !complete)} message{?s} with ETX / EOT counts other than 1 / 1 (first at byte {soh[k]}: ETX {ne[k]}, EOT {no[k]}).")
  }
  n_trunc <- sum(frag)
  frag_after <- cumsum(complete)[frag]
  frag_byte <- soh[frag]
  if (n_trunc > 0) {
    for (j in seq_len(n_trunc)) {
      if (frag_byte[j] == soh[length(soh)] && frag_after[j] == sum(complete)) {
        cli::cli_alert_warning("{basename(path)}: truncated last message (byte {frag_byte[j]} of {length(b)}, no ETX / EOT) dropped, after record {frag_after[j]}.")
      } else {
        cli::cli_alert_warning("{basename(path)}: incomplete message (byte {frag_byte[j]}, no ETX / EOT before the next SOH) dropped, after record {frag_after[j]}.")
      }
    }
    etx <- etx[complete[ie]]
    eot <- eot[complete[io]]
    soh <- soh[complete]
    cli::cli_alert_info("{basename(path)}: {length(soh)} complete message{?s} kept.")
    if (length(soh) == 0) cli::cli_abort("{.path {path}}: no complete message.")
  }
  if (any(eot != etx + 5L)) cli::cli_abort("{.path {path}}: SOH / ETX / EOT do not form complete messages (EOT not 5 bytes after ETX).")
  partial <- list(n_truncated_dropped = n_trunc,
                  fragment_pos = if (n_trunc > 0) paste0("after record ", frag_after, " (byte ", frag_byte, ")", collapse = "; ") else NA_character_)
  skyvue_with_partial({
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
    # Message 006 text lines: 1 header (ends with STX), 2 status, 3 sky condition, 4 profile header, 5 MLH, 6 profile;
    # message 004: the same without line 5 (5 text lines)
    head1 <- sub("\002$", "", vapply(lines, `[`, "", 1))
    mid <- ifelse(grepl("^CS.{4}006$", head1), "006", ifelse(grepl("^CS.{4}004$", head1), "004", NA_character_))
    if (anyNA(mid)) {
      cli::cli_abort("{.path {path}}: not message 006 or 004 ({.val {unique(head1[is.na(mid)])[1]}}).")
    }
    n_exp <- ifelse(mid == "006", 6L, 5L)
    if (any(nl != n_exp)) {
      i <- which(nl != n_exp)[1]
      cli::cli_abort("{.path {path}}: {sum(nl != n_exp)} message{?s} without the {n_exp[i]} text lines of message {mid[i]} (first record {i}).")
    }
    ln <- function(k006, k004) vapply(seq_along(lines), function(i) lines[[i]][if (mid[i] == "006") k006 else k004], "")
    L <- cbind(head1, ln(2, 2), ln(3, 3), ln(4, 4), ifelse(mid == "006", ln(5, 1), NA_character_), ln(6, 5))
    rm(lines)
    is6 <- mid == "006"
    s2 <- skyvue_fields(L[, 2], 7, "status")
    s4 <- skyvue_fields(L[, 4], 10, "profile header")
    s5 <- matrix(NA_character_, nrow(L), 6)
    if (any(is6)) s5[is6, ] <- skyvue_fields(L[is6, 5], 6, "mixing layer")
    num <- function(x) skyvue_num(x)
    ds <- substr(s2[, 1], 1, 1)
    out <- dplyr::tibble(
      message_id = mid,
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
                                         message_id = sort(unique(mid)), n_004 = sum(!is6),
                                         n_truncated_dropped = n_trunc, fragment_after = frag_after,
                                         fragment_byte = frag_byte, fragment_pos = partial$fragment_pos,
                                         record_bytes = sort(unique(len)))
    out
  }, partial)
}

# Evaluate `expr`; on an error, add the partial meta known so far (a list or an environment) to the condition as
# field `skyvue_partial` (values already there, from an inner call, take precedence) and signal it again, class
# kept. Used so that 14a can count CSV header and fragment cases on days that fail (fix-14a June 2026, item 3).
skyvue_with_partial <- function(expr, partial) {
  tryCatch(expr, error = function(e) {
    p <- if (is.environment(partial)) as.list(partial) else partial
    e$skyvue_partial <- utils::modifyList(p, if (is.null(e$skyvue_partial)) list() else e$skyvue_partial)
    stop(e)
  })
}

#' Read a SkyVUE PRO CSV export (one day)
#'
#' The CSV uses decimal commas and `/////` for missing values; the backscatter-profile field mixes decimal commas
#' with comma separators (`3,65E-06,9,7E-06,...`) and is not parsed here (spec 14), except for the first gates by
#' [skyvue_csv_profile_head()] for the decode check.
#'
#' Two CSV layouts, detected from the header (fix-14c item 9, DE-Tha cluster data): `"full"` (26 columns, from about
#' 2025-06-26: status, profile, sky condition, mixing layers and confidences) and `"short"` (7 columns until about
#' 2025-06-24, before the MLH option: `UTC Timestamp`, `Cloud Base 1..4`, `Vertical Visibility`, `Highest Signal`).
#'
#' A CSV without a header line (the first line is a record; DE-Tha 2026-06-22, 2026-06-29; close-out 2026-11-19) is
#' detected from its first line (a `dd.mm.yyyy HH:MM:SS` time stamp) and read with the header of the layout that
#' belongs to the day's message type (`message_id`: 006 = full, 004 = short, [skyvue_csv_header()]); the column count
#' must match that header, else the reader aborts. Attribute `header_missing` = `TRUE` (logged by the caller).
#'
#' @param path Path of `E1227_<ddmmyy>_Midnight.csv` [character].
#' @param message_id Message type of the day's TXT file (`"006"`, `"004"`), needed only for a CSV without a header
#'   [character].
#' @return Tibble in file order: `time_end_utc` [POSIXct UTC, end of the record interval], `cloud_base1_m` ..
#'   `cloud_base4_m`, `vertical_visibility_m`, `highest_signal_m` [m], `ceil_status`, `mlh1_m`..`mlh3_m` [m],
#'   `mlh1_q`..`mlh3_q` [1-3], `ceil_resolution_m` [m], `profile_csv` (the raw profile field, first 600 characters);
#'   in the short layout the columns it lacks are `NA`. Attributes `layout` (`"full"` or `"short"`) and
#'   `header_missing` [logical].
read_skyvue_csv <- function(path, message_id = NULL) {
  if (!file.exists(path)) cli::cli_abort("SkyVUE file {.path {path}} not found.")
  no_header <- skyvue_csv_header_missing(path)
  if (no_header) {
    if (is.null(message_id) || length(message_id) != 1 || !message_id %in% c("006", "004")) {
      cli::cli_abort("{.path {path}}: no header line, and no single message type (006 / 004) to take the layout from.")
    }
    hdr <- skyvue_csv_header(if (message_id == "006") "full" else "short")
    d <- readr::read_csv(path, col_names = FALSE, col_types = readr::cols(.default = readr::col_character()),
                         na = character(), progress = FALSE)
    if (ncol(d) != length(hdr)) {
      cli::cli_abort("{.path {path}}: no header line; {ncol(d)} columns, the {if (message_id == '006') 'full' else 'short'} layout of message {message_id} has {length(hdr)}.")
    }
    names(d) <- hdr
  } else {
    d <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), na = character(),
                         progress = FALSE, name_repair = "minimal")
  }
  col <- function(pattern) {
    k <- grep(pattern, names(d))
    if (length(k) != 1) cli::cli_abort("{.path {path}}: {length(k)} column{?s} match {.val {pattern}}.")
    d[[k]]
  }
  short_cols <- c("^UTC Timestamp", "^Cloud Base 1\\(", "^Cloud Base 2\\(", "^Cloud Base 3\\(", "^Cloud Base 4\\(",
                  "^Vertical Visibility\\(", "^Highest Signal\\(")
  has <- function(p) any(grepl(p, names(d)))
  layout <- if (any(grepl("^Mixing Layer 1\\(", names(d)))) {
    "full"
  } else if (ncol(d) == 7 && all(vapply(short_cols, has, logical(1)))) {
    "short"
  } else {
    cli::cli_abort("{.path {path}}: unknown SkyVUE CSV layout ({ncol(d)} columns: {.val {names(d)}}).")
  }
  t <- as.POSIXct(col("^UTC Timestamp"), format = "%d.%m.%Y %H:%M:%S", tz = "UTC")
  if (anyNA(t)) cli::cli_abort("{.path {path}}: {sum(is.na(t))} time stamp{?s} not in dd.mm.yyyy HH:MM:SS.")
  full <- layout == "full"
  na_num <- rep(NA_real_, nrow(d))
  out <- dplyr::tibble(
    time_end_utc = t,
    cloud_base1_m = skyvue_num(col("^Cloud Base 1\\(")), cloud_base2_m = skyvue_num(col("^Cloud Base 2\\(")),
    cloud_base3_m = skyvue_num(col("^Cloud Base 3\\(")), cloud_base4_m = skyvue_num(col("^Cloud Base 4\\(")),
    vertical_visibility_m = skyvue_num(col("^Vertical Visibility\\(")),
    highest_signal_m = skyvue_num(col("^Highest Signal\\(")),
    ceil_status = if (full) col("^Ceil Status\\(") else rep(NA_character_, nrow(d)),
    mlh1_m = if (full) skyvue_num(col("^Mixing Layer 1\\(")) else na_num,
    mlh2_m = if (full) skyvue_num(col("^Mixing Layer 2\\(")) else na_num,
    mlh3_m = if (full) skyvue_num(col("^Mixing Layer 3\\(")) else na_num,
    mlh1_q = if (full) as.integer(skyvue_num(col("^Confidence 1\\("))) else NA_integer_,
    mlh2_q = if (full) as.integer(skyvue_num(col("^Confidence 2\\("))) else NA_integer_,
    mlh3_q = if (full) as.integer(skyvue_num(col("^Confidence 3\\("))) else NA_integer_,
    ceil_resolution_m = if (full) skyvue_num(col("^Ceil Resolution\\(")) else na_num,
    profile_csv = if (full) substr(col("^Backscatter Profile\\("), 1, 600) else rep(NA_character_, nrow(d))
  )
  attr(out, "layout") <- layout
  attr(out, "header_missing") <- no_header
  out
}

#' Does a SkyVUE CSV lack its header line?
#'
#' @param path Path of the CSV [character].
#' @return `TRUE` if the first line is a record (`dd.mm.yyyy HH:MM:SS,` time stamp), else `FALSE` [logical].
skyvue_csv_header_missing <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  length(first) == 1 && grepl("^[0-9]{2}\\.[0-9]{2}\\.[0-9]{4} [0-9]{2}:[0-9]{2}:[0-9]{2},", first)
}

#' Header lines of the two SkyVUE CSV layouts (column names only)
#'
#' Full layout: the 26 columns of the DE-Tha export with the MLH option (sample of 2026-07-15); short layout: the 7
#' columns before it (`UTC Timestamp`, `Cloud Base 1..4`, `Vertical Visibility`, `Highest Signal`; units as in the full
#' header). Used for a CSV written without its header line ([read_skyvue_csv()]).
#'
#' @param layout `"full"` or `"short"`.
#' @return Column names [character].
skyvue_csv_header <- function(layout = c("full", "short")) {
  layout <- match.arg(layout)
  cb <- paste0("Cloud Base ", 1:4, "( Meters )")
  if (layout == "short") {
    return(c("UTC Timestamp", cb, "Vertical Visibility( Meters )", "Highest Signal( Meters )"))
  }
  c("UTC Timestamp", cb, "Vertical Visibility( Meters )", "Ceil Status( Status )",
    "Backscatter Profile( (100000*srad*km)^-1 )", paste0("Sky Condition ", 1:5, "( Oktas )"),
    paste0("Sky Condition Height ", 1:5, "( Distance )"), paste0("Mixing Layer ", 1:3, "( Meters )"),
    paste0("Confidence ", 1:3, "( Confidence 1-3 )"), "Highest Signal( Meters )", "Ceil Resolution( Meters )")
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

#' Per-record pairing checks of a SkyVUE CSV / TXT pair
#'
#' The checks of [read_skyvue_day()] on rows already paired by order (`a[i]` with `b[i]`); records with a CRC
#' mismatch pass every check (not compared).
#'
#' @param a Rows of [read_skyvue_csv()].
#' @param b Rows of [read_skyvue_txt()], same number as `a`.
#' @param hex Hex profile lines of `b` (at least the first `n_check` gates) [character].
#' @param layout CSV layout, `"full"` or `"short"`.
#' @param n_check Profile gates compared with the CSV [-].
#' @return List: `checks` (named list of logical vectors, `TRUE` = record passes; names = check descriptions in the
#'   order they are reported), `ok` (all checks passed) [logical], `rel` (records x gates relative difference of the
#'   decoded profile and the CSV profile where compared, else `NA`; `NULL` in the short layout) [-].
skyvue_pair_checks <- function(a, b, hex, layout, n_check = 20L) {
  same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) < 0.5)
  free <- !b$crc_ok
  ck <- list()
  cl <- b$detection_status %in% 1:4
  ck[["cloud base 1 (detection status 1-4)"]] <- free | !cl | same(round(a$cloud_base1_m), b$h1_m)
  rel <- NULL
  if (layout == "full") {
    is6 <- b$message_id == "006" # message 004 has no mixing-layer line (no MLH option): not compared
    for (k in 1:3) {
      ck[[paste("mixing layer", k, "heights")]] <- free | !is6 |
        same(a[[paste0("mlh", k, "_m")]], b[[paste0("mlh", k, "_m")]])
      qa <- dplyr::coalesce(a[[paste0("mlh", k, "_q")]], 0L)
      qb <- dplyr::coalesce(b[[paste0("mlh", k, "_q")]], 0L)
      ck[[paste("mixing layer", k, "confidences")]] <- free | !is6 | qa == qb
    }
    dec <- skyvue_decode_profile(substr(hex, 1, 5 * n_check), b$scale_pct, gates = seq_len(n_check))
    csvp <- skyvue_csv_profile_head(a$profile_csv, n_check)
    cmp <- !is.na(csvp) & !free
    rel <- abs(dec - csvp) / pmax(abs(csvp), 1e-12)
    rel[!cmp] <- NA_real_
    bad_rec <- rowSums(rel > 0.01, na.rm = TRUE) > 0
    ck[[paste0("decoded profile and CSV profile (first ", n_check, " gates)")]] <- !bad_rec
  } else {
    vv <- b$detection_status %in% 5L
    ck[["vertical visibility (detection status 5)"]] <- free | !vv | same(round(a$vertical_visibility_m), b$h1_m)
  }
  list(checks = ck, ok = Reduce(`&`, ck), rel = rel)
}

#' Find the CSV rows of dropped TXT fragments
#'
#' fix-14a June 2026: when the CSV has one row more per dropped TXT fragment ([read_skyvue_txt()]), the CSV row of
#' each fragment is searched among the rows within `reach` of its expected position (the row after the complete
#' records before it): each candidate is dropped in turn and the records within `win` of the fragment are compared
#' with [skyvue_pair_checks()] (cloud base, mixing layers and confidences, profile head; short layout: cloud base and
#' vertical visibility). The candidate with the fewest mismatching records is taken; ties (records identical in every
#' compared field) go to the candidate nearest the expected position, then the earlier one (the time stamp is then
#' uncertain by the tied rows; counted in `n_tied`). Fragments are matched in file order.
#'
#' @param a,b,hex,layout,n_check As [skyvue_pair_checks()]; `a` has `nrow(b) + length(after)` rows.
#' @param after Complete TXT records before each fragment (`fragment_after` of [read_skyvue_txt()]) [-].
#' @param reach Candidate rows on each side of the expected one [-].
#' @param win TXT records on each side of the fragment compared [-].
#' @return List: `drop` (CSV rows to drop, one per fragment), `n_tied` (candidates with the best score) [-],
#'   `n_mismatch` (mismatching records in the window with the chosen row dropped) [-].
skyvue_match_fragments <- function(a, b, hex, layout, n_check, after, reach = 3L, win = 10L) {
  nf <- length(after)
  drop <- integer(nf)
  n_tied <- integer(nf)
  n_mm <- integer(nf)
  for (f in seq_len(nf)) {
    p <- after[f]
    e <- p + f # expected CSV row: p complete records and f - 1 earlier fragment rows before it
    lo <- max(1L, p - win + 1L, if (f > 1) after[f - 1] + 1L else 1L)
    hi <- min(nrow(b), p + win, if (f < nf) after[f + 1] else nrow(b))
    i <- if (hi >= lo) lo:hi else integer()
    lo_c <- max(1L, e - reach, if (f > 1) drop[f - 1] + 1L else 1L)
    hi_c <- min(nrow(a), e + reach, if (f < nf) after[f + 1] + f else nrow(a))
    cand <- lo_c:hi_c
    score <- vapply(cand, function(c0) {
      keep <- setdiff(seq_len(nrow(a)), c(drop[seq_len(f - 1)], c0))
      if (length(i) == 0) return(0)
      sum(!skyvue_pair_checks(a[keep[i], ], b[i, ], hex[i], layout, n_check)$ok)
    }, numeric(1))
    best <- cand[score == min(score)]
    drop[f] <- best[order(abs(best - e), best)][1]
    n_tied[f] <- length(best)
    n_mm[f] <- min(score)
  }
  list(drop = drop, n_tied = n_tied, n_mismatch = n_mm)
}

#' Read one SkyVUE day (CSV + TXT) and check the pairing
#'
#' Spec 14 A1: time stamps and the three mixing layers with their confidence from the CSV; the profile header and
#' (optionally) the profile from the TXT, paired by record order. Checks ([skyvue_pair_checks()]): equal record counts
#' (a mismatch aborts with condition class `skyvue_count_mismatch`, fields `n_csv`, `n_txt`; fix-14a item 10); for
#' every record with a valid CRC, the CSV mixing layers and confidences equal the TXT line 5 (`/////` = `NA`,
#' confidence 0 when no layer), and the CSV cloud base 1 (rounded to 1 m) equals TXT `h1` where the detection status
#' is 1-4; the decoded profile (first `n_check` gates) equals the CSV values where the CSV parse is unambiguous
#' (relative difference <= 1 %, CSV values carry 3-4 significant digits). Records with a CRC mismatch are kept with
#' `crc_ok = FALSE` and their mixing layers set to `NA` (counted in the attribute).
#'
#' Short CSV layout (before the MLH option, until about 2025-06-24; fix-14c item 9): records paired by order as well;
#' checked: equal counts, cloud base 1 (detection status 1-4) and the vertical visibility (detection status 5, where
#' the TXT `h1` carries it; `CampbellSkyVUEPRO2026` Sect. 6.3.6, line 2) against the TXT; mixing layers `NA` (none in
#' the CSV, the option was not active); the CSV-profile comparison is not available (no profile column) and is
#' reported as such.
#'
#' Close-out 2026-11-19 (DE-Tha days that failed in 14a): message 004 days are read (no mixing-layer line: layers `NA`,
#' `mlh_note` "no MLH option (message 004)", the MLH comparison skipped on those records); a CSV without a header line
#' is read with the layout of the day's message type ([read_skyvue_csv()], `csv_header`).
#'
#' fix-14a June 2026 (G 2026-11-20): TXT fragments (a truncated last message or an incomplete message mid-file) are
#' dropped by [read_skyvue_txt()]. If the CSV has as many records as the TXT, the records are paired by order (the
#' fragment has no CSV row); if it has one more per fragment, the CSV row of each fragment is found by
#' [skyvue_match_fragments()], dropped and its time logged (`csv_dropped_utc`); any other difference stays a count
#' mismatch. Isolated pairing mismatches are tolerated: if at most `max_mismatch` records and at most
#' `max_mismatch_share` of the day's records fail a check, those records are dropped from both files and logged
#' (`n_pair_mismatch_dropped`, `pair_mismatch_utc`, `pair_mismatch_check`); above that the day fails as before. Each
#' case is logged with a warning. Any error after the CSV header test carries the partial meta (`csv_header`,
#' `n_truncated_dropped`, `fragment_pos`) in the condition field `skyvue_partial`.
#'
#' @param csv,txt Paths of the day's `.csv` and `.txt` files [character].
#' @param keep_profile Keep the hex profile lines [logical].
#' @param n_check Profile gates compared with the CSV [-].
#' @param max_mismatch Largest number of records failing a pairing check that are dropped instead of failing the day
#'   [-] (G 2026-11-20: 5).
#' @param max_mismatch_share Largest share of the day's records that may be dropped that way [0-1] (G 2026-11-20:
#'   0.001).
#' @return Tibble of [read_skyvue_txt()] with `time_end_utc` [POSIXct UTC] first. Attribute `skyvue_meta`: `csv`,
#'   `txt`, `csv_layout` (`full` / `short`), `csv_header` (`present` / `missing: layout of message <id>`),
#'   `message_id`, `n_004` [-], `mlh_note`, `n_csv`, `n_txt` [-, records as read; TXT after dropping fragments],
#'   `n_truncated_dropped` [-, TXT fragments], `fragment_pos` (`after record k (byte b)`, `; `-separated; `NA` if
#'   none), `n_csv_dropped` [-], `csv_dropped_utc` (time stamps of the CSV rows dropped with fragments,
#'   `YYYY-MM-DD HH:MM:SS` UTC, space-separated; `NA` if none), `n_fragment_tied` [-, fragments whose CSV row was
#'   not unique], `n_pair_mismatch_dropped` [-], `pair_mismatch_utc` (as `csv_dropped_utc`), `pair_mismatch_check`
#'   (checks failed, `; `-separated), `n_records`, `n_crc_fail`, `profile_check` (`done` / `not available (short CSV
#'   layout)`), `n_profile_records_checked`, `max_rel_diff_profile` [-], `scale_values` [%], `interval_s_median` [s].
read_skyvue_day <- function(csv, txt, keep_profile = FALSE, n_check = 20L, max_mismatch = 5L,
                            max_mismatch_share = 0.001) {
  pm <- new.env()
  pm$csv_header <- if (file.exists(csv) && skyvue_csv_header_missing(csv)) "missing" else "present"
  pm$n_truncated_dropped <- NA_integer_
  pm$fragment_pos <- NA_character_
  skyvue_with_partial({
    b <- read_skyvue_txt(txt, keep_profile = keep_profile)
    tm <- attr(b, "skyvue_txt_meta")
    pm$n_truncated_dropped <- tm$n_truncated_dropped
    pm$fragment_pos <- tm$fragment_pos
    hex <- if (keep_profile) b$profile_hex else attr(b, "profile_head_hex")
    a <- read_skyvue_csv(csv, message_id = if (length(tm$message_id) == 1) tm$message_id else NULL)
    layout <- attr(a, "layout")
    if (isTRUE(attr(a, "header_missing"))) {
      pm$csv_header <- paste("missing: layout of message", tm$message_id)
      cli::cli_alert_warning("{basename(csv)}: no header line; read with the {layout} layout of message {tm$message_id}.")
    }
    n_csv0 <- nrow(a)
    n_txt0 <- nrow(b)
    # TXT fragments: CSV rows of the fragments found by matching (fix-14a June 2026); equal counts = paired by order
    n_csv_drop <- 0L
    csv_dropped_utc <- NA_character_
    n_tied <- 0L
    if (tm$n_truncated_dropped > 0 && nrow(a) == nrow(b) + tm$n_truncated_dropped) {
      mf <- skyvue_match_fragments(a, b, hex, layout, n_check, tm$fragment_after)
      n_csv_drop <- length(mf$drop)
      csv_dropped_utc <- paste(format(a$time_end_utc[mf$drop], "%Y-%m-%d %H:%M:%S"), collapse = " ")
      n_tied <- sum(mf$n_tied > 1)
      for (f in seq_along(mf$drop)) {
        cli::cli_alert_warning("{basename(csv)}: record {mf$drop[f]} ({format(a$time_end_utc[mf$drop[f]], '%Y-%m-%d %H:%M:%S')} UTC) dropped with the TXT fragment after record {tm$fragment_after[f]}{if (mf$n_tied[f] > 1) paste0(' (', mf$n_tied[f], ' rows matched equally; nearest the fragment taken)') else ''}.")
      }
      a <- a[-mf$drop, ]
    } else if (tm$n_truncated_dropped > 0 && nrow(a) == nrow(b)) {
      cli::cli_alert_info("{basename(csv)}: as many records as the complete TXT messages: paired by order (no CSV row for the fragment{?s}).")
    }
    if (nrow(a) != nrow(b)) {
      # fix-14a item 10: a count mismatch has its own condition class (skyvue_count_mismatch, fields n_csv and n_txt)
      # so that callers can skip the day instead of counting it as a failed check; no alignment is attempted
      cli::cli_abort("SkyVUE pairing: {basename(csv)} has {nrow(a)} records, {basename(txt)} {nrow(b)}.",
                     class = "skyvue_count_mismatch", n_csv = nrow(a), n_txt = nrow(b))
    }
    pc <- skyvue_pair_checks(a, b, hex, layout, n_check)
    bad <- !pc$ok
    n_bad <- sum(bad)
    mm_utc <- NA_character_
    mm_check <- NA_character_
    if (n_bad > 0) {
      failed <- names(pc$checks)[vapply(pc$checks, function(x) !all(x), logical(1))]
      if (n_bad > max_mismatch || n_bad > max_mismatch_share * nrow(b)) {
        w <- failed[1]
        i <- which(!pc$checks[[w]])[1]
        cli::cli_abort("SkyVUE pairing ({basename(csv)} / {basename(txt)}): {w} differ on {sum(!pc$checks[[w]])} record{?s}; first record {i} ({format(a$time_end_utc[i], '%Y-%m-%d %H:%M:%S')} UTC); {n_bad} record{?s} {?fails/fail} a check (tolerated: up to {max_mismatch} and {100 * max_mismatch_share} % of {nrow(b)}).")
      }
      # fix-14a June 2026: isolated mismatches dropped from both files and logged
      mm_utc <- paste(format(a$time_end_utc[bad], "%Y-%m-%d %H:%M:%S"), collapse = " ")
      mm_check <- paste(failed, collapse = "; ")
      cli::cli_alert_warning("{basename(csv)}: {n_bad} record{?s} failing a pairing check dropped (record{?s} {which(bad)}, {mm_utc} UTC; {mm_check}).")
      a <- a[!bad, ]
      b <- b[!bad, ]
      if (!is.null(pc$rel)) pc$rel <- pc$rel[!bad, , drop = FALSE]
    }
    ok_crc <- b$crc_ok
    if (layout == "full") {
      profile_check <- "done"
      n_prof <- sum(rowSums(!is.na(pc$rel)) > 0)
      max_rel <- if (any(!is.na(pc$rel))) max(pc$rel, na.rm = TRUE) else NA_real_
    } else {
      for (k in 1:3) { # no MLH column in the short layout (option not active)
        b[[paste0("mlh", k, "_m")]] <- NA_real_
        b[[paste0("mlh", k, "_q")]] <- NA_integer_
      }
      profile_check <- "not available (short CSV layout)"
      n_prof <- NA_integer_
      max_rel <- NA_real_
    }
    for (k in 1:3) {
      b[[paste0("mlh", k, "_m")]][!ok_crc] <- NA_real_
      b[[paste0("mlh", k, "_q")]][!ok_crc] <- NA_integer_
    }
    n_004 <- sum(b$message_id == "004")
    out <- dplyr::bind_cols(dplyr::tibble(time_end_utc = a$time_end_utc), b)
    attr(out, "skyvue_meta") <- list(
      csv = basename(csv), txt = basename(txt), csv_layout = layout, csv_header = pm$csv_header,
      message_id = paste(tm$message_id, collapse = " "), n_004 = n_004,
      mlh_note = if (n_004 == nrow(out)) "no MLH option (message 004)" else if (n_004 > 0)
        paste0("no MLH option on ", n_004, " records (message 004)") else NA_character_,
      n_csv = n_csv0, n_txt = n_txt0,
      n_truncated_dropped = tm$n_truncated_dropped, fragment_pos = tm$fragment_pos,
      n_csv_dropped = n_csv_drop, csv_dropped_utc = csv_dropped_utc, n_fragment_tied = n_tied,
      n_pair_mismatch_dropped = n_bad, pair_mismatch_utc = mm_utc, pair_mismatch_check = mm_check,
      n_records = nrow(out), n_crc_fail = sum(!ok_crc),
      profile_check = profile_check, n_profile_records_checked = n_prof, max_rel_diff_profile = max_rel,
      scale_values = sort(unique(b$scale_pct)),
      interval_s_median = if (nrow(out) > 1) stats::median(diff(as.numeric(out$time_end_utc))) else NA_real_
    )
    out
  }, pm)
}

#' Find the SkyVUE day file pairs under a folder
#'
#' @param dir Folder searched recursively [character].
#' @return Tibble: `day` [Date, from `<ddmmyy>`], `csv`, `txt` [paths, `NA` if missing], `skip_reason` (`NA` for a
#'   complete pair; `"no TXT file"` / `"no CSV file"` for a day with only one of the two: kept so that the caller can
#'   list it, fix-14c item 9, e.g. DE-Tha 2025-06-25); the empty `.nc` files are ignored.
find_skyvue_days <- function(dir) {
  f <- list.files(dir, pattern = "_Midnight\\.(csv|txt)$", recursive = TRUE, full.names = TRUE)
  base <- sub("\\.(csv|txt)$", "", f)
  ext <- sub("^.*\\.", "", f)
  k <- unique(base)
  csv <- f[ext == "csv"][match(k, base[ext == "csv"])]
  txt <- f[ext == "txt"][match(k, base[ext == "txt"])]
  d <- dplyr::tibble(key = k, csv = csv, txt = txt)
  d$skip_reason <- ifelse(is.na(d$txt), "no TXT file", ifelse(is.na(d$csv), "no CSV file", NA_character_))
  half <- !is.na(d$skip_reason)
  if (any(half)) cli::cli_alert_warning("SkyVUE: {sum(half)} day{?s} with only one of .csv / .txt (skipped, listed): {.file {basename(d$key[half])}}.")
  dd <- regmatches(basename(d$key), regexpr("_[0-9]{6}_", basename(d$key)))
  d$day <- as.Date(gsub("_", "", dd), format = "%d%m%y")
  d <- d[order(d$day), c("day", "csv", "txt", "skip_reason")]
  dup <- d$day[duplicated(d$day)]
  if (length(dup)) cli::cli_abort("SkyVUE: more than one file set for {.val {format(dup)}}.")
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
    "profile decode vs CSV: max relative difference ",
    if (any(!is.na(days_meta$max_rel_diff_profile))) signif(max(days_meta$max_rel_diff_profile, na.rm = TRUE), 2) else
      "not available", " (not available on short-layout days, before the MLH option); ")
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
