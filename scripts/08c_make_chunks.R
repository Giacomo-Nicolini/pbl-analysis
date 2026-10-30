# Spec 09 A5: footprint chunk list for the cluster array, from the approved reference years only.
# Reads config/reference_years.csv (site_id, year, reason; committed after G's review of the proposal
# output/tables/08_reference_years.csv written by scripts/09_collect.R) and each site's full-period input table
# (scripts/08a_footprint_flags.R). One chunk per site x case x month of the year; a month whose estimated duration
# (half-hours x sec_per_ffp_halfhour / CPUs) exceeds chunk_target_min is split into halves (config.yml).
# Cases: every case of the input table; oper only at the extra sites (supplement) or with --with-oper.
# Output: <derived_dir>/footprint/chunk_list.csv (one row per array task) and
#         <output_dir>/tables/09_chunks_summary.csv (per site: chunks, half-hours, CPU-hours, longest chunk).
#
# Usage: Rscript scripts/08c_make_chunks.R [SITE ...] [--cpus=6] [--with-oper]
#   Without sites: every site of config/reference_years.csv. A requested site without an approved year aborts.

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
cfg <- get_cfg()
cpus <- as.integer(opt("cpus", "6"))
with_oper <- "--with-oper" %in% args
ry_file <- resolve_path(cfg$reference_years_file, must_exist = FALSE)
if (!file.exists(ry_file)) {
  cli::cli_abort(c("No approved reference years: {.path {ry_file}} not found.",
                   "i" = "Review output/tables/08_reference_years.csv (scripts/09_collect.R) and commit the approved list."))
}
ry <- readr::read_csv(ry_file, col_types = readr::cols(year = "i", .default = "c"), progress = FALSE)
need <- c("site_id", "year", "reason")
if (!all(need %in% names(ry))) cli::cli_abort("{.path {ry_file}} needs columns {.field {need}}.")
if (anyDuplicated(ry$site_id)) cli::cli_abort("{.path {ry_file}}: more than one year for {.val {ry$site_id[duplicated(ry$site_id)]}}.")
site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_ids) == 0) site_ids <- ry$site_id
miss <- setdiff(site_ids, ry$site_id[!is.na(ry$year)])
if (length(miss)) cli::cli_abort("No approved reference year for {.val {miss}} in {.path {ry_file}}.")
extra_sites <- footprint_extra_sites(resolve_path(cfg$footprint_extra_sites_file))
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)

lst <- list()
for (s in site_ids) {
  inp_file <- file.path(derived, "footprint", s, paste0(s, "_footprint_inputs.rds"))
  if (!file.exists(inp_file)) cli::cli_abort("{s}: {.path {inp_file}} not found (run step 08a).")
  inp <- readRDS(inp_file)
  cases <- sub("^h_", "", grep("^h_", names(inp$tab), value = TRUE))
  if (!(with_oper || s %in% extra_sites)) cases <- setdiff(cases, "oper")
  yr <- ry$year[ry$site_id == s]
  cl <- make_chunk_list(inp, yr, cases, sec_per_hh = cfg$sec_per_ffp_halfhour, cpus = cpus, target_min = cfg$chunk_target_min)
  if (nrow(cl) == 0) cli::cli_abort("{s}: no valid half-hour in {yr}.")
  lst[[s]] <- cl
  cli::cli_inform("{s} {yr}: {nrow(cl)} chunks ({sum(grepl('[ab]$', cl$part))} half-months), {sum(cl$n_halfhours)} FFP half-hours, {round(sum(cl$cpu_hours))} CPU-hours; cases {paste(unique(cl$case), collapse = ', ')}.")
}
cl_all <- dplyr::bind_rows(lst)
cl_all <- cl_all[order(-cl_all$est_min), ] # longest first: the array starts them early
out <- file.path(derived, "footprint", "chunk_list.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(cl_all, out)
summ <- dplyr::summarise(cl_all, .by = site_id, n_chunks = dplyr::n(), n_halfhours = sum(n_halfhours),
                         cpu_hours = round(sum(cpu_hours), 1), max_est_min = max(est_min))
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summ, file.path(tdir, "09_chunks_summary.csv"))
print(as.data.frame(summ), row.names = FALSE)
cli::cli_alert_success("{nrow(cl_all)} chunks, {round(sum(cl_all$cpu_hours))} CPU-hours ({cpus} CPUs per chunk, {cfg$sec_per_ffp_halfhour} s per half-hour): {.path {out}}")
