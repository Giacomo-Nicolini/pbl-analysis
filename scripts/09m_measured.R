# Spec 09 amendment 2: measured PBLH for one site, as run_site.R steps (sites with a primary measured source in
# config/measured_sites.csv). Runs, each in a fresh R process (output in this log), the scripts of the provider of
# the primary source (measured_site_scripts()):
#   --part=standard (step meas_std): standard files of the measured source
#     cloudnet 07c_cloudnet_mlh.R (FI-Hyy; files downloaded by 07b first), SIRTA-ReOBS 10a_sirta_reobs.R (FR-Gri),
#     MeteoSwiss 11b_meteoswiss_mlh.R (CH-Dav);
#   --part=compare (step meas_cmp, after 04): per-site comparison
#     cloudnet 07_measured_compare.R and 10_measured_compare.R, SIRTA-ReOBS 10_measured_compare.R; MeteoSwiss none
#     (compared in the multi-site scripts/11_gamma_evaluation.R, spec 11 D).
# The first failing script stops the step (exit status 1).
#
# Usage: Rscript scripts/09m_measured.R SITE --part=standard|compare
#   Profile: env var R_CONFIG_ACTIVE (default | cluster); the child processes inherit it.

for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
part <- sub("^--part=", "", grep("^--part=", args, value = TRUE))
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")
if (!isTRUE(part %in% c("standard", "compare"))) cli::cli_abort("Use {.code --part=standard} or {.code --part=compare}.")

cfg <- get_cfg()
provider <- measured_primary_provider(read_measured_sites(resolve_path(cfg$measured_sites_file)), site_id)
sc <- measured_site_scripts(provider)[[part]]
cli::cli_h1("Measured PBLH, {part}: {site_id} ({provider})")
if (length(sc) == 0) {
  cli::cli_inform("{provider}: no per-site {part} script (compared in {.file scripts/11_gamma_evaluation.R}).")
  quit(save = "no", status = 0L)
}
for (s in sc) {
  cli::cli_h2(s)
  r <- callr::rscript(here::here(s), cmdargs = site_id, stdout = "", stderr = "2>&1", fail_on_status = FALSE,
                      show = FALSE, wd = here::here())
  if (!identical(as.integer(r$status), 0L)) {
    cli::cli_alert_danger("{s} failed (exit status {r$status}).")
    quit(save = "no", status = 1L)
  }
}
cli::cli_alert_success("{site_id}: {paste(basename(sc), collapse = ', ')} ok.")
