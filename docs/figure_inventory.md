# Output inventory: tables, figures, stored footprint groups

State of `main` + `spec-11-gamma-network`, 2026-11-11. Tables in `<output_dir>/tables/`, figures in
`<output_dir>/figures/` (png). `<site>` = station code, `<per>` = period suffix of a step-08 run (`_<year>` for the
reference year, `_<YYYY-MM>` for a test month). Figures are drawn from the same in-memory data as the table listed
next to them, not by reading the table back, except where noted (*reads tables*).

Recurring column groups:
- **metrics** (`pblh_metrics()`): `n, bias, mae, rmse, r_pearson, r_spearman, n_log, n_log_excluded,
  median_log_ratio, mad_log_ratio, frac_within_factor2, fb, nrmsd, p25_log_ratio, p75_log_ratio`. Long tables store
  them as `metric, value`.
- **pairing** (spec 13): `measured_provider_site, measured_distance_km, borderline`.
- **quantiles**: `p5, p25, p50, p75, p95` (or `median, p25, p75`), plus `n`.
- **field metrics** at r = 50/70/80/90 %: `area_<r>, xmax_<r>, edge_<r>`.

## Per site: steps 01–04, 11a

| Table | Grain (one row per) | Columns | Figures |
|---|---|---|---|
| `01_inputs_<site>` | input item | `site_id, item, value` | none |
| `01_check_inputs_<site>` | input variable | `site_id, variable, n_total, n_excluded, n_valid, pct_valid, n_out_of_range, min, max` | none |
| `02_cases_<site>` | PBLHmod case | `site_id, case`, u* filter (`n_ustar, ratio_only, cap_only, both, rejected, pct, interpolated, left_na`), `gt5km, gt10km`, L check (`checked, replaced, pct_L_from_L2, sign_inconsistent`), `pct_valid, gamma_source, gamma_median, n_constant_fallback` | none |
| `02_regimes_<site>` | case × regime (all, stable, unstable_diag, convective) | `site_id, case, regime, share`, quantiles | none |
| `03_variants_<site>`, `03_reproduction_<site>` | legacy variant × stratum × statistic | `variant, stratum, statistic, value` | `03_reproduction_`, `03_waterfall[_allrows]_`, `03_one_at_a_time[_allrows]_`, `03_diurnal_<site>` |
| `03_cleanset_vs_l2_<site>` | variable | `variable, both, na_cleanset_only, na_l2_only, median_abs_diff` | none |
| `04_metrics_<site>` | pair × stratum × metric (pairs 1, 2, 6 always; 3, 7, 8 when the case exists; 4, 4b, 5 ICOS operational with `--with-oper`; strata `all`, `regime=`, `season=`, `hour_utc=`) | `pair, stratum`, metrics (long) | `04_scatter_mod_era5_`, `04_logratio_by_hour_`, `04_logratio_by_cloud_<site>` (pair 1); `09_network_logratio*` (*reads tables*, via 09) |
| `04_diurnal_<site>` | season × UTC hour × variable | `season, hour_utc, variable, median, p25, p75, n` | `04_diurnal_by_season_<site>` |
| `04_H_metrics_<site>` | H pair × stratum × metric | `pair, stratum, metric, value` | `04_H_ec_vs_era5_<site>` |
| `04_gamma_850_700_<site>` | UTC hour × month | `hour_utc, month, median, p25, p75, n` | `04_gamma_850_700_<site>` |
| (none) | ERA5 cells | none | `04_era5_spread_<site>` |
| `11a_gamma_summary_<site>` | fit layer × season | `site_id, layer, season, n_days, n_valid, share_2_levels, share_implausible, share_too_few, gamma_p25, gamma_median, gamma_p75, r2_median, gamma_orog_median, ratio_to_0.01` | none (daily γ in `11_gamma_by_month_*`) |

## Measured MLH: steps 07, 07c, 10, 10b, 13a

| Table | Grain | Columns | Figures |
|---|---|---|---|
| `07_cloudnet_inventory` | site × product × instrument × month | `site, product, instrument_name, instrument_pid, month, n_files, size_mb, n_legacy, n_volatile, mean_coverage` | none |
| `07_cloudnet_files` | Cloudnet file | `uuid, filename, site, product, measurement_date, instrument_name, instrument_pid, size_bytes, checksum_sha256, version, pid, volatile, legacy, coverage, error_level, download_url` | none |
| `07c_profile_status_<site>` | season × UTC hour × profile status | `season, hour_utc, status, Freq` | none |
| `07c_qc_by_hour_season_{UTC,local}_<site>` | season × hour × qc | `season, hour, qc, Freq` | none |
| `07c_censoring_<site>` | season × UTC hour | `season, hour_utc, n_valid, n_censored, pct_censored` | none |
| `07c_precip_rule_sensitivity_<site>` | rule × season | `rule, season, pct_profiles_precip, n_qc0…n_qc3, n_na, median_mlh_m, pct_censored` | none |
| `07c_altitude_sensitivity_<site>` | altitude variant | `variant, z_min_m, n_halfhours, pct_qc_changed, pct_censoring_changed, pct_mlh_changed, median_mlh_diff_m, pct_profiles_status_changed` | none |
| `07_metrics_<site>` | pair × stratum × metric | `pair, stratum, metric, value, comparison_class, distance_class`, pairing | `07_scatter_`, `07_logratio_by_hour_<site>` |
| `07_diurnal_<site>` | season × UTC hour × variable | `season, hour_utc, variable, median, n` | `07_diurnal_by_season_<site>` |
| `07_censoring_<site>` | season × UTC hour | `season, hour_utc, n_usable, n_below, n_lower_bound, pct_below, pct_lower_bound, frac_mod_below_zmin.<pair>` | `07_censoring_<site>` |
| `07_lower_bound_<site>` | pair × regime × season | `pair, regime, season, n, frac_mod_below_bound, median_ratio` | none |
| `07_shallow_deep_<site>` | pair × season × regime | `pair, season, regime, n_qc1, frac_le_zmin_qc1, n_qc3, frac_below_bound_qc3, median_ratio_qc3` | none |
| `07_h_timing_<site>` | season × H series | `site_id, season, series, max_H, t_max_min`, `{morning,afternoon}_{0,20}pct[_utc]`, `lag_era5_min, t_max_utc` | none |
| `10_metrics_<site>` (aerosol) | pair × stratum set × stratum × metric (sets `window_convective`, `convective`; pair 6 also `all_regimes`; strata `all`, `season=`, `hour_utc=`) | `site_id, pair, stratum_set, stratum, metric, value, no_data_reason, comparison_class, distance_class, hour_note, fb_convention`, pairing | `10_scatter_<site>` (window convective, others faded) |
| `10_compression_<site>` | pair × stratum set | `site_id, pair, stratum_set, n, n_days, r_log`, `slope_sym[_lo,_hi], intercept_sym`, `ba_mean, ba_slope[_lo,_hi], ba_intercept`, `ts_slope, ts_slope_reverse, ts_intercept, n_boot`, terciles (`ref_q33, ref_q67, n_T1–3, mlr_T1–3, f2_T1–3`), `values_note, tercile_note, reference[, comparison_class]`, pairing | `10_scatter_`, `10_bland_altman_<site>` |
| `10_bland_altman_bins_<site>` | pair × stratum set × bin | `pair, stratum_set, bin, a_lo, a_hi, gm_m, n, median_d, q25_d, q75_d, frac_within_factor2`, pairing | `10_bland_altman_<site>` |
| `10_diurnal_<site>` (aerosol) | season × UTC hour × variable | `site_id, season, hour_utc, variable, median, p25, p75, n, x_utc`, pairing | `10_diurnal_by_season_<site>` |
| `10_censored_null_<site>` (FI-Hyy) | pair × regime set × tail | `site_id, pair, regime_set, tail, n, n_days, obs_share, obs_lo, obs_hi, sim_share, sim_lo, sim_hi, sim_bound, n_exact_pool, n_fallback_regime, n_fallback_all, tail_meaning, sim_note`, pairing | none |
| `10b_cases_window_<site>` | case × hour set | `case, hours`, metrics, `reference` | none |
| `10b_growth_clear_days_<site>` | season × AST hour × variable | `season, ast_hour, variable, n, median, q25, q75, n_clear_days` | `10b_growth_clear_days_<site>` |
| `10b_h0_<site>` | day | `day_id, sunrise_utc, h0_m, h0_method, h0_lag_halfhours, n_convective, first_h_bg_m, flag_no_h0, flag_bg_stopped_gap, season, bg_run, era5_blh_m, n_win, mlr_day` | none |
| `10b_h0_summary_<site>` | season | `season, n_days, h0_p05…h0_p95, frac_h0_B1, frac_h0_lt_50m, era5_p50, era5_p25, era5_p75, rho_h0_vs_day_logratio` | none |
| `10b_stopped_days_<site>` | item | `item, value, table` | none |
| `10c_fluxes_nee_columns` | station | `site_id, l2_folder, fluxes_file, fluxnet_file, has_NEE, has_NEE_DATA_FLAG, has_NEE_VUT_REF_QC, fluxes_nee_like, fluxnet_nee_ref, n_rows, n_NEE, n_NEE_DATA_FLAG, n_NEE_DATA_FLAG_0` | none |
| `13_qc_month_hour_<site>` | month × UTC hour | `month, hour_utc, n, share_good, median_good_m, p25_good_m, p75_good_m, share_lowest_good, share_window` | `13_qc_diurnal_`, `13_qc_share_good_`, `13_qc_lowest_<site>` |
| `13_qc_monthly_<site>` | year-month | `year_month, days_with_data, days_in_month, n, share_good, share_window, share_no_temp_overlap, median_midday_m, median_night_m, share_lowest_night, segment` | `13_qc_monthly_<site>` (with `13_qc_breakpoints`) |
| `13_qc_breakpoints_<site>` | breakpoint × side | `breakpoint, segment_after, side, n, share_good, median_midday_m, median_night_m, share_lowest_night` | `13_qc_monthly_<site>` |
| `13_qc_days_<site>` | year-month | `year_month, days, missing, polar_day, polar_night, missing_polar_day, missing_polar_night` | none |
| `13_qc_sun_<site>` | sun event × season × relative hour | `event, season, h_rel, n, share_good, median_good_m` | none |
| `13_qc_events_<site>` | intervention day | `date, in_file, share_bit4` | none |
| `13_qc_temp_overlap_<site>` | firmware group | `group, months, n_night, median_night_m, p10_night_m, p90_night_m, mode_night_m, share_lowest_night, median_midday_m` | none |
| `13_qc_proposed_exclusions_<site>` | proposed month | `year_month, share_window, share_good, reason` | none |

## γ evaluation: step 11 (spec 11 D, amendment 1)

| Table | Grain | Columns | Figures |
|---|---|---|---|
| `11_gamma_cases_metrics` | site × case × reference × stratum set (`<set>:<season>` for seasons) × metric | `site_id, case, reference, stratum_set, metric, value, note`; metrics `n, bias, median_log_ratio, mad_log_ratio, frac_within_factor2, r_spearman, fb, nrmsd, n_days, mlr_lo, mlr_hi, slope_sym[_lo,_hi], ba_slope[_lo,_hi]` | `11_gamma_cases_<site>` (`default`, `clim`), `11s_ref_gamma_cases_<site>` (`gamma_med`, `gamma_site`, reference γ) |
| `11_gamma_cases_season` | site × case (incl. `era5_blh_center`) × reference × set × season | `site_id, case, reference, stratum_set, season, n, n_days, median_log_ratio, mlr_lo, mlr_hi, frac_within_factor2, fb, nrmsd, r_spearman, slope_sym, low_n, distance_class, borderline, note` | none |
| `11_gamma_decision` | site | `site_id, reference, method_class, stratum_set, distance_km, distance_class, borderline, period_start_utc, period_end_utc`; per case `<k>`: `n_<k>, n_days_<k>, mlr_<k>, mlr_lo_<k>, mlr_hi_<k>, low_n_<k>, season_max_<k>, season_max_mlr_<k>, season_max_n_days_<k>, era5_mlr_<k>, era5_n_<k>` (the last two not for `era5_blh_center`); `note` | `11_gamma_decision` |
| `11_gamma_regime_heights` | site × case × regime | `site_id, case, regime, share, n, share_B2, share_B3, p25_m, median_m, p75_m, median_gamma_bg, median_gamma_n` | none |
| `11_gamma_censored_null` | site × case × regime set × tail (FI-Hyy) | as `10_censored_null` without notes and pairing, plus `case` | none |
| `11_gamma_daily_vs_clim` | site | `site_id, n_days, median_abs_rel_dev, share_outside_25pct, sd_daily, sd_clim` | `11_gamma_by_month_<sites>` |
| `11_growth_clearday_<site>` | season × AST hour × variable | `site_id, reference, season, ast_hour, variable, n, median, q25, q75, n_clear_days` | `11_growth_clearday_<site>`, `11s_ref_growth_clearday_<site>` |

## Footprints: steps 08a, 08, 08c, 08t, 08v

| Table | Grain | Columns | Figures |
|---|---|---|---|
| `08_coverage_<site>[<per>]` | case × validity set × stratum (`all`, `daynight=`, `regime=`, `season=`) | `case, set, stratum, n_rows, n_with_h, n_inputs, n_valid, n_valid_qc, n_valid_qc_used, n_valid_qc_h, n_h_lt_zm, n_h_le_10`, `fail_code_{h10,ustar01,stab,nonneg}`, `fail_domain`, `fail_k_{rsl,el,stab}`, `fail_e_{rsl,el,tower,stab,ustar}` | none (`08_coverage_<site>` is drawn from the 08a input object, not this table) |
| `08_year_counts_<site>` | calendar year | `site_id, year, n_usable, n_months, n_measured, n_valid_qc, n_months_valid_qc` | none |
| `08_etc_multiplier_<site>` | site | `site_id, m_etc, m_loop, m_used, m_source, n_rows, zm_over_z0_median, zm_over_z0_min, n_rejected, n_rejected_loop, override_reason` | none |
| `08_nee_qc_<site>` | site (`year` NA) and year | `site_id, year, nee_qc_source, has_vut, has_cut, n_vut_qc, n_fluxnet, n_qc, n_qc0, n_qc_gapfilled, n_nee_good, n_used, status` | none |
| `08_ustar_check_<site>` | year × period | `site_id, year, period, n_good, n_ustar_filtered, n_gapfilled, pct_ustar_filtered_gapfilled, pct_gapfilled_explained_by_ustar` | none |
| `08_halfhour_metrics_<site><per>` | case × reference × set × subset × quantity (`x_peak, x_50…x_90`) × stratum × level (strata `all, daynight, hour, regime, season, zl_class, zmh_class`) | `case, reference, set, subset, quantity, stratum, level, n, median_lr, p25_lr, p75_lr, share_gt10, share_gt25, median_elasticity, median_dta_pp, p90_abs_dta_pp` | `08_dx80_by_zmh_`, `08_dx80_diurnal_<site><per>` (from the chunk `_hh` files) |
| `08_topaffected_<site><per>` | set × group (`all, top25, top10`) | `set, group, n, abs_lr_threshold, share_stable, share_unstable_diag, share_convective, share_night, share_removed_by_qc, share_removed_by_qc_or_used, share_removed_by_h_qc` | `08_topaffected_<site><per>` |
| `08_climatology_<site><per>` | case × period (`month`, `year`) × stored group (below), with `daynight` `day`, `night`, `all` | `case, period, period_type, set, subset, daynight, rowset, n, share80_in_target_area, f_domain, f_target_area`, field metrics, `jaccard80_vs_pblhmod` | `08_clim_contours_<site><per>` (80 % regions from the chunk climatologies) |
| `08_target_area_<site><per>` | case × set × subset × stratum (`all, daynight, zl_class`) × level | `case, set, subset, stratum, level, n, p10_ta_pct, median_ta_pct, p90_ta_pct, median_f_domain_pct, share_edge_90, share_f_domain_lt90, n_common, p10_dta_pp, median_dta_pp, p90_dta_pp` | none |
| `08_rebmann_<site><per>` | case × set × subset × row set | `case, set, subset, rows, n, share_ta90_case, share_ta90_pblhmod, share_cif80_case, share_cif80_pblhmod, pass_ta90_case, pass_ta90_pblhmod, verdict_differs` | none |
| `08_anchor_<site><per>` | case (±20 %) × set × z/L class | `case, set, zl_class, n, median_dxmax_pct, p90_abs_dxmax_pct, median_dfymax_pct, kljun_abs_dxmax_pct, kljun_abs_dfymax_pct` (one note row at sites without anchors) | none |
| `08_reference_years` (09) | station | `site_id, year, n_usable, n_months, n_measured, n_candidates, reason, rule, approved_year` | none |
| `08_target_area_masks[_check]` | station | `site_id, n_layers, layers, n_polygons, old_mask, old_n_layers_attr, cells_old, cells_new, area_old_m2, area_new_m2, first_layer, area_first_layer_m2, same, n_chunks, n_done, n_invalidated` | none |
| `08v_ffp_verify_<site>` | check part × quantity | `part, quantity, n, n_code_error, max_abs_diff, max_rel_diff, median_rel_diff, seconds_per_halfhour_code_A, seconds_per_halfhour_field_B` | none |
| `09_chunks_summary` (08c) | site | `site_id, n_chunks, n_halfhours, cpu_hours, max_est_min` | none |
| `08_storage_{sizes,metrics,summary}` (`scripts/dev/`) | storage option tests (spec 08, closed) | see the files | none |

`08_coverage_<site>` (figure): monthly share of valid half-hours per case, set × subset, from the 08a input object.

## Network collection: step 09 (reads the per-site tables only)

| Table | Grain | Columns | Figures |
|---|---|---|---|
| `09_network_run_status` | site × step | `site_id, step, state, start, end, exit_code, git_commit` | none |
| `09_network_inputs` | site | `site_id, file_fluxes, file_fluxnet, file_meteo, reference_start_utc, reference_end_utc, n_halfhours, n_excluded, n_fluxes, n_fluxnet, n_meteo, meteo_ta_levels, zm_min_m, zm_max_m, hc_min_m, hc_max_m, units_checked, era5_first_utc, era5_last_utc, era5_n_hours, era5_n_valid, era5_cell_dist_km` | none |
| `09_network_era5_cells` | site | `site_id, site_lat, site_lon, site_altitude_m, cell_lat, cell_lon, cell_dist_km, orog_m, dalt_m, lsm, sdor_m, orog_min_m, orog_max_m, lsm_min, lsm_max, sdor_max_m, n_cells` | none |
| `09_network_era5_metrics` | site × pair (1, 2, 6, 7, 8) × stratum (`all`, `regime=`) × metric | `site_id, pair, stratum, metric, value` | `09_network_logratio_<n>sites` (pairs 1, 6), `09_network_logratio_supp_<n>sites` (also 7, 8) |
| `09_network_regimes` | site × case (`default, clim, gamma_med, gamma_site`) × regime | as `02_regimes` | none |
| `09_network_filters` | site × case | as `02_cases` | none |
| `09_network_gamma_monthly` | site × month | `site_id, month, median, p25, p75, n_days, n_flagged` | none |
| `09_network_gamma_summary` | month | `month, n_sites, median, p25, p75, n_sites_2_levels` | none |
| `09_network_gamma_site` | site | `site_id, gamma_site, n_days` | none |
| `09_network_footprint_coverage` | site × case × set × stratum | `site_id` + `08_coverage` columns | none |
| `09_network_nee_qc` | site × year | `site_id` + `08_nee_qc` columns | none |
| `09_network_etc_multiplier` | site | `site_id, m_loop, m_used, m_source, n_rows, n_rejected, n_rejected_loop, zm_over_z0_median, zm_over_z0_min, override_reason, reference_year` | none |
| `09_network_flags` | site × flag | `site_id, flag, detail` | none |

Derived (not in `tables/`): `<derived_dir>/gamma/network_gamma.csv`, one row: `gamma_med, gamma_med_pooled, n_sites,
n_days, sites`.

## Stored footprint groups (spec 08 / 09 A5, storage c5)

Each chunk (site × case × month or half-month, `<derived_dir>/footprint/<site>/chunks/<site>_<case>_<part>`) stores:
- `_hh.csv.gz`: one row per half-hour, 1 m fields. Columns `timestamp_start_utc, case, h, x_peak, x_50, x_70, x_80,
  x_90, ffp_ok, n_nonfinite, n_retry, f_domain, f_target_area`, field metrics.
- `_sums.rds`: one summed field per group, as 5 m × 5 m block sums, float32, with the count `n` per group.
- `_clim.rds`: per `set|subset|rowset` × `daynight` (`day`, `night`, `all` = day + night): `n`, metrics, and the
  80 % region packed as a bitmask.
- `_meta.rds`: hash, `n`, seconds, workers, block.

Group key: `set|subset|daynight|rowset`.

- **Sets:** `kljun`, `kljun_norsl`, `etc`, `none`.
- **Subsets:** `valid`, `valid_qc`, `valid_qc_used`, `valid_qc_h`.
- **Day/night:** `day`, `night` (stored); `all` is built from the two.
- **Row sets:**
  - `own`: rows valid for the case. Every case.
  - `common_pe`: also valid for `pblhmod` and `era5`. Every case except `measured` and `bound`.
  - `common_pme`: also valid for `measured`. Only at sites with a measured case, and only for `pblhmod`,
    `pblhmod_clim`, `era5`, `measured`.

That makes up to 4 × 4 × 2 = 32 groups per row set: 64 per chunk for most cases, 96 for `pblhmod` / `pblhmod_clim` /
`era5` at measured sites, 32 for `measured` and `bound`. A group with no row is not stored. The `kljun` set is
always empty at the forest sites (BE-Bra, FI-Hyy), so they store 24 per row set.

Cases with chunks (reference year):
- **Every site:** `pblhmod`, `pblhmod_clim`, `era5`.
- **Extra sites** (`config/footprint_extra_sites.txt`: BE-Bra, FI-Hyy, FR-Gri, CH-Dav, SE-Svb, IT-Cp2): also
  `era5_bilinear`, `pblhmod_x0.8`, `pblhmod_x1.2`, `pblhmod_gamma001`, `pblhmod_gamma0025`, `oper`.
- **Sites with a primary measured source:** `measured`. Aerosol sources only on the window × convective rows.
- **FI-Hyy only:** `bound`, the turbulence lower bounds.

Climatologies used by step 08 (`08_climatology`): month (from the chunk, or its two halves joined) and year
(summed over the months), all from the block sums.
