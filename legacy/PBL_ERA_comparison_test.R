
library(data.table)


# Paths redacted for the public repository; in the current code, paths are set in config.local.yml.
# importo dati estratti da CT
# dati sito specifici in :
# <cluster path>/[SITO]/reanalysis-era5-single-levels__boundary_layer_height
BE_Bra_pblh_ERA5 <- fread("<local path>/BE-Bra__ERA5__reanalysis-era5-single-levels__boundary_layer_height__2024__lon+0__lat+0.csv", 
                          na.strings="-9999", sep=',', header = T, data.table = F)
# str(BE_Bra_pblh_ERA5)
# range(BE_Bra_pblh_ERA5$time)

# importo i dati calcolati da me nei L2 (scaricati da <cluster path>/[SITO]_output_FP_yyyymmdd)
BE_Bra_pblh_FFP <- fread("<local path>/PBL_heights.txt", 
                         na.strings="NA", sep='\t', header = T, data.table = F)
# str(BE_Bra_pblh_FFP)

# Trim my PBL to the ERA one
BE_Bra_pblh_FFP <- BE_Bra_pblh_FFP[which(BE_Bra_pblh_FFP$'TIMESTAMP' == BE_Bra_pblh_ERA5$time[1]):which(BE_Bra_pblh_FFP$'TIMESTAMP' == tail(BE_Bra_pblh_ERA5$time, 1)), ]


# TOCHECK
# - UTC vs local time (possible shift)
# - 30 to 60 min averaging || 60 to 30 min interpolation


windows()
plot(BE_Bra_pblh_ERA5$'time', BE_Bra_pblh_ERA5$'value', type='l')
points(BE_Bra_pblh_FFP$'TIMESTAMP', BE_Bra_pblh_FFP$'PBL_h', type='l', col=2)

windows()
plotwin <- which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-05-01 00:00:00') : which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-06-30 23:30:00')
plot(BE_Bra_pblh_FFP$'TIMESTAMP'[plotwin], BE_Bra_pblh_FFP$'PBL_h'[plotwin], type='l', col=2)
points(BE_Bra_pblh_ERA5$'time', BE_Bra_pblh_ERA5$'value', type='l')
windows()
plotwin <- which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-07-01 00:00:00') : which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-07-31 23:30:00')
plot(BE_Bra_pblh_FFP$'TIMESTAMP'[plotwin], BE_Bra_pblh_FFP$'PBL_h'[plotwin], type='l', col=2)
points(BE_Bra_pblh_ERA5$'time', BE_Bra_pblh_ERA5$'value', type='l')
windows()
plotwin <- which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-11-01 00:00:00') : which(as.character(BE_Bra_pblh_FFP$'TIMESTAMP') == '2024-12-31 23:00:00')
plot(BE_Bra_pblh_FFP$'TIMESTAMP'[plotwin], BE_Bra_pblh_FFP$'PBL_h'[plotwin], type='l', col=2, ylim=c(0,10000))
points(BE_Bra_pblh_ERA5$'time', BE_Bra_pblh_ERA5$'value', type='l')



windows()
plot(BE_Bra_pblh_ERA5$'value', BE_Bra_pblh_FFP$'PBL_h')


# downgrade (averaging) PBL from 30 to 60 minutes
BE_Bra_pblh_FFP_h <- BE_Bra_pblh_FFP %>%
  mutate(
    TIMESTAMP_h = floor_date(TIMESTAMP, unit = "60 mins")
  ) %>%
  group_by(TIMESTAMP_h) %>%
  summarise(
    PBL_h_0 = mean(PBL_h_0, na.rm = TRUE),
    PBL_h = mean(PBL_h, na.rm = TRUE),
    .groups = "drop"
  )
# head(BE_Bra_pblh_FFP_h)


windows()
plot(BE_Bra_pblh_ERA5$'value', BE_Bra_pblh_FFP_h$'PBL_h', 
     xlim=c(0,5000), ylim=c(0,5000))
abline(0,1, lty=2)


# ciclo giornaliero mensile
pblh_FFP_mon_day <- BE_Bra_pblh_FFP %>%
  mutate(
    months = floor_date(TIMESTAMP, unit = "month"),
    halfh = format(TIMESTAMP, '%H:%M'),
  ) %>%
  group_by(months, halfh) %>%
  summarise(
    median_pbl = median(PBL_h, na.rm = TRUE),
    q25 = quantile(PBL_h, 0.25, na.rm = TRUE),
    q75 = quantile(PBL_h, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

pblh_ERA_mon_day <- BE_Bra_pblh_ERA5 %>%
  mutate(
    months = floor_date(time, unit = "month"),
    halfh = format(time, '%H:%M'),
  ) %>%
  group_by(months, halfh) %>%
  summarise(
    median_pbl = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )



windows()
pblh_FFP_mon_day <- pblh_FFP_mon_day %>%
  mutate(hour = hour(Time_30min) + minute(Time_30min)/60)

# # Plot Ciclo Giornaliero

# COME AL SOLITO ggplot non va, chiedere il perchè
# pbld_plot <- ggplot(pblh_FFP_mon_day, 
#              aes(x = halfh, y = median_pbl, color = 'black')) +
#   
#   geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.22, linewidth = 0, color = NA) +
#   
#   geom_line(linewidth = 1) +
#   
#   facet_wrap(
#     ~months
#     # scales = "free_y",
#     # nrow = 1,
#     # labeller = as_labeller(units_lookup)
#   ) + 
#   scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
#   labs(x = "Hour of day", y = "PBL height") +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     aspect.ratio = 1
#   )
# pbld_plot


windows()
par(mfrow=c(3,4))
i <- 1
for(i in 1:length(unique(pblh_FFP_mon_day$'months'))) {
  
  cur.data <- which(pblh_FFP_mon_day$'months' == unique(pblh_FFP_mon_day$'months')[i])
  
  
  plot(pblh_FFP_mon_day$median_pbl[cur.data], type='l')
  polygon(c(seq_along(cur.data), rev(seq_along(cur.data))),
          c(pblh_FFP_mon_day$q25[cur.data], rev(pblh_FFP_mon_day$q75[cur.data])),
          col=adjustcolor(2, 0.5), border = NA)
  
  par(new=TRUE)
  
  
  plot(pblh_ERA_mon_day$median_pbl[cur.data], type='l')
  polygon(c(seq_along(cur.data), rev(seq_along(cur.data))),
          c(pblh_ERA_mon_day$q25[cur.data], rev(pblh_ERA_mon_day$q75[cur.data])),
          col=adjustcolor(4, 0.5), border = NA)
  
  
}


