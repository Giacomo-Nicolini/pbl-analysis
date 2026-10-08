
# ******************************************************************************
# Planetary Boundary Layer height estimation and field validation
# 
# ******************************************************************************



# TODO: ************************************************************************

# . 20190415: PBL sovrastima molto! capire se PBLe MLH sono effettivamente confrontabili o PBL è effettivamente stimato male !! DC!
#             inoltre, capire perchè rimane persistenza verso il pomeriggio (la curva è skewed a destra)
# . control and possible refine PBL formulas (VEDI PERCHè DI INVERNO SONO SEMPRE PIU ALTI)
# . find a proper way to make QC and find periods which should not be considered (e.g. rain, much cloudy, etc.)
# . find a proper way to estimate the uncertainty of my MLH (probably based on instrument specs, see the paper)
# . find a way to estimate the uncertainty of my approach
# . find a proper way to estimate the uncertainty of the temporal aggregations
# . apply Bland-Altman for each season
# . Boxplots of height differences computed at a specifc time, eg 12:00 [see Fig. 6 in Haefele 2014.]
# . find a proper way to average and compare results: see the reference paper, the considered mean over 4 hours [Hellén et al., 2019] 
# . SEE IF MEANINGFUL USING PROFILE WIND (IVAN ME LI DEVE MANDARE INSIEME AI DATI DI PBL)
# . fare notebook includendo FFP_sens_toPBL.R




# TEST DATSET:
# site: SMEAR II Hyytiala forest
# data downloaded from: AVAA open data publishing platform for research data @ https://avaa.tdata.fi/web/smart/smear/download
# period: 2014-06-01 to 2019-02-27
# SIte description:
# The Hyytiälä-SMEAR II station is located in a rather homogenous Scots pine (Pinus sylvestris L.) 
# stand on a slightly hilly terrain at Hyytiälä Forestry Field Station of the 
# University of Helsinki (61°51’N, 24°17’E, 181 m above sea level) 220 km NW from Helsinki
# vegetation map at 
# https://meta.icos-cp.eu/files/p7pE5oIyrjLNfcWVTGJ6120R/Fi-HYY_vegmap.png
# https://meta.icos-cp.eu/files/Dg3k9UW2-ao8O-IEYQwc6fAm/Fi-Hyy_VEGMAPinfo.txt


# libraries --------------------------------------------------------------------
kpacks <- c(
  'readtext',
  'R.matlab',
  'Biobase',
  'lubridate',
  'xts',
  'zoo',
  'data.table',
  'bigleaf',
  'kehra',
  'sirad',
  'BlandAltmanLeh',
  'viridis',
  'RColorBrewer'
)   
new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
if(length(new.packs)){
  install.packages(new.packs)
}
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)


# set data path
data.path <- 'C:/WORK/DATA/PBL_height/'



# Site metadata and roughness parameters ---------------------------------------

LAT <- 61.84741    # latitude [decimal degrees]
LON <- 24.29477    # longitude [decimal degrees]
ALT <- 181.0       # altitude
zm. <- 24          # ACTUAL measurement height [m]
hc. <- 18          # canopy height (height 15-20 m)  # NOTE: VERIFICARE PRECISIONE DI QUESTA MISURA

d. <- 2/3 * hc.    # displacement height [m]
z0. <- 0.15 * hc.  # roughness length [m]



# Set the analyses period ------------------------------------------------------

# available LiDar mixed Layer are from 2016 to 2018.
# currently only 2016 have been acquired.
# f it's worth, ask Ivan Mammarella for the other 2 years

t.start <- strptime('20160101000000','%Y%m%d%H%M%S', tz = 'GMT')
t.end <- strptime('2016123123300','%Y%m%d%H%M%S', tz = 'GMT')




# MIXING LAYER HEIGHT (MLH) field data -----------------------------------------
# TKE profile and MLH from Doppler lidar measurements at Hyytiälä
# Minimum MLH that can be determined is 60m above instrument. 
# MLH = 0 indicates MLH < 60 m.
# MLH is not defined in rain, but rain filtering may not be perfect.

# PBL height and Mixing layer references:
# - small hint
file.show("C:/WORK/DATA/PBL_height/FI-Hyy_LiDar/readme_tke_mlh.txt")
# # reference paper
# openPDF("C:/WORK/new_references/Hellén_2018_ACP_Long-term measurements of volatile organic compounds highlight.pdf")
# # other references in Mendeley

# Data are daily files in Matlab *.MAT ﬁle structure 
mat.path.name <- list.files('C:/WORK/DATA/PBL_height/FI-Hyy_LiDar/', ".mat", full.names = TRUE) 
MLH.ls <- list(NULL)
TKE.ls <- list(NULL)
TKE.elev <- c(NULL)

i <- 1
for(i in 1:length(mat.path.name)){
  
  cur.mat <- readMat(mat.path.name[i]) 
  
  # Convert timestamp *************
  # "TKE_time is time stamp for both MLH and TKE."
  # Matlab uses "January 01, 0000", a fictional reference date, to calculate its date number. 
  # The origin of time for the "POSIXct" class in R is, ‘1970-01-01 00:00.00 UTC’.
  # Where 719529 is ‘1970-01-01 00:00.00 UTC’ in Matlab's datenum and 86400 the number of seconds in an standard UTC day
  cur.ts <- as.POSIXct((cur.mat$TKE.time - 719529)*86400, origin = "1970-01-01", tz = "GMT")
  # round actual timestamp to 30 minutes sharp
  cur.ts.hh <- round_date(cur.ts, '30 min')
  
  # Create MLH daily dataframe list *************
  # "MLH is determined with constant threshold. Applying lower thresholds (column 2) give a little higher MLH."
  MLH.ls[[i]] <- data.frame(Timestamp.act = cur.ts, 
                            Timestamp = cur.ts.hh, 
                            MHL_1 = cur.mat$MLH[,1],
                            MHL_2 = cur.mat$MLH[,2])
  
  # Create TKE daily dataframe list *************
  # NOTE: TKE are currently not used, but are interesting to explore
  TKE.ls[[i]] <- data.frame(Timestamp.act = cur.ts, 
                            Timestamp = cur.ts.hh, 
                            cur.mat$TKE)
  # TKE_elevation is m above instrument. (all the same)
  if(i == 1){TKE.elev <- cur.mat$TKE.elevation}

}

# convert and merge lists to dataframe
MLH.df <- do.call(rbind, MLH.ls)
TKE <- do.call(rbind, TKE.ls)
# head(MLH.df)
# str(MLH.df)
# dim(MLH.df)

# Time specification ******************************
# "Note that some days are missing (no empty files created for gaps)"

# summary(diff(as.double(MLH.df$Timestamp)))
# length(MLH.df$Timestamp)
# windows()
# plot(as.double(MLH.df$Timestamp))

# remove duplicated times
# NOTE: data are not the same, but are only 12 reocrds so are removed
# MLH.df$Timestamp[which(duplicated(MLH.df$Timestamp))]
# MLH.df[which(duplicated(MLH.df$Timestamp)), ]
# cbind(MLH.df[which(duplicated(MLH.df$Timestamp)), ], MLH.df[which(duplicated(MLH.df$Timestamp))-1, ])
MLH.df <- MLH.df[-which(duplicated(MLH.df$Timestamp)), ]

# fill the gaps with empty rows
# Align TA and master timestamps and in case of missing values, fill with -9999
MLH.df.xts <- xts(MLH.df, order.by = MLH.df$Timestamp, tzone = 'GMT') # covert to xts object
MLH.df.xts <- merge(MLH.df.xts, 
                 xts(rep(1, length(seq(first(MLH.df$Timestamp), last(MLH.df$Timestamp), 1800))), 
                     order.by = seq(first(MLH.df$Timestamp), last(MLH.df$Timestamp), 1800), tzone = 'GMT')) # align
MLH.df <- as.data.frame(coredata(MLH.df.xts)[, -5]) # and get back to vector (no xts object)
MLH.df$Timestamp <- time(MLH.df.xts)

# timestamp [%Y-%m-%d %H:%M:%S]
# > use this timestamp as the MASTER for all other variables !!
# MLH.df.TS <- strptime(MLH.df$Timestamp, '%Y%m%d%H%M%S', tz = 'GMT')
MLH.df.TS <- MLH.df$Timestamp
# range(MLH.df.TS)
# summary(diff(as.double(MLH.df.TS)))


# Time factors for period analyses *******************
# NOTE: CB timestamps below must be replaced by these
MLH.TS.Y <- format(MLH.TS, '%Y')
MLH.TS.seas <- getSeason(MLH.TS)
MLH.TS.mon <- months(MLH.TS)
MLH.TS.HHMM <- format(MLH.TS, '%H:%M')
# define day/night time
# see EP Help @ https://www.licor.com/env/support/EddyPro/faqs.html
MLH.TS.daytime <- (hour(MLH.TS) > 6 & hour(MLH.TS) < 18)
# potential radiation :
#  TODO: la radpot è a risoluzione oraria: o trovo il modo di interpolare a 30 minuti o trovo un altra funzione
# rad.pot <- extrat(unique(dayOfYear(MLH.TS)), # day number in the year (julian day)
#                   lat = radians(LAT)) # latitude in radians
# windows()
# plot(rad.pot$ExtraTerrestrialSolarRadiationHourly[1:50], type='l')


# Mixing Layer height 
# "MLH is determined with constant threshold (MHL_1). Applying lower thresholds (MHL_2) give a little higher MLH."
MLH.1 <- MLH.df$'MHL_1'
MLH.2 <- MLH.df$'MHL_2'


summary(MLH.1)
summary(MLH.2)


windows()
plot(MLH.TS, MLH.2, col=4)
points(MLH.TS, MLH.1, col=1)

windows()
plot(tapply(MLH.2, MLH.TS.HHMM, mean, na.rm=T), type = 'l', col=4, ylim=c(0,800), lwd=2)
points(tapply(MLH.1, MLH.TS.HHMM, mean, na.rm=T), type = 'l', lwd=2)
points(tapply(MLH.2, MLH.TS.HHMM, median, na.rm=T), type = 'l', col=4)
points(tapply(MLH.1, MLH.TS.HHMM, median, na.rm=T), type = 'l')
abline(h=0)


# NOTE: POSSIBLY ADD A QC FILTER HERE
#       option1: remove data for which the differece between the two masures is above a certain threshold
#       option2: remove data for which the linear fit residuals are above a certain threshold
#



# CLOUD BASE HEIGHT ------------------------------------------------------------
# downloaded independently from the meteorology data, NO QUALITY CHECKED

# NOTA: decidere se mantenere questi dati e confrontarli o NON considerarli più per l'analisi

# Cloud base height layer 1 and Cloud base height layer 2, Ceilometer
# A ceilometer is a device that uses a laser or other light source to determine the height of a cloud ceiling or cloud base. 
# When based on laser, it is a type of atmospheric lidar.


# PBL height and cloud cover references:
# # - small hint
# file.show("C:/WORK/DATA/PBL_height/PBL-clouds_hint1.docx")
# # reference paper
# openPDF("C:/WORK/new_references/Ungureanu_2010_RRP_Investigation of the cloud cover and planetary.pdf")
# other references in Mendeley

# data
CB.H <- fread("C:/WORK/DATA/PBL_height/FI-Hyy_smeardata_2014-2019_Cloud-base-height-layer.txt", header = TRUE, sep = ',', na.strings = 'NaN', data.table = FALSE)

# timestamp [%Y-%m-%d %H:%M:%S]
# > use this timestamp as the master for all other variables !!
CB.H.TS <- strptime(paste(CB.H$Year, 
                          formatC(CB.H$Month, digits = 0, width = 2, format = "f", flag = "0"),
                          formatC(CB.H$Day, digits = 0, width = 2, format = "f", flag = "0"), 
                          formatC(CB.H$Hour, digits = 0, width = 2, format = "f", flag = "0"),
                          formatC(CB.H$Minute, digits = 0, width = 2, format = "f", flag = "0"),
                          formatC(CB.H$Second, digits = 0, width = 2, format = "f", flag = "0"),
                          sep = ''), 
                    '%Y%m%d%H%M%S', tz = 'GMT')
# range(CB.H.TS)
# summary(diff(as.double(CB.H.TS)))

# Cut the dataset to the desidered time period
# Mixed Layer data should be for the whole 2018 
CB.H <- CB.H[which(CB.H.TS == t.start):which(CB.H.TS == t.end), ]
CB.H.TS <- CB.H.TS[which(CB.H.TS == t.start):which(CB.H.TS == t.end)]

# # Time factors for period analyses *******************
# CB.H.TS.Y <- format(CB.H.TS, '%Y')
# CB.H.TS.seas <- getSeason(CB.H.TS)
# CB.H.TS.mon <- months(CB.H.TS)
# CB.H.TS.HHMM <- format(CB.H.TS, '%H:%M')
# CB.H.TS.daytime <- (hour(CB.H.TS) > 6 & hour(CB.H.TS) < 18)


# Cloud base height layer 1
CB.H.1 <- CB.H$'HYY_META.Cloud_layer_1'
# Cloud base height layer 2
CB.H.2 <- CB.H$'HYY_META.Cloud_layer_2'

# windows()
# plot(CB.H.TS, CB.H.1)
# points(CB.H.TS, CB.H.2, col=4)
# windows()
# plot(CB.H.1, CB.H.2)


# NOTE: POSSIBLY ADD A QC FILTER HERE
#       option1: remove data for which the differece between the two masures is above a certain threshold
#       option2: remove data for which the linear fit residuals are above a certain threshold
# BUT take cae because dataset #2 is sensibly shorter than #1
# length(na.omit(CB.H.1)); length(na.omit(CB.H.2))




# FLUXES and TURBULENCE data ---------------------------------------------------

# # see variables list here:
# file.show("C:/WORK/DATA/PBL_height/FI-Hyy_smeardata_2014-2019_fluxes-turbulence_variables-list.txt")
# # daata are split in different files


flux.turb.files <- list.files(data.path, pattern = 'flux', full.names = TRUE)[-length(list.files(data.path, pattern = 'flux', full.names = TRUE))] # to remove the variable list file
flux.turb.list <- list(NULL)
for(i in 1:length(flux.turb.files)){
  
  flux.turb.list[[i]] <- fread(flux.turb.files[i], header = TRUE, sep = ',', na.strings = 'NaN', data.table = FALSE)
  
  flux.turb.TS <- strptime(paste(flux.turb.list[[i]]$Year, 
                                 formatC(flux.turb.list[[i]]$Month, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(flux.turb.list[[i]]$Day, digits = 0, width = 2, format = "f", flag = "0"), 
                                 formatC(flux.turb.list[[i]]$Hour, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(flux.turb.list[[i]]$Minute, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(flux.turb.list[[i]]$Second, digits = 0, width = 2, format = "f", flag = "0"),
                                 sep = ''), 
                           '%Y%m%d%H%M%S', tz = 'GMT')
  
  # covert to xts object
  flux.turb.list[[i]] <- xts(flux.turb.list[[i]], order.by = flux.turb.TS, tzone = 'GMT')
  
}
# lapply(flux.turb.list, dim)

# Align the datasets in the list according with timestamps and in case of missing values, fill with -9999 
# Go for a coffe, it takes a while...
flux.turb.xts <- do.call(merge, flux.turb.list)
# head(flux.turb.xts)
# range(time(flux.turb.xts))

# Align fluxes and CB-H master timestamps and in case of missing values, fill with -9999
flux.turb.xts <- merge(flux.turb.xts, 
                       xts(seq_along(CB.H.TS), order.by = CB.H.TS, tzone = 'GMT'))


# Cut the dataset to the desidered time period
# Mixed Layer data should be for the whole 2018 
flux.turb.xts <- window(flux.turb.xts, start = t.start, end = t.end)

# fluxes timestamp
flux.turb.TS <- time(flux.turb.xts)
# range(flux.turb.TS)

# and get back to normal dataset (no more xts object)
flux.turb <- data.frame(flux.turb.xts) 
# str(flux.turb)


# set sensible variables *******************************************************
# NOTE: considering the 24 m tower as primary flux measuring setup

# Monin-Obukov length @ 24 m [m] 
MOL.24 <- flux.turb$'HYY_EDDY233.MO_length'
# summary(MOL.24)
# Monin-Obukov length @ 27 m [m] 
MOL.27 <- flux.turb$'HYY_EDDYMAST.MO_length_270'
# summary(MOL.27)
# Monin-Obukov length @ 3m (subcanopy) [m] 
MOL.subm <- flux.turb$'HYY_EDDYSUB.MO_length_subm'
# summary(MOL.subm)

# friction velocity @ 24 m [ms-2]
ustar.24 <- flux.turb$'HYY_EDDY233.u_star'
# friction velocity @ 27 m [ms-2]
ustar.27 <- flux.turb$'HYY_EDDYMAST.u_star_270'
# friction velocity @ 3 m (subcanopy) [ms-2]
ustar.subm <- flux.turb$'HYY_EDDYSUB.u_star_subm'
# windows()
# plot(ustar.24, ustar.27)

# sensible heat @ 24 m [Wm-2]
H.24 <- flux.turb$'HYY_EDDY233.H'
# sensible heat @ 27 m [Wm-2]
H.27 <- flux.turb$'HYY_EDDYMAST.H_270'
# sensible heat @ 3 m (subcanopy) [Wm-2]
H.subm <- flux.turb$'HYY_EDDYSUB.H_subm'
# windows()
# plot(H.24, H.27)



# METEOROLOGY data -------------------------------------------------------------

# # see variables list here:
# file.show("C:/WORK/DATA/PBL_height/FI-Hyy_smeardata_2014-2019_meteorology_variables-list.txt")

# NOTE: PROFILE DATA ARE HERE!


meteo.files <- list.files(data.path, pattern = 'meteo', full.names = TRUE)[-length(list.files(data.path, pattern = 'meteo', full.names = TRUE))] # to remove the variable list file
meteo <- fread(meteo.files, header = TRUE, sep = ',', na.strings = 'NaN', data.table = FALSE)
meteo.TS <- strptime(paste(meteo$Year, 
                           formatC(meteo$Month, digits = 0, width = 2, format = "f", flag = "0"),
                           formatC(meteo$Day, digits = 0, width = 2, format = "f", flag = "0"), 
                           formatC(meteo$Hour, digits = 0, width = 2, format = "f", flag = "0"),
                           formatC(meteo$Minute, digits = 0, width = 2, format = "f", flag = "0"),
                           formatC(meteo$Second, digits = 0, width = 2, format = "f", flag = "0"),
                           sep = ''), 
                     '%Y%m%d%H%M%S', tz = 'GMT')

# covert to xts object
meteo.xts <- xts(meteo, order.by = meteo.TS, tzone = 'GMT')

# Align fluxes and CB-H master timestamps and in case of missing values, fill with -9999
meteo.xts <- merge(meteo.xts,
                   xts(seq_along(CB.H.TS), order.by = CB.H.TS, tzone = 'GMT'))
# head(meteo.xts)
# range(time(meteo.xts))

# Cut the dataset to the desidered time period
# Mixed Layer data should be for the whole 2018 
meteo.xts <- window(meteo.xts, start = t.start, end = t.end)


# and get back to normal dataset (no more xts object)
meteo <- data.frame(meteo.xts) 


# set sensible variables *******************************************************
# NOTE: considering the 24 m tower as primary flux measuring setup


# Atmospheric pressure at ground level (180 m a.s.l.) [from hPa to kPa]
AP <- meteo$'HYY_META.Pamb0' * 0.1
# summary(AP)

# Air temperature at 16.8 m height [°C] 
AT.16 <- meteo$'HYY_META.T168'

# Air temperature at 33.6 m height [°C] 
AT.33 <- meteo$'HYY_META.T336'

# Derive air temperature at 24 m height [°C] 
# wheighted mean, wheights are considered as the reciprocal of the difference in height 
AT.24 <- apply(cbind(AT.16, AT.33), 1, function(x) weighted.mean(x, w = c(1/7, 1/10), na.rm = TRUE))

# windows()
# plot(meteo.TS, AT.16)
# points(meteo.TS, AT.33, col=4)
# points(meteo.TS, AT.24, col=2)

# Air density of moist air [kg m-3]
rho <- bigleaf::air.density(Tair = AT.24, 
                            pressure = AP)

# Specific heat of air for constant pressure (J K-1 kg-1)
cp <- bigleaf.constants()$'cp'


# Precipitation (liquid water equivalent) accumulated during previous 1 min
# Vaisala FD12P weather sensor at 18 m height
precipitation <- meteo$'HYY_META.Precipacc'
# windows()
# plot(meteo.TS, precipitation, col=4)


# Precipitation intensity during previous 1 min
# Vaisala FD12P weather sensor at 18 m height
precipitation.int <- meteo$'HYY_META.Precipint'
# par(new=T)
# plot(meteo.TS, precipitation.int, col=1, type = 'h')



# ANCILLARY data ---------------------------------------------------------------

# # see variables list here:
# file.show("C:/WORK/DATA/PBL_height/FI-Hyy_smeardata_2014-2019_ancillary_variables-list.txt")
# daata are split in different files


ancillary.files <- list.files(data.path, pattern = 'ancillary', full.names = TRUE)[-length(list.files(data.path, pattern = 'ancillary', full.names = TRUE))] # to remove the variable list file
ancillary.list <- list(NULL)
for(i in 1:length(ancillary.files)){
  
  ancillary.list[[i]] <- fread(ancillary.files[i], header = TRUE, sep = ',', na.strings = 'NaN', data.table = FALSE)
  
  ancillary.TS <- strptime(paste(ancillary.list[[i]]$Year, 
                                 formatC(ancillary.list[[i]]$Month, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(ancillary.list[[i]]$Day, digits = 0, width = 2, format = "f", flag = "0"), 
                                 formatC(ancillary.list[[i]]$Hour, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(ancillary.list[[i]]$Minute, digits = 0, width = 2, format = "f", flag = "0"),
                                 formatC(ancillary.list[[i]]$Second, digits = 0, width = 2, format = "f", flag = "0"),
                                 sep = ''), 
                           '%Y%m%d%H%M%S', tz = 'GMT')
  
  # covert to xts object
  ancillary.list[[i]] <- xts(ancillary.list[[i]], order.by = ancillary.TS, tzone = 'GMT')
  
}
# lapply(ancillary.list, dim)

# Align the datasets in the list according with timestamps and in case of missing values, fill with -9999 (it takes awhile)
ancillary.xts <- do.call(merge, ancillary.list)
# head(ancillary.xts)
# range(time(ancillary.xts))

# Align fluxes and CB-H master timestamps and in case of missing values, fill with -9999
ancillary.xts <- merge(ancillary.xts, 
                       xts(seq_along(CB.H.TS), order.by = CB.H.TS, tzone = 'GMT'))

# Cut the dataset to the desidered time period
# Mixed Layer data should be for the whole 2018 
meteo.xts <- window(meteo.xts, start = t.start, end = t.end)

# fluxes timestamp
ancillary.TS <- time(ancillary.xts)

# and get back to normal dataset (no more xts object)
ancillary <- data.frame(ancillary.xts) 
# str(ancillary)



# PBL HEIGHT ESTIMATION ========================================================

source('calc_PBL_height.R')

PBL.H <- calc_PBL_height(TIMESTAMP = MLH.TS,   # timestamp [%Y-%m-%d %H:%M:%S]
                         LAT = LAT,             # point latitude [decimal degrees]
                         zm = zm.,              # measurement height [m]
                         d = d.,                # displacement height [m]
                         MOL = MOL.24,          # Monin-Obukov length [m]
                         ustar = ustar.24,      # friction velocity [ms-2]
                         H = H.24,              # sensible heat flux [Wm-2]
                         AT = AT.24,            # air tmperature [K]
                         AP = AP,               # air pressure [Pa]
                         rho = rho,             # density of dry air [kg m-3]
                         cp = cp,               # air_heat_capacity [J kg-1 K-1]
                         altit = ALT,           # site altitude [m]
                         # processing options
                         GF.PBL = TRUE,         # fill less than 1 hour gaps in PBL heights by linear interpolation 
                         smooth.PBL = TRUE,     # Smooth PBL heights ONLY IF temporal resolution is less than 1 hour
                         plot.atmc = FALSE,     # set TRUE to plot atmospheric conditions frequency
                         plot.PBL = FALSE,      # set TRUE to plot PBL height estimations
                         save.PBL = FALSE)
# head(PBL.H)
# summary(PBL.H)

# convert from 'ts' object to numeric
PBL.H <- c(PBL.H)
# class(PBL.H)

# PBL statistics
tapply(PBL.H, MLH.TS.Y, summary)
tapply(MLH.1, MLH.TS.Y, summary)
quantile(na.omit(PBL.H), seq(0,1,0.1))
quantile(na.omit(MLH.2), seq(0,1,0.1))
windows(); hist(PBL.H, 50)
table(cut(PBL.H, 20))
table(PBL.H < 60) # 60 m is the plausibility limit of MLH estimation (Hellén et al., 2019)
# data below should be marked as 0

# POSSIBLE PBL QC HERE #################


# PBL HEIGHT VALIDATION ########################################################


# DATA PREPARATION *************************************************************

MLH.1.f <- MLH.1

# 1. Remove CB data during rain 
MLH.1.f[precipitation > 0] <- NA

windows() 
plot(MLH.TS, MLH.1, col=2)
points(MLH.TS, MLH.1.f, col=8)
points(MLH.TS, PBL.H)

windows()
plot(MLH.1, PBL.H)
points(MLH.1.f, PBL.H, col=4)


Nhh <- 1:(48*7)
windows()
plot(PBL.H[Nhh], type='o')
points(MLH.1[Nhh], col=2, type='o')




# Yearly daily cycles -----------------------------------------------------

MLH.1.y.dc <- aggregate(MLH.1, list(MLH.TS.Y, MLH.TS.HHMM), mean, na.rm=T)
MLH.1.y.dc <- MLH.1.y.dc[order(MLH.1.y.dc$Group.1), ]
PBL.H.y.dc <- aggregate(PBL.H, list(MLH.TS.Y, MLH.TS.HHMM), mean, na.rm=T)
PBL.H.y.dc <- PBL.H.y.dc[order(PBL.H.y.dc$Group.1), ]

windows()
# par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(1,1,.1,.1))
i <- 0
for(i in 1:length(unique(MLH.1.y.dc$Group.1))){
  plot(MLH.1.y.dc$x[MLH.1.y.dc$Group.1 %in% unique(MLH.1.y.dc$Group.1)[i]], type='o', ylim=c(0,3000), 
       col = 1, xlab = 'DoW', ylab = 'Mixing height [m]')
  points(PBL.H.y.dc$x[PBL.H.y.dc$Group.1 %in% unique(PBL.H.y.dc$Group.1)[i]], type='o', col=4)
  abline(h=0, v=24, lty=c(1,32))
  if(i==1)legend('topleft', legend=c('MLH', 'PBL'), col=c(1,4), lty=1, pch=1, bty = 'n')
}



# Yearly/seasonal daily cycles -------------------------------------------------

MLH.1.y.s.dc <- aggregate(MLH.1, list(MLH.TS.Y, MLH.TS.seas, MLH.TS.HHMM), mean, na.rm=T)
# order it according to year and season
MLH.1.y.s.dc <- MLH.1.y.s.dc[order(MLH.1.y.s.dc$Group.1, MLH.1.y.s.dc$Group.2), ]

PBL.H.y.s.dc <- aggregate(PBL.H, list(MLH.TS.Y, MLH.TS.seas, MLH.TS.HHMM), mean, na.rm=T)
# order it according to year and season
PBL.H.y.s.dc <- PBL.H.y.s.dc[order(PBL.H.y.s.dc$Group.1, PBL.H.y.s.dc$Group.2), ]


# PLOT on YEARLY base ********************
# set season color (in this order: unique(MLH.1.y.s.dc$Group.2))
# seas.col <- cividis(4, begin = 0.3)
# seas.col <- viridis(4, begin = 0.3)
seas.col <- brewer.pal(11, 'RdYlBu')[c(9,4,2,11)]

windows()
# par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(1,1,.1,.1), pty='s')

i <- 0
for(i in 1:length(unique(MLH.1.y.s.dc$Group.1))){
  
  y.indx <- MLH.1.y.s.dc$Group.1 %in% unique(MLH.1.y.s.dc$Group.1)[i]
  
  ii <- 0
  for(ii in 1:length(unique(MLH.1.y.s.dc$Group.2))){
    
    s.indx <- Reduce('&', list(y.indx, MLH.1.y.s.dc$Group.2 %in% unique(MLH.1.y.s.dc$Group.2)[ii]))
    
    if(ii == 1){
      plot(MLH.1.y.s.dc$x[s.indx], type='o', pch = 19, ylim=c(0,1800), col = seas.col[ii], cex = 1.5, xlab = 'Time [# halfhours]', ylab = 'Height [m]')
    }else{
      points(MLH.1.y.s.dc$x[s.indx], type='o', pch = 19, cex = 1.5, col = seas.col[ii])
    }
    points(PBL.H.y.s.dc$x[s.indx], type='o', pch = 1, cex = 1.5, col = seas.col[ii])
  }
  mtext(unique(MLH.1.y.s.dc$Group.1)[i], 3, 0)
  abline(h=0, v=24, lty=c(1,32), col=c(1,8))
  
  
  if(i==1){
    legend('topleft', legend=unique(MLH.1.y.s.dc$Group.2), col=seas.col, lty=1, lwd=2, bty = 'n')
    legend('topright', legend=c('MLH','PBL'), col='darkgrey', pch=c(19,1), pt.cex=1.5, bty = 'n')
  }
  
}


# # PLOT on SEASONAL base ********************
# MEANINGFUL ONLY IN CASE OF MULTIPLE YEARS DATASET
# 
# # set season color (in this order: unique(MLH.1.y.s.dc$Group.1))
# # year.col <- cividis(4, end = 0.7)
# year.col <- viridis(length(unique(MLH.1.y.s.dc$Group.1))-1, end = 0.7)
# 
# windows()
# par(mfrow=c(2,2), mar=c(1,1,1,1), oma=c(1,1,.1,.1), pty='s')
# 
# i <- 1
# for(i in 1:length(unique(MLH.1.y.s.dc$Group.2))){
#   
#   s.indx <- MLH.1.y.s.dc$Group.2 %in% unique(MLH.1.y.s.dc$Group.2)[i]
#   
#   ii <- 0
#   for(ii in 1:(length(unique(MLH.1.y.s.dc$Group.1))-1)){
#     
#     y.indx <- Reduce('&', list(s.indx, MLH.1.y.s.dc$Group.1 %in% unique(MLH.1.y.s.dc$Group.1)[ii]))
#     
#     if(ii == 1){
#       plot(MLH.1.y.s.dc$x[y.indx], type='o', pch = 19, ylim=c(0,3000), col = year.col[ii], cex = 1.5, xlab = 'Time [# halfhours]', ylab = 'Height [m]')
#     }else{
#       points(MLH.1.y.s.dc$x[y.indx], type='o', pch = 19, cex = 1.5, col = year.col[ii])
#     }
#     points(PBL.H.y.s.dc$x[y.indx], type='o', pch = 1, cex = 1.5, col = year.col[ii])
#   }
#   mtext(unique(MLH.1.y.s.dc$Group.2)[i], 3, -1)
#   abline(h=0, v=24, lty=c(1,32), col=c(1,8))
#   
#   
#   if(i==1){
#     legend('topleft', legend=unique(MLH.1.y.s.dc$Group.1)[1:5], col = year.col, lty=1, lwd=2, bty = 'n')
#     legend('topright', legend=c('CB','PBL'), col='darkgrey', pch=c(19,1), pt.cex=1.5, bty = 'n')
#     }
#   
# }

# Yearly/monthly daily cycles -------------------------------------------------


MLH.1.y.m.dc <- aggregate(MLH.1, list(MLH.TS.Y, MLH.TS.mon, MLH.TS.HHMM), median, na.rm=T)
# order it according to year and month
# TODO: GROUP MONTHS ON SEASONAL BASIS, (CREATE A PROPER RANKING COLUMN FOR MONTHS, E.G. starting on NOV-FEB is winter and so on)
MLH.1.y.m.dc <- MLH.1.y.m.dc[order(MLH.1.y.m.dc$Group.1, MLH.1.y.m.dc$Group.2), ]
# unique(MLH.1.y.m.dc$Group.2)
# mon.rank <- c(6,10,2,4,3,9,8,5,7,1,12,11) # starting on NOV
mon.rank <- c(5,9,1,3,2,8,7,4,6,12,11,10) # starting on DEC
names(mon.rank) <- unique(MLH.1.y.m.dc$Group.2)
MLH.1.y.m.dc$mon.rank <- rep(mon.rank, each=48)
MLH.1.y.m.dc <- MLH.1.y.m.dc[order(MLH.1.y.m.dc$mon.rank, MLH.1.y.m.dc$Group.3), ]

PBL.H.y.m.dc <- aggregate(PBL.H, list(MLH.TS.Y, MLH.TS.mon, MLH.TS.HHMM), median, na.rm=T)
# order it according to year and season
PBL.H.y.m.dc <- PBL.H.y.m.dc[order(PBL.H.y.m.dc$Group.1, PBL.H.y.m.dc$Group.2), ]
PBL.H.y.m.dc$mon.rank <- rep(mon.rank, each=48)
PBL.H.y.m.dc <- PBL.H.y.m.dc[order(PBL.H.y.m.dc$mon.rank, PBL.H.y.m.dc$Group.3), ]


# PLOT on YEARLY base ********************
# set season color (in this order: unique(MLH.1.y.m.dc$Group.2))
month.col <- cividis(12, begin = 0.3)
# month.col <- viridis(12, begin = 0.3)
# month.col <- brewer.pal(11, 'RdYlBu')[c(9,4,2,11)]

windows()
par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(1,1,.1,.1))#, pty='s'

i <- 0
for(i in 1:length(unique(MLH.1.y.m.dc$Group.1))){
  
  y.indx <- MLH.1.y.m.dc$Group.1 %in% unique(MLH.1.y.m.dc$Group.1)[i]
  
  ii <- 0
  for(ii in 1:length(unique(MLH.1.y.m.dc$Group.2))){
    
    s.indx <- Reduce('&', list(y.indx, MLH.1.y.m.dc$Group.2 %in% unique(MLH.1.y.m.dc$Group.2)[ii]))
    
    
    plot(MLH.1.y.m.dc$x[s.indx], type='o', pch = '', ylim=c(0,1600), col = month.col[ii], cex = 1.5, xlab = 'Time [# halfhours]', ylab = 'Height [m]')
    points(PBL.H.y.m.dc$x[s.indx], type='o', pch = '', cex = 1.5, lwd=2, col = month.col[ii])
    
    if(ii == 1){
      legend('topleft', legend=c('MLH','PBLH'), col=c('darkgrey','darkgrey'), lty=1, lwd=c(1,2), bty = 'n')
    }
    
    mtext(unique(MLH.1.y.m.dc$Group.2)[ii], 3, 0)
    abline(h=0, v=24, lty=c(1,32), col=c(1,8))
  }
  
  # mtext(unique(MLH.1.y.m.dc$Group.1)[i], 3, 0)
  # abline(h=0, v=24, lty=c(1,32), col=c(1,8))
  
  
  # if(i==1){
  #   legend('topleft', legend=unique(MLH.1.y.m.dc$Group.2), col=month.col, lty=1, lwd=2, bty = 'n')
  #   legend('topright', legend=c('CB','PBL'), col='darkgrey', pch=c(19,1), pt.cex=1.5, bty = 'n')
  # }
  
}



# Bland-Altman ------------------------------------------------------------

# from https://en.wikipedia.org/wiki/Bland%E2%80%93Altman_plot:
# Bland and Altman plots are extensively used to evaluate the agreement among two different instruments
# or two measurements techniques. Bland and Altman plots allow us to investigate the existence of any 
# systematic difference between the measurements (i.e., fixed bias) and to identify possible outliers.
# The mean difference is the estimated bias, and the SD of the differences measures the random fluctuations around this mean.
# If the mean value of the difference differs significantly from 0 on the basis of a 1-sample t-test, this indicates the presence of fixed bias. 
# If there is a consistent bias, it can be adjusted for by subtracting the mean difference from the new method. 
# It is common to compute 95% limits of agreement for each comparison (average difference ± 1.96 standard deviation of the difference),
# which tell us how far apart measurements by 2 methods were more likely to be for most individuals. 
# If the differences within mean ± 1.96 SD are not clinically important, the two methods may be used interchangeably. 
# Bland and Altman plots were also used to investigate any possible relationship of the discrepancies between 
# the measurements and the true value (i.e., proportional bias). The existence of proportional bias indicates 
# that the methods do not agree equally through the range of measurements (i.e., the limits of agreement will depend on the actual measurement). 
# To evaluate this relationship formally, the difference between the methods should be regressed on the average of the 2 methods. 
# When a relationship between the differences and the true value was identified (i.e., a significant slope of the regression line), 
# regression-based 95% limits of agreement should be provided
# vedi anche:
# Giavarina D (2015) Understanding Bland Altman analysis. Biochem Medica 25:141–151. doi: 10.11613/BM.2013.003

BA.stats <- bland.altman.stats(MLH.1, PBL.H, two = 1.96, mode = 1, conf.int = 0.95)
# str(BA.stats)
# summary(BA.stats$means)
# summary(BA.stats$diffs)

# plot embed
windows()
bland.altman.plot(MLH.1, PBL.H,
                  two = 1.96, mode = 1,
                  graph.sys = "base", silent = TRUE, sunflower = FALSE,
                  xlab="mean measurement", ylab="differences", conf.int = 0.95) 

# per plot custom vedi PAPER_FiguresLayout.R in /STORAGE_advex



# Boxplots ---------------------------------------------------------------------
# of height differences between computed at 12:00

PBL.H.abs.err <- PBL.H - MLH.1
PBL.H.rel.err <- PBL.H.abs.err/PBL.H *100
summary(PBL.H.rel.err)

windows()
# absolute error
boxplot(PBL.H.abs.err[MLH.TS.daytime] ~ MLH.TS.HHMM[MLH.TS.daytime])
# relativ eerror
boxplot(PBL.H.rel.err[MLH.TS.daytime] ~ MLH.TS.HHMM[MLH.TS.daytime],  outline = FALSE)


windows()
boxplot(PBL.H.abs.err[MLH.TS.daytime] ~ MLH.TS.HHMM[MLH.TS.daytime] * MLH.TS.seas[MLH.TS.daytime], 
        col=rep(seas.col, each = 22), ylim=c(-4000, 3000), outline = FALSE)
legend('bottomleft', legend=unique(MLH.1.y.s.dc$Group.2), col=seas.col, pch=15, bty = 'n')
abline(h=0)


