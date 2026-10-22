


# # Parte usata per debugging: quando ok via***************************************
# 
# TIMESTAMP <- set.cleaned.timestamp         # timestamp [%Y-%m-%d %H:%M:%S]
# LAT <- EC.coord$'lat'              # point latitude [decimal degrees]
# zm <- FFP.input$'zm.'                # measurement height [m]
# d <- FFP.input$'d'                 # displacement height [m]
# MOL <- FFP.input$'ol'               # Monin-Obukov length [m]
# ustar <- FFP.input$'ustar'             # friction velocity [ms-2]
# H <- set.cleaned$'H_cleaned'                  # sensible heat flux [Wm-2]
# AT <- set.cleaned$'air_temperature'                # air tmperature [K]
# AP <- set.cleaned$'air_pressure'         # air pressure [Pa]
# rho <- set.cleaned$'air_density'        # density of dry air [kg m-3]
# cp <- set.cleaned$'air_heat_capacity'
# 
# # rm(TIMESTAMP, LAT,zm,d,MOL,ustar,H,AT,AP,rho,cp)
# 
#********************************************************************************


calc_PBL_height <- function(PBL.dir = NULL,    # directory of PBL estimations I/O
                            TIMESTAMP,         # timestamp [%Y-%m-%d %H:%M:%S], either POSIXlt or character
                            LAT,               # point latitude [decimal degrees]
                            zm,                # measurement height [m]
                            d,                 # displacement height [m]
                            MOL,               # Monin-Obukov length [m]
                            ustar,             # friction velocity [ms-2]
                            H,                 # sensible heat flux [Wm-2]
                            AT,                # air tmperature [K]
                            AP = NULL,         # air pressure [Pa]
                            rho = NULL,        # density of dry air [kg m-3]
                            cp = NULL,         # air_heat_capacity [J kg-1 K-1]
                            altit = NULL,      # site altitude [m]
                            plot.atmc = NULL,  # set TRUE to plot atmospheric conditions frequency
                            plot.PBL = NULL,   # set TRUE to plot PBL height estimations
                            save.PBL = NULL    # set TRUE to save result (tabbed ASCII)
){
  
  
  # possibly load key packages
  kpacks <- c('zoo')   # for FFP Climatology
  new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
  if(length(new.packs)){install.packages(new.packs)}
  lapply(kpacks, require, character.only=T)
  remove(kpacks, new.packs)
  
  # Timestamp formatting
  if (class(TIMESTAMP)[1] != "POSIXlt") {
    TIMESTAMP <- strptime(TIMESTAMP, '%Y-%m-%d %H:%M:%S', tz='GMT')
  }
  
  # check for time continuity
  if (any(diff(summary(diff(as.double(TIMESTAMP)))) != 0)) {
    cat('\n> summary(diff(as.double(TIMESTAMP))):\n')
    cat(summary(diff(as.double(TIMESTAMP))))
    stop('[ETPPS::calc_PBL_height]: There are jumps in the time-line! \n\n', call. = FALSE)
  }
  
  
  ### Derivation of the boundary layer height
  # Based on Kljun et al 2015 (Appendix B) and other references
  
  ## Constants and variables
  omega. <- 2 * pi / (24*60*60)           # equatorial angular speed of Earth's rotation [rad s-1]
  f.cor <- 2 * omega. * sin(LAT * pi/180) # Coriolis parameter [s-1]
  zdl <- (zm - d) / MOL                   # stability paramenter
  
  ## B&G model parameter 
  gamma. <- 0.0025                        # gradient of potential temperature above the convective boundary layer (normally approximated to ≈ 0.01 Km−1)
  k. <- 0.41                              # von Karman constant
  g. <- 9.80665                           # earth-surface gravitational acceleration [m s-1]
  dt. <- median(diff(as.double(TIMESTAMP))) # time steps [s]
  A. <- 0.2
  B. <- 2.5
  C. <- 8
  R. <- 287.058                           # gas constant for dry air [J kg-1 K-1]
  
  if(is.null(AP)){                        # Air pressure. If not provided as input, calculate it according to Campbell and Norman (1998)
    if(!is.null(altit)){                  # but if neither site altitude is not provided
      AP <- 101300 * exp(-(altit+zm)/8200)
    }else{
      warning('One among air pressure or site altitude must be provided')
      break()
    }
  }
  
  if(is.null(rho)){                        # Density of dry air [kg m-3]. If not provided as input, calculate it according to the ideal gas law (https://en.wikipedia.org/wiki/Density_of_air)
    rho <- AP / (R. * AT)
  }
  
  if(is.null(cp)){                         # Dry air heat capacity at constant pressure [J kg-1 K-1]. If not provided as input, calculate it according to (https://www.licor.com/env/help/eddypro/topics_eddypro/Calculating_Micromet_Variables.html)
    cp <- 1005 + (((AT-273.15)+23.12)^2 / 3364) 
  }
  
  wTheta <- H / (rho * cp)                 # Surface (kinematic) turbulent flux of sensible heat [K m s-1]
  
  
  # Atmospheric conditions and PBL height formula calculation
  
  # Receiving data structures
  atmos.cond <- c(NULL)                  # atmospheric conditions
  PBL.h.0 <- c(NULL)                     # original boundary layer height
  
  i <- 0
  for(i in 1:length(TIMESTAMP)){
    
    if(!is.na(1/MOL[i]) | !is.na(1/ustar[i])){ # if MOL is NOT NA, calculate it, otherwise set it to NA (missing data)
      
      
      ## Current Atmospheric condition
      # Pasquill-Gifford Stability Class Correlated with Monin-Obukhov Length: [Essa et al., 2006]
      # A: strongly unstable (convective)    1/L < -0.056
      # B: unstable                          -0.056 <= 1/L < −0.016
      # C: slightly unstable (near neutral)  −0.016 <= 1/L < −0.004
      # D: neutral                           −0.004 <= 1/L < 0.002
      # E: stable                            0.002 <= 1/L < 0.006
      # F: strongly stable                   0.006 <= 1/L < 0.022
      # G: nighttime stable conditions       1/L >= 0.022
      
      if(1/MOL[i] < -0.056){
        atmos.cond[i] <- 'strongly-unstable'
      }else if(1/MOL[i] >= -0.056 & 1/MOL[i] < -0.003){ # Note that this value has been adjusted to -0.003 from -0.016
        atmos.cond[i] <- 'unstable'
      }else if(1/MOL[i] >= -0.016 & 1/MOL[i] < -0.004){
        atmos.cond[i] <- 'near-neutral'
      }else if(1/MOL[i] >= -0.004 & 1/MOL[i] < 0.002){
        atmos.cond[i] <- 'neutral'
      }else if(1/MOL[i] >= 0.002 & 1/MOL[i] < 0.006){
        atmos.cond[i] <- 'stable'
      }else if(1/MOL[i] >= 0.006 & 1/MOL[i] < 0.022){
        atmos.cond[i] <- 'strongly-stable'
      }else{
        atmos.cond[i] <- 'nighttime-stable'
      }
      
      
      ## Strongly unstable and unstable conditions (convective, daytime)
      if(atmos.cond[i] == 'strongly-unstable' || atmos.cond[i] == 'unstable'){ 
        PBL.h.0[i] <- NA # will be estimated later using the Batchvarova and Gryning (1991)
      }
      
      
      ## Near-neutral conditions [Hanna and Chang (1993); Sugiyama and Nasstrom (1999); Hessa et al., 2006]
      if(atmos.cond[i] == 'near-neutral'){ # & abs(ustar[i] / abs(f.cor) * MOL[i]) < 4
        PBL.h.0[i] <- 0.3 * ustar[i] / abs(f.cor)
      }
      
      ## Neutral conditions
      if(atmos.cond[i] == 'neutral'){
        PBL.h.0[i] <- 0.3 * ustar[i] / abs(f.cor) # Hanna and Chang (1993); Hessa et al., 2006
        # If the above measurement is not applicable [Nieuwstadt (1981)]
        if(is.na(PBL.h.0[i])){
          PBL.h.0[i] <- MOL[i]/3.8 * (-1 + (1 + 2.28 * ustar[i] / (f.cor * MOL[i]))^(1/2))
        }
      }
      
      ## Stable conditions [Zilitinkevich et al. (2012); Sugiyama and Nasstrom (1999); Hessa et al., 2006]
      if(atmos.cond[i] == 'stable'){
        PBL.h.0[i] <- 0.4 * (ustar[i] * MOL[i] / abs(f.cor))^(1/2)
        # If the above measurement is not applicable [Nieuwstadt (1981)]
        if(is.na(PBL.h.0[i])){
          PBL.h.0[i] <- MOL[i]/3.8 * (-1 + (1 + 2.28 * ustar[i] / (f.cor * MOL[i]))^(1/2))
        }
      }
      
      ## Strongly stable [Zilitinkevich et al (2012)]
      if(atmos.cond[i] == 'strongly-stable' || atmos.cond[i] == 'nighttime-stable'){
        PBL.h.0[i] <- 0.4 * (ustar[i] * MOL[i] / abs(f.cor))^(1/2)
      }
      
    }else{
      
      PBL.h.0[i] <- NA
      
    }
  }
  # windows()
  # plot(PBL.h.0[1:500], type='o')
  
  # interpolate NA values (MAX 1 hour gap)
  PBL.h.0.int <- na.approx(PBL.h.0, maxgap=2, na.rm = FALSE)
  # points(PBL.h.0.int, type='o', col='darkgrey', pch=19)
  # points(PBL.h.0, type='o')
  
  
  ## SLAB model: Batchvarova and Gryning (1991)
  # For convective conditions, the boundary layer height cannot be diagnosed due 
  # to the nearly symmetric diurnal cycle of the surface heat flux. It must therefore 
  # be integrated employing a prognostic expression, and starting at sunrise 
  # (before the surface heat flux first becomes positive) when the initial height 
  # is diagnosed using one of the above expressions (Kljun et al., 2015)
  
  # use all input variables from the actual time i and the pbl height of time i-1 
  # to calculate the dh for time i, then add it to h(i-1)
  # and you need to take the dt into account.
  
  # determine the PBL.h rate of change and consequently the PBL height
  dh <- c(NULL)
  PBL.h.BG <- c(NULL)
  i <- 0
  for(i in 1:length(TIMESTAMP)){
    
    if(!is.na(H[i]) & H[i] > 0 & i != 1){
      
      if(!is.na(PBL.h.0.int[i-1])){
        
        # B&G model
        dh[i] <- dt. * (wTheta[i]/gamma.) / 
          ((PBL.h.0.int[i-1]^2 / ((1+2*A.)*PBL.h.0.int[i-1] - 2*B.*k.*MOL[i])) + ((C.*AT[i]*ustar[i]^2) / (gamma.*g.*((1+A.)*PBL.h.0.int[i-1] - B.*k.*MOL[i]))))
        
        PBL.h.BG[i] <- PBL.h.0.int[i-1] + dh[i]
        
      }else{
        
        dh[i] <- dt. * (wTheta[i]/gamma.) /
          ((PBL.h.BG[i-1]^2 / ((1+2*A.)*PBL.h.BG[i-1] - 2*B.*k.*MOL[i])) + ((C.*AT[i]*ustar[i]^2) / (gamma.*g.*((1+A.)*PBL.h.BG[i-1] - B.*k.*MOL[i]))))
        
        PBL.h.BG[i] <- PBL.h.BG[i-1] + dh[i]
      }
      
    }else{
      
      dh[i] <- NA
      PBL.h.BG[i] <- NA
      
    }
    
  }
  # points(PBL.h.BG, type='o', col=2, pch=19)
  # SEMBRA SHIFTATO IN AVANTI DI MEZZ'ORA, CAPIRE PERCHE
  
  # interpolate NA values (MAX 1 hour gap)
  PBL.h.BG.int <- na.approx(PBL.h.BG, maxgap=2, na.rm = FALSE)
  # points(PBL.h.BG.int, type='o', col=2)
  
  ### Replace convective conditions values of PBL.h.0 (i.e. NA) with modeled values
  PBL.h.1 <- replace(PBL.h.0.int, is.na(PBL.h.0.int), PBL.h.BG.int[is.na(PBL.h.0.int)])
  # points(PBL.h.1, type='o', col=4)
  
  ### Smooth PBL heights ONLY IF temporal resolution is less than 1 hour
  # Smoothed symmetrically: average of current sample, 4 future samples, and 4 past samples, i.e. average over 2 hours windows
  PBL.h <- stats::filter(PBL.h.1, 
                         filter = rep(1/5, 5),    # a vector of filter coefficients in reverse time order
                         method = 'convolution',  # if "convolution" a moving average is used: if "recursive" an autoregression is used
                         sides = 2)               # If sides = 1 the filter coefficients are for past values only; if sides = 2 they are centred around lag 0. In this case the length of the filter should be odd
  # points(PBL.h, type='l', col=4, lwd=2)
  
  
  # rendering
  
  # atmospheric conditions
  if (isTRUE(plot.atmc)) {
    if (is.null(PBL.dir)) {
      windows()
      barplot(prop.table(table(atmos.cond)), ylim=c(0,1), las=1)
    } else {
      png(paste(PBL.dir,'Atmospheric_conditions.png', sep=''), width = 14, height = 14, unit = 'cm', res=300)  
      barplot(prop.table(table(atmos.cond)), ylim=c(0,1), las=1)
      dev.off()
    }
  }
  
  # PBL heigth
  if (isTRUE(plot.PBL)) {
    
    if (is.null(PBL.dir)) {
      windows()
    } else {
      png(paste(PBL.dir,'PBL_heights.png', sep=''), width = 14, height = 14, unit = 'cm', res=300)  
    }
    
    if (mean(diff(as.double(TIMESTAMP)), na.rm=T) <= 3600) {
      plot(PBL.h.0, pch=19, type='o', ylab='PBL height [m]', xlab='time\n [progressive half-hours]', ylim=c(0,max(PBL.h.0, na.rm=T)))
      points(PBL.h, col = 4, type='l', pch=1)
      legend('topleft', c('original PBL H', 'modelled/interpolated/smoothed PBL H'), col=c(1,4), pch=c(19,NA), bty='n')
    } else {
      plot(PBL.h.0, pch=19, type='p', ylab='PBL height [m]', xlab='time\n [progressive half-hours]')
      legend('topleft', c('original PBL H'), col=1, pch=19, bty='n')
    }
    
    if (!is.null(PBL.dir)) { dev.off() }
    
  }
  
  if (isTRUE(save.PBL)) {
    
    # save output data
    write.table(data.frame(TIMESTAMP = as.character(TIMESTAMP),
                           atmos_cond = atmos.cond,
                           PBL_h_0 = round(PBL.h.0, 2),
                           PBL_h = round(PBL.h, 2)),
                paste0(PBL.dir, 'PBL_heights.txt'), 
                quote=F, row.names = F, sep='\t')
  }
  
  return(PBL.h)
  
}

# rm(TIMESTAMP,LAT,zm, d, MOL, H, AT, AP, rho, cp)
