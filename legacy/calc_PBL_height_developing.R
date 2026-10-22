

# TODO: plan to use a flagging vector 
# TODO: findo a better smoothing approach (the current is fine but it removes quite some data in gaps proximity)



calc_PBL_height <- function(TIMESTAMP,         # timestamp [%Y-%m-%d %H:%M:%S]
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
                            # processing options
                            GF.PBL = TRUE,    # fill 1 hour gaps in PBL heights by linear interpolation 
                            smooth.PBL = TRUE,# Smooth PBL heights ONLY IF temporal resolution is less than 1 hour
                            plot.atmc = NULL, # set TRUE to plot atmospheric conditions frequency
                            plot.PBL = NULL,  # set TRUE to plot PBL height estimations
                            save.PBL = NULL)  # set TRUE to save result (tabbed ASCII)
{
  
  # Parte usata per debugging: quando o via***************************************
  
  # TIMESTAMP <- flux.turb.TS         # timestamp [%Y-%m-%d %H:%M:%S]
  # LAT <- LAT               # point latitude [decimal degrees]
  # zm <- zm.                # measurement height [m]
  # d <- d.                 # displacement height [m]
  # MOL <- MOL.24               # Monin-Obukov length [m]
  # ustar <- ustar.24             # friction velocity [ms-2]
  # H <- H.24                 # sensible heat flux [Wm-2]
  # AT <- AT.24                # air tmperature [K]
  # AP <- AP.24         # air pressure [Pa]
  # rho <- rho        # density of dry air [kg m-3]
  # cp <- cp
  
  #********************************************************************************
  
  
  # Libraries ..................................................................
  kpacks <- c(
    'zoo',
    'data.table',
    'bigleaf',
    'crayon'
  )   
  new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
  if(length(new.packs)){
    install.packages(new.packs)
  }
  lapply(kpacks, require, character.only=T)
  remove(kpacks, new.packs)
  # ............................................................................
  
  
  # Timestamp formatting
  TIMESTAMP <- strptime(TIMESTAMP, '%Y-%m-%d %H:%M:%S', tz='GMT')
  # check for time continuity: there must be a continuous time-series
  if(any(diff(summary(diff(as.double(TIMESTAMP)))) != 0)){
    cat(yellow(bold('WARNING:') %+% ' There are jumps in the time-line! \n you have to provide a continuous timeseries (unse NAs for missing values)\n\n'))
    stop(call. = FALSE, domain = NA)
  }

  
  
  # Derivation of the boundary layer height ************************************
  # Based on formulation reported in Kljun et al 2015 (Appendix B) and other references
  # Kljun, N., Calanca, P., Rotach, M. W. and Schmid, H. P.: A simple two-dimensional parameterisation for Flux Footprint Predictions FFP, Geosci. Model Dev., 8, 3695–3713, doi:10.5194/gmdd-8-6757-2015, 2015.
  
  
  # Constants and variables ....................................................
  
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
  
  # Air pressure [kPa]
  if(is.null(AP)){
    if(!is.null(altit)){ # If not provided as input, calculate it according to Campbell and Norman (1998)
      AP <- 101300 * exp(-(altit+zm)/8200)
    }else{ # but if neither site altitude is not provided
      cat(yellow(bold('WARNING:') %+% ' one among air pressure and site altitude must be provided\n\n'))
      stop()
    }
  }
  
  # Density of dry air [kg m-3]
  # If not provided as input, calculate it according to the ideal gas law (https://en.wikipedia.org/wiki/Density_of_air)
  if(is.null(rho)){                        
    rho <- AP / (R. * AT)
  }
  
  # Dry air heat capacity at constant pressure [J kg-1 K-1]
  # If not provided as input, calculate it according to (https://www.licor.com/env/help/eddypro/topics_eddypro/Calculating_Micromet_Variables.html)
  if(is.null(cp)){                         
    cp <- 1005 + (((AT-273.15)+23.12)^2 / 3364) 
  }
  
  # Surface (kinematic) turbulent flux of sensible heat [K m s-1]
  wTheta <- H / (rho * cp)                
  # ............................................................................
  
  
  # Atmospheric conditions and PBL height calculation **************************
  
  # Receiving data structures
  atmos.cond <- c(NULL)                  # atmospheric conditions
  PBL.h.0 <- c(NULL)                     # original boundary layer heights
  
  i <- 0
  for(i in 1:length(TIMESTAMP)){
    
    if(!is.na(MOL[i]) & !is.na(ustar[i])){ # if neither MOL or u* is NOT NA, estimate PBL, otherwise set it to NA (missing data)
      
      ## Atmospheric conditions ................................................
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
      }else if(1/MOL[i] >= -0.056 & 1/MOL[i] < -0.016){ 
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
      
      
      # PBL height .............................................................
      
      ## Strongly unstable and unstable conditions (convective, daytime)
      # set to NA as will be estimated later using the Batchvarova and Gryning (1991)
      if(atmos.cond[i] == 'strongly-unstable' || atmos.cond[i] == 'unstable'){ 
        PBL.h.0[i] <- NA 
      }
      
      
      ## Near-neutral conditions 
      # (Hanna and Chang, 1993; Sugiyama and Nasstrom, 1999; Hessa et al., 2006)
      if(atmos.cond[i] == 'near-neutral'){ # & abs(ustar[i] / abs(f.cor) * MOL[i]) < 4
        PBL.h.0[i] <- 0.3 * ustar[i] / abs(f.cor)
      }
      
      ## Neutral conditions
      # (Hanna and Chang, 1993; Hessa et al., 2006)
      if(atmos.cond[i] == 'neutral'){
        PBL.h.0[i] <- 0.3 * ustar[i] / abs(f.cor) 
        # # If the above measurement is not applicable [Nieuwstadt (1981)]
        # if(is.na(PBL.h.0[i])){
        #   PBL.h.0[i] <- MOL[i]/3.8 * (-1 + (1 + 2.28 * ustar[i] / (f.cor * MOL[i]))^(1/2))
        # }
      }
      
      ## Stable conditions 
      # (Zilitinkevich et al., 2012; Sugiyama and Nasstrom, 1999; Hessa et al., 2006)
      if(atmos.cond[i] == 'stable'){
        PBL.h.0[i] <- 0.4 * (ustar[i] * MOL[i] / abs(f.cor))^(1/2)
        # # If the above measurement is not applicable [Nieuwstadt (1981)]
        # if(is.na(PBL.h.0[i])){
        #   PBL.h.0[i] <- MOL[i]/3.8 * (-1 + (1 + 2.28 * ustar[i] / (f.cor * MOL[i]))^(1/2))
        # }
      }
      
      ## Strongly stable
      # Zilitinkevich et al., 2012)
      if(atmos.cond[i] == 'strongly-stable' || atmos.cond[i] == 'nighttime-stable'){
        PBL.h.0[i] <- 0.4 * (ustar[i] * MOL[i] / abs(f.cor))^(1/2)
      }
      
    }else{
      
      PBL.h.0[i] <- NA
      
    }
  }
  
  
  ## interpolate NA values (MAX 1 hour gap)
  if(GF.PBL){
    PBL.h.0 <- na.approx(PBL.h.0, maxgap = 2, na.rm = FALSE)
  }
  
  
  ## SLAB model: Batchvarova and Gryning (1991) ********************************
  # For convective conditions, the boundary layer height cannot be diagnosed due 
  # to the nearly symmetric diurnal cycle of the surface heat flux. It must therefore 
  # be integrated employing a prognostic expression, and starting at sunrise 
  # (before the surface heat flux first becomes positive) when the initial height 
  # is diagnosed using one of the above expressions (Kljun et al., 2015)
  
  # use all input variables from the actual time i and the PBL height of time i-1 
  # to calculate the dH for time i, then add it to PBLH(i-1)
  # and you need to take the dt into account.
  
  dh <- c(NULL) #  rate of change
  PBL.h.BG <- c(NULL) # modelled PBLH
  i <- 1
  for(i in 1:length(TIMESTAMP)){
    
    # if(!is.na(H[i]) & H[i] > 5){ # modell it in case of daytime conditions
    if(!is.na(H[i]) & H[i] > 0){ # modell it in case of daytime conditions
      
      if(!is.na(PBL.h.0[i-1])){ # if the parameterized PBLH in the previous half-hour is NOT NA, use it
        
        # B&G model
        dh[i] <- dt. * (wTheta[i]/gamma.) /
          ((PBL.h.0[i-1]^2 / ((1+2*A.)*PBL.h.0[i-1] - 2*B.*k.*MOL[i])) + ((C.*AT[i]*ustar[i]^2) / (gamma.*g.*((1+A.)*PBL.h.0[i-1] - B.*k.*MOL[i]))))
        
        PBL.h.BG[i] <- PBL.h.0[i-1] + dh[i]
        
      }else{ # otherwise, use the modelled value
        
        dh[i] <- dt. * (wTheta[i]/gamma.) /
          ((PBL.h.BG[i-1]^2 / ((1+2*A.)*PBL.h.BG[i-1] - 2*B.*k.*MOL[i])) + ((C.*AT[i]*ustar[i]^2) / (gamma.*g.*((1+A.)*PBL.h.BG[i-1] - B.*k.*MOL[i]))))
        
        PBL.h.BG[i] <- PBL.h.BG[i-1] + dh[i]
      }
      
    }else{
      
      dh[i] <- NA
      PBL.h.BG[i] <- NA
      
    }
  }
  
  # interpolate modelled PBL NA values (MAX 1 hour gap)
  if(GF.PBL){
    PBL.h.BG <- na.approx(PBL.h.BG, maxgap=2, na.rm = FALSE)
  }
  
  
  ### Replace convective conditions values of PBL.h.0 (i.e. NA) with modeled values
  PBL.h.1 <- replace(PBL.h.0, is.na(PBL.h.0), PBL.h.BG[is.na(PBL.h.0)])
  
  
  ### Smooth PBL heights ONLY IF temporal resolution is less than 1 hour
  # Smoothed symmetrically: average of current sample, 4 future samples, and 4 past samples, i.e. average over 2 hours windows
  if(smooth.PBL){
    PBL.h <- filter(PBL.h.1, filter=rep(1/5, 5), sides = 2)
  }
  
  
  # Rendering ..................................................................
  
  ## atmospheric conditions
  if(isTRUE(plot.atmc)){
    windows()
    barplot(prop.table(table(atmos.cond)), ylim=c(0,1), las=1)
  }
  
  ## PBL heigth
  if(isTRUE(plot.PBL)){
    if(mean(diff(as.double(TIMESTAMP)), na.rm=T) <= 3600){
      windows()
      plot(PBL.h.0, pch=19, type='o', ylab='PBL height [m]', xlab='time\n [progressive half-hours]')
      # abline(v=seq_along(PBL.h.0), col=8, lty=32)
      points(PBL.h, col = 2, type='l', pch=1)
      legend('topleft', c('original PBL H', 'smoothed/interpolated PBL H'), col=c(1,4), pch=c(19,1), bty='n')
    }else{
      windows()
      plot(PBL.h.0, pch=19, type='p', ylab='PBL height [m]', xlab='time\n [progressive half-hours]')
      legend('topleft', c('original PBL H'), col=1, pch=19, bty='n')
    }
    
  }
  
  if(isTRUE(save.PBL)){
    # save output data
    write.table(data.frame(TIMESTAMP = as.character(TIMESTAMP),
                           atmos_cond = atmos.cond,
                           PBL_h_0 = round(PBL.h.0, 2),
                           PBL_h = round(PBL.h, 2)),
                'PBL_heights.txt', quote=F, row.names = F, sep='\t')
  }
  
  return(PBL.h)
  
}

# rm(TIMESTAMP,LAT,zm, d, MOL, H, AT, AP, rho, cp)
