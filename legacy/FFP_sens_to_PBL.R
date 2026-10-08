
# ******************************************************************************
# Footprint model sensitivity to Planetary Boundary Layer height 
# 
# ******************************************************************************

# TODO:
# fare analisi per matrice tipo 4x4 con 4 scenari di MOL (stabilità) e 4 di z0 (o vedere quali sono i parametri più influenti)


# Libraries and Functions ......................................................

kpacks <- c(
  'leaflet',
  'itsadug'
)   
new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
if(length(new.packs)){
  install.packages(new.packs)
}
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)


# Single footprint
source("C:/WORK/R/PROJECTS/FOOTPRINT/calc_footprint_FFP.R")
# Footprint climatology
source("C:/WORK/R/PROJECTS/FOOTPRINT/calc_footprint_FFP_climatology.R")

# ..............................................................................



# Analyses settings ------------------------------------------------------------

PBL.set <- seq(50, 500, 20) # Vector of boundary layer height set [m]
sap.sens <- c(80) # Percentage of source area
zm.sens <- 10 # Measurement height above displacement height (i.e. z‐d) [m]
z0.sens <- 0.02 # Roughness length [m]
umean.sens <- NA # Mean wind speed at zm [ms‐1]
ol.sens <- -500 # Obukhov length [m]
sigmav.sens <- 0.3 # standard deviation of lateral velocity fluctuations [ms‐1]
ustar.sens <- 0.3 # friction velocity [ms‐1]
wind_dir.sens <- 0 # wind direction in degrees (of 360) for rotation of the footprint

# SET the PBL heights on which base the sensitivity estimation
# TODO: aggiustare il valore mediano, scegliendo un punto che discrimina la sensitività al meglio
PBL.checkpoints <- c(min(PBL.set), 150, max(PBL.set)) 

# plot settings
# leaflet::previewColors(colorNumeric("Blues", domain = PBL.set), PBL.set)
# leaflet::previewColors(colorNumeric("RdYlBu", domain = PBL.set), PBL.set)
FFP.col <- colorNumeric("RdYlBu", domain = PBL.set)


# Sensitivity analyses on 2D Footprint climatology =============================
FFP.pbl.sens.Fmax <- c(NULL)
i <- 1
for(i in 1:length(PBL.set)){
  
  FFP.pbl.sens <- calc_footprint_FFP_climatology(zm = 10,
                                                 z0 = z0.sens, 
                                                 umean = umean.sens, 
                                                 h = PBL.set[i],
                                                 ol = ol.sens,
                                                 sigmav = sigmav.sens, 
                                                 ustar = ustar.sens, 
                                                 wind_dir = wind_dir.sens,
                                                 r = sap.sens,
                                                 domain = c(-100,1000,-100,1000), nx = 1100, smooth_data=1
                                                 )
  
  # store max distance (Fetch) at 2 given PBL heights
  isop.indx <- which(sap.sens == 80) # internal index to spot out the 80% contribution isoplethe
  FFP.pbl.sens.Fmax[i] <- max(FFP.pbl.sens$yr[[isop.indx]])
  
  
  # plot two‐dimensional footprint with contour lines of 80%
  if(i==1){
    windows()
    plot(FFP.pbl.sens$xr[[isop.indx]], FFP.pbl.sens$yr[[isop.indx]], type="l", col=FFP.col(PBL.set[i]), las=1, xlab = 'X [m]', ylab = 'Y [m]')
    points(0,0,pch=3)
  }else{
    lines(FFP.pbl.sens$xr[[isop.indx]], FFP.pbl.sens$yr[[isop.indx]], type="l", col=FFP.col(PBL.set[i]))
  }
  
  # add a reference line @ PBL checkpoints and its max H
  if(PBL.set[i] %in% PBL.checkpoints){
    abline(h = FFP.pbl.sens.Fmax[i], lty=32, col=FFP.col(PBL.set[i]))
    # TODO: add fetch value and corresponding PBL above lines (?)
  }
  
}
# legend
gradientLegend(valRange=PBL.set,
               color = FFP.col(PBL.set),
               border.col = NA,
               pos=c(.02,.75,.045,.98), # xleft, ybottom, xright, ytop of a rectangle
               side=4)


# footprint sensitivity to PBL height (fetch gradient) -------------------------

PBL.fetch <- as.data.frame(cbind(PBL.set, FFP.pbl.sens.Fmax))

windows()
plot(PBL.fetch$PBL.set, PBL.fetch$FFP.pbl.sens.Fmax, col=FFP.col(PBL.set), type='p', pch=19, las=1, xlab = 'PBL height [m]', ylab = '80% contribution fetch [m]')

# sensitivity for shallow and higher PBL (basing on heights decided for the checkpoints)

sens.indx <- which(PBL.fetch$PBL.set %in% c(PBL.checkpoints, max(PBL.set)))
# average fetch decrease (m) for 1 meter of PBL heigth increase 
FFP.sens.2PBL <- abs(diff(PBL.fetch$FFP.pbl.sens.Fmax[sens.indx]) / diff(PBL.fetch$PBL.set[sens.indx]))
FFP.sens.2PBL





