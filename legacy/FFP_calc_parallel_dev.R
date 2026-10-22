
# NOTE_20251104:
# Codice aggiornato per rimpiazzare rgdal (dismesso) con funzioni "spatial" più moderne e efficienti
# 
# Modifiche principali:
# 
# rgdal::ogrListLayers() → sf::st_layers()
# rgdal::readOGR() / readOGR() → sf::st_read()
# sp::spTransform() → sf::st_transform()
# SpatialPoints(...) → sf::st_as_sf(...) (con coords=)
# rgeos::gBuffer(..., width=0) → sf::st_buffer(..., dist = 0) o sf::st_make_valid() per problemi topologici
# rgeos::gUnion() → sf::st_union()
# controllo validità: rgeos::gIsValid() → sf::st_is_valid() / sf::st_make_valid()
# plotting aggiornato per oggetti sf (uso plot(st_geometry(...), ...) e plot(..., add = TRUE))



FFPCalcParallel <- function(
    
  site.ID,
  input.data,  # a vector containing the following variables: TIMESTAMP [YYYY-mm-dd HH:MM:SS]
  #                                                           daytime (1 yes, 2 no)
  #                                                           qcf: QC flag of the record (FALSE = OK, TRUE = record to discard)
  #                                                           hm: measurement height [m]
  #                                                           d: displacement height [m]
  #                                                           zm: dynamic measurement height (hm - d) [m]
  #                                                           hc: canopy height [m]
  #                                                           z0: roughness length [m] 
  #                                                           umean: wind speed [ms-1]
  #                                                           ol: Monin-Obukov length [m]
  #                                                           sigmav: standard deviation ov lateral wind velocity
  #                                                           ustar: friction velocity [ms-1]
  #                                                           wind_dir: wind direction [°]
  #                                                           PBL: planetary boundary layer height [m]
  FFPSP.output.dir,         # Path to the output directory
  FFP.func.dir,             # Path to the directory containing the FFP function by Kljun et al. 2015
  FFP.R,                    # Percentage of source area, i.e. a value between 10% and 90%. * 2
  ffp.domain = 1000,        # FFP domain size as an array of (xmin xmax ymin ymax) [m]
  plot.FFP.cwi = FALSE,     # Plot the crosswind‐integrated footprint visualization
  SnW.thr = 2,            # Axis ratio threshold of the ellipses fit to the 70% contribution. Ratios above it imply data discard
  plot.ellipse.fit = FALSE, # Plot the ellipses fit to the 70% contribution
  plot.FFP.clim = FALSE,    # Plot the footprint climatology
  EC.tower.utm,             # EC tower coordinates as SpatialPoints object, in UTM
  plot.spat.FFP = FALSE,    # plot the 70% and 80% FFP contributions over site's spatial features
  TA.utm,                   # Target Area polygon as SpatialPolygonsDataFrame object, in UTM
  plot.TA.FFP = FALSE,      # plot the TA FFP intersection
  DO.EA.analysis = FALSE,   # whether run the analysis in EA or not
  EA.utm,                   # Exclusion Area(s) polygon as a single SpatialPolygonsDataFrame object, in UTM
  plot.EA.FFP = FALSE,      # plot the EA FFP intersection
  DO.LCT.analysis = TRUE,   # whether run the analysis on LCTs or not
  LCT.utm,                  # Land Cover Typologies(s) polygon as a list of SpatialPolygonsDataFrame objects, in UTM
  plot.LCT.FFP = FALSE      # plot the LCT FFP intersection
  
){
  
  # Libraries ....................................................................
  kpacks <- c(
    "sf",
    "terra",
    "ggplot2",
    "viridisLite",
    "EBImage",
    "crayon")
  new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
  if (length(new.packs)){install.packages(new.packs)} # install packages not already installed
  loaded.packs <- lapply(kpacks, require, character.only=T) # load them
  if(!all(unlist(loaded.packs))){cat(yellow $ bold(paste('WARNING: not able to load the following libraries:\n', paste(kpacks[which(unlist(loaded.packs)==F)], collapse=', '))))}
  remove(kpacks, new.packs)
  
  
  # define message themes
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  # Source FFP functions 
  # FFP (Kljun et al., 2015)
  source(paste0(FFP.func.dir, 'calc_footprint_FFP.R'))
  source(paste0(FFP.func.dir, 'repmat.R'))
  source(paste0(FFP.func.dir, 'calc_footprint_FFP_climatology.R'))
  # Fit ellipse to data points
  source(paste0(FFP.func.dir, 'EllipseFit.R'))
  
  # set colors for FFP CWI isoplethes
  FFP.cwi.col <- viridisLite::cividis(5)
  
  # set/create a directory for the individual FFP plot figures
  # Output directory with exemplary individual FFP maps
  FFP.dir.FFP <- paste0(FFPSP.output.dir, "individual_FFP/"); suppressWarnings(dir.create(FFP.dir.FFP))
  
  # set the current half-hour timestamp in ISO format
  cur.ts.iso <- format(strptime(input.data$'TIMESTAMP', '%Y-%m-%d %H:%M:%S', tz='GMT'), '%Y%m%d%H%M%S')
  if (is.na(cur.ts.iso)) { # midnight cases
    cur.ts.iso <- format(strptime(input.data$'TIMESTAMP', '%Y-%m-%d', tz='GMT'), '%Y%m%d%H%M%S')
  }
  
  
  if (input.data$'qcf' | any(is.na(input.data[c('umean', 'ol', 'sigmav', 'ustar', 'wind_dir')]))) { # IF FC IS QC FILTERED OR A FUNDAMENTAL VAR IS NA, DO NOT COMPUTE THE FFP
    
    # Footprint 
    FFP.cwi.flag <- NA # FFP CWI error (0/1)
    FFP.cwi.peak <- NA # FFP CWI peak location [m]
    FFP.cwi.offset <- NA # FFP CWI offset distance [m]
    FFP.cwi.cum10.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 10% cumulative contribution
    FFP.cwi.cum20.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 20% cumulative contribution
    FFP.cwi.cum30.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 30% cumulative contribution
    FFP.cwi.cum40.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 40% cumulative contribution
    FFP.cwi.cum50.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 50% cumulative contribution
    FFP.cwi.cum60.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 60% cumulative contribution
    FFP.cwi.cum70.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 70% cumulative contribution
    FFP.cwi.cum80.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 80% cumulative contribution
    FFP.cwi.cum90.dist <- NA # upwind distance at which the crosswind integrated FFP reaches the 90% cumulative contribution
    FFP.el.m <- NA # Minor semi-axis of the fitted ellepses to the FFP CWI 70% contribution isoplethe
    FFP.el.M <- NA # Major semi-axis of the fitted ellepses to the FFP CWI 70% contribution isoplethe
    FFP.el.M_m <- NA # Ellipse axis ratio
    SnW.FFP <- NA # Unrealistic FFP shape (named Short and Wide (SnW) FFP): 1 yes <- NA 0 no
    # Footprint climatology and spatial analysis
    FFP.flag.err <- NA # FFP CLIM error (0/1)
    FFP.70.iso.surf <- NA # Area of the 70% contribution isoplethe geometry in the units of the current projection >> [m2]
    FFP.80.iso.surf <- NA # Area of the 80% contribution isoplethe geometry in the units of the current projection >> [m2]
    FFP.TA.intersect.surf <- NA # Area of intersection between FFP and the Target Area [m2] (FFP_TA_SURF)
    FFP.TA.intersect.surf.perc <- NA # Area of intersection between FFP and the Target Area as percentage of the 80% contribution isoplethe area (FFP_TA_INTER_%SURF)
    FFP.TA.intersect.contr <- NA # FFP cumulative contribution of the Target Area (intersection) [m2 m-2 <- NA %] (FFP_TA_INTER_CONTIRIB) 
    FFP.EA.intersect.surf <- NA # Area of intersection between FFP and the Exclusion Area [m2] (FFP_EA_SURF)
    FFP.EA.intersect.surf.perc <- NA # # Area of intersection between FFP and the Exclusion Area as percentage of the 80% contribution isoplethe area (FFP_EA_INTER_%SURF)
    FFP.EA.intersect.contr <- NA # FFP cumulative contribution of the Exclusion Area (intersection) [m2 m-2 <- NA %] (FFP_EA_INTER_CONTIRIB)
    
    if (DO.LCT.analysis) {
      FFP.LCT.intersect.surf <- rep(NA, length(LCT.utm))       # Areas of intersections between FFP and the Land Cover Typology surfaces [m2]
      FFP.LCT.intersect.surf.perc <- rep(NA, length(LCT.utm))  # Area of intersection between FFP and the Land Cover Typology surfaces as percentage of the 80% contribution isoplethe area (FFP_LCTxx_INTER_%SURF)
      FFP.LCT.intersect.contr <- rep(NA, length(LCT.utm))    # FFP cumulative contribution of the Land Cover Typologies (intersections) [m2 m-2, %] (FFP_LCTxx_INTER_CONTIRIB)
    }
    
    
  } else { # IF FC QC IS GOOD AND ALL THE FUNDAMENTAL VARIABLES ARE PRESENT, COMPUTE FFP AND INTERSECTIONS
    
    cat(note(paste0('\nComputing footprint analysis for the: ', as.character(input.data$TIMESTAMP),' half-hour.')))
    
    
    # *************************************************************************
    # Cross-wind integrated footprint -----------------------------------------
    # NOTE: calculate single FFP flux footprint, in addition to the single footprint climatology below, 
    #       because 'calc_footprint_FFP' provide the following  3 additional output:
    #       . x_ci_max = x location of footprint peak (distance from measurement) [m]
    #       . x_ci = x array of crosswind integrated footprint [m]
    #       . f_ci = Footprint function values of crosswind integrated footprint [m‐1]
    
    FFP.cwi <- calc_footprint_FFP(
      
      # mandatory input 
      zm = input.data$'zm',                                          # Measurement height above displacement height (i.e. z‐d) [m]
      z0 = input.data$'z0',                                          # Roughness length [m] ‐ enter [NaN] if not known
      umean = input.data$'umean',                                    # Mean wind speed at zm [ms‐1] ‐ enter [NaN] if not known (Either z0 or umean is required.
      h = ifelse(is.na(input.data$'PBL'), 1000, input.data$'PBL'),   # Boundary layer height [m]
      ol = input.data$'ol',                                          # >> outputbukhov length [m]
      sigmav = input.data$'sigmav',                                  # Standard deviation of lateral velocity fluctuations [ms‐1]
      ustar = input.data$'ustar',                                   # Friction velocity [ms‐1]
      # outputptional input
      wind_dir = input.data$'wind_dir',                              # Wind direction in degrees (of 360) for rotation of the footprint
      rslayer = 1,                                                   # Calculate footprint even if zm within roughness sublayer: set rslayer = 1. Note that this only gives a rough estimate of the footprint as the model is not valid within the roughness sublayer. Default is 0 (i.e. no footprint for within RS). z0 is needed for estimation of the RS.
      r = FFP.R                                                      # Percentage of source area, i.e. a value between 10% and 90%.
    )
    
    # store the CWI error 
    FFP.cwi.flag <- FFP.cwi$'flag_err'    # >> output
    
    if (FFP.cwi$'flag_err' == 0) { # in case of NO errors (flag_err == 0) in FFP
      
      # Footprint peak location [m]: the distance from the anemometer in the direction from which the largest relative individual contribution to the flux originates (upwind distance from measurement point) 
      # FFP.cwi.peak <- FFP.cwi$'x_ci_max' # >> output: This output from the K model R function is wrong. The peak is now calculated 
      FFP.cwi.peak <- FFP.cwi$'x_ci'[birk::which.closest(FFP.cwi$'f_ci', max(FFP.cwi$'f_ci'))] # Correct computation [20221116]
      
      # >> offset distance [m]: the area between the anemometer and the offset provides up to 1% of the total flux [EP manual].
      # FFP.cwi.offset <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), 0.01)] # >> output
      FFP.cwi.offset <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.01)] # Correct computation [20221116]
      
      # Store the upwind distance at which the crosswind integrated FFP reaches the 10% to 90% cumulative contribution
      FFP.cwi.cum10.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.1)] # >> output
      FFP.cwi.cum20.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.2)] # >> output
      FFP.cwi.cum30.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.3)] # >> output
      FFP.cwi.cum40.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.4)] # >> output
      FFP.cwi.cum50.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.5)] # >> output
      FFP.cwi.cum60.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.6)] # >> output
      FFP.cwi.cum70.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.7)] # >> output
      FFP.cwi.cum80.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.8)] # >> output
      FFP.cwi.cum90.dist <- FFP.cwi$'x_ci'[birk::which.closest(cumsum(FFP.cwi$'f_ci'), sum(FFP.cwi$'f_ci') * 0.9)] # >> output
      
      # Calculate the area under the FFPcwi curve (integral) following the trapezoid rule: should be 1
      FFP.cwi.integral <- flux::auc(FFP.cwi$'x_ci', FFP.cwi$'f_ci', thresh = NULL, dens = 100)  # >> output
      
      
      # Fit (try to) an ellipse to the FFP 70% isoplethe and calculate its axes
      ffp.el.fit <- NULL
      FFP.cwi.2DEllFit.went.wrong <- ''
      FFP.cwi.2DEllFit.went.wrong <- tryCatch(
        
        # Fit an ellipses to the 70% cumulative contribution
        ffp.el.fit <- fit.ellipse(x = FFP.cwi$'xr'[[which(FFP.R == 70)]], 
                                  y = FFP.cwi$'yr'[[which(FFP.R == 70)]]),
        
        error = function(e) {e}
      )
      
      if (!inherits(FFP.cwi.2DEllFit.went.wrong, "error")) { # If the ellipse fit function converges
        
        ell.fit <- TRUE
        
        # minor semi-axis length (alongwind axis)
        FFP.el.m <- ffp.el.fit$'minor'    # >> output
        # major semi-axis length (cross-wind axis)
        FFP.el.M <- ffp.el.fit$'major'    # >> output
        # ellipse axis ratio
        FFP.el.M_m <- FFP.el.M / FFP.el.m    # >> output
        
        if (plot.ellipse.fit) {
          # utility function which takes a fitted ellipse and returns a matrix of vertices for plotting
          ffp.el <- NULL
          ffp.el <- get.ellipse(ffp.el.fit)
          
          png(paste(FFP.dir.FFP, site.ID, '_','Ellispses_fit_', cur.ts.iso, '.png', sep=''),
              width = 14, height = 14, unit = 'cm', res=300)
          # windows()
          i.ffpcwi <- 1
          for (i.ffpcwi in length(FFP.R):1) {
            if (i.ffpcwi == length(FFP.R)){
              plot(FFP.cwi$'xr'[[i.ffpcwi]], FFP.cwi$'yr'[[i.ffpcwi]], type="l", col=FFP.cwi.col[i.ffpcwi], lwd=2)
              points(0,0,pch=3)
            }else{
              lines(FFP.cwi$'xr'[[i.ffpcwi]], FFP.cwi$'yr'[[i.ffpcwi]], col=FFP.cwi.col[i.ffpcwi], lwd=2)
            }
          }
          # plot fit
          lines(ffp.el, col="red") 
          mtext(input.data$TIMESTAMP, 3) 
          legend('bottomleft', legend=paste(FFP.R, '%'), lty=1, col=FFP.cwi.col[length(FFP.R)], bty='n')
          legend('bottomright', legend='ellipses fit to 70%', lty=1, col=2, bty='n')
          
          dev.off()
        }
        
        # Footprints (data) with unrealistic shape, named Short and Wide (SnW) FFP
        SnW.FFP <- ifelse(FFP.el.M_m > SnW.thr, 1, 0)     # >> output
        # USARE QUESTA FLAG, UNA VOLTA SETTATA LA THRESHOLD (in base a analisi dei dati), 
        # COME CONDIZIONE PER ESEGUIRE O MENO I CALCOLI DI SEGUITO INSIEME ALLA
        # CONDIZIONE SETTATAT CON remove.SnW (T/F)
        
        rm(FFP.cwi.2DEllFit.went.wrong)
        
      } else {
        
        ell.fit <- FALSE
        
        FFP.el.m <- NA    # minor semi-axis length (alongwind axis)
        FFP.el.M <- NA    # major semi-axis length (cross-wind axis)
        FFP.el.M_m <- NA  # ellipse axis ratio
        SnW.FFP <- NA     # set to NA in case of unrealistic distances
        
        cat(note(paste('\nWARN: the crosswind integrated footprint estimation for the', as.character(input.data$TIMESTAMP),'half-hour')))
        cat(note(paste(' is likely wrongly estimated (unrealistic cumulative footprint distances). \nIt was not possible to fit an ellipse to the 70% contribution isoplethe.')))
        cat(note(paste('\nError message from the fit function:', FFP.cwi.2DEllFit.went.wrong$'message')))
        
      }
      
      
      ## Crosswind‐integrated footprint visualization
      if (plot.FFP.cwi & ell.fit) { 
        
        png(paste(FFP.dir.FFP, site.ID, '_','FFP_CWI_', cur.ts.iso, '.png', sep=''),
            width = 14, height = 14, unit = 'cm', res=300)
        # windows()
        par(mfrow=c(2,1))
        plot(FFP.cwi$'x_ci', FFP.cwi$'f_ci',
             axes = FALSE, frame = TRUE, type="l", col=1, xlim=c(0, 600), ylab='', xlab='')
        axis(2, las=1); mtext('cross-wind integrated FFP', 2, 3)
        # mark the peak location
        points(FFP.cwi$'x_ci_max', 0.05, col=adjustcolor(3, 0.5), pch=16, type='h')
        # points(FFP.cwi$'x_ci_max' - FFP.cwi.offset, 0.05, col=adjustcolor(3, 0.5), pch=16, type='h')
        legend('topright', 'peak distances', pch='|', col=3, bty='n')
        # cumulative
        plot(cumsum(FFP.cwi$'f_ci'),
             type="l", col=1, xlim=c(0, 600), ylim=c(0, 2), ylab='', xlab='', las=1)
        mtext('upwind distance [m]', 1, 2.5)
        mtext('cumulative\ncross-wind integrated FFP', 2, 3)
        # mark the distance at which FFP reaches the 70% cumulative distribution
        abline(v = FFP.cwi.cum70.dist, lwd=3, col=adjustcolor(2, 0.3))
        
        dev.off()
      }
      
      
    } else { # in case of error in FFP 1D
      
      FFP.cwi.peak <- NA
      FFP.cwi.offset <- NA
      FFP.cwi.cum10.dist <- NA
      FFP.cwi.cum20.dist <- NA
      FFP.cwi.cum30.dist <- NA
      FFP.cwi.cum40.dist <- NA
      FFP.cwi.cum50.dist <- NA
      FFP.cwi.cum60.dist <- NA
      FFP.cwi.cum70.dist <- NA
      FFP.cwi.cum80.dist <- NA
      FFP.cwi.cum90.dist <- NA
      FFP.cwi.integral <- NA
      
      FFP.el.m <- NA
      FFP.el.M <- NA
      FFP.el.M_m <- NA
      SnW.FFP <- NA
    }
    
    # Free memory from FFP.cwi
    rm(FFP.cwi)
    
    
    # *************************************************************************
    # Footprint climatology ---------------------------------------------------
    # *1 Footprint will be calculated for a measurement at [0 0 zm] m. Default is smallest area including the r% footprint or (‐1000, 1000, ‐1000, 1000) m,
    #    whichever smallest (80% footprint if r not given)
    # *2 Can be either a single value (e.g., "80") or an array of increasing percentage values (e.g., seq(10, 80, 10))
    #    Expressed either in percentages ("80") or in fractions of 1 ("0.8") Default is [10:10:80]. Set to "NaN" for no output of percentages
    # *3 Note that this only gives a rough estimate of the footprint as the model is not valid within the roughness sublayer. 
    #    Default is 0 (i.e. no footprint for within RS). z0 is needed for estimation of the RS.
    
    # FFP domain array
    ffp.domain <- ffp.domain * c(-1,1,-1,1)
    
    FFP <- NULL
    FFP <- calc_footprint_FFP_climatology(
      
      # mandatory input 
      zm = input.data$'zm',                  # Measurement height above displacement height (i.e. z‐d) [m]
      z0 = input.data$'z0',                  # Roughness length [m] ‐ enter [NaN] if not known
      umean = input.data$'umean',            # Mean wind speed at zm [ms‐1] ‐ enter [NaN] if not known (Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint)
      h = input.data$'PBL',                  # Vector of boundary layer height [m]
      ol = input.data$'ol',                  # Vector of Obukhov length [m]
      sigmav = input.data$'sigmav',          # Vector of standard deviation of lateral velocity fluctuations [ms‐1]
      ustar = input.data$'ustar',            # Vector of friction velocity [ms‐1]
      wind_dir = input.data$'wind_dir',      # Vector of wind direction in degrees (of 360) for rotation of the footprint
      
      # optional input
      domain = ffp.domain,                   # Domain size as an array of (xmin xmax ymin ymax) [m] *1.
      dx = 1, dy = 1,                        # Cell size of domain [m] (default is dx = dy = 2 m)
      nx = 1000, ny = 1000,                  # Two integer scalars defining the number of grid elements in x and y (default is nx = ny = 1000)
      r = FFP.R,                             # Percentage of source area, i.e. a value between 10% and 90%. *2
      rslayer =  1,                          # Calculate footprint even if zm within roughness sublayer: set rslayer = 1 *3
      smooth_data = 1,                       # Apply convolution filter to smooth footprint climatology if smooth_data = 1 (default)
      # crop = 1,                            # Crop output area to size of the 80% footprint or the largest r given if crop = 1
      pulse = 1                              # Display progress of footprint calculations every pulse‐th footprint (e.g., "100")
      
    )
    
    
    # ***************************************************************************
    # Spatial analysis ---------------------------------------------------------
    
    # record FFP model error flag
    FFP.flag.err <- FFP$'flag_err'    # >> output
    
    # sometimes x-grid has NA -> treat as error
    
    # sometimes the x‐grid of footprint climatology is NA (NaN) also in case of no errors
    # or the distances estimation is unrealistic (see FFP CWI).
    # in this case force the flag to be 1 (as an error)
    if (any(is.na(FFP$'x_2d'))) FFP.flag.err <- 1
    if (!is.na(SnW.FFP) && SnW.FFP == 1) FFP.flag.err <- 1
    
    
    # Perform the analysis only in case of NO errors (flag) in the FFP estimation
    # 1 in case of error, 2 if not all contour plots (r%) within specified domain, 0 otherwise
    
    if (FFP.flag.err == 0) {
      
      # FFP isoplethes on which the analysis is based: this is a position index of the FFP.R vector
      X.perc.indx <- which(FFP.R %in% c(70, 80))  # ADATTARE CODICE IN MODO DA PREVEDERE QUALSIASI PERC. ATTUALMENTE 70 E 80 CI DEOVNO ESSERE PER FORZA
      
      # Plot the footprint climatology
      # OR_code
      # if (plot.FFP.clim) { 
      #   
      #   png(paste(FFP.dir.FFP, site.ID, '_','Ellispses_fit_', cur.ts.iso, '.png', sep=''),
      #       width = 14, height = 14, unit = 'cm', res=300)
      #   # windows()
      #   # Two‐dimensional footprint with contour lines of R% (using the fields package)
      #   # ffp2Dcol <- designer.colors(100, c("transparent", viridis(20, alpha = 0.8)) ) 
      #   ffp2Dcol <- c('transparent', gplots::rich.colors(100, alpha = 0.5))
      #   # breaks with a gap of 10 to 17 assigned the white color
      #   brks <- quantile(c(FFP$'fclim_2d'), na.rm = T)
      #   
      #   if (FFP$'flag_err' == 0 & !all(is.na(FFP$'x_2d'))) {
      #     image(x = FFP$'x_2d'[1, ], y = FFP$'y_2d'[ , 1], z = FFP$'fclim_2d', col = ffp2Dcol, xlab='', ylab='')
      #     # image(x = FFP$'x_2d'[1, ], y = FFP$'y_2d'[ , 1], z = FFP$'fclim_2d', col = ffp2Dcol, xlab='', ylab='', xlim=c(-200,200), ylim=c(-200,200))
      #     # windows()
      #     # EC tower
      #     points(0,0, pch=3)
      #     # isoplethes
      #     j <- 0
      #     for (j in 1:length(FFP.R)){
      #       lines(FFP$'xr'[[j]], FFP$'yr'[[j]], type="l", col="black")
      #     }
      #     lines(FFP$'xr'[[X.perc.indx[1]]], FFP$'yr'[[X.perc.indx[1]]], type="l", col="red", lwd=1)
      #     mtext(input.data$TIMESTAMP, 3, line=-2)
      #     mtext('Northing [m]', 2, line=2); mtext('Easting [m]', 1, line=2.5)
      #     legend('topleft', c('50%', '70%', '80%'), lty=1, col=c(1,2,1), bty='n')
      #     #
      #   }
      # }
      # NEW_code
      # Plot footprint climatology (use ggplot2 if preferred) --------------------
      if (plot.FFP.clim) {
        png(file = paste0(FFP.dir.FFP, site.ID, '_Ellispses_fit_', cur.ts.iso, '.png'),
            width = 14, height = 14, units = 'cm', res = 300)
        # windows()
        # palette
        ffp2d_col <- c('transparent', viridisLite::viridis(100, alpha = 0.5))
        if (FFP$flag_err == 0 && !all(is.na(FFP$x_2d))) {
          # base raster-like image plotting using terra for stability
          x_coords <- as.numeric(FFP$x_2d[1, ])
          y_coords <- as.numeric(FFP$y_2d[, 1])
          z_mat <- FFP$fclim_2d
          # simple base image (keeps behaviour similar to previous code)
          graphics::image(x = x_coords, y = y_coords, z = z_mat, col = ffp2d_col, xlab = '', ylab = '')
          points(0, 0, pch = 3) # EC tower
          for (j in seq_along(FFP.R)) {
            lines(as.numeric(FFP$xr[[j]]), as.numeric(FFP$yr[[j]]), col = "black")
          }
          # highlight requested percents
          lines(as.numeric(FFP$xr[[X.perc.indx[1]]]), as.numeric(FFP$yr[[X.perc.indx[1]]]), col = "red", lwd = 1)
          mtext(input.data$TIMESTAMP, 3, line = -2)
          mtext('Northing [m]', 2, line = 2); mtext('Easting [m]', 1, line = 2.5)
          legend('topleft', legend = c('50%', '70%', '80%'), lty = 1, col = c(1, 2, 1), bty = 'n')
        }
        dev.off()
      }
      
      
      
      ## FFP areas --------------------------------------------------------------
      
     
      # ### FFP 70% : ISOPLETHE (polygon) .........................................
      #
      # OR_code
      # # Isoplethe relative coordinates (northward/eastward)
      # # NOTA: CI POSSONO ESSERE NA NELLE COORDINATE, PER ORA SKIPPATI MA CAPIRE COME GESTIRLI !!!
      # FFP.70.iso.df <- data.frame(NULL)
      # FFP.70.iso.df <- data.frame(LON = as.numeric(FFP$'xr'[[X.perc.indx[1]]]),
      #                             LAT = as.numeric(FFP$'yr'[[X.perc.indx[1]]]),
      #                             row.names = NULL)
      # if (any(is.na(FFP.70.iso.df))) {
      #   FFP.70.iso.df$'LON' <- zoo::na.approx(FFP.70.iso.df$'LON')
      #   FFP.70.iso.df$'LAT' <- zoo::na.approx(FFP.70.iso.df$'LAT')
      # }
      # 
      # # NORTH/EAST to UTM and convert to object of class Polygons
      # if ('FFP.70.iso.utm.p' %in% ls()) {rm(FFP.70.iso.utm.p)}
      # FFP.70.iso.utm.p <- sp::Polygon(cbind(FFP.70.iso.df$'LON' + sp::coordinates(EC.tower.utm)[,1],
      #                                       FFP.70.iso.df$'LAT' + sp::coordinates(EC.tower.utm)[,2]))
      # 
      # # convert to object of class SpatialPolygons (assign projection)
      # FFP.70.iso <- sp::SpatialPolygons(list(Polygons(list(FFP.70.iso.utm.p), 'FFP.70.iso')),
      #                                   proj4string = EC.tower.utm@proj4string)
      # 
      # # check for validity of the FFP.70 Spatial polygon (self-intersection error)
      # if (!rgeos::gIsValid(FFP.70.iso, reason = F)) { # if gIsValid returns FALSE (NO OK):
      #   FFP.70.iso <- rgeos::gBuffer(FFP.70.iso, byid=TRUE, width=0) # use a zero-width buffer to clean up the topology problem
      # }
      # 
      # # Area of the geometry in the units of the current projection >> [m2]
      # FFP.70.iso.surf <- rgeos::gArea(FFP.70.iso, byid = FALSE) # >> output
      
      # Build 70% isopleth polygon (sf)
      # build data.frame of coords and interpolate NAs with approx
      iso70_lon <- as.numeric(FFP$xr[[X.perc.indx[1]]])
      iso70_lat <- as.numeric(FFP$yr[[X.perc.indx[1]]])
      if (any(is.na(iso70_lon)) || any(is.na(iso70_lat))) {
        # linear interpolation for NA (approx). preserves length ordering.
        not_na_lon <- !is.na(iso70_lon)
        if (sum(not_na_lon, na.rm = TRUE) >= 2) {
          iso70_lon <- approx(seq_along(iso70_lon)[not_na_lon], iso70_lon[not_na_lon],
                              xout = seq_along(iso70_lon), rule = 2)$y
        }
        not_na_lat <- !is.na(iso70_lat)
        if (sum(not_na_lat, na.rm = TRUE) >= 2) {
          iso70_lat <- approx(seq_along(iso70_lat)[not_na_lat], iso70_lat[not_na_lat],
                              xout = seq_along(iso70_lat), rule = 2)$y
        }
      }
      
      # add EC tower offset (EC.tower.utm assumed sf POINT)
      ec_coords <- sf::st_coordinates(EC.tower.utm)[1, ]
      iso70_coords <- cbind(iso70_lon + ec_coords[1], iso70_lat + ec_coords[2])
      # ensure polygon closed
      if (!all(iso70_coords[1, ] == iso70_coords[nrow(iso70_coords), ])) {
        iso70_coords <- rbind(iso70_coords, iso70_coords[1, ])
      }
      FFP.70.iso_sf <- sf::st_sfc(sf::st_polygon(list(iso70_coords)), crs = sf::st_crs(EC.tower.utm))
      # validity check and fix
      if (!all(sf::st_is_valid(FFP.70.iso_sf))) {
        FFP.70.iso_sf <- sf::st_make_valid(FFP.70.iso_sf)
        if (!all(sf::st_is_valid(FFP.70.iso_sf))) {
          # as last resort buffer zero
          FFP.70.iso_sf <- sf::st_buffer(FFP.70.iso_sf, 0)
        }
      }
      # area
      FFP.70.iso.surf <- as.numeric(sf::st_area(FFP.70.iso_sf))
      
      
      ### FFP 80% : ISOPLETHE (polygon) .........................................
      
      # OR_code
      # # Isoplethe relative coordinates (northward/eastward)
      # FFP.80.iso.df <- data.frame(NULL)
      # FFP.80.iso.df <- data.frame(LON = as.numeric(FFP$'xr'[[X.perc.indx[2]]]),
      #                             LAT = as.numeric(FFP$'yr'[[X.perc.indx[2]]]),
      #                             row.names = NULL)
      # if (any(is.na(FFP.80.iso.df))) { 
      #   FFP.80.iso.df$'LON' <- zoo::na.approx(FFP.80.iso.df$'LON')
      #   FFP.80.iso.df$'LAT' <- zoo::na.approx(FFP.80.iso.df$'LAT')
      # }
      # 
      # # NORTH/EAST to UTM and convert to object of class Polygons
      # if ('FFP.80.iso.utm.p' %in% ls()) {rm(FFP.80.iso.utm.p)}
      # FFP.80.iso.utm.p <- sp::Polygon(cbind(FFP.80.iso.df$'LON' + sp::coordinates(EC.tower.utm)[,1],
      #                                       FFP.80.iso.df$'LAT' + sp::coordinates(EC.tower.utm)[,2]))
      # 
      # # convert to object of class SpatialPolygons (assign projection)
      # FFP.80.iso <- NULL
      # FFP.80.iso <- sp::SpatialPolygons(list(Polygons(list(FFP.80.iso.utm.p), 'FFP.80.iso')),
      #                                   proj4string = EC.tower.utm@proj4string)
      # 
      # # check for validity of the FF80 Spatial polygon (self-intersection error)
      # if (!rgeos::gIsValid(FFP.80.iso, reason = F)) { # if gIsValid returns FALSE (NO OK):
      #   FFP.80.iso <- rgeos::gBuffer(FFP.80.iso, byid=TRUE, width=0) # use a zero-width buffer to clean up the topology problem
      # }
      # # repeat it (sometimes 2 iterations are needed)
      # if (!rgeos::gIsValid(FFP.80.iso, reason = F)) { # if gIsValid returns FALSE (NO OK):
      #   FFP.80.iso <- rgeos::gBuffer(FFP.80.iso, byid=TRUE, width=0) # use a zero-width buffer to clean up the topology problem
      # }
      
      # Build 80% isopleth polygon (sf)
      iso80_lon <- as.numeric(FFP$xr[[X.perc.indx[2]]])
      iso80_lat <- as.numeric(FFP$yr[[X.perc.indx[2]]])
      if (any(is.na(iso80_lon)) || any(is.na(iso80_lat))) {
        not_na_lon <- !is.na(iso80_lon)
        if (sum(not_na_lon, na.rm = TRUE) >= 2) {
          iso80_lon <- approx(seq_along(iso80_lon)[not_na_lon], iso80_lon[not_na_lon],
                              xout = seq_along(iso80_lon), rule = 2)$y
        }
        not_na_lat <- !is.na(iso80_lat)
        if (sum(not_na_lat, na.rm = TRUE) >= 2) {
          iso80_lat <- approx(seq_along(iso80_lat)[not_na_lat], iso80_lat[not_na_lat],
                              xout = seq_along(iso80_lat), rule = 2)$y
        }
      }
      iso80_coords <- cbind(iso80_lon + ec_coords[1], iso80_lat + ec_coords[2])
      if (!all(iso80_coords[1, ] == iso80_coords[nrow(iso80_coords), ])) {
        iso80_coords <- rbind(iso80_coords, iso80_coords[1, ])
      }
      FFP.80.iso_sf <- sf::st_sfc(sf::st_polygon(list(iso80_coords)), crs = sf::st_crs(EC.tower.utm))
      if (!all(sf::st_is_valid(FFP.80.iso_sf))) {
        FFP.80.iso_sf <- sf::st_make_valid(FFP.80.iso_sf)
        if (!all(sf::st_is_valid(FFP.80.iso_sf))) FFP.80.iso_sf <- sf::st_buffer(FFP.80.iso_sf, 0)
      }
      FFP.80.iso.surf <- as.numeric(sf::st_area(FFP.80.iso_sf))
      
      
      ## FFP 80% : FOOTPRINT VALUES matrix ......................................
      
      # OR_code
      # # Exctract the 80% contribution footprint values (matrix)
      # FFP.80.mx <- matrix(NA)
      # FFP.80.mx <- FFP$'fclim_2d'
      # FFP.80.mx[FFP.80.mx < FFP$'fr'[X.perc.indx[2]]] <- NA # set to Nas FFP values outside the 80% bounds
      # # sum(FFP.80.mx, na.rm = TRUE)
      # # windows()
      # 
      # # NORTH/EAST to UTM and promote to object of class RasterLayer
      # FFP.80.r <- raster::raster(list(x = FFP$x_2d[1, ] + sp::coordinates(EC.tower.utm)[,1],
      #                                 y = FFP$y_2d[, 1] + sp::coordinates(EC.tower.utm)[,2],
      #                                 z = FFP.80.mx),
      #                            crs = EC.tower.utm@proj4string)
      # proj4string(FFP.80.r) <- EC.tower.utm@proj4string
      # 
      # 
      # # Compute the reference FFP contribution to calculate the individual features weighted contributions
      # FFP.80.ref.contr <- sum(FFP.80.r@data@values, na.rm = TRUE) # should be 0.8
      # 
      # # Area of the geometry in the units of the current projection >> [m2]
      # FFP.80.iso.surf <- rgeos::gArea(FFP.80.iso, byid = FALSE) # >> output
      
      # FFP 80% matrix -> terra SpatRaster ...................................
      # take fclim_2d and mask values outside 80% threshold
      FFP.80.mx <- FFP$fclim_2d
      thr80 <- FFP$fr[X.perc.indx[2]]
      FFP.80.mx[FFP.80.mx < thr80] <- NA
      
      x_coords <- as.numeric(FFP$x_2d[1, ])
      y_coords <- as.numeric(FFP$y_2d[, 1])
      
      # Build an xyz data.frame for terra::rast(type="xyz")
      # expand.grid ensures consistent ordering: x varying fastest to match original matrix layout
      grid_df <- expand.grid(x = x_coords + ec_coords[1], y = y_coords + ec_coords[2])
      # matrix is oriented as [length(y) x length(x)]; use as.vector(t(...)) to get row-major order
      z_vals <- as.vector(t(FFP.80.mx))
      xyz <- data.frame(x = grid_df$x, y = grid_df$y, z = z_vals)
      
      # create SpatRaster; terra ignores NA cells automatically
      FFP.80.r <- tryCatch({
        terra::rast(xyz, type = "xyz", crs = sf::st_crs(EC.tower.utm)$wkt)
      }, error = function(e) {
        # fallback: create empty raster with extent and fill values
        r <- terra::rast(nrows = length(y_coords), ncols = length(x_coords),
                         xmin = min(x_coords + ec_coords[1], na.rm = TRUE),
                         xmax = max(x_coords + ec_coords[1], na.rm = TRUE),
                         ymin = min(y_coords + ec_coords[2], na.rm = TRUE),
                         ymax = max(y_coords + ec_coords[2], na.rm = TRUE),
                         crs = sf::st_crs(EC.tower.utm)$wkt)
        terra::values(r) <- z_vals
        r
      })
      
      # compute reference contribution
      FFP.80.ref.contr <- sum(terra::values(FFP.80.r), na.rm = TRUE)
      
      
      
      
      # Spatial data visualization and FFP .....................................
      # OR_code
      # if (plot.spat.FFP){
      #   
      #   png(paste(FFP.dir.FFP, site.ID, '_','SPfeat-FFP_', cur.ts.iso, '.png', sep=''),
      #       width = 14, height = 14, unit = 'cm', res=300)
      #   # windows()
      #   # TA + EC
      #   plot(TA.utm, axes=T,
      #        xlim=c((sp::coordinates(EC.tower.utm)[,1]-400), (sp::coordinates(EC.tower.utm)[,1]+400)),
      #        ylim=c((sp::coordinates(EC.tower.utm)[,2]-400), (sp::coordinates(EC.tower.utm)[,2]+400)))
      #   points(EC.tower.utm, pch=3)
      #   # EA
      #   if (length(EA.utm) > 0) {
      #     plot(EA.utm, add=T, col=adjustcolor(2, alpha.f = 0.2), border=NA)
      #     legend('topleft', legend=c('EC','TA','EA'), col = c(1,1, adjustcolor(2, alpha.f = 0.2)), pch=c(3,0,15), bty='n')
      #   } else {
      #     legend('topleft', legend=c('EC','TA'), col = c(1,1), pch=c(3,0), bty='n')
      #   }
      #   # LCT
      #   LCT.color <- terrain.colors(length(LCT.utm)+3)[-c((length(LCT.utm)+1):(length(LCT.utm)+3))]
      #   if (length(LCT.utm) > 0) {
      #     j <- 1
      #     for(j in 1:length(LCT.utm)){
      #       plot(LCT.utm[[j]], add=T, border = LCT.color[j], col=adjustcolor(LCT.color[j], alpha.f = 0.2))
      #       text(x=sp::coordinates(LCT.utm[[j]])[1], y=sp::coordinates(LCT.utm[[j]])[2], LCT.utm[[j]]@data$Description, pos=1, cex=.5)
      #     }
      #     legend('topright', legend=paste('LCT_0', seq_along(LCT.utm), sep=''), col = adjustcolor(LCT.color, alpha.f = 0.4), pch=15, bty='n')
      #     plot(TA.utm, lwd=3, add=T)
      #   }
      #   
      #   # FFP2D 80% footprint
      #   # plot(FFP.80.r, col=adjustcolor(ffp2Dcol, alpha.f = 0.6), add=T, legend=FALSE)
      #   # FFP2D isoplethes
      #   plot(FFP.70.iso, border=2, add=T)
      #   plot(FFP.80.iso, border=4, add=T)
      #   # axis(1)
      #   legend('bottomleft', c('FFP 70% contribution','FFP 80% contribution'), pch='0', col=c(2,4), bty='n')
      #   mtext(as.character(input.data$'TIMESTAMP'), side=3)
      #   
      #   dev.off()
      # }
      
      # Plot spatial features + FFP 
      if (plot.spat.FFP) {
        png(paste0(FFP.dir.FFP, site.ID, '_SPfeat-FFP_', cur.ts.iso, '.png'),
            width = 14, height = 14, units = 'cm', res = 300)
        # plot TA extent centered on EC tower
        # TA.utm assumed sf
        plot(sf::st_geometry(TA.utm), reset = FALSE,
                 xlim = c(ec_coords[1] - 400, ec_coords[1] + 400),
                 ylim = c(ec_coords[2] - 400, ec_coords[2] + 400))
        plot(sf::st_geometry(EC.tower.utm), add = TRUE, pch = 3)
        if (length(EA.utm) > 0) {
          plot(sf::st_geometry(EA.utm), add = TRUE, col = adjustcolor("red", alpha.f = 0.2), border = NA)
          graphics::legend('topleft', legend = c('EC', 'TA', 'EA'), col = c(1, 1, adjustcolor("red", alpha.f = 0.2)), pch = c(3, 0, 15), bty = 'n')
        } else {
          graphics::legend('topleft', legend = c('EC', 'TA'), col = c(1, 1), pch = c(3, 0), bty = 'n')
        }
        # LCT polygons if present (assumed list of sf)
        if (length(LCT.utm) > 0) {
          LCT_colors <- grDevices::terrain.colors(length(LCT.utm) + 3)[-c((length(LCT.utm) + 1):(length(LCT.utm) + 3))]
          for (j in seq_along(LCT.utm)) {
            plot(sf::st_geometry(LCT.utm[[j]]), add = TRUE, border = LCT_colors[j], col = adjustcolor(LCT_colors[j], alpha.f = 0.2))
          }
          plot(sf::st_geometry(TA.utm), add = TRUE, lwd = 2)
        }
        # plot 70/80 isopleths
        plot(FFP.70.iso_sf, add = TRUE, border = 2)
        plot(FFP.80.iso_sf, add = TRUE, border = 4)
        graphics::legend('bottomleft', c('FFP 70% contribution', 'FFP 80% contribution'), pch = '0', col = c(2, 4), bty = 'n')
        dev.off()
      }
      
      
      # Target Area - FFP intersections -----------------------------------------
      
      # # OR_code
      # # area of the geometry in the units of the current projection >> [m2]
      # TA.utm.surf <- rgeos::gArea(TA.utm, byid = FALSE) 
      # # areaPolygon(TA) - TA.utm.surf # difference with calculating the area using LAT/LONG
      # 
      # ### Check FFP/TA spatial intersection (contain or intersect) and possibly perform the analyis
      # if (isTRUE(rgeos::gContains(TA.utm, FFP.80.iso)) || isTRUE(rgeos::gIntersects(TA.utm, FFP.80.iso))) {
      #   
      #   
      #   ### Intersection bounds and its PHYSICAL area (using the isoplethe polygon)
      #   # Intersection spatial feature (SpatialPolygons). If the geometries do not intersect then an empty geometry is returned.
      #   FFP.TA.intersect <- NULL
      #   FFP.TA.intersect <- rgeos::gIntersection(spgeom1=TA.utm, spgeom2=FFP.80.iso, 
      #                                            byid = T, checkValidity=T) 
      #   # Area of intersection [m2]
      #   FFP.TA.intersect.surf <- rgeos::gArea(FFP.TA.intersect)    # >> output (FFP_TA_SURF)
      #   # Area of intersection [%] (of the whole FFP area) 
      #   FFP.TA.intersect.surf.perc <- FFP.TA.intersect.surf * 100 / FFP.80.iso.surf # >> output (FFP_TA_INTER_%SURF)
      #   
      #   ### Intersection FFP-WEIGHTED area (using the FFP80 raster)
      #   # Intersection spatial feature (RasterLayer)
      #   FFP.TA.intersect.r <- NULL
      #   FFP.TA.intersect.r <- raster::mask(FFP.80.r, FFP.TA.intersect, byid = TRUE)
      #   # plot(FFP.TA.intersect.r, col=adjustcolor(rich.colors(100), 0.5), add=F)
      #   
      #   # Cumulative FFP contribution of the intersection [m2 m-2, %]
      #   FFP.TA.intersect.contr <- sum(FFP.TA.intersect.r@data@values, na.rm = TRUE)  # >> output (FFP_TA_INTER_CONTIRIB) (must be as much as close to 0.8)
      #   
      # } else {
      #   
      #   FFP.TA.intersect.surf <- NA       # Area of FFP/TA intersection [m2] >> OUTPUT (FFP70_TA_INTER_SURF)
      #   FFP.TA.intersect.surf.perc <- NA  # Area of FFP/TA intersection [%] >> OUTPUT (FFP70_TA_INTER_%SURF)
      #   FFP.TA.intersect.contr <- NA      # FFP-weighted area of FFP/TA intersection [m2 m-2, %] >> OUTPUT (FFP70_TA_INTER_%SURF_W)
      #   
      # }
      # 
      # # visualize FFP/TA interaction
      # if (plot.TA.FFP){
      #   
      #   png(paste(FFP.dir.FFP, site.ID, '_','TA-FFP_', cur.ts.iso, '.png', sep=''),
      #       width = 14, height = 14, unit = 'cm', res=300, bg = 'transparent')
      #   # windows()
      #   plot(TA.utm, axes=F, border='white')
      #   # box()
      #   points(EC.tower.utm, pch=3, col='white')
      #   # ffp into TA
      #   plot(FFP.TA.intersect, add=T, col=adjustcolor('grey', 0.3), border=NA)
      #   # isoplethes
      #   plot(FFP.70.iso, add=T, border='red')
      #   plot(FFP.80.iso, add=T, border='orange')
      #   # mtext(as.character(input.data$'TIMESTAMP'), cex=0.8)
      #   # legend('topleft', c('TA', 'EC', 'FFP70', 'FFP80', 'FFP/TA inters'), pch=c(0,3,1,1,19), pt.cex=c(1,1,1,1,2), 
      #   # col=c(1,1,2,4,adjustcolor(3,0.5)), bty='n')
      #   dev.off()
      #   
      # }  
      
      # NEW_code
      # TA intersection with FFP ...............................................
      # area TA (in m2) using sf
      TA.utm.surf <- as.numeric(sf::st_area(TA.utm))
      
      # check intersection/contain with sf predicates
      intersects_TA <- as.logical(sf::st_intersects(TA.utm, FFP.80.iso_sf, sparse = FALSE)[1, 1])
      contains_TA <- as.logical(sf::st_contains(TA.utm, FFP.80.iso_sf, sparse = FALSE)[1, 1])
      
      if (intersects_TA || contains_TA) {
        
        # Do the intersection
        FFP.TA.intersect_sf <- sf::st_intersection(TA.utm, FFP.80.iso_sf)
        # And calculate the surfaces
        FFP.TA.intersect.surf <- as.numeric(sf::st_area(FFP.TA.intersect_sf))
        FFP.TA.intersect.surf.perc <- FFP.TA.intersect.surf * 100 / FFP.80.iso.surf
        # .........................
        # 2) Rendi valide le geometrie (per sicurezza)
        FFP.TA.intersect_sf <- sf::st_make_valid(FFP.TA.intersect_sf)
        # 3) Estrai solo POLYGON / MULTIPOLYGON dalla collection
        FFP.TA.intersect_poly <- sf::st_collection_extract(FFP.TA.intersect_sf, "POLYGON")
        # 4) Se ci sono più poligoni li fondi, altrimenti li lasci così
        # caso 0: nessuna intersezione → oggetto vuoto
        if (nrow(FFP.TA.intersect_poly) == 0) {
          FFP.TA.intersect_union <- FFP.TA.intersect_poly  # rimane vuoto
          # caso 1: più poligoni → unisci in uno solo
        } else if (nrow(FFP.TA.intersect_poly) > 1 ||
                   any(sf::st_geometry_type(FFP.TA.intersect_poly) == "MULTIPOLYGON")) {
          FFP.TA.intersect_union <- sf::st_union(FFP.TA.intersect_poly)
          # caso 2: già un solo poligono semplice → non serve union
        } else {
          FFP.TA.intersect_union <- FFP.TA.intersect_poly
        }
        
        # 5) mask raster by intersection geometry, no attributes (terra requires SpatVector)
        iv <- terra::vect(FFP.TA.intersect_union) # converting the sf object to a SpatVector
        masked <- try(terra::mask(FFP.80.r, iv), silent = TRUE)
        if (inherits(masked, "try-error")) {
          FFP.TA.intersect.contr <- NA
        } else {
          FFP.TA.intersect.contr <- sum(terra::values(masked), na.rm = TRUE)
        }
      } else {
        FFP.TA.intersect.surf <- NA
        FFP.TA.intersect.surf.perc <- NA
        FFP.TA.intersect.contr <- NA
      }
      
        
      if (plot.TA.FFP) {
        png(paste0(FFP.dir.FFP, site.ID, '_TA-FFP_', cur.ts.iso, '.png'),
            width = 14, height = 14, units = 'cm', res = 300, bg = 'transparent')
        # windows()
        plot(sf::st_geometry(TA.utm), axes = FALSE, border = 'black')
        plot(sf::st_geometry(EC.tower.utm), add = TRUE, pch = 3, col = 'black')
        if (!is.na(FFP.TA.intersect.surf) && length(FFP.TA.intersect_sf) > 0) {
          plot(FFP.TA.intersect_sf, add = TRUE, col = adjustcolor('grey', 0.3), border = NA)
        }
        plot(FFP.70.iso_sf, add = TRUE, border = 'red')
        plot(FFP.80.iso_sf, add = TRUE, border = 'orange')
        dev.off()
      }
      
      
      # ****************************************************************************
      # Exclusion Area - FFP intersections --------------------------------------
      
      # NOTE: This is NOT reported on the official description document
      
      # NOTE 20220317: currently the EA is treated as a unique, single object derived from the joining of the different EA
      #                in case the analysis must be performed on each individual EA, the script must be updated
      
      # OR_code
      # if (DO.EA.analysis & length(EA.utm) > 0) {
      #   
      #   ### EA areas [m2]
      #   EA.surf <- rgeos::gArea(EA.utm)
      #   
      #   ### Check FFP/EA spatial interaction (contain or intersect) and possibly perform the analyis
      #   FFP.contain.EA <- NULL
      #   FFP.contain.EA <- rgeos::gContains(FFP.80.iso, EA.utm, byid = TRUE)
      #   FFP.intersect.EA <- NULL
      #   FFP.intersect.EA <- rgeos::gIntersects(FFP.80.iso, EA.utm, byid = TRUE)
      #   
      #   if (FFP.contain.EA || FFP.intersect.EA) {
      #     
      #     ### Calculate possible intersection bounds and its PHYSICAL area (using the isoplethe polygon)
      #     # Intersection spatial feature (SpatialPolygons)
      #     FFP.EA.intersect <- NULL
      #     FFP.EA.intersect <- rgeos::gIntersection(FFP.80.iso, EA.utm, byid = TRUE) 
      #     
      #     # Areas of intersection [m2]
      #     FFP.EA.intersect.surf <- rgeos::gArea(FFP.EA.intersect)  # >> output  (FFP_EA_SURF)
      #     # Areas of FFP/EA intersection [%]  
      #     FFP.EA.intersect.surf.perc <- FFP.EA.intersect.surf * 100 / FFP.80.iso.surf # >> output (FFP_EA_INTER_%SURF)
      #     
      #     ### Intersection FFP-WEIGHTED area (using the FFP80 raster)
      #     # Intersection spatial features (RasterLayer)
      #     FFP.EA.intersect.r <- NULL 
      #     FFP.EA.intersect.r <- raster::mask(FFP.80.r, FFP.EA.intersect, byid = T)
      #     
      #     # Cumulative FFP contribution of the EA intersections
      #     FFP.EA.intersect.contr <- sum(FFP.EA.intersect.r@data@values, na.rm = TRUE) # >> output FFP_EA_INTER_CONTIRIB
      #     
      #     # BUGfix: Force the intersection contribution to be 80% maximum
      #     #         the possible error results from an error in the masking (from 'mask' function)
      #     if (FFP.EA.intersect.contr > 0.8) { FFP.EA.intersect.contr <- 0.8 }
      #     
      #     
      #     # Visualize FFP/EA interaction 
      #     if (plot.EA.FFP) {
      #       
      #       png(paste(FFP.dir.FFP, site.ID, '_','EA-FFP_', cur.ts.iso, '.png', sep=''),
      #           width = 14, height = 14, unit = 'cm', res=300)
      #       # windows()
      #       # TA + EC
      #       plot(TA.utm, axes=T)
      #       points(EC.tower.utm, pch=3); box()
      #       mtext(as.character(input.data$'TIMESTAMP'), cex=0.8)
      #       # EA
      #       plot(EA.utm, add=T, border = 8, col = adjustcolor(8, alpha.f = 0.3))
      #       # EA intersections
      #       if(!is.null(FFP.EA.intersect)){
      #         plot(FFP.EA.intersect, add=T, col=adjustcolor('darkgrey', alpha.f = 0.9), border=NA)
      #       }
      #       # isoplethes
      #       if(!is.null(FFP.80.iso)){
      #         plot(FFP.70.iso, border=2, lwd=1, add=T)
      #         plot(FFP.80.iso, border=4, add=T)
      #       }
      #       legend('topleft', c('TA','EC','FFP70','FFP80','FFP/EA inters'), pch=c(0,3,1,1,19), pt.cex=c(1,1,1,1,2), col=c(1,1,2,4,adjustcolor(8,0.9)), bty='n')
      #       
      #       dev.off()
      #     }
      #     
      #     
      #   } else {
      #     
      #     FFP.EA.intersect.surf <- 0
      #     FFP.EA.intersect.surf.perc <- 0
      #     FFP.EA.intersect.contr <- 0
      #     
      #   }
      #   
      #   FFP.EA.intersect.surf <- NA
      #   FFP.EA.intersect.surf.perc <- NA
      #   FFP.EA.intersect.contr <- NA
      #   
      # } else {
      #   
      #   FFP.EA.intersect.surf <- NA         # Area of intersection between FFP and the Exclusion Area [m2]
      #   FFP.EA.intersect.surf.perc <- NA    # Area of intersection between FFP and the Exclusion Area as percentage of the 80% contribution isoplethe area
      #   FFP.EA.intersect.contr <- NA        # FFP cumulative contribution of the Exclusion Area (intersection) [m2 m-2, %]
      # }
      
      # EA intersections (if requested) ........................................
      if (DO.EA.analysis && length(EA.utm) > 0) {
        EA_surf <- as.numeric(sf::st_area(EA.utm))
        contain_EA <- as.logical(sf::st_contains(FFP.80.iso_sf, EA.utm, sparse = FALSE)[1, 1])
        intersect_EA <- any(as.logical(sf::st_intersects(FFP.80.iso_sf, EA.utm, sparse = FALSE)))
        
        if (contain_EA || intersect_EA) {
          FFP.EA.intersect_sf <- sf::st_intersection(FFP.80.iso_sf, EA.utm)
          FFP.EA.intersect.surf <- as.numeric(sf::st_area(FFP.EA.intersect_sf))
          FFP.EA.intersect.surf.perc <- FFP.EA.intersect.surf * 100 / FFP.80.iso.surf
          
          iv <- terra::vect(FFP.EA.intersect_sf)
          masked_ea <- try(terra::mask(FFP.80.r, iv), silent = TRUE)
          if (inherits(masked_ea, "try-error")) {
            FFP.EA.intersect.contr <- NA
          } else {
            FFP.EA.intersect.contr <- sum(terra::values(masked_ea), na.rm = TRUE)
            if (!is.na(FFP.EA.intersect.contr) && FFP.EA.intersect.contr > 0.8) FFP.EA.intersect.contr <- 0.8
          }
        } else {
          FFP.EA.intersect.surf <- 0
          FFP.EA.intersect.surf.perc <- 0
          FFP.EA.intersect.contr <- 0
        }
        
        if (plot.EA.FFP) {
          png(paste0(FFP.dir.FFP, site.ID, '_EA-FFP_', cur.ts.iso, '.png'),
              width = 14, height = 14, units = 'cm', res = 300)
          windows()
          plot(sf::st_geometry(TA.utm))
          plot(sf::st_geometry(EC.tower.utm), add = TRUE, pch = 3)
          plot(sf::st_geometry(EA.utm), add = TRUE, border = 8, col = adjustcolor(8, alpha.f = 0.3))
          if (exists("FFP.EA.intersect_sf") && length(FFP.EA.intersect_sf) > 0) {
            plot(FFP.EA.intersect_sf, add = TRUE, col = adjustcolor('darkgrey', alpha.f = 0.6), border = NA)
          }
          if (exists("FFP.70.iso_sf")) plot(FFP.70.iso_sf, add = TRUE, border = 2)
          if (exists("FFP.80.iso_sf")) plot(FFP.80.iso_sf, add = TRUE, border = 4)
          dev.off()
        }
      } else {
        FFP.EA.intersect.surf <- NA
        FFP.EA.intersect.surf.perc <- NA
        FFP.EA.intersect.contr <- NA
      }
      
      
      # *************************************************************************
      # Land Cover Types - FFP intersections ------------------------------------
      
      # OR_code
      # if (DO.LCT.analysis & length(LCT.utm) > 0) {
      #   
      #   ### LCT areas [m2]
      #   LCT.surf <- lapply(LCT.utm, function(x) rgeos::gArea(x))
      #   
      #   
      #   ### Check FFP/LCT spatial interation (contain or intersect) and possibly perform the analyis
      #   if (any(unlist(lapply(LCT.utm, function(x) rgeos::gContains(FFP.80.iso, x, byid = TRUE)))) ||
      #       any(unlist(lapply(LCT.utm, function(x) rgeos::gIntersects(FFP.80.iso, x, byid = TRUE))))) {
      #     
      #     ### Calculate possible intersection bounds and its PHYSICAL area (using the isoplethe polygon)
      #     # Intersection spatial feature (SpatialPolygons)
      #     FFP.LCT.intersect <- list(NULL)
      #     FFP.LCT.intersect <- lapply(LCT.utm, function(x) rgeos::gIntersection(FFP.80.iso, x, byid = TRUE)) 
      #     
      #     # Areas of intersection [m2]
      #     # Note: a loop must be used because gArea crashes with NULL features
      #     FFP.LCT.intersect.surf <- list(NULL)
      #     ii <- 0
      #     for (ii in 1 : length(FFP.LCT.intersect)) {
      #       if (!is.null(FFP.LCT.intersect[[ii]])) {
      #         FFP.LCT.intersect.surf[[ii]] <- rgeos::gArea(FFP.LCT.intersect[[ii]])
      #       } else {
      #         FFP.LCT.intersect.surf[[ii]] <- NA
      #       }
      #     }
      #     # simplify the list to vector for easier data management
      #     FFP.LCT.intersect.surf <- unlist(FFP.LCT.intersect.surf) # >> output (FFP_LCTxx_SURF)
      #     # gArea(FFP.80.iso); sum(FFP.LCT.intersect.surf)
      #     
      #     # Areas of FFP/LCT intersection [%]
      #     FFP.LCT.intersect.surf.perc <- unlist(
      #       lapply(FFP.LCT.intersect.surf, function(x) x * 100 / FFP.80.iso.surf)   # >> output (FFP_LCTxx_INTER_%SURF)
      #     )
      #     
      #     ### Intersection FFP-WEIGHTED areas (using the FFP80 raster)
      #     
      #     # Intersection spatial features (RasterLayer)
      #     FFP.LCT.intersect.r <- list(NULL) 
      #     ii <- 1
      #     for (ii in 1 : length(FFP.LCT.intersect)) {
      #       if (!is.null(FFP.LCT.intersect[[ii]])) {
      #         FFP.LCT.intersect.r[[ii]] <- raster::mask(FFP.80.r, FFP.LCT.intersect[[ii]], byid = T)
      #       }else{
      #         FFP.LCT.intersect.r[[ii]] <- NA
      #       }
      #     }
      #     
      #     # Cumulative FFP contribution of the (each) LCT intersections
      #     FFP.LCT.intersect.contr <- list(NULL) 
      #     ii <- 0
      #     for (ii in 1 : length(FFP.LCT.intersect)) {
      #       if (!is.null(FFP.LCT.intersect[[ii]])) {
      #         FFP.LCT.intersect.contr[ii] <- sum(FFP.LCT.intersect.r[[ii]]@data@values, na.rm = TRUE)
      #       }else{
      #         FFP.LCT.intersect.contr[ii] <- NA
      #       }
      #     }
      #     # simplify the list to vector for easier data management
      #     FFP.LCT.intersect.contr <- unlist(FFP.LCT.intersect.contr)  # >> output FFP_LCTxx_INTER_CONTIRIB
      #     
      #     
      #     # BUGfix: Force the intersection contribution to be 80% maximum by rescaling the LCT contributions sum to 0.8
      #     #         the possible error results from an error in the masking (from 'mask' function)
      #     if (sum(FFP.LCT.intersect.contr, na.rm = TRUE) > 0.8) {
      #       FFP.LCT.intersect.contr <- FFP.LCT.intersect.contr * (0.8/sum(FFP.LCT.intersect.contr))
      #     }
      #     
      #     
      #     # Visualize FFP/LCT interaction .......................
      #     
      #     if (plot.LCT.FFP) {
      #       
      #       LCT.color <- terrain.colors(length(LCT.utm)+3)[-c((length(LCT.utm)+1):(length(LCT.utm)+3))]
      #       
      #       png(paste(FFP.dir.FFP, site.ID, '_','LCT-FFP_', cur.ts.iso, '.png', sep=''),
      #           width = 14, height = 14, unit = 'cm', res=300)
      #       # windows()
      #       # TA + EC
      #       plot(TA.utm, axes=T)
      #       points(EC.tower.utm, pch=3); box()
      #       mtext(as.character(input.data$'TIMESTAMP'), cex=0.8)
      #       # EA
      #       if (length(EA.utm) > 0) {
      #         plot(EA.utm, add=T, border = 8, col = adjustcolor(8, alpha.f = 0.3))
      #       }
      #       # LCT
      #       i.lct <- 0
      #       for (i.lct in 1:length(LCT.utm)) {
      #         plot(LCT.utm[[i.lct]], add=T, border = LCT.color[i.lct], col=adjustcolor(LCT.color[i.lct], alpha.f = 0.2))
      #         text(x=sp::coordinates(LCT.utm[[i.lct]])[1], y=sp::coordinates(LCT.utm[[i.lct]])[2], LCT.utm[[i.lct]]@data$Description, pos=1, cex=.8)
      #       }
      #       plot(TA.utm, lwd=2, add=T)
      #       
      #       # LCT intersections
      #       i.lct <- 0
      #       for (i.lct in 1:length(LCT.utm)) {
      #         if (!is.null(FFP.LCT.intersect[[i.lct]])) {
      #           plot(FFP.LCT.intersect[[i.lct]], add=T, col=adjustcolor(LCT.color[i.lct], alpha.f = 0.6), border=NA)
      #         } else {
      #           next
      #         }
      #       }
      #       # FFP2D isoplethes
      #       if (!is.null(FFP.80.iso)) {
      #         plot(FFP.70.iso, border=2, lwd=1, add=T)
      #         plot(FFP.80.iso, border=4, add=T)
      #       }
      #       
      #       legend('bottomleft', legend=paste('LCT_0', seq_along(LCT.utm), sep=''), col = adjustcolor(LCT.color, alpha.f = 0.8), pch=15, bty='n')
      #       legend('topleft', c('TA','EC','FFP70','FFP80',paste0('FFP/LCT_0', 1:length(LCT.utm), ' inters')), pch=c(0,3,1,1,rep(19, length(LCT.utm))), 
      #              pt.cex=c(1,1,1,1,rep(2, length(LCT.utm))), col=c(1,1,2,4,adjustcolor(LCT.color, alpha.f = 0.8)), bty='n')
      #       
      #       # # - visualizzare paramteri turbolence su plot esempi ffp
      #       # # could be a not possible combination of roughness length and turbulence input (N. Kljun pers. comm.)
      #       # legend('bottomright', legend=c(bquote(bar(u) == .(round(input.data$umean, 2))),
      #       #                                bquote(z[0] == .(round(input.data$z0, 2))),
      #       #                                bquote(z[m]/L == .(round(input.data$zm/input.data$ol, 2))),
      #       #                                bquote('L' == .(round(input.data$ol, 2))),
      #       #                                bquote(sigma[v] == .(round(input.data$sigmav, 2))),
      #       #                                bquote(u^'*' == .(round(input.data$ustar, 2))),
      #       #                                as.expression(bquote(PBL == .(round(input.data$PBL, 2)))),
      #       #                                bquote(SnW[ratio] == .(round(FFP.el.M_m, 2)))),
      #       #                                # bquote(HWH[test] == .(set.cleaned$HWHTest))),
      #       #        bty='n')
      #       
      #       dev.off()
      #       
      #     }
      #     # .....................................................
      #     
      #   } else {
      #     
      #     FFP.LCT.intersect.surf <- rep(NA, length(LCT.utm))       # list of list with areas of FFP/LCTs intersection [m2] >> OUTPUT (FFP70_LCTxx_INTER_SURF)
      #     FFP.LCT.intersect.surf.perc <- rep(NA, length(LCT.utm))  # list of list with areas of FFP/LCTs intersection [%] >> OUTPUT (FFP70_LCTxx_INTER_%SURF)
      #     FFP.LCT.intersect.contr <- rep(NA, length(LCT.utm))      # list of list with cumulative FFP contribution of the intersections with LCT
      #     
      #   }
      #   
      # }
      # 
      # NEW_code
      # LCT intersections (if requested) .......................................
      if (DO.LCT.analysis && length(LCT.utm) > 0) {
        
        # areas
        LCT_surf <- lapply(LCT.utm, function(x) as.numeric(sf::st_area(x)))
        
        # check any intersects/contains
        intersects_any <- any(vapply(LCT.utm, function(x) {
          any(sf::st_intersects(FFP.80.iso_sf, x, sparse = FALSE))
        }, logical(1)))
        
        if (intersects_any) {
          FFP.LCT.intersect <- lapply(LCT.utm, function(x) {
            intr <- try(sf::st_intersection(FFP.80.iso_sf, x), silent = TRUE)
            if (inherits(intr, "try-error") || length(intr) == 0) return(NULL)
            intr
          })
          
          # areas per LCT
          FFP.LCT.intersect.surf <- vapply(FFP.LCT.intersect, function(g) {
            if (is.null(g)) return(NA_real_)
            sum(as.numeric(sf::st_area(g)), na.rm = TRUE)
          }, numeric(1))
          
          FFP.LCT.intersect.surf.perc <- vapply(FFP.LCT.intersect.surf, function(a) {
            if (is.na(a)) return(NA_real_)
            a * 100 / FFP.80.iso.surf
          }, numeric(1))
          
          # masked rasters + contributions
          FFP.LCT.intersect.contr <- vapply(FFP.LCT.intersect, function(g) {
            if (is.null(g)) return(NA_real_)
            iv <- terra::vect(g)
            masked <- try(terra::mask(FFP.80.r, iv), silent = TRUE)
            if (inherits(masked, "try-error")) return(NA_real_)
            sum(terra::values(masked), na.rm = TRUE)
          }, numeric(1))
          
          # bugfix: scale to max 0.8 if sum > 0.8
          total_contr <- sum(FFP.LCT.intersect.contr, na.rm = TRUE)
          if (!is.na(total_contr) && total_contr > 0.8) {
            FFP.LCT.intersect.contr <- FFP.LCT.intersect.contr * (0.8 / total_contr)
          }
          
          # plotting if requested
          if (plot.LCT.FFP) {
            LCT_colors <- grDevices::terrain.colors(length(LCT.utm) + 3)[-c((length(LCT.utm) + 1):(length(LCT.utm) + 3))]
            png(paste0(FFP.dir.FFP, site.ID, '_LCT-FFP_', cur.ts.iso, '.png'),
                width = 14, height = 14, units = 'cm', res = 300)
            sf::plot(sf::st_geometry(TA.utm))
            sf::plot(sf::st_geometry(EC.tower.utm), add = TRUE, pch = 3)
            if (length(EA.utm) > 0) sf::plot(sf::st_geometry(EA.utm), add = TRUE, border = 8, col = adjustcolor(8, alpha.f = 0.3))
            for (i_lct in seq_along(LCT.utm)) {
              sf::plot(sf::st_geometry(LCT.utm[[i_lct]]), add = TRUE, border = LCT_colors[i_lct], col = adjustcolor(LCT_colors[i_lct], alpha.f = 0.2))
            }
            sf::plot(sf::st_geometry(TA.utm), add = TRUE, lwd = 2)
            # intersections overlays
            for (i_lct in seq_along(FFP.LCT.intersect)) {
              if (!is.null(FFP.LCT.intersect[[i_lct]])) {
                sf::plot(FFP.LCT.intersect[[i_lct]], add = TRUE, col = adjustcolor(LCT_colors[i_lct], alpha.f = 0.6), border = NA)
              }
            }
            if (exists("FFP.70.iso_sf")) sf::plot(FFP.70.iso_sf, add = TRUE, border = 2)
            if (exists("FFP.80.iso_sf")) sf::plot(FFP.80.iso_sf, add = TRUE, border = 4)
            dev.off()
          }
          
        } else {
          FFP.LCT.intersect.surf <- rep(NA_real_, length(LCT.utm))
          FFP.LCT.intersect.surf.perc <- rep(NA_real_, length(LCT.utm))
          FFP.LCT.intersect.contr <- rep(NA_real_, length(LCT.utm))
        }
        
      } else {
        # default NA outputs if LCT not enabled
        FFP.LCT.intersect.surf <- rep(NA_real_, ifelse(exists("LCT.utm"), length(LCT.utm), 0))
        FFP.LCT.intersect.surf.perc <- rep(NA_real_, ifelse(exists("LCT.utm"), length(LCT.utm), 0))
        FFP.LCT.intersect.contr <- rep(NA_real_, ifelse(exists("LCT.utm"), length(LCT.utm), 0))
      }
      
      
    } else { # In case of ERROR in the FFP climatology computation
      
      # FFP CWI (actually not needed as they should be already set to NA)
      FFP.cwi.peak <- NA
      FFP.cwi.offset <- NA
      FFP.cwi.cum10.dist <- NA
      FFP.cwi.cum20.dist <- NA
      FFP.cwi.cum30.dist <- NA
      FFP.cwi.cum40.dist <- NA
      FFP.cwi.cum50.dist <- NA
      FFP.cwi.cum60.dist <- NA
      FFP.cwi.cum70.dist <- NA
      FFP.cwi.cum80.dist <- NA
      FFP.cwi.cum90.dist <- NA
      FFP.cwi.integral <- NA
      FFP.el.m <- NA
      FFP.el.M <- NA
      FFP.el.M_m <- NA
      SnW.FFP <- SnW.FFP
      
      ### Single footprint 2D
      FFP.70.iso.surf <- NA             # Area of the FFP geometry in the units of the current projection >> UTM, [m2]
      FFP.80.iso.surf <- NA             # Area of the FFP geometry in the units of the current projection >> UTM, [m2]
      FFP.80.ref.contr <- NA            # Reference FFP contribution used to calculate the individual features weighted contributions
      
      FFP.TA.intersect.surf <- NA       # Area of FFP/TA intersection [m2] >> OUTPUT (FFP70_TA_INTER_SURF)
      FFP.TA.intersect.surf.perc <- NA  # Area of FFP/TA intersection [%] >> OUTPUT (FFP70_TA_INTER_%SURF)
      FFP.TA.intersect.contr <- NA      # FFP-weighted area of FFP/TA intersection [m2 m-2, %] >> OUTPUT (FFP70_TA_INTER_%SURF_W)
      
      # if (DO.EA.analysis) {
      FFP.EA.intersect.surf <- NA        # Areas of FFP/EA intersection [m2] >> OUTPUT (FFP70_EA_INTER_SURF)
      FFP.EA.intersect.surf.perc <- NA   # Areas of FFP/EA intersection [%] >> OUTPUT (FFP70_EA_INTER_%SURF)
      FFP.EA.intersect.contr <- NA       # FFP-weighted area of FFP/TA intersection [m2 m-2, %] >> OUTPUT (FFP70_TA_INTER_%SURF_W)
      # }
      
      if (DO.LCT.analysis) {
        FFP.LCT.intersect.surf <- rep(NA, length(LCT.utm))       # Areas of FFP/LCTs intersection [m2] >> OUTPUT (FFP70_LCTxx_INTER_SURF)
        FFP.LCT.intersect.surf.perc <- rep(NA, length(LCT.utm))  # Areas of FFP/LCTs intersection [%] >> OUTPUT (FFP70_LCTxx_INTER_%SURF)
        FFP.LCT.intersect.contr <- rep(NA, length(LCT.utm))      # Cumulative FFP contribution of the intersections with LCT
      }
      
    }
    
  } # COMPUTE FFP AND INTERSECTIONS CHUNK END
  
  
  # Output result ---------------------------------------------------------------
  FFP.result <- c(NULL)
  FFP.result <- c(
    TIMESTAMP_END = cur.ts.iso,
    # Single footprint 
    FFP_cwi_FLAG = FFP.flag.err,                  # FFP CWI error (0/1)
    FETCH_MAX = FFP.cwi.peak,                  # FFP CWI peak location [m]
    FETCH_OFFSET = FFP.cwi.offset,                # FFP CWI offset distance [m]
    FETCH_10 = FFP.cwi.cum10.dist,            # upwind distance at which the crosswind integrated FFP reaches the 10% cumulative contribution
    FETCH_20 = FFP.cwi.cum20.dist,            # upwind distance at which the crosswind integrated FFP reaches the 20% cumulative contribution
    FETCH_30 = FFP.cwi.cum30.dist,            # upwind distance at which the crosswind integrated FFP reaches the 30% cumulative contribution
    FETCH_40 = FFP.cwi.cum40.dist,            # upwind distance at which the crosswind integrated FFP reaches the 40% cumulative contribution
    FETCH_50 = FFP.cwi.cum50.dist,            # upwind distance at which the crosswind integrated FFP reaches the 50% cumulative contribution
    FETCH_60 = FFP.cwi.cum60.dist,            # upwind distance at which the crosswind integrated FFP reaches the 60% cumulative contribution
    FETCH_70 = FFP.cwi.cum70.dist,            # upwind distance at which the crosswind integrated FFP reaches the 70% cumulative contribution
    FETCH_80 = FFP.cwi.cum80.dist,            # upwind distance at which the crosswind integrated FFP reaches the 80% cumulative contribution
    FETCH_90 = FFP.cwi.cum90.dist,            # upwind distance at which the crosswind integrated FFP reaches the 90% cumulative contribution
    FFP_cwi_ellmsa = FFP.el.m,                      # Minor semi-axis of the fitted ellepses to the FFP CWI 70% contribution isoplethe
    FFP_cwi_ellMsa = FFP.el.M,                      # Major semi-axis of the fitted ellepses to the FFP CWI 70% contribution isoplethe
    FFP_cwi_ellMmR = FFP.el.M_m,                    # Ellipse axis ratio
    SnW_FFP_FLAG = SnW.FFP,                       # Unrealistic FFP shape (named Short and Wide (SnW) FFP): 1 yes, 0 no
    # Footprint climatology and spatial analysis
    FFP_FLAG = FFP.flag.err,                  # FFP CLIM error (0/1)
    FOOTPRINT_70_SURF = FFP.70.iso.surf,               # Area of the 70% contribution isoplethe geometry in the units of the current projection [m2]
    FOOTPRINT_80_SURF = FFP.80.iso.surf,               # Area of the 80% contribution isoplethe geometry in the units of the current projection [m2]
    FOOTPRINT_TA_SURF = FFP.TA.intersect.surf,         # Area of intersection between FFP and the Target Area [m2]
    FOOTPRINT_TA_PSURF = FFP.TA.intersect.surf.perc,    # Area of intersection between FFP and the Target Area as percentage of the 80% contribution isoplethe area
    FOOTPRINT_TA_CONTR = FFP.TA.intersect.contr,        # FFP cumulative contribution of the Target Area (intersection) [m2 m-2, %] 
    FOOTPRINT_EA_SURF = FFP.EA.intersect.surf,         # Area of intersection between FFP and the Exclusion Area [m2]
    FOOTPRINT_EA_PSURF = FFP.EA.intersect.surf.perc,    # Area of intersection between FFP and the Exclusion Area as percentage of the 80% contribution isoplethe area
    FOOTPRINT_EA_CONTR = FFP.EA.intersect.contr        # FFP cumulative contribution of the Exclusion Area (intersection) [m2 m-2, %] 
    
  )
  
  
  # Add the LCT analysis results
  if (DO.LCT.analysis & length(LCT.utm) > 0) {
    
    names(FFP.LCT.intersect.surf) <- paste0('FOOTPRINT_LCT_0', 1:length(LCT.utm), '_SURF')       # Areas of intersections between FFP and the Land Cover Typology surfaces [m2] (FFP_LCTxx_SURF)
    names(FFP.LCT.intersect.surf.perc) <- paste0('FOOTPRINT_LCT_0', 1:length(LCT.utm), '_PSURF')   # Area of intersection between FFP and the Land Cover Typology surfaces as percentage of the 80% contribution isoplethe area (FFP_LCTxx_INTER_%SURF)
    names(FFP.LCT.intersect.contr) <- paste0('FOOTPRINT_LCT_0', 1:length(LCT.utm), '_CONTR')        # FFP cumulative contribution of the Land Cover Typologies (intersections) [m2 m-2, %] (FFP_LCTxx_INTER_CONTIRIB)
    FFP.result <- c(FFP.result, FFP.LCT.intersect.surf, FFP.LCT.intersect.surf.perc, FFP.LCT.intersect.contr)
  }
  
  return(FFP.result)
  
}

