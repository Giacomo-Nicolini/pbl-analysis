

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



FFPanalysis <- function(
    
  site.ID,
  proc.mode = 'L2',                   # processing mode: can be 'L2' for standard processing or 'TEST' for ICOS T2/T3 ouptuts
  FFP.func.dir,                       # Directory of the functions (custom and from Kljun) needed fort the analysis
  FFP.dir,                            # Directory of the footprint analysis: will contain the input spatial data and the results
  MM.data.dir,                        # Directory of the micro-meteorological (QCed) input dataset file. In ICOS is the "cleanset" file created by RFlux.
  MM.data.url = NULL,                 # Alternative to [MM.data.dir]: link to a GoogleDrive folder containing the micro-meteorological (QCed) input dataset file. In ICOS is the "cleanset" file created by RFlux.
  MM.data.patern = 'clean',           # Name pattern of the micro-meteorological (QCed) input dataset file
  SP.data.dir,                        # Directory containing the spatial data features (Target Area, possibly Exclusion areass and Land cover typologies polygons). If NULL, it is tried to download them from the ICOS database.
  QC.done = TRUE,                     # The QC on half-hourly data has been done and results are available to use
  ecmd.available = TRUE,              # Station metadata are available
  ecmd.path = NULL,                   # EC station metadata (table) path: it can be a direct path or NULL In the latter case it could be downloaded by means of the ETC web-service (if available)
  Nrec.min = 50,                      # Minimum percentage of N records (after QC filter) to perform the analysis 
  DO.EA.analysis = FALSE,             # Perform the analysis on FFP/EA interaction
  save.FFPEA.inter = FALSE,           # In case the the analysis on FFP/EA interaction is done, save the FFP/EA intersection polygons?
  DO.LCT.analysis = FALSE,             # whether run the analysis on LCTs or not
  DO.ustar.filter = FALSE,            # Perform USTAR filtering for FFP data. Note that u* > 0.2 ms-1 is already a conditon for FFP model
  ustar.thr = 0.2,                    # In case USTAR filtering for FFP data is used, set ustar threshold (ms-1) 
  DO.wind.sect.excl = FALSE,          # Perform the wind sector exclusion
  remove.SnW = FALSE,                  # Remove data flagged for Short and Wide footprints # RIVEDERE E ADATTAER
  SnW.thr = 5,                        # In case Short and Wide footprints are removed, set the threshold to flag a record as the ratio between major and minor ellipse axes of the 70% isoplethe)# ??? # DEVE ESSERE SETTATO IN BASE A ANALISI DATI
  # FFP settings 
  FFP.R = c(50,70,80),                # FFP cumulative % to compute (can be multiple)
  ffp.domain = 1000,                  # FFP domain size as an array of (xmin xmax ymin ymax) [m]
  cell.size = 1,                      # Cell size of domain [m] (default is dx = dy = 2 m)
  DO.climatology = TRUE,             # Compute and plot the footprint climatology
  DO.climatology.NT = FALSE,          # Compute and plot the footprint climatology during nightime
  # Plotting 
  plot.spatial.feats = TRUE,         # plot site spatial features
  plot.FFP.cwi = FALSE,               # plot cross-wind integrated FFP 
  plot.ellipse.fit = FALSE,           # plot fitted ellipses
  show.clim.raster = FALSE            # Whether to over-impose the climatology raster on the figure
){
  
  # ******************************************************************************
  # ........................ FOOTPRINT ANALYSIS ...............................
  # ******************************************************************************
  
  # FARE HEADER
  
  # ******************************************************************************
  
  # Libraries ....................................................................
  kpacks <- c("birk",
              "crayon",
              "data.table",
              "dismo",
              "EBImage",      # for FFP Climatology
              "fields",
              "flux",         # Library for Estimate gas flux rates using non-steady-state closed chamber data >> CHECK IT
              "geosphere",
              "gplots",
              "gstat",
              "maptools",     # for reading shape files  along with it come rgeos and sp
              "openair",
              "plotKML",
              "qdapTools",
              "raster",
              "RColorBrewer",
              "reshape2",
              # "rgdal",       # for reading kml files > RIMOSSA
              "rgeos",         # opzionale: non più necessaria per le parti riscritte con 'sf'
              "sp",            # # mantenuto per compatibilità minima, ma il flusso spaziale ora usa 'sf'
              "sf",            # NUOVO: gestione spaziale moderna
              "lwgeom",        # NUOVO: per funzioni topologiche avanzate e st_make_valid
              "spcosa",        
              "xts",
              "zoo",
              'foreach',
              'doSNOW',
              'doParallel')
  new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
  if (length(new.packs)){install.packages(new.packs)} # install packages not already installed
  loaded.packs <- lapply(kpacks, require, character.only=T) # load them
  if (!all(unlist(loaded.packs))){cat(yellow $ bold(paste('WARNING: not able to load the following libraries:\n', 
                                                          paste(kpacks[which(unlist(loaded.packs)==F)], collapse=', '))))
  }
  remove(kpacks, new.packs)
  
  # # in case of EBImage problems:
  # # With R version 3.5 or greater, install Bioconductor packages using BiocManager
  # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # BiocManager::install()
  # BiocManager::install("EBImage")
  # library(EBImage)
  # ##  previous version:
  # source("https://bioconductor.org/biocLite.R")
  # biocLite()
  # biocLite("EBImage")
  # library(EBImage)
  # install.packages("")
  
  # Set global options
  options(scipen = 999)
  
  # Define message themes
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  # *****************************************************************************
  cat(prog.mes(paste0('\n\n[',Sys.time(),'] Computing FOOTPRINTS and perform the spatial analysis\n')))
  cat(prog.mes('***************************************************************************\n'))
  
  
  
  ## Functions ...................................................................
  
  # Custom
  source(paste0(FFP.func.dir, 'GDrive_files_dwld.R'))         # download files from GoogleDrive
  source(paste0(FFP.func.dir, 'SpatialData_ICOSDB_dwnld.R'))  # download spatial data files from ICOS database: it contains 2 functions,one to download polygond one for points (EC, CP, SPI and SPII)
  source(paste0(FFP.func.dir, 'long2UTM.R'))                  # longlat to UTM conversion
  source(paste0(FFP.func.dir, 'calc_PBL_height.R'))           # Planetary boundary layer (Mixing) height estimation (GN) ATTUALMENTE USO QUESTA VERSIONE MA LA VERSIONE IN ("<local path>") E' PIU ATTUALE E CI STO LAVORANDO
  source(paste0(FFP.func.dir, 'proc_PBL_height.R'))           # Planetary boundary layer (Mixing) height quality check and gapfilling
  source(paste0(FFP.func.dir, 'EllipseFit.R'))                # Fit ellipse to data points
  # source(paste0(FFP.func.dir, 'FFP_calc_parallel.R'))       # Parallel execution of FFP calculation and analysis of the interactions with spatial data
  source(paste0(FFP.func.dir, 'FFP_calc_parallel_dev.R'))     # Parallel execution of FFP calculation, libraries/functions UPDATED to replace 'rgdal'
  # TODO 20221123: FARE FUNZIONE CHE LAVORA IN SERIALE DA USARE IN CASO QUELLA PARALLELA FALLISCE
  
  # Helper: funzione per leggere KML (se più layer, prendi il primo se non trovi il layer con 'TA')
  read_kml_first_layer <- function(kml_path) {
    layers <- sf::st_layers(kml_path)$name
    chosen <- layers[1]
    # prefer layer con 'TA' o 'Target' se presente
    idx <- grep("TA|Target", layers, ignore.case = TRUE)
    if (length(idx) > 0) chosen <- layers[idx[1]]
    sf::st_read(dsn = kml_path, layer = chosen, quiet = TRUE)
  }
  
  # FFP 2D function [Kljun et al. 2015]
  source(paste0(FFP.func.dir, 'calc_footprint_FFP.R'))
  source(paste0(FFP.func.dir, 'repmat.R'))
  source(paste0(FFP.func.dir, 'calc_footprint_FFP_climatology.R'))
  
  
  # Directories ..................................................................
  # Store the current project directory
  
  
  if (is.null(FFP.dir)) {
    proc.date <- format(Sys.Date(), '%Y%m%d')
    # Footprint directory: will contain the input spatial data footprint analysis results
    proj.dir <- getwd()
    FFP.dir <- paste0(proj.dir, '/', site.ID, '_output_FP_', proc.date, '/'); suppressWarnings(dir.create(FFP.dir))
    # Micrometeorological dataset input file
    MM.data.dir <- paste0(proj.dir, '/', site.ID, '_output_RF_', proc.date, '/'); suppressWarnings(dir.create(MM.data.dir))
  }
  
  # Spatial features directory (as sub-folder of the main directory)
  if (is.null(SP.data.dir)) {
    SP.data.dir <- paste0(FFP.dir, 'spatial_data/'); suppressWarnings(dir.create(SP.data.dir))
  }
  
  
  # Local settings
  # FFP domain array
  ffp.domain <- ffp.domain * c(-1,1,-1,1)
  # FFP isoplethes on whic the analysis is based: this is a position index of the FFP.R vector
  X.perc.indx <- which(FFP.R %in% c(70, 80))  # ADATTARE CODICE IN MODO DA PREVEDERE QUALSIASI PERC. ATTUALMENTE 70 E 80 CI DEOVNO ESSERE PER FORZA
  
  # color palette 
  ffp2Dcol <- c('transparent', rich.colors(100, alpha = 0.5))
  
  # Create a vector which will contain the record numbers and stats for data to use for the FFP
  FFP.data.summary <- c(NULL)
  
  
  # STATIC metadata  ============================================================
  # Load station metadata from an ad-hoc web-service. HTML table from database.
  # This is the same table used for QC metadata (yet dynamically read online)
  
  
  if (ecmd.available) { # if station metadata are available
    
    if (!is.null(ecmd.path)) {
      # Read it
      ecmd.file <- grep('ecmd', list.files(ecmd.path, full.names = TRUE), value = TRUE)
      EC.MD <- fread(ecmd.file, sep=',', header=T, data.table = F)
    } else {
      # Download the ecmd table from the server
      EC.MD <- fread(paste0(scan('ETCWebService', what = 'character')[1], tolower(site.ID)), sep=',', header=T, data.table = F)
    }
    
    if (nrow(EC.MD) > 0) {
      
      # station coordinates
      # NOTE: currently only the last position is considered
      st.lat <- EC.MD$'LATITUDE'[nrow(EC.MD)]
      st.lon <- EC.MD$'LONGITUDE'[nrow(EC.MD)]
      st.alt <- EC.MD$'ALTITUDE'[nrow(EC.MD)]
      
      ## set site-specific spatial reference in UTM
      site.UTM.proj <- paste("+proj=utm +zone=", long2UTM(st.lon), " +ellps=WGS84 +units=m +datum=WGS84 +no_defs", sep='')
      
    } else {
      
      cat(error('\nError: Station metadata are NOT available. '))
      stop()
    }
    
  } else {
    
    cat(error('\nError: Station metadata must be supplied. '))
    stop()
    
  }
  
  
  # MICROMET data ===================================================================
  # Import half-hourly data (fluxes and micrometeorological variables) processed in EP
  # and quality controlled by mean of the ICOS QC procedure
  # This file has the same structure as the fulloutput di EddyPro, except for the first 4 columns
  # replaced by TIMESTAMP ISO (yyyymmddHHMM) referred to the start of the averaging period
  # X_cleaned: X flux, filtered by the whole error flags
  # X_QC_flag: 0 retained record, 1 rejected record (for whatever reason) - USE THIS FOR FFP DATA FILETRING
  # X_Er_flag: 0 No Er, 1 Mod Er, 2 Sev Er
  # X_Outlying: 0 no outlier, 1 outlier
  
  
  if (QC.done) { # If the data Quality Control has been done
    
    # In case an URL is provided, download the QC file from Drive
    # TODO: TO BE CHECKED ON CLUSTER
    if (!is.null(MM.data.url)) {
      if (grepl('https', MM.data.url)) {
        GD_files_dwld(GD_folder_url = MM.data.url,
                      site.ID = site.ID,
                      output_folder = MM.data.dir,
                      pattern = 'clean',
                      return.folder.path = FALSE,
                      download = TRUE)  
      } else {
        cat(warn('\nWARNING: The URL to the MM dataset must be properly provided. PLease check it out and retry.\n'))
      }
    }
    
    ## read it
    if (any(grepl(MM.data.patern, list.files(MM.data.dir)))) {
      set.cleaned <- data.table::fread(paste0(MM.data.dir, grep('clean', list.files(MM.data.dir), value = TRUE)),
                                       sep=',', header = TRUE, na.strings = c('NA','-9999','-9999.0'), data.table = FALSE)
    } else {
      
      cat(error('\nError: a dataset containing micromet variables must be supplied.\n'))
      stop()
    }
    
    # Time specification
    # NOTE: Refer time to the following half-hour
    set.cleaned.timestamp <- strptime(set.cleaned$'TIMESTAMP_END', format = '%Y%m%d%H%M', tz='GMT')
    
    # start filling the data summary log
    FFP.data.summary$'date_start' <- format(first(set.cleaned.timestamp), '%Y%m%d%H%M%S') # first day of data
    FFP.data.summary$'date_end' <- format(last(set.cleaned.timestamp), '%Y%m%d%H%M%S') # last day of data
    FFP.data.summary$'n_days' <- ndays(set.cleaned.timestamp) # number of days of data
    FFP.data.summary$'N_raw_data' <- nrow(set.cleaned) # number of half-hours
    
    
    # . NEE data filtering mask (QC flag) ----
    # from QC policy:
    # 0 = good data
    # 1 = moderate flag + outlier;
    # 2 = severe flag
    # -9999 = originally missing data
    NEE.QC.flag <- set.cleaned$'NEE_DATA_FLAG' %in% c(1,2)  # TRUE: remove data point, FALSE: retain data point
    # table(NEE.QC.flag)/nrow(set.cleaned)*100 # percentage of data removed by QC
    FFP.data.summary$'N_QC_flag' <- unname(table(NEE.QC.flag)['TRUE']) # data flagged for bad NEE
    
    # . Horizontal wind component homogeneity test (QC flag), TRUE = disomogeneity ----
    # TODO: ATTUALMENTE NON PIU RIPORTATO IN cleanset > sentire domenico
    HWHT.flag <- set.cleaned$'FootprintFlag' == 1 # new flag name
    # table(HWHT.flag)/nrow(set.cleaned) * 100 # percentage of data flagged for hunomoganeities in horizontal wind direction (Domenico's custom test)
    # data flagged for wind inhomogeneiies
    FFP.data.summary$'N_HWHT_flag' <- unname(table(HWHT.flag)['TRUE'])
    # data flagged for wind inhomogeneiies & NOT yet excluded by the NEE QC
    FFP.data.summary$N_HWHTf_over_QCf <- ifelse(is.na(unname(table(HWHT.flag[!NEE.QC.flag])['TRUE'])), 0, unname(table(HWHT.flag[!NEE.QC.flag])['TRUE']))
    
    # . USTAR data filtering mask ----
    if (DO.ustar.filter){
      # NOTE: This was the original agreement, now implemented with a custom threshold
      # TRUE: remove data point, FALSE: retain data point
      USTAR.QC.flag <- set.cleaned$'USTAR' < ustar.thr
      # table(USTAR.QC.flag)/nrow(set.cleaned) # percentage of data flagged low turbulence
      FFP.data.summary$'N_USTAR_QC_flag' <- unname(table(USTAR.QC.flag)['TRUE']) # data flagged for u* below the threshold
      FFP.data.summary$'N_USTAR_QC_flag_over_NEEQC' <- ifelse(is.na(unname(table(USTAR.QC.flag[!NEE.QC.flag])['TRUE'])), 0, unname(table(USTAR.QC.flag[!NEE.QC.flag])['TRUE'])) # data flagged for u* below the threshold & NOT yet excluded by the NEE QC
    }
    
    # . Potential number of data over which perform FFP computation ----
    if (DO.ustar.filter) {
      FFP.data.summary$'N_data_FFP_pot' <- FFP.data.summary$'N_raw_data' - (FFP.data.summary$'N_QC_flag' + FFP.data.summary$'N_USTAR_QC_flag_over_NEEQC')
    }else{
      FFP.data.summary$'N_data_FFP_pot' <- FFP.data.summary$'N_raw_data' - FFP.data.summary$'N_QC_flag'
    }
    # If percentage of original data, and if less than Nrec.min [%], STOP
    FFP.data.summary$'N_data_FFP_pot_perc' <- round(FFP.data.summary$'N_data_FFP_pot' / FFP.data.summary$'N_raw_data' * 100, 1)
    if(FFP.data.summary$'N_data_FFP_pot_perc' < Nrec.min){
      cat(error('\nWARN: there are not enough data to get a robust footprint estimation\n'))
      stop()
    } else {
      cat(note('\nOK, there are enough data to get a robust footprint estimation\n'))
    }
    
  } else { # in case the QC is not available
    
    cat(error('\nError: A micrometeorological dataset must be provided. The analysis cannot be done.'))
    stop()
    
  }
  
  
  # DYNAMIC metadata ============================================================
  # NOTE: these data are currently taken directly from the input dataset (set.cleaned, QC controllod)
  
  # Average vegetation height [m], may change during the period
  H.canopy <- set.cleaned[, grep('canopy', names(set.cleaned), ignore.case = TRUE)] # QC dataset
  if (median(H.canopy, na.rm=T) != median(as.numeric(EC.MD$'CANOPY_HEIGHT'), na.rm=T)) {
    cat(warn('\nThe median canopy heights used for flux computation and QC do not match the one in the metadata table.'))
    cat(note('\nThe formaer will be used.\n'))
  }
  
  # Height of the eddy covariance system [m], may change during the period
  H.EC <- set.cleaned[, grep('sa_hei', names(set.cleaned), ignore.case = TRUE)] # QC dataset
  if (median(H.EC, na.rm=T) != median(as.numeric(EC.MD$'SA_HEIGHT'))) {
    cat(warn('\nEC system heigh used for QC is different from that in the metadata table.'))
  }
  
  # Displacement height
  d. <- H.canopy * 2/3  # EddyPro computes the displacement height as: d = 0.67 x canopyheight
  
  # Roughness length
  z0. <- H.canopy * 0.15  # EddyPro computes the roughness length as: z0 = 0.15 x canopyheight
  
  # Dynamic measurement height (i.e. z‐d) [m]
  zm. <- H.EC - d.
  
  
  # SPATIAL DATA =================================================================
  # TA, EA, LCT* and EC tower coordinates, downloaded from ICOS database
  # * Land Cover Typologies: Vegetation characteristics as homogeneous patches in *.kml files made by ETC-BE basing on vegetation analysis results
  
  # In case spatial features data are not provided, try download them from the ICOS-ETC server
  if (length(list.files(SP.data.dir)) == 0) {
    
    cat(prog.mes('\nNOTE: spatial data were not provided. Trying to download them from the ICOS database ...\n'))
    cat(note("Trying fetching the spatial features as *kml files\n"))
    
    # Start trying with *.kml files first ...................
    SPd.dwnld.went.wrong <- NULL
    SPd.dwnld.went.wrong <- tryCatch(
      SpDataDL_poly(site.ID,
                    data.dir = SP.data.dir,
                    file.type = 'kml',
                    data.path = FALSE,
                    download = TRUE),
      error = function(e) {e}
    )
    
    # If the kml query succeeded AND the TA object was actually downloaded  
    if (!inherits(SPd.dwnld.went.wrong, "error") & any(grepl('_TA', list.files(SP.data.dir), ignore.case = F))) { 
      spd.ok <- TRUE
      spd.type <- 'klm'
      cat(prog.mes('done.\n'))
      rm(SPd.dwnld.went.wrong)
      
    } else { 
      
      # Then try with *.shp files ...................
      cat(note("WARNIING [SpDataDL_poly]: the spatial features *kml files were not fetched. Trying fetching them as shape files ..."))
      
      # delete all the temporary objects possibly downloaded but the LCTs
      unlink(grep('LCT', list.files(SP.data.dir, full.names = T), value = T, invert = T)) 
      
      SPd.dwnld.went.wrong <- ''
      SPd.dwnld.went.wrong <- tryCatch(
        SpDataDL_poly(site.ID,
                      data.dir = SP.data.dir,
                      file.type = 'shape',
                      data.path = FALSE,
                      download = TRUE),
        error = function(e) {e}
      )
      
      # If the SHP query succeeded AND the TA object was actually downloaded
      if (!inherits(SPd.dwnld.went.wrong, "error")  & any(grepl('_TA', list.files(SP.data.dir), ignore.case = F))) {
        spd.ok <- TRUE
        spd.type <- 'shp'
        cat(prog.mes('done.\n'))
        rm(SPd.dwnld.went.wrong)
        
      } else {
        
        # If an error persist in the download, even with shapes objects
        spd.ok <- FALSE
        cat(error('\nERROR [ETCPS::SpDataDL_poly]: an error occurred while downloading the spatial data! Data processing cannot be completed.'))
        cat(note(paste('\nETCPS::SpDataDL_poly Error:', SPd.dwnld.went.wrong$'message\n')))
        stop('\n')
      }
      
    }
    
  } # Spatial data download chunk END
  
  
  # . Target Area (TA) ---------------------------------
  
  TA.filelist <- list.files(SP.data.dir, pattern = 'TA')
  if (length(TA.filelist) == 0) {
    stop("No Target Area (TA) files found in SP.data.dir")
  }
  
  # # OR_code
  # if (any(grepl('*.kml', TA.filelist))) { # If input file is a KML 
  #   # set kml name
  #   TA.kmlname <- TA.filelist[grepl('*.kml', TA.filelist)]
  #   # read kml
  #   TA.layername <- rgdal::ogrListLayers(dsn = paste0(SP.data.dir, TA.kmlname))[1]
  #   if (!grepl('TA|*Target*', TA.layername, ignore.case = TRUE)) {TA.layername <- substr(TA.kmlname, 1, 9) }
  #   # open kml 
  #   TA <- suppressWarnings(rgdal::readOGR(dsn = paste0(SP.data.dir, TA.kmlname), layer = TA.layername, verbose = T))
  # }
  
  # KML 
  if (any(grepl('\\.kml$', TA.filelist, ignore.case = TRUE))) {
    TA.kmlname <- TA.filelist[grepl('\\.kml$', TA.filelist, ignore.case = TRUE)]
    # se ci sono più file prendo il primo: adattare se necessario
    TA.kmlpath <- file.path(SP.data.dir, TA.kmlname[1])
    TA <- suppressWarnings(read_kml_first_layer(TA.kmlpath))
    # assicurati che abbia un CRS; se no usa EPSG:4326 (KML usualmente lon/lat WGS84)
    if (is.na(sf::st_crs(TA))) sf::st_crs(TA) <- 4326
  }
  
  # # # OR_code
  # if (any(grepl('*.shp', TA.filelist))) {  # If input file is a SHP
  #   if(!any(grepl('*.shx', TA.filelist))){warning('There is no shape index format file for Target Area'); break()}
  #   if(!any(grepl('*.dbf', TA.filelist))){warning('There is no attribute format file for Target Area'); break()}
  #   if(!any(grepl('*.prj', TA.filelist))){warning('There is no projection format file for Target Area'); break()}
  #   # set shape name
  #   TA.shpname <- TA.filelist[grepl('*.shp', TA.filelist)]
  #   # open the shape file 
  #   TA <- readOGR(dsn = paste0(SP.data.dir, TA.shpname))
  # } 
  
  # Shapefile
  if (any(grepl('\\.shp$', TA.filelist, ignore.case = TRUE))) {
    TA.shpname <- TA.filelist[grepl('\\.shp$', TA.filelist, ignore.case = TRUE)]
    TA.shppath <- file.path(SP.data.dir, TA.shpname[1])
    TA <- sf::st_read(dsn = TA.shppath, quiet = TRUE)
    if (is.na(sf::st_crs(TA))) {
      # se manca .prj provo ad attribuire la proiezione originale salvata in proj.orig (se definita)
      if (exists("proj.orig") && !is.null(proj.orig)) {
        sf::st_crs(TA) <- proj.orig
      } else {
        stop("Shapefile TA senza .prj: specificare un CRS.")
      }
    }
  }
  
  
  # # Site Coordinate Reference Systems
  # proj.orig <- TA@proj4string
  # Save original proj (come nel codice precedente)
  proj.orig <- sf::st_crs(TA)
  
  
  # # Re-project to long/lat and add DATUM projection to comply with KML requirements
  # TA <- sp::spTransform(TA, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
  # Reproject to long/lat WGS84 for KML compatibility
  TA <- sf::st_transform(TA, crs = 4326)
  
  # # TA coordinates in UTM
  # TA.utm <- sp::spTransform(TA, CRS(site.UTM.proj))
  # TA coordinates in UTM (usa la stringa site.UTM.proj creata sopra)
  TA.utm <- sf::st_transform(TA, crs = site.UTM.proj)
  
  
  # # . EC tower ------------------------------------.
  # if (exists('st.lon')) {
  #   # EC tower spatial object from ecmd data
  #   EC.tower <- SpatialPoints(list(lon = st.lon, lat = st.lat), proj4string = TA@proj4string)
  # } else {
  #   # If not present, try catching EC tower coordinates from DB
  #   # TODO: CHECK THIS ON CLUSTER
  #   EC.tower <- SpDataDL_points(site.ID, p4s = as.character(TA@proj4string))$'EC'
  # }
  # 
  # # EC tower coordinates in UTM: USE THIS AS POINT 0 FOR FFP
  # EC.tower.utm  <- sp::spTransform(EC.tower, CRS(site.UTM.proj))
  
  # . EC tower ------------------------------------
  if (exists('st.lon')) {
    # Creazione punto EC come sf
    EC.tower <- sf::st_as_sf(data.frame(lon = st.lon, lat = st.lat),
                             coords = c("lon","lat"), crs = sf::st_crs(TA))
  } else {
    # fallback: scarica dal DB (assumi che SpDataDL_points ritorni oggetti compatibili)
    # se la funzione ritorna sp objects, convertili in sf:
    pts_list <- SpDataDL_points(site.ID, p4s = as.character(sf::st_crs(TA)))
    if (!is.null(pts_list$EC)) {
      EC.tower <- pts_list$EC
      if (!inherits(EC.tower, "sf")) EC.tower <- sf::st_as_sf(EC.tower)
      sf::st_crs(EC.tower) <- sf::st_crs(TA)
    } else {
      stop("EC tower coordinates not available")
    }
  }
  
  # EC tower in UTM
  EC.tower.utm <- sf::st_transform(EC.tower, crs = site.UTM.proj)
  
  
  # . Exclusion Are(s) -----------------------------------------------------------
  
  EA.filelist <- list.files(SP.data.dir, pattern = 'EA')
  
  # # OR_code
  # if (length(EA.filelist) > 0) { # check if there are EA, if yes, load them
  #   
  #   if (any(grepl('*.kml', EA.filelist))) { # If input file is a KML
  #     
  #     # EA number (may be multiple files)
  #     EA.file.num <- length(which(grepl('*.kml', EA.filelist)))
  #     
  #     if (EA.file.num == 1) {  # If only one file exist for EA
  #       
  #       EA.kmlname <- EA.filelist[grepl('*.kml', EA.filelist)]
  #       EA.layername <- rgdal::ogrListLayers(dsn = paste(SP.data.dir, EA.kmlname, sep=""))[1]
  #       EA <- suppressWarnings(rgdal::readOGR(dsn = paste(SP.data.dir, EA.kmlname, sep=""), layer = EA.layername, verbose = F))
  #       # Assign to EA the same Target area projections
  #       EA <- sp::spTransform(EA, TA@proj4string)
  #       # get the actual EA number (can be merged during spatial sampling processing)
  #       EA.num <- length(EA)
  #       
  #     } else {    # If multiple  files exist for EA
  #       
  #       EA.kmlname <- EA.filelist[grepl('*.kml', EA.filelist)]
  #       EA.list <- list(NULL)
  #       i.ea <- 1
  #       for (i.ea in 1:EA.file.num) {
  #         cur.EA.kmlname <- rgdal::ogrListLayers(dsn = paste(SP.data.dir, EA.kmlname[i.ea], sep=""))[1]
  #         EA.list[[i.ea]] <- suppressWarnings(rgdal::readOGR(dsn = paste(SP.data.dir, EA.kmlname[i.ea], sep=""), layer = cur.EA.kmlname, verbose = F))
  #         # Assign to EA the same Target area projections
  #         EA.list[[i.ea]] <- sp::spTransform(EA.list[[i.ea]], TA@proj4string)
  #         # Using a zero-width buffer cleans up many topology problems in R.
  #         if(!gIsValid(EA.list[[i.ea]])){
  #           EA.list[[i.ea]] <- rgeos::gBuffer(EA.list[[i.ea]], byid=TRUE, width=0)
  #         }
  #       }
  #       # Merge all the excluded areas (Note: the loop is needed as union accept 2 polygons only)
  #       i.ea <- 1
  #       for (i.ea in 1:length(EA.list)) {
  #         if (i.ea == 1) {
  #           EA <- EA.list[[i.ea]]
  #         } else {
  #           # if (site.ID == 'CH-Dav') { # NOTA: union non funziona in ch-dav perchè la nuov EA è un track: faccio così ma perdo individualità
  #           # EA <- rbind(EA, EA.list[[i.ea]])
  #           # } else {
  #           # EA <- raster::union(EA, EA.list[[i.ea]])
  #           EA <- rgeos::gUnion(EA, EA.list[[i.ea]])
  #           # }
  #         }
  #       }
  #       if(site.ID == 'CH-Dav'){
  #         # merge possible multiple polygons (downgrading to SpatialPolygon)
  #         EA <- aggregate(EA, dissolve=T)
  #         # Reconvert data frame to SpatialPolygons
  #         EA <- SpatialPolygonsDataFrame(EA, data = data.frame(1), match.ID = F)
  #       }
  #       # Assign to EA the same Target area projections
  #       EA <- sp::spTransform(EA, TA@proj4string)
  #       # get the actual EA number (can be merged during spatial sampling processing)
  #       EA.num <- length(EA)
  #     } # multiple EA loop end
  #     
  #   } # EA kml loop end
  #   
  #   
  #   if (any(grepl('*.shp', EA.filelist))) { # If input file is a SHP
  #     
  #     # EA number (may be multiple files)
  #     EA.num <- length(which(grepl('*prj', EA.filelist)))
  #     
  #     if (EA.num == 1) { # If only one file exist for EA
  #       # To be properly read a shape file must be accompained by a minimum set of mandatory files:
  #       # .shp — shape format, the feature geometry itself
  #       # .shx — shape index format, a positional index of the feature geometry to allow seeking forwards and backwards quickly
  #       # .dbf — attribute format, columnar attributes for each shape, in dBase IV format
  #       # .prj — projection format, the coordinate system and projection information, a plain text file describing the projection using well-known text format
  #       if(!any(grepl('*.shx', EA.filelist))){warning('There is no shape index format file for Exclusion Area'); break()}
  #       if(!any(grepl('*.dbf', EA.filelist))){warning('There is no attribute format file for Exclusion Area'); break()}
  #       if(!any(grepl('*.prj', EA.filelist))){warning('There is no projection format file for Exclusion Area'); break()}
  #       # set shape name
  #       EA.shpname <- EA.filelist[grepl('*.shp', EA.filelist)]
  #       EA <- readOGR(dsn = paste(dir.input, EA.shpname, sep=""))
  #       # Assign to EA the same Target area projections 
  #       # EA@proj4string
  #       EA <- spTransform(EA, TA@proj4string)
  #       
  #     } else {    # If multiple  files exist for EA
  #       
  #       # set kml names
  #       EA.shpname <- grep('.xml', EA.filelist[grepl('*.shp', EA.filelist)], invert = T, value = T)
  #       # read all of them
  #       EA.list <- list(NULL)
  #       i <- 0
  #       for (i in 1:EA.num) {
  #         # set actual kml names and read it
  #         EA.list[[i]] <- readOGR(dsn = paste(SP.data.dir, EA.shpname[i], sep=""))
  #         # Assign to EA the same Target area projections
  #         # EA.list[[i]]@proj4string
  #         if(!is.na(EA.list[[i]]@proj4string)){     # se i file hanno le proiezioni, trasformo semplicemente...
  #           EA.list[[i]] <- spTransform(EA.list[[i]], TA@proj4string)
  #         } else {                                  # ...ALTRIMENTI, gli forzo quelle della TA e trasformo
  #           EA.list[[i]]@proj4string <- proj.orig
  #           EA.list[[i]] <- spTransform(EA.list[[i]], TA@proj4string)
  #         }
  #       }
  #       # Merge all the excluded areas Note: the loop is needed as union accept 2 polygons only
  #       i <- 1
  #       for(i in 1:length(EA.list)){
  #         if(i == 1){
  #           EA <- EA.list[[i]]
  #         }else{
  #           if(all(names(EA.list[[i]]) %in% names(EA))){
  #             cat(yellow(black$bgYellow$bold('WARNING:') %+% "'union' cannot modify data containing duplicate names. Some names of EA n.", i, " were thus modified adding '_1' to the original name"))
  #             name.change.id <- which(names(EA.list[[i]]) %in% names(EA))
  #             names(EA.list[[i]])[name.change.id] <- paste(names(EA.list[[i]])[name.change.id], '_', i, sep='')
  #           }
  #           EA <- union(EA, EA.list[[i]])
  #         }
  #       }
  #     }
  #   } # EA shp loop end
  #   
  #   EA.color <- 'red'
  #   
  # } else {
  #   
  #   cat(note('\nExclusion Areas were NOT submitted.'))
  #   EA.num <- NA # set it to NA to use it as pointer
  #   
  # }
  # 
  # # EA coordinates in UTM
  # if (!is.na(EA.num)) {
  #   EA.utm <- sp::spTransform(EA, CRS(site.UTM.proj))
  #   if (length(EA.utm) != 1) {
  #     # Merge the possible EA patches into one (NEEDED IN CASE OF EA ANALYSIS)
  #     EA.utm.uni <- EA.utm[1, ]
  #     for (i.ea in 2:length(EA.utm)) {
  #       EA.utm.uni <- rgeos::gUnion(EA.utm.uni, EA.utm[i.ea, ])
  #     }
  #     EA.utm <- EA.utm.uni
  #   }
  # }
  
  if (length(EA.filelist) > 0) {
    
    # KML EAs
    if (any(grepl('\\.kml$', EA.filelist, ignore.case = TRUE))) {
      
      EA.kmls <- EA.filelist[grepl('\\.kml$', EA.filelist, ignore.case = TRUE)]
      EA.list <- lapply(EA.kmls, function(fname) {
        p <- file.path(SP.data.dir, fname)
        # legge il primo layer o quello con 'EA' se presente
        layers <- sf::st_layers(p)$name
        chosen <- layers[1]
        idx <- grep("EA|Exclusion", layers, ignore.case = TRUE)
        if (length(idx) > 0) chosen <- layers[idx[1]]
        e <- suppressWarnings(sf::st_read(dsn = p, layer = chosen, quiet = TRUE))
        # assegna CRS della TA se mancante
        if (is.na(sf::st_crs(e))) sf::st_crs(e) <- sf::st_crs(TA)
        e <- sf::st_transform(e, crs = sf::st_crs(TA))
        # clean invalid geometries
        if (!all(sf::st_is_valid(e))) e <- sf::st_make_valid(e)
        e
      })
    } # kml chunck end
    
    # SHP EAs
    if (any(grepl('\\.shp$', EA.filelist, ignore.case = TRUE))) {
      EA.shps <- EA.filelist[grepl('\\.shp$', EA.filelist, ignore.case = TRUE)]
      EA.list <- lapply(EA.shps, function(fname){
        p <- file.path(SP.data.dir, fname)
        e <- sf::st_read(dsn = p, quiet = TRUE)
        if (is.na(sf::st_crs(e))) sf::st_crs(e) <- proj.orig
        e <- sf::st_transform(e, crs = sf::st_crs(TA))
        if (!all(sf::st_is_valid(e))) e <- sf::st_make_valid(e)
        e
      })
    } # shape chunk end
    
    
    # EA manipulation to combine the elements 
    
    # 1. verifico che tutti gli elementi abbiano lo stesso CRS, altrimenti lo forzo
    if (length(unique(lapply(EA.list, st_crs))) > 1) {
      EA.list <- lapply(EA.list, st_transform, crs = 4326)  # WGS84
    }
    
    # 2. porto tutto in 2D e stesso tipo geometria
    #    scarto tutto il resto
    
    # EA.list.clean <- lapply(EA.list, function(x) {
    #   x <- st_zm(x, drop = TRUE, what = "ZM")            # togli Z se presente (diventa XY)
    #   x <- st_cast(x, "MULTIPOLYGON")                    # porta tutto a MULTIPOLYGON
    #   x
    # })
    EA.list.clean <- lapply(EA.list, 
                            function(x) {
                              # 1) porta tutto in 2D
                              x <- st_zm(x, drop = TRUE, what = "ZM")
                              # 2) tieni solo POLYGON / MULTIPOLYGON
                              geom_type <- unique(st_geometry_type(x))
                              if (!geom_type %in% c("POLYGON", "MULTIPOLYGON")) {
                                return(NULL)   # scarta LINESTRING ecc.
                              }
                              # 3) cast sicuro
                              x <- st_cast(x, "MULTIPOLYGON")
                              x}
                            )
    # rimuovo i NULL (LINESTRING eliminati)
    EA.list.clean <- Filter(Negate(is.null), EA.list.clean)
    
    # 3. estraggo le geometrie
    EA_geoms <- lapply(EA.list.clean, st_geometry)
    
    # 4. estraggo attributi come data.frame, rimuovendo la colonna geometry
    EA_attrs <- lapply(EA.list.clean, function(x) st_drop_geometry(x))
    
    # 5. unisco gli attributi con full-join delle colonne (bind_rows fa già questo)
    EA_attrs_combined <- dplyr::bind_rows(EA_attrs)
    
    # 6. combino le geometrie in una sola sfc
    EA_geoms_combined <- do.call(c, EA_geoms)
    
    # 7. ricostruisco EA sf (assumendo stesso CRS per tutti)
    crs_ea <- st_crs(EA.list.clean[[1]])
    EA <- st_sf(EA_attrs_combined, geometry = EA_geoms_combined, crs = crs_ea)
    
    # store the new EA number
    EA.num <- length(EA)
    
    EA.color <- 'red'
    
  } else {
    cat(note('\nExclusion Areas were NOT submitted.'))
    EA.num <- NA
  }
  
  # EA coordinates in UTM (se presenti)
  if (!is.na(EA.num)) {
    EA.utm <- sf::st_transform(EA, crs = site.UTM.proj)
    # se più poligoni, unisco in uno solo
    if (length(EA.utm) > 1) {
      EA.utm <- sf::st_union(EA.utm)
      EA.utm <- sf::st_sf(geometry = EA.utm)
    }
  }
  
  
  
  # . Land Cover Typologies -----------------------
  # These are the LCT categories that will be used for the analysis of the
  # "representativeness of different ecosystems included in the Target Area". 
  # The polygons are drawn based on the site characterization results as achieved and reported by the PI.
  # Vegetation typology is relative to each CP or SPI and the concerning variable
  # differs according to the actual ecosystem (data to be downloaded from http://www.europe-fluxdata.eu/icosdataview.aspx):
  # The polygons are made by the Belgian ETC Team and send to Italian ETC
  
  # List (possible) LCT polygons
  LCT.filelist <- list.files(SP.data.dir, pattern = 'LCT')
  
  if (length(LCT.filelist) != 0) {
    LCT.present <- TRUE
    LCT.resultant <- FALSE # modify this in case needed
  } else {
    LCT.present <- FALSE
  }
  
  # OR_code
  # # LCT import and manipulation (polygons are made in GoogleEarth by ETC-BE)
  # if (LCT.present) {
  #   
  #   # LCT number (may be multiple files)
  #   LCT.num <- length(which(grepl('*.kml', LCT.filelist)))
  #   LCT.kmlname <- LCT.filelist[grepl('*.kml', LCT.filelist)]
  #   
  #   LCT.list <- list(NULL)
  #   i.lct <- 1
  #   for (i.lct in 1:LCT.num) {
  #     cur.LCT.kmlname <- rgdal::ogrListLayers(dsn = paste(SP.data.dir, LCT.kmlname[i.lct], sep=""))[1]
  #     LCT.list[[i.lct]] <- suppressWarnings(rgdal::readOGR(dsn = paste(SP.data.dir, LCT.kmlname[i.lct], sep=""), layer = cur.LCT.kmlname))
  #     # Project LCT to the same Target area projections
  #     LCT.list[[i.lct]] <- sp::spTransform(LCT.list[[i.lct]], TA@proj4string)
  #   }
  #   
  #   # Possibly define a further LCT as resultant of the difference between TA and the actual LCT
  #   if (LCT.resultant) {
  #     
  #     # Merge all the LCT areas just to compute their sum (Note: the loop is needed as union accept 2 polygons only)
  #     ii <- 0
  #     for(ii in 1:length(LCT.list)){
  #       if(ii == 1){
  #         LCT.tmp <- LCT.list[[ii]]
  #       } else {
  #         LCT.tmp <- raster::union(LCT.tmp, LCT.list[[ii]])
  #       }
  #     }
  #     LCT.list[[length(LCT.list)+1]] <- TA - LCT.tmp
  #     
  #     # update the LCT number according to the resultant LCT
  #     LCT.num <- length(LCT.list)
  #     
  #     # remove the tmp object
  #     rm(LCT.tmp)
  #   }
  #   
  #   # set LCT colors
  #   LCT.color <- terrain.colors(LCT.num+3)[-c((LCT.num+1):(LCT.num+3))]
  #   
  #   ## check and possibly correct for topology problems
  #   if (any(!unlist(lapply(LCT.list, function(x) rgeos::gIsValid(x, byid=TRUE))))) {
  #     LCT.list <- lapply(LCT.list, function(x) rgeos::gBuffer(x, byid=TRUE, width=0))
  #   }
  #   # check correction
  #   # lapply(LCT.list, function(x) gIsValid(x, byid=TRUE))
  #   
  #   
  #   ### LCT coordinates in UTM
  #   if (!is.na(LCT.num)) { 
  #     LCT.list.utm <- lapply(LCT.list, function(x) sp::spTransform(x, CRS(site.UTM.proj)))
  #   }
  #   
  #   ### check and possibly correct for topology problems
  #   if (any(!unlist(lapply(LCT.list.utm, function(x) rgeos::gIsValid(x, byid=TRUE))))) {
  #     LCT.list.utm <- lapply(LCT.list.utm, function(x) rgeos::gBuffer(x, byid=TRUE, width=0))
  #   }
  #   # check correction
  #   # lapply(LCT.list.utm, function(x) gIsValid(x, byid=TRUE))
  #   
  #   
  # } else {
  #   
  #   cat(note('\nNO Vegetation Polygons (LCT) have been identified (or at least reported)'))
  #   LCT.num <- NA # set it to NA to use it as pointer
  #   LCT.list.utm <- NULL
  #   
  # }
  
  # NEW_code
  if (LCT.present) {
    LCT.kmls <- LCT.filelist[grepl('\\.kml$', LCT.filelist, ignore.case = TRUE)]
    LCT.list <- lapply(LCT.kmls, function(fname){
      p <- file.path(SP.data.dir, fname)
      layers <- sf::st_layers(p)$name
      chosen <- layers[1]
      if (any(grepl("LCT", layers, ignore.case = TRUE))) {
        idx <- grep("LCT", layers, ignore.case = TRUE)
        chosen <- layers[idx[1]]
      }
      l <- suppressWarnings(sf::st_read(dsn = p, layer = chosen, quiet = TRUE))
      if (is.na(sf::st_crs(l))) sf::st_crs(l) <- sf::st_crs(TA)
      l <- sf::st_transform(l, crs = sf::st_crs(TA))
      if (!all(sf::st_is_valid(l))) l <- sf::st_make_valid(l)
      l
    })
    if (LCT.resultant) {
      # Merge all LCT to compute their union
      LCT.tmp <- Reduce(function(a,b) sf::st_union(a,b), LCT.list)
      # resultant patch = TA - union(LCT)
      res_geom <- sf::st_difference(TA, LCT.tmp)
      LCT.list[[length(LCT.list)+1]] <- sf::st_sf(geometry = res_geom)
    }
    
    # couleurs
    LCT.num <- length(LCT.list)
    LCT.color <- terrain.colors(LCT.num+3)[-c((LCT.num+1):(LCT.num+3))]
    
    # LCT in UTM
    LCT.list.utm <- lapply(LCT.list, function(x) sf::st_transform(x, crs = site.UTM.proj))
    # validazione
    LCT.list.utm <- lapply(LCT.list.utm, function(x) {
      if (!all(sf::st_is_valid(x))) x <- sf::st_make_valid(x)
      x
    })
  } else {
    cat(note('\nNO Vegetation Polygons (LCT) have been identified (or at least reported)'))
    LCT.num <- NA
    LCT.list.utm <- NULL
  }
  
  
  # # OR_code
  # # Site spatial features visualization
  # if (plot.spatial.feats) {
  #   
  #   png(paste(FFP.dir, site.ID, '_','SpatialFeatures_', format(Sys.Date(), '%Y%m%d'), '.png', sep=''),
  #       width = 14, height = 14, unit = 'cm', res=300)
  #   # windows()
  #   
  #   # TA + EC
  #   plot(TA.utm, axes=T)
  #   points(EC.tower.utm, pch=3)
  #   
  #   # EA
  #   if(!is.na(EA.num)){
  #     if(EA.num == 1){
  #       plot(EA.utm, add=T, col=adjustcolor(EA.color, alpha.f = 0.6), border=EA.color)
  #     }else{
  #       j <- 0
  #       for(j in 1:length(EA)){
  #         plot(EA.utm[j, ], add=T, col=adjustcolor(EA.color[j], alpha.f = 0.6), border=EA.color[j])
  #         text(x=coordinates(EA.utm[j, ])[1], y=coordinates(EA.utm[j, ])[2], EA.utm[j, ]@data$Description, pos=1, cex=.5)
  #       }
  #     }
  #     legend('topleft', legend=c('EC','TA',paste0('EA_0',seq_along(EA))), col = c(1,1, adjustcolor(EA.color, alpha.f = 0.6)), pch=c(3,0, rep(15,length(EA))), cex=0.8, bty='n')
  #   }else{
  #     legend('topleft', legend=c('EC','TA'), col = 1, pch=c(3,0), cex=0.8, bty='n')
  #   }
  #   
  #   # LCT
  #   if(LCT.present){
  #     ii <- 0
  #     for(ii in 1:length(LCT.list)){
  #       plot(LCT.list.utm[[ii]], add=T, border = LCT.color[ii], col=adjustcolor(LCT.color[ii], alpha.f = 0.6))
  #       # text(x=coordinates(LCT.list.utm[[ii]])[1], y=coordinates(LCT.list.utm[[ii]])[2], LCT.list[[ii]]@data$Description, pos=1, cex=.5)
  #     }
  #     plot(TA.utm, lwd=2, add=T)
  #     legend('topright', legend=paste('LCT_0', seq_along(LCT.list), sep=''), col = adjustcolor(LCT.color, alpha.f = 0.6), cex=0.8, pch=15, bty='n')
  #   }
  #   
  #   dev.off()
  #   
  # }
  
  # Site spatial features visualization (plot adattato per sf)
  if (plot.spatial.feats) {
    png(paste(FFP.dir, site.ID, '_','SpatialFeatures_', format(Sys.Date(), '%Y%m%d'), '.png', sep=''),
        width = 14, height = 14, unit = 'cm', res=300)
    # windows()
    plot(sf::st_geometry(TA.utm), reset = FALSE)
    plot(sf::st_geometry(EC.tower.utm), add = TRUE, pch = 3)
    
    if (LCT.present){
      for(ii in seq_along(LCT.list.utm)){
        plot(sf::st_geometry(LCT.list.utm[[ii]]), add = TRUE, col = adjustcolor(LCT.color[ii], 0.5))
      }
      plot(sf::st_geometry(TA.utm), add = TRUE, lwd = 2)
    }
    
    if (!is.na(EA.num)){
      plot(sf::st_geometry(EA.utm), add = TRUE, col = adjustcolor(2, alpha.f = 0.5))
    }
    
    dev.off()
  }
  
  # # remove spatial data files once imported
  # unlink(SP.data.dir, recursive = TRUE)
  
  
  # PBL height estimation =====================================================
  
  # Compute Planetary boundary layer heights 
  PBL.h.0 <- calc_PBL_height(TIMESTAMP = set.cleaned.timestamp,   # timestamp [%Y-%m-%d %H:%M:%S], either POSIXlt or character
                             LAT = st.lat,                        # point latitude [decimal degrees]
                             zm = H.EC,                           # measurement height [m]
                             d = d.,                              # displacement height [m]
                             MOL = data.table::nafill(set.cleaned$'MO_LENGTH', "locf"),               # Monin-Obukov length [m]
                             ustar = set.cleaned$'USTAR',         # friction velocity [ms-2]
                             H = set.cleaned$'H',                 # sensible heat flux [Wm-2]
                             AT = set.cleaned$'AT',                # air temperature [K]
                             AP = set.cleaned$'AP',                # air pressure [Pa]
                             rho = set.cleaned$'rho',              # density of dry air [kg m-3]
                             cp = set.cleaned$'cp',                # air heat capacity [J kg-1 K-1]
                             altit = st.alt,                       # site altitude [m]
                             plot.atmc = FALSE,                     # set TRUE to plot atmospheric conditions frequency
                             plot.PBL = TRUE,                      # set TRUE to plot PBL height estimations
                             save.PBL = TRUE,                      # set TRUE to save result (tabbed ASCII)
                             PBL.dir = FFP.dir)                      # directory of PBL estimations output
  
  
  # Do some adjustement of Planetary boundary layer heights 
  PBL.h <- proc_PBL_height(TIMESTAMP = set.cleaned.timestamp,               # timestamp [%Y-%m-%d %H:%M:%S]
                           PBL.h = PBL.h.0,                   # PBL heights 
                           plot.PBL.h.freq = FALSE,  # plot the PBL heights frequency histogram
                           PBL.dir = FFP.dir           # directory of PBL estimations
  )
  
  
  # FFP invalid condition data filtering mask ------------------------------------
  
  # The footprint parameterization is valid only in certain ranges of micrometeorological conditions.
  # In particular, the model is claimed to be valid if the following conditions hold:
  # - The measurement height is lower than the boundary layer height;
  # - The terrain is dynamically homogeneous;
  # - The stability parameter is in the range of: -200<Z<1 ;
  # - The friction velocity is larger than a specific threshold: u. ≥ 0.2 m s-1;
  # - The measurement height is larger than 1 m: hm≥ 1m
  # Additional limiting conditions are taken according to Kljun et al., 2015, Sect 2 and Eq. 27
  
  # TRUE: remove data point, FALSE: retain data point
  
  ## Model simulations ranges:
  
  # . Receptor heights at zm/h:  0.0005 - 0.5 (Tab. 1)
  if (min(zm./PBL.h, na.rm = TRUE) < 0.0005 | max(zm./PBL.h, na.rm = TRUE) > 0.5) {
    cat(warn('\nWARNING: Some PBL normalized Receptor heights (zm/h) values are out of FFP model simulations ranges'))
  }
  
  # . Roughness length: 0.01 - 3.0 (Tab. 1)
  if (min(z0., na.rm = TRUE) < 0.01 | max(z0., na.rm = TRUE) > 3) {
    cat(warn('\nWARNING: Some roughness lengths are out of FFP model simulations ranges'))
  }
  
  # . Roughness Layer height (Kljun 2015: z* = 2.75 hrs, with hrs is approximated by 10 * z0)
  zstar <- 2.75 * 10 * z0.
  
  # . Entertainment Layer height [Holtslag and Nieuwstadt, 1986, in Kljun et al., 2015]
  EL.h <- PBL.h * 0.8
  if (any(isTRUE(table(zm. > EL.h)))) {
    cat(warn('\nWARNING: there are some half-hours in which the dynamic measurement height is haigher that the concurrent entertainment layer height'))
  }
  
  
  # FFP MODEL invalid conditions (thresholds from N. Kljun 2015)
  # TRUE:  remove data point, FALSE: retain data point
  
  # Iterative process to scale for the "zm. < 5 * z0" requirement:
  # As it has been noted that 5 often imply a large data rejection, the multiplier is progressively reduced down 
  # until being acceptable, i.e. less than 50% of FALSE (good data) records, or 2 at least.
  z0.multi <- 5
  i.z0 <- 1
  req.met <- FALSE
  for (i.z0 in 1:30) {
    # frequency table of the current setting
    req.tab <- table(zm. < z0.multi * z0.)
    if (length(req.tab) == 1 && grepl('FALSE', names(req.tab))) { # If all data are ok
      req.met <- TRUE # exit
    } else {
      if (length(req.tab) == 1 && grepl('TRUE', names(req.tab))) { # If all data ore bad
        z0.multi <- z0.multi - 0.1 # decrease the multiplier and continue
      } else {
        if (req.tab['TRUE'] <= 0.5 * req.tab['FALSE']) { # If the percentage is met
          req.met <- TRUE # exit
        } else {
          z0.multi <- z0.multi - 0.1 # otherwise decrease the multiplier and continue
        }
      }
    }
    if (req.met) {break}
  }
  
  # If z0.multi is not 5, save a log file with the exception
  if (z0.multi < 5) {data.table::fwrite(list(paste('The footprint model requirement of measurement height vs roughness sublayer limit was set to',
                                                   'zm. <', z0.multi, '* z0, instead of zm. < 5 * z0.')),
                                        paste0(FFP.dir, 'FFP_z0-exceptions_', format(Sys.Date(), '%Y%m%d'), '.txt'),
                                        quote = F, row.names = F)}
  
  # Set invalid records positions for FFP model function
  # TRUE if violates model assumptions
  FFP.invalid.cond.flag <- Reduce('|',
                                  list(
                                    ## . measurement height vs roughness sublayer limit [Eq. 27]
                                    # RL height: z* = n hrs , where: 2 < n < 5 (Kljun used 2.75), hrs = mean roughness element height = 10 z0
                                    # Requirement: minimum RL height: zm > 20 z0 [Eq. 27] # NOTE, RELAXED CONDITION: should be 20, but is too much restrictive
                                    zm. < z0.multi * z0.,  # table(zm. < 1 * z0.); summary(5 * z0.); summary(zm.)
                                    ## . measurement height vs entertainment layer height [Eq. 27]
                                    zm. > EL.h,  # table(zm. > EL.h)
                                    ## . measurement height vs planetary boundary layer
                                    H.EC > PBL.h, # table(PBL.h < H.EC)
                                    ## . stability [Eq. 27]
                                    set.cleaned$'ZL' <= -15.5 | set.cleaned$'ZL' > 1, # table(set.cleaned$'ZL' <= -15.5); table(set.cleaned$'ZL' > 1)
                                    ## . friction velocity
                                    set.cleaned$'USTAR' < 0.2 # table(set.cleaned$USTAR < 0.18); table(set.cleaned$USTAR[!NEE.QC.flag] < 0.2)
                                  )
  )
  # Simplify the object [TRUE == drop]
  FFP.invalid.cond.flag <- c(FFP.invalid.cond.flag)
  
  ## Data effectively removed because of invalid condition:
  # table(FFP.invalid.cond.flag)
  # table(FFP.invalid.cond.flag[!NEE.QC.flag])
  # table(FFP.invalid.cond.flag)/nrow(set.cleaned)
  
  ## Summary statistics of invalid condition:
  FFP.invcond.df <- merge(merge(merge(table(PBL.h > H.EC),                                    
                          as.data.frame(table(set.cleaned$'USTAR' < 0.2)), by='Var1', all=T),                       
                          merge(table(set.cleaned$'ZL' < -15.5 | set.cleaned$'ZL' > 1), 
                          table(zm. < z0.multi * z0.), by='Var1', all=T), by='Var1', all=T),                            
                          table (zm. > EL.h), by='Var1', all=T)
  colnames(FFP.invcond.df) <- c('FLAG','PBL.h>H.EC', 'USTAR<0.2', '-15.5<ZL<1', 'zm<*z0', 'zm>EL.h')
  data.table::fwrite(as.data.frame(FFP.invcond.df),
                     paste0(FFP.dir, 'FFP_conditions_stats_', format(Sys.Date(), '%Y%m%d'), '.txt'),
                     quote = F, row.names = T, sep='\t')
  
  # populate the FFP stats vector
  
  # data flagged for invalid FFP conditions
  FFP.data.summary$'N_INVCOND_FFP_flag' <- unname(table(FFP.invalid.cond.flag)['TRUE'])
  
  # data flagged for invalid FFP conditions & NOT already excluded by the NEE QC (or possibly USTAR)
  if (DO.ustar.filter) {
    
    FFP.data.summary$N_INVCOND_FFP_flag_over_QCUstf <- ifelse(is.na(unname(table(FFP.invalid.cond.flag[!NEE.QC.flag & !USTAR.QC.flag])['TRUE'])), 0, unname(table(FFP.invalid.cond.flag[!NEE.QC.flag & !USTAR.QC.flag])['TRUE'])) # data flagged for invalid FFP conditions & NOT excluded by the NEE and USTAR QC
    
    FFP.data.summary$'N_data_FFP' <- FFP.data.summary$N_raw_data -
      FFP.data.summary$N_QC_flag -
      FFP.data.summary$N_USTAR_QC_flag_over_QCf -
      FFP.data.summary$N_INVCOND_FFP_flag_over_QCUstf
    
  } else {
    
    FFP.data.summary$'N_INVCOND_FFP_flag_over_QCf' <- ifelse(is.na(unname(table(FFP.invalid.cond.flag[!NEE.QC.flag])['TRUE'])),
                                                             0,
                                                             unname(table(FFP.invalid.cond.flag[!NEE.QC.flag])['TRUE']))
    
    FFP.data.summary$'N_data_FFP' <- FFP.data.summary$N_raw_data -
      FFP.data.summary$N_QC_flag -
      FFP.data.summary$N_INVCOND_FFP_flag_over_QCf
  }
  
  
  ## percentage of data actually used for FFP
  FFP.data.summary$'N_data_FFP_perc' <- round(FFP.data.summary$'N_data_FFP' / FFP.data.summary$'N_raw_data' * 100, 1)
  
  # if percentage of original data, and if less than 50%, STOP
  if (FFP.data.summary$N_data_FFP_perc < 50) {
    if (FFP.data.summary$N_data_FFP_perc >= 30) {
      cat(warn('\nWARNING: less than 50% of data, but more than 30%, survived for the footprint estimation.\n'))
    } else {
      cat(error('\nWARNING: there are a few data (less than 30%) to get a robust footprint estimation.\n'))
      # stop()
    }
  } else {
    cat(note('\nOK: there are enough data to get a robust footprint estimation.'))
  }
  
  
  ### overall filtering mask ###
  
  # set the proper filtering policy
  if (DO.ustar.filter & DO.wind.sect.excl) {
    qcf <- Reduce("|", list(NEE.QC.flag, USTAR.QC.flag, WD.excl.QC.flag, FFP.invalid.cond.flag))
  } else if (DO.ustar.filter & !DO.wind.sect.excl) {
    qcf <- Reduce("|", list(NEE.QC.flag, USTAR.QC.flag, FFP.invalid.cond.flag))
  } else {
    qcf <- Reduce("|", list(NEE.QC.flag, FFP.invalid.cond.flag))
  }
  
  # QCF check: FALSE length must be == to FFP.data.summary$N_data_FFP
  if (!all(qcf)) {
    if (unname(table(qcf)['FALSE']) != FFP.data.summary$'N_data_FFP') {
      cat(warn('\nWARNING: QCF good flags are not the same as the number of FFP valid record reported in FFP.data.summary$N_data_FFP.'))
    } else {
      cat(note('\nOK: QCF good flags are the same as the number of FFP valid record reported in FFP.data.summary$N_data_FFP'))
    }
  } else {
    cat(error('\nError: no data left for FFP analysis after the FFP quality check. The analysis cannot proceed.'))
    stop()
  }
  
  # possibly set NAs to TRUE
  if (any(is.na(qcf))) {
    qcf[is.na(qcf)] <- TRUE # Add input missing data to invalid conditions: SEE IF IT WORTH TO SPLIT THE TWO
  }
  
  
  # FFP input -----------------------------------------------------------------
  
  # If daytime is not reported in the dat, calculate it
  if (!is.null(set.cleaned$'daytime')) {
    daytime <- set.cleaned$'daytime'
  } else {
    # TODO: currently based on time. Use the potential radiation instead (see e.g. https://www.licor.com/env/support/EddyPro/faqs.html)
    daytime <- c(as.numeric(format(set.cleaned.timestamp, '%H')) > 6 & as.numeric(format(set.cleaned.timestamp, '%H')) < 18) * 1
  }
  
  FFP.input <- data.frame(NULL)
  FFP.input <- data.frame(
    TIMESTAMP = set.cleaned.timestamp,
    daytime = daytime,
    # custom flag for filtering out half-hourly records based on NEE and possibly ustar data
    qcf = qcf,
    # site/setup characteristics
    hm = H.EC,                    # Height of the eddy covariance system [m], may change during the period
    hc = H.canopy,                # Average vegetation height, may change during the period
    d = d.,                       # Displacement height
    z0 = z0.,                     # Roughness length
    zm = zm.,                     # Dynamic measurement height (i.e. z‐d) [m]
    # turbulence
    umean = set.cleaned$'WS',     # Mean wind speed at zm [ms‐1] ‐ enter [NaN] if not known (Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint)
    ol = set.cleaned$'MO_LENGTH', # Vector of Obukhov length [m]
    sigmav = set.cleaned$'V_SIGMA', # Vector of standard deviation of lateral velocity fluctuations [ms‐1]
    ustar = set.cleaned$'USTAR',  # Vector of friction velocity [ms‐1]
    wind_dir = set.cleaned$'WD',  # Vector of wind direction in degrees (of 360) for rotation of the footprint
    PBL = coredata(PBL.h)         # Vector of boundary layer height [m]
  )
  
  # Save FFP input data file
  data.table::fwrite(FFP.input, 
                     paste0(FFP.dir, site.ID, '_FFPinput_', format(Sys.Date(), '%Y%m%d'), '.csv'), 
                     sep=',', quote = FALSE, row.names = FALSE, na = '-9999')
  
  
  #****************************************************************************
  # FFP - SP INTERACTION ======================================================
  # Parallel computing of FFPs and their interactions with spatial data 
  
  cat(prog.mes('\n::Running the FFP computation and spatial analysis in parallel::\n'))
  
  # Split the whole period in years 
  # range(FFP.input$'TIMESTAMP')
  # range(FFP.input$'TIMESTAMP' - 1800)
  FFP.input.split <- list(NULL)
  FFP.input.split <- split(FFP.input, 
                           format(FFP.input$'TIMESTAMP' - 1800, "%Y")) # shift 30 minutes back to avoid passing the year with the last half-hour
  # length(FFP.input.split)
  # lapply(FFP.input.split, dim)
  
  # Split the cleanset too, just to match the NEE flag (see the output section)
  set.cleaned.split <- list(NULL)
  set.cleaned.split <- split(set.cleaned, 
                             format(set.cleaned.timestamp - 1800, "%Y")) 
  # lapply(set.cleaned.split, dim)
  
  # Split the invalid condition flag list too, just to match the NEE flag (see the output section)
  FFP.invalid.cond.flag.split <- list(NULL)
  FFP.invalid.cond.flag.split <- split(FFP.invalid.cond.flag, 
                                       format(set.cleaned.timestamp - 1800, "%Y")) 
  # lapply(FFP.invalid.cond.flag.split, dim)
  
  
  i.ffp <- 1
  for (i.ffp in 1:length(FFP.input.split)) { # loop over years
    
    cat(prog.mes('\nComputing FFP for the year ' %+% names(FFP.input.split)[i.ffp] %+% ' ...\n'))
    
    # Register the doParallel parallel back-end specifying a cluster with the total number of cores
    cl <- switch((.Platform$OS.type == 'windows')+1, 
                 parallel::makeCluster(6),                         # @NARTEN (or unix systems) >> use 6 nodes and launch max 4 job in parallel 
                 parallel::makeCluster(parallel::detectCores()-2)) # @windowsPC >> keep 2 nodes free not to overload the computer
    
    # doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl)
    
    iterations <- nrow(FFP.input.split[[i.ffp]])
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    pbopts <- list(progress = progress)
    
    FFP.res <- data.frame(NULL)
    FFP.res <- as.data.frame(
      
      foreach::foreach(i = 1:nrow(FFP.input.split[[i.ffp]]), 
                       .combine = dplyr::bind_rows,
                       # .packages = c('birk','flux','sp','rgeos','raster'),
                       .packages = c('birk','flux', 'sf', 'terra', 'ggplot2', 'viridisLite'),
                       .options.snow = pbopts,
                       # .verbose = TRUE, # for debugging
                       .errorhandling = 'stop'  # execution will be stopped via the stop function if an error occurs
                       # .errorhandling = 'pass'  #  the error object generated by task evaluation will be included with the rest of the results
                       # .errorhandling = 'remove'  #  the result for that task will not be returned, or passed to the .combine function
                       
      ) %dopar% {
        # ) %do% {
        
        # source(paste0(FFP.func.dir, 'FFP_calc_parallel.R'))       # Parallel execution of FFP calculation and analysis of the interactions with spatial data
        source(paste0(FFP.func.dir, 'FFP_calc_parallel_dev.R'))     # Parallel execution of FFP calculation and analysis of the interactions with spatial data
        # i <- 5790
        tmp <- FFP.input.split[[i.ffp]][i, ]
        # print(tmp)
        FFPCalcParallel(
          site.ID = site.ID,
          input.data = tmp,                # a vector FFP.input[1,]
          FFPSP.output.dir = FFP.dir,      # Path to the output directory
          FFP.func.dir = FFP.func.dir,     # Path to the directory containing the FFP function by Kljun et al. 2015
          FFP.R = c(50,70,80),             # Percentage of source area, i.e. a value between 10% and 90%. * 2
          plot.FFP.cwi = FALSE,            # Plot the crosswind‐integrated footprint visualization
          plot.ellipse.fit = FALSE,        # plot the ellipses fit to the 70% contribution
          plot.FFP.clim = FALSE,           # Plot the footprint climatology
          EC.tower.utm = EC.tower.utm,     # EC tower coordinates as SpatialPoints object, in UTM
          plot.spat.FFP = FALSE,           # plot the 70% and 80% FFP contributions over site's spatial features
          TA.utm = TA.utm,                 # Target Area polygon as SpatialPolygonsDataFrame object, in UTM
          plot.TA.FFP = FALSE,             # plot the TA - FFP intersection
          DO.EA.analysis = DO.EA.analysis, # whether run the analysis in EA or not
          EA.utm = switch(is.na(EA.num)+1 , EA.utm, NULL),                   # Exclusion Area(s) polygon as a single SpatialPolygonsDataFrame object, in UTM
          plot.EA.FFP = FALSE,             # plot the EA - FFP intersection
          DO.LCT.analysis = DO.LCT.analysis,  # whether run the analysis on LCTs or not
          LCT.utm = switch(LCT.present+1 , NULL, LCT.list.utm),  # Land Cover Typologies(s) polygon as a list of SpatialPolygonsDataFrame objects, in UTM
          plot.LCT.FFP = FALSE             # plot the LCT FFP intersection
        )
        
      }
      
    )
    
    close(pb)
    parallel::stopCluster(cl)
    # t1 <- Sys.time()
    # t1-t0
    cat(prog.mes('done.'))
    
    # class(FFP.res)
    # drop row names from the results dataframe
    rownames(FFP.res) <- NULL
    
    # character to numeric
    numcolumns <- names(FFP.res)[-c(1)]
    FFP.res[, numcolumns] <- lapply(numcolumns, function(x) as.numeric(FFP.res[[x]]))
    
    
    
    # OUTPUT ================================================================
    
    ## . Build footprint FLAG ----------------------------------------------------------
    # combine and format footprint flags 
    # from 2D model error, invalid conditions, flagged by QC NEE
    # 0 = data used for footprint 
    # 1 = 2D footprint model flag error 
    # 2 = data removed because of invalid conditions for footprint model 
    # 3 = NEE removed by RFlux::QC
    
    # summary(FFP.flag.err) # 2D model error: 1 in case of error, 2 if not all contour plots (r%) within specified domain, 0 otherwise
    # summary(FFP.invalid.cond.flag) # model invalid conditions
    # summary(NEE.QC.flag) # NEE QC (ATTUALMENTE NON CONSIDERA USTAR, WIND EXCLUSION E HWHT FLAG, ma i primi due non sono fatti generalmente e la terza ci sono pochissimi punti)
    
    # base vector (use the standard format for originally missed values)
    FOOTPRINT_FLAG <- rep(NA, nrow(FFP.input.split[[i.ffp]]))
    
    # 3 = flagged by QC (0 = ok, 1 = sev.err, 2 = mod.err+out, -9999 originally missed)
    # FOOTPRINT_FLAG[NEE.QC.flag == TRUE] <- 3
    FOOTPRINT_FLAG <- set.cleaned.split[[i.ffp]]$'NEE_DATA_FLAG' 
    FOOTPRINT_FLAG[set.cleaned.split[[i.ffp]]$'NEE_DATA_FLAG' != 0] <- 3
    
    # 0 = data used for footprint & 1 = 2D footprint model flag error
    FOOTPRINT_FLAG[!is.na(FFP.res$'FFP_FLAG')] <- replace(FFP.res$'FFP_FLAG', FFP.res$'FFP_FLAG' == 2, 1)[!is.na(FFP.res$'FFP_FLAG')]
    
    # 2 = data removed because of invalid conditions for footprint model
    FOOTPRINT_FLAG[which(FFP.invalid.cond.flag.split[[i.ffp]] == TRUE)] <- 2
    
    # table(FOOTPRINT_FLAG)
    
    
    ## . Summary table ----------------------------------------------------------
    # Make result the table (partial, in case of more than one year)
    FFP.res.part <- data.frame(NULL)
    
    if (proc.mode == 'TEST') {
      
      # 1D FFP and 2D FFP/TA intersection
      FFP.res.part <- data.frame(
        TIMESTAMP = format(FFP.input.split[[i.ffp]]$'TIMESTAMP', '%Y%m%d%H%M'),                      # TO BE REMOVED FOR L2 PRODUCT
        daytime = FFP.input.split[[i.ffp]]$'daytime',                          # Daytime condition, 0 = night 1 = day
        # QC_REMOVED = NEE.QC.flag,                             # data removed by QC (removed = 1)
        # FETCH_FLAG = FFP.input$qcf*1,                         # Footprint flag 0 / 1: 0 = data used for footprint, 1 = data removed by specific FFP filtering (e.g. invalid conditions for model)
        PBLH = round(FFP.input.split[[i.ffp]]$'PBL', 2),                                 	# Planetary boundary layer height
        
        # FFP 1D OUTPUT
        FETCH_MAX = round(FFP.res$'FETCH_MAX', 2),                               # Distance at which footprint contribution is maximum	[m] * 
        # the following variables are currently NOT reported
        FETCH_OFFSET = round(FFP.res$'FETCH_OFFSET', 2),                        # Distance at which footprint cumulative probability is 1%	[m] CHIEDERE CONFERMA A GER
        FETCH_10 = round(FFP.res$'FETCH_10', 2),                        # Distance at which footprint cumulative probability is 10%	[m]
        FETCH_20 = round(FFP.res$'FETCH_20', 2),                        # Distance at which footprint cumulative probability is 20%	[m]
        FETCH_30 = round(FFP.res$'FETCH_30', 2),                        # Distance at which footprint cumulative probability is 30%	[m]
        FETCH_40 = round(FFP.res$'FETCH_40', 2),                        # Distance at which footprint cumulative probability is 40%	[m]
        FETCH_50 = round(FFP.res$'FETCH_50', 2),                        # Distance at which footprint cumulative probability is 50%	[m]
        FETCH_60 = round(FFP.res$'FETCH_60', 2),                        # Distance at which footprint cumulative probability is 60%	[m]
        FETCH_70 = round(FFP.res$'FETCH_70', 2),                          # Distance at which footprint cumulative probability is 70%	[m] *
        FETCH_80 = round(FFP.res$'FETCH_80', 2),                          # Distance at which footprint cumulative probability is 80%	[m] *
        FETCH_90 = round(FFP.res$'FETCH_90', 2),                          # Distance at which footprint cumulative probability is 90%	[m] *
        
        # FFP 2D - TA OUTPUT
        FFP_80_SURF = round(FFP.res$'FOOTPRINT_80_SURF', 2),                # 2D FFP 80% isoplethe area [m2]
        # FFP - TA INTERSECTION
        FFP_TA_SURF = round(FFP.res$'FOOTPRINT_TA_SURF', 2),          # Area of FFP/TA intersection [m2]
        FFP_TA_PSURF = round(FFP.res$'FOOTPRINT_TA_PSURF', 2),        # Area of FFP/TA intersection as percentage of the whol FFP area [%]
        FFP_TA_CONTR = round(FFP.res$'FOOTPRINT_TA_CONTR', 2),        # Cumulative FFP contribution of the FFP/TA intersection [m2 m-2, %]
        # FOOTPRINT_FLAG	= FFP.res$'FFP_FLAG',                         # Footprint flag (0, 1, 2, 3, -9999). 0 = data used for footprint, 1 = 2D footprint model flag error,  2 = not all 2D footprint model isoplethes within specified domain, 3 = data removed because of invalid conditions for footprint model, -9999 missing because of QC
        FOOTPRINT_FLAG = FOOTPRINT_FLAG,                         # Footprint flag (0, 1, 2, 3, -9999). 0 = data used for footprint, 1 = 2D footprint model flag error,  2 = not all 2D footprint model isoplethes within specified domain, 3 = data removed because of invalid conditions for footprint model, -9999 missing because of QC
        # FFP - EA INTERSECTION
        FFP_EA_SURF = round(FFP.res$'FOOTPRINT_EA_SURF', 2),          # Area of FFP/EA intersection [m2]
        FFP_EA_PSURF = round(FFP.res$'FOOTPRINT_EA_PSURF', 2),        # Area of FFP/EA intersection as percentage of the whol FFP area [%]
        FFP_EA_CONTR = round(FFP.res$'FOOTPRINT_EA_CONTR', 2)         # Cumulative FFP contribution of the FFP/EA intersection [m2 m-2, %]
      )
      
      # 2D FFP/LCT intersection (if meant)
      if (LCT.present) {
        # FFP - LCT INTERSECTION (for each LCT)
        
        # Area of FFP/LCT intersections [m2]
        # lct.surf <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.surf), 2, function(x) round(x,2)))     
        lct.surf <- round(FFP.res[, grepl('LCT.*_SURF', names(FFP.res))], 2)
        names(lct.surf) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_SURF', sep='')
        FFP.res.part <- cbind(FFP.res.part, lct.surf)
        
        # Area of FFP/LCT intersection as percentage of the whole FFP area [%]
        # lct.surf.p <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.surf.perc), 2, function(x) round(x,2)))  
        lct.surf.p <- round(FFP.res[, grepl('LCT.*PSURF', names(FFP.res))], 2)
        names(lct.surf.p) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_PSURF', sep='')
        FFP.res.part <- cbind(FFP.res.part, lct.surf.p)
        
        # Cumulative FFP contribution of the FFP/LCT intersection [m2 m-2, %]
        # lct.contr <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.contr), 2, function(x) round(x,2)))   
        lct.contr <- lct.surf.p <- round(FFP.res[, grepl('LCT.*CONTR', names(FFP.res))], 2)
        names(lct.contr) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_CONTR', sep='')
        FFP.res.part <- cbind(FFP.res.part, lct.contr)
        
      }
    } # proc.mode == 'TEST' condition END
    
    
    if (grepl('L2|NRT', proc.mode)) { # output for L2(I) PRODUCT
      
      FFP.res.part <- data.frame()
      
      # 1D FFP and 2D FFP/TA intersection
      FFP.res.part <- data.frame(
        TIMESTAMP_END = format(FFP.input.split[[i.ffp]]$'TIMESTAMP', '%Y%m%d%H%M'),                      # TO BE REMOVED FOR L2 PRODUCT
        daytime = FFP.input.split[[i.ffp]]$'daytime',                          # Daytime condition, 0 = night 1 = day
        # QC_REMOVED = NEE.QC.flag,                             # data removed by QC (removed = 1)
        # FETCH_FLAG = FFP.input$qcf*1,                         # Footprint flag 0 / 1: 0 = data used for footprint, 1 = data removed by specific FFP filtering (e.g. invalid conditions for model)
        PBLH = round(c(FFP.input.split[[i.ffp]]$'PBL'), 2),                                 	# Planetary boundary layer height
        
        # FFP 1D OUTPUT
        FETCH_MAX = round(FFP.res$'FETCH_MAX', 2),                               # Distance at which footprint contribution is maximum	[m] * 
        # the following variables are currently NOT reported
        FETCH_OFFSET = round(FFP.res$'FETCH_OFFSET', 2),                        # Distance at which footprint cumulative probability is 1%	[m] CHIEDERE CONFERMA A GER
        # FETCH_10 = round(FFP.res$'FETCH_10', 2),                        # Distance at which footprint cumulative probability is 10%	[m]
        # FETCH_20 = round(FFP.res$'FETCH_20', 2),                        # Distance at which footprint cumulative probability is 20%	[m]
        # FETCH_30 = round(FFP.res$'FETCH_30', 2),                        # Distance at which footprint cumulative probability is 30%	[m]
        # FETCH_40 = round(FFP.res$'FETCH_40', 2),                        # Distance at which footprint cumulative probability is 40%	[m]
        FETCH_50 = round(FFP.res$'FETCH_50', 2),                        # Distance at which footprint cumulative probability is 50%	[m]
        # FETCH_60 = round(FFP.res$'FETCH_60', 2),                        # Distance at which footprint cumulative probability is 60%	[m]
        FETCH_70 = round(FFP.res$'FETCH_70', 2),                          # Distance at which footprint cumulative probability is 70%	[m] *
        FETCH_80 = round(FFP.res$'FETCH_80', 2),                          # Distance at which footprint cumulative probability is 80%	[m] *
        FETCH_90 = round(FFP.res$'FETCH_90', 2),                          # Distance at which footprint cumulative probability is 90%	[m] *
        
        # FFP 2D - TA OUTPUT
        FFP_80_SURF = round(FFP.res$'FOOTPRINT_80_SURF', 2),                # 2D FFP 80% isoplethe area [m2]
        # FFP - TA INTERSECTION
        FFP_TA_SURF = round(FFP.res$'FOOTPRINT_TA_SURF', 2),          # Area of FFP/TA intersection [m2]
        # FFP_TA_PSURF = round(FFP.res$'FOOTPRINT_TA_PSURF', 2),        # Area of FFP/TA intersection as percentage of the whol FFP area [%]
        FFP_TA_CONTR = round(FFP.res$'FOOTPRINT_TA_CONTR', 2),        # Cumulative FFP contribution of the FFP/TA intersection [m2 m-2, %]
        FOOTPRINT_FLAG = FOOTPRINT_FLAG                         # Footprint flag (0, 1, 2, 3, -9999). 0 = data used for footprint, 1 = 2D footprint model flag error,  2 = not all 2D footprint model isoplethes within specified domain, 3 = data removed because of invalid conditions for footprint model, -9999 missing because of QC
        # # FFP - EA INTERSECTION
        # FFP_EA_SURF = round(FFP.res$'FOOTPRINT_EA_SURF', 2),          # Area of FFP/EA intersection [m2]
        # FFP_EA_PSURF = round(FFP.res$'FOOTPRINT_EA_PSURF', 2),        # Area of FFP/EA intersection as percentage of the whol FFP area [%]
        # FFP_EA_CONTR = round(FFP.res$'FOOTPRINT_EA_CONTR', 2)         # Cumulative FFP contribution of the FFP/EA intersection [m2 m-2, %]
      )
      
      # 2D FFP/LCT intersection (if meant)
      # CURRENTLY NOT REPORTED FOR L2 PRODUCT (because LCT are not official features)
      # if (LCT.present) {
      #   # FFP - LCT INTERSECTION (for each LCT)
      #   
      #   # Area of FFP/LCT intersections [m2]
      #   # lct.surf <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.surf), 2, function(x) round(x,2)))     
      #   lct.surf <- round(FFP.res[, grepl('LCT.*_SURF', names(FFP.res))], 2)
      #   names(lct.surf) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_SURF', sep='')
      #   FFP.res.part <- cbind(FFP.res.part, lct.surf)
      #   
      #   # Area of FFP/LCT intersection as percentage of the whole FFP area [%]
      #   # lct.surf.p <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.surf.perc), 2, function(x) round(x,2)))  
      #   lct.surf.p <- round(FFP.res[, grepl('LCT.*PSURF', names(FFP.res))], 2)
      #   names(lct.surf.p) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_PSURF', sep='')
      #   FFP.res.part <- cbind(FFP.res.part, lct.surf.p)
      #   
      #   # Cumulative FFP contribution of the FFP/LCT intersection [m2 m-2, %]
      #   # lct.contr <- as.data.frame(apply(do.call('rbind', FFP.LCT.intersect.contr), 2, function(x) round(x,2)))   
      #   lct.contr <- lct.surf.p <- round(FFP.res[, grepl('LCT.*CONTR', names(FFP.res))], 2)
      #   names(lct.contr) <- paste('FFP_LCT', sprintf("%02.f", 1:LCT.num), '_CONTR', sep='')
      #   FFP.res.part <- cbind(FFP.res.part, lct.contr)
      #   
      
    }  # proc.mode == 'L2||NRT' condition END
    
    # Save yearly results
    write.table(FFP.res.part, 
                paste(FFP.dir, site.ID, '_FFP_partres_', names(FFP.input.split)[[i.ffp]], '.csv', sep=''), 
                sep=',', quote = F, row.names = F)
    
    
  } # FFP computation loop end
  
  
  # Merge yearly datasets and SAVE definitive RESULTS
  # this will be used in ICOS
  # FFP.analysis.results <- 
  FFP.y.files <- grep('partres', list.files(FFP.dir, full.names = TRUE), value = TRUE)
  FFP.analysis.results <- data.frame(NULL)
  FFP.analysis.results <- as.data.frame(
    data.table::rbindlist(lapply(FFP.y.files, function(x) fread(x, na.strings = '-9999', header = TRUE, data.table = FALSE)),
                          use.names = TRUE, fill = TRUE))
  
  # check FFP flag and FFP data summary consistency

  if (all(FFP.analysis.results$'FOOTPRINT_FLAG' == '3')) { # All NEE records discarded by RFlux 
    cat(warn('\nFFP_anaysis::WARNING - NO records have been processed by FFP because all NEE data were flagged by RFlux::QC.'))
    cat(note('\nAll footptint variables will be set to NA.\n'))
  } 
  
  if (!any(FFP.analysis.results$'FOOTPRINT_FLAG' == '0')) { # FFP not computed because of any reason (FOOTPRINT_FLAG == 1|2|3)
    cat(warn('\nFFP_anaysis::WARNING - NO records have been processed by FFP because all data were bad-flagged (either by RFlux::QC or the FFP flagging policy).'))
    cat(note('\nAll footptint variables will be set to NA.\n'))
  } 
  
  # if (FFP.data.summary$'N_data_FFP' == table(FFP.analysis.results$'FOOTPRINT_FLAG')['0']) {
  if (FFP.data.summary$'N_data_FFP' == sum(table(FFP.analysis.results$'FOOTPRINT_FLAG')[c('0','1')], na.rm=T)) { # "0"=OK, "1"=FFP_error
    # cat(note('\nOK, the number of data according to the FFP flag and the one reported in the FFP data summary are consistent'))
    cat(note('\nThe number of processed FFP records (FFP_FLAG == 0|1) matches the available FFP record count reported in the data summary ("N_data_FFP").'))
  } else {
    # if (FFP.data.summary$'N_data_FFP' == table(FFP.analysis.results$'FOOTPRINT_FLAG')['0'] + table(FFP.analysis.results$'FOOTPRINT_FLAG')['1']) {
    # cat(note('\nOK, the number of data according to the FFP flag and the one reported in the FFP data summary are consistent'))
    # } else {
    cat(note('\nFFP_anaysis::WARNING - the number of processed FFP records (FFP_FLAG == 0|1) DOES NOT match the available FFP record count reported in the data summary ("N_data_FFP").'))
    # }
  }
  
  # Save FFP results
  write.table(FFP.analysis.results, 
              paste(FFP.dir, site.ID, '_FFP_results_', format(Sys.Date(), '%Y%m%d'), '.csv', sep=''), 
              sep=',', quote = F, row.names = F)
  
  
  
  # free some memeory
  rm(FFP.input.split, set.cleaned.split)
  
  
  # -------------------------------------------------------------------.
  # IN 'L2' MODE THE ANALYSIS TERMINATES HERE
  
  if (grepl('L2|NRT', proc.mode)) {
    cat(prog.mes(paste0('\nOK! The FFP output dataset for ', bold(proc.mode),' production has been created.\n\n')))
    # cat('\n')
    # stop('\rIn L2 mode the FFP processing terminates here.') # putting \r in the beginning of the message, overwrite Error: in the message
  }
  
  
  # -------------------------------------------------------------------.
  # IN 'TEST' MODE THE ANALYSIS CONTINUES
  
  if (grepl('TEST', proc.mode)) {
    
    # . SPATIAL ANALYSIS RESULTS ###################################################
    
    FFP.an.res.4spatial <- FFP.analysis.results
    
    # Filter for SnW footprint ***********
    if (remove.SnW) {
      FFP.an.res.4spatial[which(SnW.FFP == 1), ] <- NA
    }
    
    # Remove FFP NAs (FFP_FLAG != 0) form the FFP results dataframe ***********
    # before going through the spatial analysis 
    FFP.an.res.4spatial <- FFP.an.res.4spatial[-which(FFP.an.res.4spatial$FOOTPRINT_FLAG != 0), ]
    
    
    ## . . FFP data summary -----------------------------------------------------
    
    # Total number of data (half-hours) used for spatial analysis 
    FFP.data.summary$'N_data_SpatAnal' <- nrow(FFP.an.res.4spatial) # !!!!!! VERIFICARE, sembra troppo alto
    FFP.data.summary$'N_data_SpatAnal_perc' <- round(FFP.data.summary$'N_data_SpatAnal' / FFP.data.summary$'N_dataset' * 100, 0)
    # FFP.data.summary$'N_data_SpatAnal_daytime' <- length(which(FFP.an.res.4spatial$'daytime' == 1))
    # FFP.data.summary$'N_data_SpatAnal_nighttime' <- length(which(FFP.an.res.4spatial$'daytime' == 0))
    
    FFP.data.summary.out <- as.data.frame(do.call(rbind, FFP.data.summary))
    # save results
    write.table(FFP.data.summary.out, paste(FFP.dir, site.ID, '_FFP_data_summary_', format(Sys.Date(), '%Y%m%d'), '.csv', sep=''), sep='\t', quote = F, row.names = TRUE, col.names = FALSE)
    
    
    # ***************************************************************************
    ## . . TEST 2 - Percentage of data coming from the Target Area --------------
    # ***************************************************************************
    #
    # TEST
    # During Step1 and also in the spatial sampling design a Target Area is defined by the PI. This is the area that
    # we aim to monitor with the ICOS station (see definition in the Sampling Design Instruction). The test will
    # analyze using a footprint model (Klijun et al. 2015) the estimated contribution area for each half hour and
    # check how many records have a contribution coming mainly from the target area. The footprint model has
    # clearly uncertainty associated to the estimations and this has been considered in the definition of the target
    # value that are relatively low (see below). However the fact that the same model is applied to all the stations
    # and following the same scheme makes the inter-station comparison robust.
    
    # TARGET
    # It will be calculated the number of half hours that has a cumulated contribution of at least 70% from the
    # Target Area and the target values is that they are at least 70% of the “good data” (data that already passed the QAQC).
    #
    # ACTION if the target value is not reached
    # An in depth analysis of the conditions that lead to an high percentage of data outside the target area is
    # conducted and a discussion with the PI is started in order to find solutions and alternatives, in particular
    # changing the measurement height or wind sectors to exclude.
    #
    # ******************************************************************************
    
    # Split the whole dataset into N subperiods
    # current splitting strategy: divide the whole period into 3 subperiods of approximately the same length
    Nsp <- 4 # number of SubPeriods
    subperiods <- cut(1:nrow(FFP.an.res.4spatial), Nsp, labels=paste0('SP#', 1:Nsp))
    FFP.res.list <- split.data.frame(FFP.an.res.4spatial, subperiods) # periods
    
    # SET thresholds (according to the MSA document)
    TA.contr.target <- 0.7   # minimum cumulated contribution of 70% from the Target Area
    TA.record.target <- 0.7  # minimum 70% of the “good data” must have a cumulated contribution of at least 70% from the Target Area
    
    
    # set structures
    FFP.TA.N <- c(NULL)
    FFPinTA.N <- c(NULL)
    FFPinTA.perc <- c(NULL)
    FFP.TA.N.DT <- c(NULL)
    FFPinTA.N.DT <- c(NULL)
    FFPinTA.perc.DT <- c(NULL)
    FFP.TA.N.NT <- c(NULL)
    FFPinTA.N.NT <- c(NULL)
    FFPinTA.perc.NT <- c(NULL)
    
    # create empty structure to fold the stats to plot
    FFP.TA.mx.SP <- list(NULL)
    FFP.TA.mx.SP.p <- list(NULL)
    
    i.span <- 1
    for(i.span in 1:(length(FFP.res.list)+1)){ # first loop for the whole dataset, subsequent loops for the periods
      
      if (i.span == 1){ # WHOLE PERIOD ANALYSIS  ********************************
        
        # 48H dataset analysis ..................................................
        
        # Reference number of records (ignore NAs)
        FFP.TA.N[i.span] <- length(which(!is.na(FFP.an.res.4spatial[, 'FFP_TA_CONTR'])))
        
        # Percentage of data having more than 70% contribution from the TA
        FFPinTA.N[i.span] <- length(which(FFP.an.res.4spatial[, 'FFP_TA_CONTR'] > TA.contr.target)) # number of records in TA
        FFPinTA.perc[i.span] <- FFPinTA.N[i.span]/FFP.TA.N[i.span]                                             # percentage
        if (FFPinTA.perc[i.span] < TA.record.target) {
          warning(paste('Wait! Only the', round(FFPinTA.perc[i.span] * 100, 2), '% of the actual whole dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        } else {
          cat(paste('OK! The', round(FFPinTA.perc[i.span] * 100, 2), '% of the actual whole dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        }
        
        # Daytime dataset analysis ..............................................
        daytime <- FFP.an.res.4spatial[, 'daytime'] == 1
        
        # Reference number of records
        FFP.TA.N.DT[i.span] <- length(which(!is.na(FFP.an.res.4spatial[, 'FFP_TA_CONTR'][daytime])))
        
        # Percentage of data having more than 70% contribution from the TA
        FFPinTA.N.DT[i.span] <- length(which(FFP.an.res.4spatial[, 'FFP_TA_CONTR'][daytime] > TA.contr.target)) # number of records
        FFPinTA.perc.DT[i.span] <- FFPinTA.N.DT[i.span]/FFP.TA.N.DT[i.span]
        if (FFPinTA.perc.DT[i.span] < TA.record.target) {
          warning(paste('Wait! Only the', round(FFPinTA.perc.DT[i.span] * 100, 2), '% of the actual daytime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        } else {
          cat(paste('OK! The', round(FFPinTA.perc.DT[i.span] * 100, 2), '% of the actual daytime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        }
        
        
        # Nightime dataset analysis .............................................
        nightime <- !daytime
        
        # Reference number of records
        FFP.TA.N.NT[i.span] <- length(which(!is.na(FFP.an.res.4spatial[, 'FFP_TA_CONTR'][nightime])))
        
        # Percentage of data having more than 70% contribution from the TA
        FFPinTA.N.NT[i.span] <- length(which(FFP.an.res.4spatial[, 'FFP_TA_CONTR'][nightime] > TA.contr.target)) # number 
        FFPinTA.perc.NT[i.span] <- FFPinTA.N.NT[i.span]/FFP.TA.N.NT[i.span]                                                                       # percentage
        if (FFPinTA.perc.NT[i.span] < TA.record.target) {
          warning(paste('Wait! Only the', round(FFPinTA.perc.NT[i.span] * 100, 2), '% of the actual nightime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        } else {
          cat(paste('OK! The', round(FFPinTA.perc.NT[i.span] * 100, 2), '% of the actual nightime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        }
        
        # Create data for FFP TA Whole Period analysis visualization
        FFP.TA.mx.WP <- matrix(rep(NA, 6) , nrow=2)
        colnames(FFP.TA.mx.WP) <- c("whole data","daytime","nightime")
        rownames(FFP.TA.mx.WP) <- c("IN","OUT")
        FFP.TA.mx.WP[, "whole data"] <- c(FFPinTA.N[i.span], FFP.TA.N[i.span]-FFPinTA.N[i.span])
        FFP.TA.mx.WP[, "daytime"] <- c(FFPinTA.N.DT[i.span], FFP.TA.N.DT[i.span]-FFPinTA.N.DT[i.span])
        FFP.TA.mx.WP[, "nightime"] <- c(FFPinTA.N.NT[i.span], FFP.TA.N.NT[i.span]-FFPinTA.N.NT[i.span])
        #Transform this data in %
        FFP.TA.mx.WP.p <- apply(FFP.TA.mx.WP, 2, function(x) x * 100 / sum(x, na.rm=T))
        
        
      } else { # SUB-PERIODS ANALYSIS *******************************************
        
        j <- i.span-1
        
        # Reference number of records
        FFP.TA.N[i.span] <- length(which(!is.na(FFP.res.list[[j]][, 'FFP_TA_CONTR'])))
        
        # Percentage of data having more than 70% contribution from the TA
        FFPinTA.N[i.span] <- length(which(FFP.res.list[[j]][, 'FFP_TA_CONTR'] > TA.contr.target))
        FFPinTA.perc[i.span] <- FFPinTA.N[i.span]/FFP.TA.N[i.span]
        if (FFPinTA.perc[i.span]  < TA.record.target) {
          warning(paste('Wait! Only the', round(FFPinTA.perc[i.span] * 100, 2), ' of the period', levels(subperiods)[j], 'have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        } else {
          cat(paste('OK! The', round(FFPinTA.perc[i.span] * 100, 2), '% of the', levels(subperiods)[j], 'dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        }
        
        # Daytime analysis ......................................................
        daytime <- FFP.res.list[[j]][, 'daytime'] == 1
        
        # Reference number of records
        FFP.TA.N.DT[i.span] <- length(which(!is.na(FFP.res.list[[j]][, 'FFP_TA_CONTR'][daytime])))
        
        # Percentage of data having more than 70% contribution from the TA
        FFPinTA.N.DT[i.span] <- length(which(FFP.res.list[[j]][, 'FFP_TA_CONTR'][daytime] > TA.contr.target)) # number of records
        FFPinTA.perc.DT[i.span] <- FFPinTA.N.DT[i.span]/FFP.TA.N.DT[i.span]
        if (FFPinTA.perc.DT[i.span] < TA.record.target) {
          warning(paste('Wait! Only the', round(FFPinTA.perc.DT[i.span] * 100, 2), '% of the period', levels(subperiods)[j], ' during daytime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        } else {
          cat(paste('OK! The', round(FFPinTA.perc.DT[i.span] * 100, 2), '% of the period', levels(subperiods)[j], ' during daytime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
        }
        
        # Nightime dataset analysis .............................................
        nightime <- !daytime
        
        # Reference number of records
        FFP.TA.N.NT[i.span] <- length(which(!is.na(FFP.res.list[[j]][, 'FFP_TA_CONTR'][nightime])))
        
        if (FFP.TA.N.NT[i.span] != 0) {
          # Percentage of data having more than 70% contribution from the TA
          FFPinTA.N.NT[i.span] <- length(which(FFP.res.list[[j]][, 'FFP_TA_CONTR'][nightime] > TA.contr.target))
          FFPinTA.perc.NT[i.span] <- FFPinTA.N.NT[i.span]/FFP.TA.N.NT[i.span]
          if (FFPinTA.perc.NT[i.span] < TA.record.target) {
            warning(paste('Wait! Only the', round(FFPinTA.perc.NT[i.span] * 100, 2), '% of the period', levels(subperiods)[j], ' during nightime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
          } else {
            cat(paste('OK! The', round(FFPinTA.perc.NT[i.span] * 100, 2), '% of the period', levels(subperiods)[j], ' during nightime dataset have a cumulated contribution of at least', TA.contr.target*100, '% from the Target Area\n'))
          }
        } else {
          FFPinTA.N.NT[i.span] <- 0
          FFPinTA.perc.NT[i.span] <- 0
          warning(paste('Wait! There are NO data during nightime in this subperiod with any contribution from the Target Area\n'))
        }
        
        # Create data for sub-periods FFP TA analysis visualization
        FFP.TA.mx.SP[[j]] <-  matrix(rep(NA, 6) , nrow=2)
        colnames(FFP.TA.mx.SP[[j]]) <- c("whole data","daytime","nightime")
        rownames(FFP.TA.mx.SP[[j]]) <- c("IN","OUT")
        # whole dataset column
        FFP.TA.mx.SP[[j]][, "whole data"] <- c(FFPinTA.N[i.span], FFP.TA.N[i.span] - FFPinTA.N[i.span])
        # daytime column
        FFP.TA.mx.SP[[j]][, "daytime"] <- c(FFPinTA.N.DT[i.span], FFP.TA.N.DT[i.span]-FFPinTA.N.DT[i.span])
        # nightime column
        FFP.TA.mx.SP[[j]][, "nightime"] <- c(FFPinTA.N.NT[i.span], FFP.TA.N.NT[i.span]-FFPinTA.N.NT[i.span])
        # Transform this data in %
        FFP.TA.mx.SP.p[[j]] <-  matrix(rep(NA, 6) , nrow=2)
        FFP.TA.mx.SP.p[[j]] <- apply(FFP.TA.mx.SP[[j]], 2, function(x) x * 100 / sum(x, na.rm=T))
        
      }
      
    } # period loops end
    
    ## Stacked bar-plot 
    # cnovert stats list to df (and merge whole and sub-periods results)
    FFP.TA.mx.SP.p <- do.call('rbind', FFP.TA.mx.SP.p)
    FFP.TA.mx <- rbind(FFP.TA.mx.WP.p, FFP.TA.mx.SP.p)
    FFP.TA.df <- as.data.frame(FFP.TA.mx, row.names = F)
    # create a column for periods
    FFP.TA.df$period <- rep(c("whole period", paste0("sub-period #", seq_along(FFP.res.list))), each=2)
    # melt it and rename the automatic name "variable" into cat (category)
    FFP.TA.df <- melt(FFP.TA.df, 'period', value.name = "FFP_perc")
    names(FFP.TA.df)[2] <- 'cat'
    # properly sort it
    FFP.TA.df <- plyr::arrange(FFP.TA.df, period, cat) 
    # add the boolean variable for 'IN','OUT' cases
    FFP.TA.df$FFP_TA <- c('IN','OUT')
    # create a period ordering factor
    FFP.TA.df$period_f = factor(FFP.TA.df$period,  levels = c("whole period", paste0("sub-period #", seq_along(FFP.res.list))))
    # rename the levels in cat factor
    FFP.TA.df$cat <- plyr::revalue(FFP.TA.df$cat, c("whole data"="24H", "daytime"="D", 'nightime'='N')) 
    
    
    FFP_perc_label <- round(FFP.TA.df$FFP_perc, 1)
    FFP_perc_label[seq_along(FFP_perc_label)[c(F,T)]] <- ''
    
    # Do the barplot
    ffpbp <-  ggplot(data=FFP.TA.df, aes(x=cat, y=FFP_perc, fill=FFP_TA)) +
      geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
      # according to periods
      facet_grid(.~ period_f) +
      # add percentages
      geom_text(aes(y=FFP_perc, label=FFP_perc_label), vjust=1.6, 
                color="white", size=2.8)+
      # scale_fill_brewer(palette="Paired")+
      scale_fill_brewer(palette="Dark2")+
      theme_minimal() +
      # add threshold line
      geom_hline(yintercept=70, linetype="dashed", color = "black") +
      xlab(" ") +
      scale_y_continuous("% of records") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio = 1.5,
            legend.position="right",
            legend.title = element_text(size = 8),
            legend.text = element_text(size =8),
            legend.key.size = unit(.5, "cm")
      ) + 
      # Add titles
      labs(
        # title = "TEST 2 results",
        # subtitle = paste0('N of observations in each sub-period = ', table(subperiods)[1]),
        subtitle = paste0('TEST 2 results - ', site.ID),
        caption = paste0("Total N of observations = ", FFP.data.summary$N_data_SpatAnal, ' , ', table(subperiods)[1], ' in each subperiod.')) +
      theme(
        # plot.title = element_text(hjust = 0, size = 14),
        plot.subtitle = element_text(hjust = 0, size=10, vjust=-2),    
        plot.caption = element_text(hjust = 1, size=7, vjust=10)
      )
    
    # # visualize it or
    # windows()
    # ffpbp
    
    # save it
    png(paste(FFP.dir, site.ID, '_','FFP-TA_results_', format(Sys.Date(), '%Y%m%d'), '.png', sep=''), width = 18, height = 6, unit = 'cm', res=300)
    ffpbp
    dev.off()
    
    
    # save stats table
    FFP.TA.stats <- data.frame(
      ref_period = c('whole period', substr(levels(subperiods), 1, 7)),
      N_records = FFP.TA.N, 
      N_days = round(FFP.TA.N/48, 1),
      N_rec_inTA = FFPinTA.N,
      N_records_DT = FFP.TA.N.DT, 
      N_rec_inTA_DT = FFPinTA.N.DT,
      N_records_NT = FFP.TA.N.NT, 
      N_rec_inTA_NT = FFPinTA.N.NT
    )
    
    
    write.table(FFP.TA.stats, paste(FFP.dir, site.ID, '_FFP_TA_Nrec_', format(Sys.Date(),'%Y%m%d'), '.csv', sep=''), sep='\t', quote = F, row.names = F)
    
    
    
    # ******************************************************************************
    ## . . TEST 3 - LCT representativeness in the TA -------------------------------
    # ******************************************************************************
    #
    # TEST
    # During Step 1 a map of the different vegetation classes present in the Target Area is submitted. 
    # These information are heterogeneous (in terms of quality and quantity of the information) across the sites and for
    # this reason will be complemented using the data collected during the station characterization. 
    # The aim is to identify areas that are characterized by different species composition or different management (and consequently biomass and density).
    # It is assumed that other differences (e.g. soil) will be reflected in the vegetation characteristics. 
    # The test will analyze, using the same footprint model (Kljun et al. 2015) the amount of records 
    # coming from the different ecosystems and check their representativeness in terms of day-night conditions and in the period analyzed.
    #
    # TARGET
    # For the ecosystem patches in the footprint that contribute with at least 70% of the fluxes in at least 20% of
    # the data (good data after filtering for QAQC) it will be calculated the number of halfhourly measurements
    # collected during daytime, nighttime and for each of two periods obtained dividing the dataset in two parts.
    # The target values is that each group includes at least 20% of data.
    #
    # ACTION (if the target value is not reached)
    # An in depth analysis of the conditions that lead to have an ecosystem patch contributing to the
    # measurements but with under-representativeness of specific conditions is conducted and a discussion with
    # the PI is started in order to find solutions and alternatives, in particular changing the measurement height or
    # wind sectors to exclude.
    # ******************************************************************************
    
    
    # TRHRESHOLDS (according to the MSA document)
    # for the main ecosystem types, at least 20% of the data during night and during day and also distributed along the period analysed
    LCT.record.target <- 0.2 # FREQUENCY THRESHOLD %
    LCT.contr.target <- 0.7  # FLUX CONTRIBUTION THRESHOLD %
    
    # Reference number of records (ignore NAs): use the records for which the TA contribution id more than 70%
    LCT.ref.rec <- FFPinTA.N[1] # The whole period record from the T2 analysis
    
    # Create an ad-hoc dataframe for LCT results
    FFP.LCT.df <- FFP.an.res.4spatial[, c(1,2, grep('LCT', names(FFP.an.res.4spatial)))]
    
    # # summary of the contribution for each LCT
    lct.contr.summary <- apply(FFP.LCT.df[, grep('CONTR', names(FFP.LCT.df))], 2, summary)
    
    
    ## Overall contribution from all areas ......................................
    
    # locate the columns with the weighted FFP area contribution
    LCT.contr.indx <- grep('CONTR', names(FFP.LCT.df))
    
    ## Flag for good halfhours (1 ok, 0 no ok):
    # *** half hours with a cumulated contribution of at least 70% from any LCT Area ***
    FFP.LCT.df <- cbind(FFP.LCT.df, 
                        matrix(0, 
                               nrow=nrow(FFP.LCT.df), ncol=length(LCT.contr.indx), 
                               dimnames = list(NULL, paste('OK_FLAG', sprintf("%02.f", 1:length(LCT.contr.indx)), sep=''))))
    
    # frequency of contribution split into 10% classes 
    FFP.LCT.Num <- apply(FFP.LCT.df[, LCT.contr.indx], 2, function(x) tapply(x, cut(x, seq(0,0.8,0.1)), length)) # numerosity
    FFP.LCT.Freq <- round(apply(FFP.LCT.df[, LCT.contr.indx], 2, function(x) tapply(x, cut(x, seq(0,0.8,0.1)), length) / LCT.ref.rec * 100), 2) # percentage
    FFP.LCT.Freq <- as.data.frame(FFP.LCT.Freq)
    FFP.LCT.Freq$'Contribution' <- row.names(FFP.LCT.Freq)
    FFP.LCT.Freq <- FFP.LCT.Freq[, c(ncol(FFP.LCT.Freq),1:(ncol(FFP.LCT.Freq)-1))]
    
    write.table(FFP.LCT.Freq, paste(FFP.dir, site.ID, '_FFP-LCT_Freq_', format(Sys.Date(),'%Y%m%d'), '.csv', sep=''), sep='\t', quote = F, row.names = FALSE)
    
    
    # Find out the LCT which is/are sampled the most *******************************
    # Current rule: 70% contribution in 20% of data at least
    
    # histograms of LCT contributions
    png(paste(FFP.dir, site.ID, '_','FFP-LCT_contributions_', format(Sys.Date(), '%Y%m%d'), '.png', sep=''), width = 13, height = 13, unit = 'cm', res=300)
    # windows()
    
    lct.avgcont <- c()
    ylimup <- 50
    for(i in 1:length(LCT.contr.indx)){
      
      cur.hist <- hist(FFP.LCT.df[, LCT.contr.indx[i]], plot=FALSE, breaks=seq(0,1,0.02))
      cur.hist$counts <- cur.hist$counts / length(FFP.an.res.4spatial$FFP_TA_CONTR >= 0.7) * 100 # calculate frequency respective to the tota number of observations
      # lct.cumcont[i] <- round(sum(cur.hist$counts * cur.hist$mids), 1) # calculate the cumulative relative frequency of the LCTs
      lct.avgcont[i] <- round(mean(FFP.LCT.df[, LCT.contr.indx[i]], na.rm=T)*100, 1) # average cumulative LCTs contribution
      if (i == 1) {plot(cur.hist, 
                        freq = TRUE,
                        col = adjustcolor(LCT.color[i], 0.5), border=NA, 
                        xaxs='i', yaxs='i', xlab='LCTs cumulative contribution', ylab='Frequency [% ]', main='', las=1,
                        xlim=c(0,.82), ylim=c(0,ylimup)
      )
      } else {
        plot(cur.hist, add=TRUE, freq = TRUE, col = adjustcolor(LCT.color[i], 0.5), border=NA)
      }
    }
    abline(v=0.7, lty=32, col=2); text.default(0.7, ylimup, pos=2,'LCT contribution target', cex=.8, srt = 90)
    # rect(xleft=0.7, ybottom=0, xright=.8, ytop=100, density = NULL, angle = NULL,col = NA, border = 2, lwd=1, lty=32)
    legend('topleft', legend=paste('LCT_0', seq_along(LCT.list), ' (', lct.avgcont, '%)', sep=''), 
           col = adjustcolor(LCT.color, alpha.f = 0.6), pch=15, pt.cex=1.5, bty='n', cex=.8)
    
    dev.off()
    
    # ## percentage of LCT contributions > 70%
    # apply(FFP.LCT.df[, LCT.contr.indx], 2, function(x) length(which(x > LCT.contr.target))) / LCT.ref.rec
    
    # Main LCT:
    LCT.sampled <- which(apply(FFP.LCT.df[, LCT.contr.indx], 2, function(x) 
      length(which(x > LCT.contr.target))) / LCT.ref.rec > LCT.record.target)
    
    
    # LCT contribution analysis ****************************************************
    
    if (length(LCT.sampled) == 1) {
      
      cat(note('YO! There is only one LCT which substantially contributes to the fluxes. The TEST 3 is consequently not needed.\n '))
      
    } else {
      
      # index list for the 2 subperiods
      FFP.LCT.df.SG.indx <- list(SG1 = c(1:trunc((nrow(FFP.LCT.df)/2))),
                                 SG2 = c((trunc((nrow(FFP.LCT.df)/2))+1)):nrow(FFP.LCT.df))
      
      LCT.dtime <- list(NULL)
      LCT.ntime <- list(NULL)
      FFPinLCT.N.dt <- list(NULL)
      FFPinLCT.N.nt <- list(NULL)
      FFPinLCT.perc.dtdt <- list(NULL)
      FFPinLCT.perc.ntnt <- list(NULL)
      FFPinLCT.perc.dt <- list(NULL)
      FFPinLCT.perc.nt <- list(NULL)
      
      i <- 1
      for (i in 1:length(FFP.LCT.df.SG.indx)) {
        
        # current daytime index
        LCT.dtime[[i]] <- c(NULL)
        LCT.dtime[[i]] <- FFP.LCT.df$daytime[FFP.LCT.df.SG.indx[[i]]] == 1
        # current nighttime index
        LCT.ntime[[i]] <- c(NULL)
        LCT.ntime[[i]] <- FFP.LCT.df$daytime[FFP.LCT.df.SG.indx[[i]]] == 0
        
        # how many records with at least 70% of the fluxes
        
        ### During daytime
        cur.dt.df <- data.frame()
        cur.dt.df <- FFP.LCT.df[FFP.LCT.df.SG.indx[[i]], LCT.contr.indx[LCT.sampled]][LCT.dtime[[i]], ]
        
        # absolute records
        FFPinLCT.N.dt[[i]] <- apply(cur.dt.df, 2, function(x) which(x > LCT.contr.target))
        # percentage OVER DAYTIME DATA
        FFPinLCT.perc.dtdt[[i]] <- lapply(FFPinLCT.N.dt[[i]], function(x) length(which(x > LCT.contr.target)) / length(which(LCT.dtime[[i]])) * 100) 
        # percentage OVER WHOLE SUBPERIOD DATA
        FFPinLCT.perc.dt[[i]] <- lapply(FFPinLCT.N.dt[[i]], function(x) length(which(x > LCT.contr.target)) / length(FFP.LCT.df.SG.indx[[i]]) * 100) 
        
        ### During nighttime
        cur.nt.df <- data.frame()
        cur.nt.df <- FFP.LCT.df[FFP.LCT.df.SG.indx[[i]], LCT.contr.indx[LCT.sampled]][LCT.ntime[[i]], ]
        
        # absolute records
        FFPinLCT.N.nt[[i]] <- apply(cur.nt.df, 2, function(x) which(x > LCT.contr.target))
        # percentage OVER DAYTIME DATA
        FFPinLCT.perc.ntnt[[i]] <- lapply(FFPinLCT.N.nt[[i]], function(x) length(which(x > LCT.contr.target)) / length(which(LCT.ntime[[i]])) * 100)
        # percentage OVER WHOLE SUBPERIOD DATA
        FFPinLCT.perc.nt[[i]] <- lapply(FFPinLCT.N.nt[[i]], function(x) length(which(x > LCT.contr.target)) / length(FFP.LCT.df.SG.indx[[i]]) * 100) 
        
      }
      
      FFPinLCT.T3.out <- data.frame(NULL)
      FFPinLCT.T3.out <- data.frame(
        subperiod = rep(rep(c(1,2), each=length(FFPinLCT.N.dt)), 2),
        LCT = rep(substr(names(LCT.sampled), 5,9), length(FFPinLCT.N.dt)),
        ToD = rep(c('daytime','nighttime'), each=length(FFPinLCT.N.dt)*2),
        N_record = c(rep(unlist(lapply(LCT.dtime, function(x) length(which(x)))), each=length(LCT.sampled)),
                     rep(unlist(lapply(LCT.ntime, function(x) length(which(x)))), each=length(LCT.sampled))),
        N_record_inLCT = c(unlist(lapply(FFPinLCT.N.dt, function(x) lapply(x, length))),
                           unlist(lapply(FFPinLCT.N.nt, function(x) lapply(x, length)))),
        Perc_record_inLCT = round(c(unlist(FFPinLCT.perc.dtdt), unlist(FFPinLCT.perc.ntnt)), 1)
      )
      
      FFPinLCT.T3.out <- FFPinLCT.T3.out[order(FFPinLCT.T3.out$'subperiod'), ]
      
      
      # Save output data for FFP/LCT interaction
      write.table(FFPinLCT.T3.out, paste(FFP.dir, site.ID, '_TEST3_results.txt', sep = ''), quote = F, row.names=F, sep='\t')
      
      # Histograms of LCT DAY-NIGHT percentage contributions
      FFPLCT.mtx <- matrix(FFPinLCT.T3.out$'Perc_record_inLCT', nrow=length(FFPinLCT.N.dt)*2, byrow = T)
      colnames(FFPLCT.mtx) <- substr(names(LCT.sampled), 5,9)
      rownames(FFPLCT.mtx) <- c('DT_SP1','NT_SP1','DT_SP2','NT_SP2')
      FFPLCT.mtx <- t(FFPLCT.mtx)
      
      png(paste(dir.output, site.ID, '_','FFP-LCT_percentages_', format(Sys.Date(), '%Y%m%d'), '.png', sep=''), width = 13, height = 13, unit = 'cm', res=300)
      # windows()
      
      FFPLCT.bp <- barplot(FFPLCT.mtx, beside = TRUE,
                           col = LCT.color,
                           legend = rownames(FFPLCT.mtx), ylim = c(0,100), 
                           main = "TEST 3 RESULTS", font.main = 1, sub = " ", cex.names = 1, las=1)
      mtext('Frequency', 2, 2.5)
      # segments(FFPLCT.bp, FFPLCT.mtx, FFPLCT.bp, FFPLCT.mtx + 2*sqrt(1000*FFPLCT.mtx/100), col = mybarcol, lwd = 1.5) # Faked upper 2*sigma error bars
      stopifnot(dim(FFPLCT.bp) == dim(FFPLCT.mtx))  # corresponding matrices
      abline(h = 20, col=2)
      # mtext(side = 1, at = FFPLCT.bp, line = -2, text = paste('SP', FFPinLCT.T3.out$'subperiod', FFPinLCT.T3.out$'ToD'), cex=.8)
      text(x = FFPLCT.bp, y = FFPinLCT.T3.out$Perc_record_inLCT, formatC(FFPLCT.mtx), cex=.8, pos=1)
      
      dev.off()
      
    }
    
    
    # ***************************************************************************
    ## . . TEST X - EA contibution analysis -------------------------------------
    # ***************************************************************************
    # TEST
    # This is not an ICOS standard: text to do
    # TARGET
    # to do 
    # ACTION (if the target value is not reached)
    # to do
    # ******************************************************************************
    
    # TODO: import from desktop code and update 
    
  }
  
  
  
  
  # FOOTPRINT CLIMATOLOGY =====================================================
  
  ## . 24 H climatology -------------------------------------------------------
  if (DO.climatology) {
    
    cat(note('\n............................................\n'))
    cat(prog.mes('Compute the footprint climatology:\n'))
    
    FFP.input.4clim <- FFP.input
    
    # Remove bad flagged data after having created an ad-hoc dataframe for FFP climatology
    # . Remove rocords for which NEE  was discarded from QC procedure
    FFP.input.4clim[which(FFP.input$'qcf'), ] <- NA
    # . Remove records with not possible combination of roughness length and turbulence input (N. Kljun pers. comm.)
    # which cause the footprint ot be extremely short and wide
    if (remove.SnW) { FFP.input.4clim[which(FFP.res$'SnW_FFP_FLAG' == 1), ] <- NA}
    
    # Possibly further remove rows if even 1 record of any variable is NA
    FFP.input.4clim <- na.omit(FFP.input.4clim)
    
    # Compute footprint climatology
    FFP.clim <- calc_footprint_FFP_climatology(
      
      # mandatory input
      zm = FFP.input.4clim$'zm',                      # Measurement height above displacement height (i.e. z‐d) [m]
      z0 = FFP.input.4clim$'z0',                      # Roughness length [m] ‐ enter [NaN] if not known
      umean = FFP.input.4clim$'umean',                # Mean wind speed at zm [ms‐1] ‐ enter [NaN] if not known (Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint)
      h = FFP.input.4clim$'PBL',                      # Vector of boundary layer height [m]
      ol = FFP.input.4clim$'ol',                      # Vector of Obukhov length [m]
      sigmav = FFP.input.4clim$'sigmav',              # Vector of standard deviation of lateral velocity fluctuations [ms‐1]
      ustar = FFP.input.4clim$'ustar',                # Vector of friction velocity [ms‐1]
      wind_dir = FFP.input.4clim$'wind_dir',          # Vector of wind direction in degrees (of 360) for rotation of the footprint
      # optional input
      domain = ffp.domain,           # Domain size as an array of (xmin xmax ymin ymax) [m] *1.
      dx = cell.size, dy = cell.size,# Cell size of domain [m] (default is dx = dy = 2 m)
      nx = 1000, ny = 1000,          # Two integer scalars defining the number of grid elements in x and y (default is nx = ny = 1000)
      r = FFP.R,                     # Percentage of source area, i.e. a value between 10% and 90%. *2
      rslayer =  1,                  # Calculate footprint even if zm within roughness sublayer: set rslayer = 1 *3
      smooth_data = 1,               # Apply convolution filter to smooth footprint climatology if smooth_data = 1 (default)
      pulse = 10,                    # Display progress of footprint calculations every pulse‐th footprint (e.g., "100")
      fig = 0
    )
    
    # Number of footprints actually calculated and included in footprint climatology:
    FFP.data.summary$'N_data_FFPclim' <- FFP.clim$'n'
    # as % of original data
    FFP.data.summary$'N_data_FFPclim_perc' <- round(FFP.data.summary$'N_data_FFPclim' / nrow(FFP.input.4clim) * 100, 1)
    # save.image(file = paste0(FFP.dir, site.ID,  '_T2_data_', format(Sys.Date(), '%Y%m%d'), '.Rdata'))
    
    # # Two‐dimensional average footprint climatology - RASTER
    # windows()
    # x.lim <- y.lim <- c(-1000,1000)
    # # plot(TA, axes=T)
    # image.plot(x = FFP.clim$x_2d[1, ], y = FFP.clim$y_2d[ , 1], z = FFP.clim$fclim_2d,
    #            col = ffp2Dcol, ylab='Northing [m]', xlab='Easting [m]',
    #            xlim=x.lim, ylim=y.lim)
    # # EP FFP mean peak
    # abline(v=c(-median(set.cleaned$'x_peak', na.rm=T), median(set.cleaned$'x_peak', na.rm=T)), lty=32, col='white')
    # # EP FFP 70% cumulative contribution distance
    # abline(v=c(-median(set.cleaned$'x_70.', na.rm=T), median(set.cleaned$'x_70.', na.rm=T)), lty=32, col='red')
    # # EC tower
    # points(0,0, pch=3)
    # # isoplethes
    # j<-0
    # for (j in 1:length(FFP.R)){
    #   lines(FFP.clim$xr[[j]], FFP.clim$yr[[j]], type="l", col="black")
    # }
    # lines(FFP.clim$xr[[X.perc.indx[1]]], FFP.clim$yr[[X.perc.indx[1]]], type="l", col="red", lwd=1)
    
    
    # Footprint climatology spatial visualization .............................
    
    # Isoplethe relative coordinates (northward/eastward)
    FFP.CLIM.iso.df <- data.frame(NULL)
    FFP.CLIM.iso.df <- data.frame(LON = as.numeric(FFP.clim$xr[[X.perc.indx[1]]]),
                                  LAT = as.numeric(FFP.clim$yr[[X.perc.indx[1]]]),
                                  row.names = NULL)
    FFP.CLIM.iso.df <- na.omit(FFP.CLIM.iso.df)
    # NORTH/EAST to UTM and convert to object of class Polygons
    FFP.CLIM.iso.utm.p <- Polygon(cbind(FFP.CLIM.iso.df$LON + coordinates(EC.tower.utm)[,1],
                                        FFP.CLIM.iso.df$LAT + coordinates(EC.tower.utm)[,2]))
    # convert to object of class SpatialPolygons (assign projection)
    FFP.CLIM.iso <- NULL
    FFP.CLIM.iso <- SpatialPolygons(list(Polygons(list(FFP.CLIM.iso.utm.p), 'FFP.CLIM.iso')),
                                    proj4string = EC.tower.utm@proj4string)
    ## check for validity of the FF80 Spatial polygon (self-intersection error)
    if (!gIsValid(FFP.CLIM.iso, reason = F)) { # if gIsValid returns FALSE (NO OK):
      # use a zero-width buffer to clean up the topology problem
      FFP.CLIM.iso <- gBuffer(FFP.CLIM.iso, byid=TRUE, width=0)
    }
    
    # Exctract the 80% contribution footprint values (matrix)
    FFP.CLIM.mx <- matrix(NA)
    FFP.CLIM.mx <- FFP.clim$fclim_2d
    FFP.CLIM.mx[FFP.CLIM.mx < FFP.clim$fr[X.perc.indx[2]]] <- NA
    # NORTH/EAST to UTM and promote to object of class RasterLayer
    FFP.CLIM.r <- raster(list(x = FFP.clim$x_2d[1, ] + coordinates(EC.tower.utm)[,1],
                              y = FFP.clim$y_2d[, 1] + coordinates(EC.tower.utm)[,2],
                              z = FFP.CLIM.mx),
                         crs = EC.tower.utm@proj4string)
    proj4string(FFP.CLIM.r) <- EC.tower.utm@proj4string
    
    
    # Whole spatial data visualization
    
    png(paste(FFP.dir, site.ID, '_', 'FFPclim_', format(Sys.Date(), '%Y%m%d'), '_r.png', sep=''), width = 13, height = 13, unit = 'cm', res=300)
    # windows()
    
    # TA + EC
    plot(TA.utm, axes=T, ylab='Northing [m]', xlab='Easting [m]')
    points(EC.tower.utm, pch=3)
    # EA
    if (!is.na(EA.num)){
      plot(EA.utm, add=T, col=adjustcolor(8, alpha.f = 0.5), border=NA)
      legend('topleft', legend=c('EC','TA','EA'), col = c(1,1, adjustcolor(8, alpha.f = 0.2)), pch=c(3,0,15), bty='n')
    } else {
      legend('topleft', legend=c('EC','TA'), col = 1, pch=c(3,0), cex=0.8, bty='n')
    }
    
    # LCT
    if (LCT.present) {
      j <- 0
      for(j in 1:length(LCT.list)){
        plot(LCT.list.utm[[j]], add=T, border = LCT.color[j], col=adjustcolor(LCT.color[j], alpha.f = 0.5))
        text(x=coordinates(LCT.list.utm[[j]])[1], y=coordinates(LCT.list.utm[[j]])[2], LCT.list[[j]]@data$Description, pos=1, cex=.5)
      }
      legend('topright', legend=paste('LCT_0', seq_along(LCT.list), sep=''), col = adjustcolor(LCT.color, alpha.f = 0.4), pch=15, bty='n')
    }
    plot(TA.utm, lwd=2, add=T)
    
    # FFP2D isoplethes
    j<-0
    for (j in 1:length(FFP.R)){
      lines(FFP.clim$xr[[j]] + coordinates(EC.tower.utm)[,1], FFP.clim$yr[[j]] + coordinates(EC.tower.utm)[,2],
            type="l", col=rich.colors(9)[2:4][j], lwd=2)
    }
    # lines(FFP.clim$xr[[X.perc.indx[1]]] + coordinates(EC.tower.utm)[,1], FFP.clim$yr[[X.perc.indx[1]]] + coordinates(EC.tower.utm)[,2], type="l", col="red", lwd=2)
    
    if (show.clim.raster){
      # FFP2D CLIM footprint
      plot(FFP.CLIM.r, col=adjustcolor(ffp2Dcol, alpha.f = 0.8), legend = FALSE, add = TRUE)
      # FFP2D custom legend ex: https://stackoverflow.com/questions/9436947/legend-properties-when-legend-only-t-raster-package
      # plot(r, legend.only=TRUE, col=topo.colors(100),
      #      legend.width=1, legend.shrink=0.75,
      #      axis.args=list(at=seq(r.range[1], r.range[2], 25),
      #                     labels=seq(r.range[1], r.range[2], 25),
      #                     cex.axis=0.6),
      #      legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))
      # legend('bottomleft', c('FFP 70% contribution'), pch='0', col=c(2), bty='n')
    } else {
      legend('bottomleft', c('FFP 50%', 'FFP 70%','FFP 80%'), lty=1, col=rich.colors(9)[2:4], lwd=2, bty='n')
    }
    
    dev.off()
  }
  cat(prog.mes('done!\n'))
  
  ## . Nighttime climatology ---------------------------------------------------
  if (DO.climatology.NT) {
    
    cat(note('\n............................................\n'))
    cat(prog.mes('Compute the nightime footprint climatology:\n'))
    
    # extract nightime time data
    FFP.input.4clim.NT <- FFP.input.4clim
    FFP.input.4clim.NT <- FFP.input.4clim.NT[FFP.input.4clim.NT$'daytime' == 0, ]
    
    # Compute footprint climatology
    FFP.clim.NT <- calc_footprint_FFP_climatology(
      
      # mandatory input
      zm = FFP.input.4clim.NT$'zm',                      # Measurement height above displacement height (i.e. z‐d) [m]
      z0 = FFP.input.4clim.NT$'z0',                      # Roughness length [m] ‐ enter [NaN] if not known
      umean = FFP.input.4clim.NT$'umean',                # Mean wind speed at zm [ms‐1] ‐ enter [NaN] if not known (Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint)
      h = FFP.input.4clim.NT$'PBL',                      # Vector of boundary layer height [m]
      ol = FFP.input.4clim.NT$'ol',                      # Vector of Obukhov length [m]
      sigmav = FFP.input.4clim.NT$'sigmav',              # Vector of standard deviation of lateral velocity fluctuations [ms‐1]
      ustar = FFP.input.4clim.NT$'ustar',                # Vector of friction velocity [ms‐1]
      wind_dir = FFP.input.4clim.NT$'wind_dir',          # Vector of wind direction in degrees (of 360) for rotation of the footprint
      # optional input
      domain = ffp.domain,           # Domain size as an array of (xmin xmax ymin ymax) [m] *1.
      dx = cell.size, dy = cell.size,# Cell size of domain [m] (default is dx = dy = 2 m)
      nx = 1000, ny = 1000,          # Two integer scalars defining the number of grid elements in x and y (default is nx = ny = 1000)
      r = FFP.R,                     # Percentage of source area, i.e. a value between 10% and 90%. *2
      rslayer =  1,                  # Calculate footprint even if zm within roughness sublayer: set rslayer = 1 *3
      smooth_data = 1,               # Apply convolution filter to smooth footprint climatology if smooth_data = 1 (default)
      pulse = 10,                    # Display progress of footprint calculations every pulse‐th footprint (e.g., "100")
      fig = 1
    )
    
    # Nightime spatial data visualization 
    png(paste(FFP.dir, site.ID, '_', 'FFPclim_NT_', format(Sys.Date(), '%Y%m%d'), '_r.png', sep=''), width = 13, height = 13, unit = 'cm', res=300)
    
    # TA + EC
    plot(TA.utm, axes=T, ylab='Northing [m]', xlab='Easting [m]')
    points(EC.tower.utm, pch=3)
    # EA
    if (!is.na(EA.num)){
      plot(EA.utm, add=T, col=adjustcolor(8, alpha.f = 0.5), border=NA)
      legend('topleft', legend=c('EC','TA','EA'), col = c(1,1, adjustcolor(8, alpha.f = 0.2)), pch=c(3,0,15), bty='n')
    } else {
      legend('topleft', legend=c('EC','TA'), col = 1, pch=c(3,0), cex=0.8, bty='n')
    }
    
    # LCT
    if (LCT.present) {
      j <- 0
      for(j in 1:length(LCT.list)){
        plot(LCT.list.utm[[j]], add=T, border = LCT.color[j], col=adjustcolor(LCT.color[j], alpha.f = 0.5))
        text(x=coordinates(LCT.list.utm[[j]])[1], y=coordinates(LCT.list.utm[[j]])[2], LCT.list[[j]]@data$Description, pos=1, cex=.5)
      }
      legend('topright', legend=paste('LCT_0', seq_along(LCT.list), sep=''), col = adjustcolor(LCT.color, alpha.f = 0.4), pch=15, bty='n')
    }
    plot(TA.utm, lwd=2, add=T)
    
    # FFP2D isoplethes
    j<-0
    for (j in 1:length(FFP.R)){
      lines(FFP.clim.NT$xr[[j]] + coordinates(EC.tower.utm)[,1], FFP.clim.NT$yr[[j]] + coordinates(EC.tower.utm)[,2],
            type="l", col=rich.colors(9)[2:4][j], lwd=2)
    }
    # lines(FFP.clim$xr[[X.perc.indx[1]]] + coordinates(EC.tower.utm)[,1], FFP.clim$yr[[X.perc.indx[1]]] + coordinates(EC.tower.utm)[,2], type="l", col="red", lwd=2)
    
    if (show.clim.raster) {
      # FFP2D CLIM footprint
      plot(FFP.CLIM.r, col=adjustcolor(ffp2Dcol, alpha.f = 0.8), legend = FALSE, add = TRUE)
      # plot(r, legend.only=TRUE, col=topo.colors(100),
      #      legend.width=1, legend.shrink=0.75,
      #      axis.args=list(at=seq(r.range[1], r.range[2], 25),
      #                     labels=seq(r.range[1], r.range[2], 25),
      #                     cex.axis=0.6),
      #      legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))
      # legend('bottomleft', c('FFP 70% contribution'), pch='0', col=c(2), bty='n')
    } else {
      legend('bottomleft', c('FFP 50%', 'FFP 70%','FFP 80%'), lty=1, col=rich.colors(9)[2:4], lwd=2, bty='n')
    }
    dev.off()
  }
  cat(prog.mes('done!\n'))
  
} # END OF FUNCTION
