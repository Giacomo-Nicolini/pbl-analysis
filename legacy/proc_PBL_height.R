
proc_PBL_height <- function(
  
  TIMESTAMP,               # timestamp [%Y-%m-%d %H:%M:%S]
  PBL.h,                   # PBL heights 
  plot.PBL.h.freq = TRUE,  # plot the PBL heights frequency histogram
  PBL.dir = NULL           # directory of PBL estimations
  
) {
  
  
  # summary(PBL.h)
  if (plot.PBL.h.freq) {
    
    if (is.null(PBL.dir)) {
      windows() 
    } else {
      png(paste(PBL.dir, 'PBL_heights_histogram.png', sep=''), width = 14, height = 14, unit = 'cm', res=300)  
    }
    # hist(PBL.h, breaks = quantile(PBL.h, seq(0,1,.05), na.rm=T))
    hist(PBL.h, breaks = pretty(PBL.h, 100))
    abline(v=H.EC, col=2, lwd=2) # measurement height
    abline(v=100, col=2)
    
    if (!is.null(PBL.dir)) {dev.off()}
  }
  
  
  # manipultae PBL height
  PBL.h.qc <- PBL.h
  
  # filtro < 0
  if (any(PBL.h.qc < 0, na.rm = TRUE)) {
    PBL.h.qc[PBL.h.qc < 0] <- NA
  }
  

    # # check for PBL heights smaller than measurement height 
  # # summary(PBL.h.qc < H.EC)
  # hm.lessthan.PBL <- PBL.h.qc < H.EC # mask for conditions in which PBL is less than the EC measurement height
  # # NOTE: make a possible warning giving the % of these conditions (once the PBL estimate is robust enough)
  # # possible filter
  # PBL.h.qc[hm.lessthan.PBL] <- 50
  # #(above displacement height)
  # # if(any(PBL.h.qc < (H.EC - H.canopy * 2/3), na.rm = T)){
  # #   PBL.h.qc[PBL.h.qc < (H.EC - H.canopy * 2/3)] <- 100
  # # }
  
  # filtro PBL < 50
  # summary(PBL.h.qc[!hm.lessthan.PBL] < 50) # other than the previously forced values
  if (any(PBL.h.qc < 50, na.rm = TRUE)) {
    PBL.h.qc[PBL.h.qc < 50] <- 50
  }
  
  
  # Approssimazione per rimpiazzare i NAs, gestirla meglio
  # summary(PBL.h)
  # Replace NA values with average day and night values
  dtime <- c(as.numeric(format(TIMESTAMP, '%H')) >= 6 & as.numeric(format(TIMESTAMP, '%H')) <= 18)
  # table(dtime)
  i <- 0
  for (i in seq_along(PBL.h.qc)) {
    if (is.na(PBL.h.qc[i])) {
      if (!is.na(dtime[i]) & dtime[i]) {
        PBL.h.qc[i] <- mean(PBL.h.qc[dtime], na.rm=T)
        # PBL.h.qc[i] <- 1000
      } else {
        PBL.h.qc[i] <- mean(PBL.h.qc[!dtime], na.rm=T)
        # PBL.h.qc[i] <- 100
      }
    }
  }
  
  # summary(PBL.h)
  # summary(PBL.h.qc)
  
  
  if (is.null(PBL.dir)) {
    windows() 
  } else {
    png(paste(PBL.dir, 'PBL_heights_QCed.png', sep=''), width = 14, height = 14, unit = 'cm', res=300)  
  }
  plot(tapply(PBL.h, format(TIMESTAMP, '%H:%M', tz='GMT'), mean, na.rm=T), type='o', ylim=c(0,max(PBL.h, na.rm = T)), col=1, pch=19)
  points(tapply(PBL.h, format(TIMESTAMP, '%H:%M', tz='GMT'), median, na.rm=T), type='o', col=adjustcolor(1, 0.5))
  points(tapply(PBL.h.qc, format(TIMESTAMP, '%H:%M', tz='GMT'), mean, na.rm=T), type='o', col=4, pch=19)
  points(tapply(PBL.h.qc, format(TIMESTAMP, '%H:%M', tz='GMT'), median, na.rm=T), type='o', col=adjustcolor(4, 0.5))
  abline(h=0, v=24)
  legend('topleft', c('PBLH_avg','PBLH_med','PBLH_QC_avg','PBLH_QC_med'), lty=1, pch=c(19,1,19,1), col=c(1,1,4,4))
  
  if (!is.null(PBL.dir)) {dev.off()}
  
  
  
  # # IF OK
  # PBL.h <- PBL.h.qc
  return(PBL.h.qc)
  
  
}

