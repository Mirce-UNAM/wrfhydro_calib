#.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
#library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)

#########################################################
# SETUP
#########################################################

source("calib_utils.R")

# Multi-core
parallelFlag <- TRUE
ncores <- 8
if (parallelFlag && ncores>1) {
        library(doParallel)
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
}


# Metrics
metrics <- c("cor", "rmse", "bias", "nse", "nselog", "nsewt", "kge", "msof")

#########################################################
# MAIN CODE
#########################################################

# First loop check
if (file.exists("proj_data.Rdata")) { 
   load("proj_data.Rdata")
} else {
   # First run so need to initialize
   ReadNamelist("namelist.calib")
   cyclecount <- 0

   # Setup plot directory
   writePlotDir <- paste0(runDir, "/plots")
   dir.create(writePlotDir)

   # Load obs so we have them for next iteration
   load(obsFile)
   obsDT$q_cms <- NULL

   # Find the index of the gage
#   rtLink <- ReadRouteLink(rtlinkFile)
#   rtLink <- data.table(rtLink)
#   linkId <- which(trimws(rtLink$gages) %in% siteId)

   # Setup value lists from paramBnds
   xnames <- paramBnds$param
   x0 <- paramBnds$ini
   names(x0) <- xnames
   x_min <- paramBnds$min
   names(x_min) <- xnames
   x_max <- paramBnds$max
   names(x_max) <- xnames

   # Initialize parameter archive DF
   message("Initialize parameter archive")
   x_archive <- as.data.frame(matrix(, nrow=1, ncol=length(xnames)+2+length(metrics)))
   names(x_archive) <- c("iter", xnames, "obj", metrics)

   # Output parameter set
   x_new <- x0
   cyclecount <- 1

   x_new_out <- c(cyclecount, x_new)
   names(x_new_out)[1] <- "id"
   write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")

   dir.create("archive")

   # Save and exit
   save.image("proj_data.Rdata")
   quit("no")
}

if (cyclecount > 0) {

   file.rename("params_new.txt", paste0("CALIB_RESULTS/OUTPUT", cyclecount, "/params_new.txt"))

   # Read model out and calculate performance metric
   outPath <- paste0(runDir, "/CALIB_RESULTS/OUTPUT", cyclecount)
   print(outPath)

   # Read files
#   message("Reading model out files.")
#   system.time({
#   filesList <- list.files(path = outPath,
#                          pattern = glob2rx("*.CHRTOUT_DOMAIN*"),
#                          full.names = TRUE)
#   filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
#   whFiles <- which(filesListDate >= startDate)
#   filesList <- filesList[whFiles]
#   if (length(filesList) == 0) stop("No matching files in specified directory.")
#   chrt <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
#   })

   # Read files
  message("Reading model out files.")
  frxst <- as.data.table(read.table (paste0(outPath, "/frxst_pts_out.txt"),sep=",",stringsAsFactors=FALSE,
                       col.names=c("time_sec","POSIXct","stn","lon","lat","q_cms","q_cfs","head")))
  frxst$POSIXct <- as.POSIXct(frxst$POSIXct,format = "%Y-%m-%d %H:%M:%S",tz="UTC")

  chrt <- subset(frxst, stn == siteId)

   # Convert to daily
   chrt.d <- Convert2Daily(chrt)
   chrt.d[, site_no := siteId]
   assign(paste0("chrt.d.", cyclecount), chrt.d)
   save(list=c(paste0("chrt.d.", cyclecount)), file=paste0("archive/", paste0("chrt.d.", cyclecount), ".Rdata"))

   # Merge
   setkey(chrt.d, "site_no", "POSIXct")
   setkey(obsDT, "site_no", "POSIXct")
   chrt.d <- merge(chrt.d, obsDT, by.x = c("site_no", "Date"), by.y=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
   F_new <- objFn(chrt.d$q_cms, chrt.d$obs)
   print(F_new)

 # Calc stats
   statCor <- cor(chrt.d$q_cms, chrt.d$obs)
   statRmse <- Rmse(chrt.d$q_cms, chrt.d$obs, na.rm=TRUE)
   statBias <- PBias(chrt.d$q_cms, chrt.d$obs, na.rm=TRUE)
   statNse <- Nse(chrt.d$q_cms, chrt.d$obs, na.rm=TRUE)
   statNseLog <- NseLog(chrt.d$q_cms, chrt.d$obs, na.rm=TRUE)
   statNseWt <- NseWt(chrt.d$q_cms, chrt.d$obs)
   statKge <- Kge(chrt.d$q_cms, chrt.d$obs, na.rm=TRUE)
   statMsof <- Msof(chrt.d$q_cms, chrt.d$obs)

   # Archive results
  # x_archive[cyclecount,] <- c(cyclecount, x_new, F_new)
   x_archive[cyclecount,] <- c(cyclecount, x_new, F_new, statCor, statRmse, statBias, statNse, statNseLog, statNseWt, statKge, statMsof)

   # Evaluate performance metric
   if (cyclecount == 1) {
      x_best <- x_new
      F_best <- F_new
   } else if (F_new <= F_best) {
      x_best <- x_new
      F_best <- F_new
   }

   if (cyclecount < m) {
      # Select next parameter set
      x_new <- DDS.sel(i=cyclecount, m=m, r=r, xnames=xnames, x_min=x_min, x_max=x_max, x_best=x_best)
      cyclecount <- cyclecount+1  

      # Output next parameter set
      x_new_out <- c(cyclecount, x_new)
      names(x_new_out)[1] <- "id"
      write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")
   }

   # Stop cluster
   if (parallelFlag) stopCluster(cl)

   # Update plot
   gg <- ggplot(data=x_archive, aes(x=iter, y=obj)) + 
              geom_point() + theme_bw() + 
              labs(x="run", y="objective function")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_calib_run_obj.png"),
              plot=gg, units="in", width=6, height=5, dpi=300)

   # Archive output
   if (!archiveOutput) system(paste0("rm -r ", outPath), intern=FALSE)

   # Archive model run dir
   modFromPath <- paste0(runDir, "/RUN.CALTMP")
   modToPath <- paste0(runDir, "/CALIB_RUNS/RUN.CALTMP", cyclecount)
   if (archiveRun) {
      system(paste0("mv ", modFromPath, " ", modToPath), intern=FALSE)
   } else {
      system(paste0("rm -r ", modFromPath), intern=FALSE)
   }

   # Save and exit
   save.image("proj_data.Rdata")
   quit("no")

}



