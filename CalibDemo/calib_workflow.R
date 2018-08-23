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
   lastcycle <- FALSE

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
  # Extra (Time_CST), change only for this particular case where daily is an average from 8:00 to 8:00 Local time
 frxst$Time_CST <- as.POSIXct(format(frxst$POSIXct, tz="America/Mexico_City",usetz=TRUE))

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
      iter_best <- cyclecount
   } else if (F_new <= F_best) {
      x_best <- x_new
      F_best <- F_new
     iter_best <- cyclecount
   }

   if (cyclecount < m) {
      # Select next parameter set
      x_new <- DDS.sel(i=cyclecount, m=m, r=r, xnames=xnames, x_min=x_min, x_max=x_max, x_best=x_best)
      cyclecount <- cyclecount+1  

      # Output next parameter set
      x_new_out <- c(cyclecount, x_new)
      names(x_new_out)[1] <- "id"
      write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")
   } else {
      lastcycle <- TRUE
   }

   # Stop cluster
   if (parallelFlag) stopCluster(cl)


#########################################################
# PLOTS
#########################################################
# First we check if all the objective function values are less than the threshold (here 5), define it as no outlier in the iterations
# If there are objFun values greater than the threshold in the objFun, then calulate the 90% of the objFun
# Any iteration with objFun values above the 90% would be flagged as outlier. And then two plots will be created
# one with all iteration including the outliers, two only 90% of the data if there was an outlier in the model.

objFunThreshold <- 5
objFunQuantile <- quantile(x_archive$obj, 0.9)

if (any(x_archive$obj > objFunThreshold)) {
   write("Outliers found!", stdout())

   # Check which outlier threshold to use
   if (any(x_archive$obj <= objFunThreshold)) {
     x_archive_plot <- subset(x_archive, x_archive$obj <= objFunThreshold)
     x_archive_plot_count <- nrow(x_archive) - nrow(x_archive_plot)
     x_archive_plot_threshold <- objFunThreshold
   } else {
     x_archive_plot <- subset(x_archive, x_archive$obj <= objFunQuantile)
     x_archive_plot_count <- nrow(x_archive) - nrow(x_archive_plot)
     x_archive_plot_threshold <- objFunQuantile
   }

   if (!exists("x_archive_plot_count_track")) x_archive_plot_count_track <- data.frame()
   x_archive_plot_count_track <- rbind(x_archive_plot_count_track, data.frame(iter=ifelse(lastcycle, cyclecount, cyclecount-1), outliers=nrow(x_archive)-nrow(x_archive_plot)))

   # Outlier count
   if (nrow(x_archive_plot_count_track) > 0) {
       write("Outlier count plot...", stdout())
       gg <- ggplot(data=x_archive_plot_count_track, aes(x=iter, y=outliers)) +
            geom_point() + theme_bw() +
            labs(x="run", y="count of outlier cycles")
       ggsave(filename=paste0(writePlotDir, "/", siteId, "_calib_outliers.png"),
            plot=gg, units="in", width=6, height=5, dpi=300)
   }

} else {
  write("No outliers found.", stdout())
  # All the objFun vlaues are less than the threshold defined above, therefore, there will not be any outliers specified
   x_archive_plot <- x_archive
   x_archive_plot_count <- 0
   x_archive_plot_threshold <- objFunThreshold
}


#**************************************************************************************************************************************
#                                   Create the plots with outlier
#**************************************************************************************************************************************

   # Update basic objective function plot
   write("Basin objective function plot...", stdout())
   gg <- ggplot(data=x_archive, aes(x=iter, y=obj)) +
              geom_point() + theme_bw() +
              labs(x="run", y="objective function")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_calib_run_obj_outlier.png"),
              plot=gg, units="in", width=6, height=5, dpi=300)

   # Update the Objective function versus the parameter variable
   write("Obj function vs. params...", stdout())
   DT.m1 = melt(x_archive[, setdiff(names(x_archive), metrics)], id.vars = c("obj"), measure.vars = setdiff( names(x_archive), c(metrics, "iter", "obj")))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(value, obj))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free_x")
   gg <- gg + ggplot2::ggtitle(paste0("Scatter Plot of Obj. function versus parameters: ", siteId))
   gg <- gg + ggplot2::xlab("Parameter Values")+theme_bw()+ggplot2::ylab("Objective Function")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_obj_vs_parameters_calib_run_outlier.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)


   # Plot the variables as a function of calibration runs
   write("Params over runs...", stdout())
   DT.m1 = melt(x_archive[, setdiff(names(x_archive), metrics)], id.vars = c("iter"), measure.vars = setdiff(names(x_archive), c("iter", metrics)))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(iter, value))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free")
   gg <- gg + ggplot2::ggtitle(paste0("Parameter change with iteration: ", siteId))
   gg <- gg + ggplot2::xlab("Calibration Iteration")+theme_bw()
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_parameters_calib_run_outlier.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)

   # Plot all the stats
   write("Metrics plot...", stdout())
   DT.m1 = melt(x_archive[,which(names(x_archive) %in% c("iter", "obj", "cor", "rmse", "bias", "nse", "nselog", "nsewt", "kge", "msof"))],
               iter.vars = c("iter"), measure.vars = c("obj", "cor", "rmse", "bias", "nse", "nselog", "nsewt", "kge", "msof"))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(iter, value))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free")
   gg <- gg + ggplot2::ggtitle(paste0("Metric Sensitivity: ", siteId))
   gg <- gg + ggplot2::xlab("Calibration Iteration No.")+theme_bw()+ylab("Value")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_metric_calib_run_outlier.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)

#############################################################################################################################################################################
#                      Create the plots without outliers
############################################################################################################################################################################3

  # Update basic objective function plot
   write("Basin objective function plot...", stdout())
   gg <- ggplot(data=x_archive_plot, aes(x=iter, y=obj)) +
              geom_point() + theme_bw() +
              labs(x="run", y="objective function") +
              ggtitle(paste0("ObjFun: ", siteId,  ", No. outliers = ", x_archive_plot_count, ", Threshold = ",  formatC(x_archive_plot_threshold, digits  = 4)))

   ggsave(filename=paste0(writePlotDir, "/", siteId, "_calib_run_obj.png"),
              plot=gg, units="in", width=6, height=5, dpi=300)

   # Update the Objective function versus the parameter variable
   write("Obj function vs. params...", stdout())
   DT.m1 = melt(x_archive_plot[, setdiff(names(x_archive_plot), metrics)], id.vars = c("obj"), measure.vars = setdiff( names(x_archive_plot), c(metrics, "iter", "obj")))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(value, obj))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free_x")
   gg <- gg + ggplot2::ggtitle(paste0("ObjFun vs. Params: ", siteId,  ", No. outliers = ", x_archive_plot_count, ", Threshold = ",  formatC(x_archive_plot_threshold, digits  = 4)))
   gg <- gg + ggplot2::xlab("Parameter Values")+theme_bw()+ggplot2::ylab("Objective Function")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_obj_vs_parameters_calib_run.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)


   # Plot the variables as a function of calibration runs
   write("Params over runs...", stdout())
   DT.m1 = melt(x_archive_plot[, setdiff(names(x_archive_plot), metrics)], id.vars = c("iter"), measure.vars = setdiff(names(x_archive_plot), c("iter", metrics)))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(iter, value))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free")
   gg <- gg + ggplot2::ggtitle(paste0("Parameter vs. iteration: ", siteId,  ", No. outliers = ", x_archive_plot_count, ", Threshold = ",  formatC(x_archive_plot_threshold, digits  = 4)))
   gg <- gg + ggplot2::xlab("Calibration Iteration")+theme_bw()
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_parameters_calib_run.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)

   # Plot all the stats
   write("Metrics plot...", stdout())
   DT.m1 = melt(x_archive_plot[,which(names(x_archive_plot) %in% c("iter", "obj", "cor", "rmse", "bias", "nse", "nselog", "nsewt", "kge", "msof"))],
               iter.vars = c("iter"), measure.vars = c("obj", "cor", "rmse", "bias", "nse", "nselog", "nsewt", "kge", "msof"))
   DT.m1 <- subset(DT.m1, !is.na(DT.m1$value))
   gg <- ggplot2::ggplot(DT.m1, ggplot2::aes(iter, value))
   gg <- gg + ggplot2::geom_point(size = 1, color = "red", alpha = 0.3)+facet_wrap(~variable, scales="free")
   gg <- gg + ggplot2::ggtitle(paste0("Metric Sensitivity: ", siteId, ", No. outliers = ", x_archive_plot_count, ", Threshold = ",  formatC(x_archive_plot_threshold, digits  = 4)))
   gg <- gg + ggplot2::xlab("Calibration Iteration No.")+theme_bw()+ylab("Value")
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_metric_calib_run.png"),
         plot=gg, units="in", width=8, height=6, dpi=300)


   # Plot the time series of the observed, control, best calibration result and last calibration iteration
   write("Hydrograph...", stdout())
   # The first iteration is the control run  called chrt.d.1
   controlRun <- copy(chrt.d.1)
   controlRun [, run := "Control Run"]
   # We have already advanced the cyclescount, so subtract 1 to get last complete
   lastRun <- copy(get(paste0("chrt.d.", ifelse(lastcycle, cyclecount, cyclecount-1))))
   lastRun [ , run := "Last Run"]
   # the best iteration should be find
   bestRun <- copy(get(paste0("chrt.d.", iter_best)))
   bestRun [ , run := "Best Run"]

   obsStrDataPlot <- copy(obsDT)
   setnames(obsStrDataPlot, "obs", "q_cms")
   obsStrDataPlot <- obsStrDataPlot[, c("Date", "q_cms", "POSIXct", "site_no"), with=FALSE]
   obsStrDataPlot$POSIXct <- as.POSIXct(paste0(as.character(obsStrDataPlot$POSIXct), "_08"), format= "%Y-%m-%d_%H", tz = "UTC")

   obsStrDataPlot <- obsStrDataPlot[as.integer(POSIXct) >= min(as.integer(controlRun$POSIXct)) & as.integer(POSIXct) <= max(as.integer(controlRun$POSIXct)),]
   obsStrDataPlot[ , run := "Observation"]

   chrt.d_plot <- rbindlist(list(controlRun, lastRun, bestRun, obsStrDataPlot), use.names = TRUE, fill=TRUE)

   gg <- ggplot2::ggplot(chrt.d_plot, ggplot2::aes(POSIXct, q_cms, color = run))
   gg <- gg + ggplot2::geom_line(size = 0.3, alpha = 0.7)
   gg <- gg + ggplot2::ggtitle(paste0("Streamflow time series for ", siteId))
   #gg <- gg + scale_x_datetime(limits = c(as.POSIXct("2008-10-01"), as.POSIXct("2013-10-01")))
   gg <- gg + ggplot2::xlab("Date")+theme_bw( base_size = 15) + ylab ("Streamflow (cms)")
   gg <- gg + scale_color_manual(name="", values=c('black', 'dodgerblue', 'orange' , "dark green"),
                                 limits=c('Observation','Control Run', "Best Run", "Last Run"),
                                  label=c('Observation','Control Run', "Best Run", "Last Run"))

   ggsave(filename=paste0(writePlotDir, "/", siteId, "_hydrograph.png"),
           plot=gg, units="in", width=8, height=4, dpi=300)


# Plot the scatter plot of the best, last and control run.
   write("Scatterplot...", stdout())
   maxval <- max(max(chrt.d_plot$q_cms, na.rm = TRUE), max(obsDT$obs, na.rm=TRUE))
   gg <- ggplot()+ geom_point(data = merge(chrt.d_plot [run %in% c("Control Run", "Last Run", "Best Run")], obsDT, by.x = c("site_no", "Date"), by.y=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE),
                              aes (obs, q_cms, color = run), alpha = 0.5)
   gg <- gg + scale_color_manual(name="", values=c('dodgerblue', 'orange' , "dark green"),
                                 limits=c('Control Run', "Best Run", "Last Run"),
                                 label=c('Control Run', "Best Run", "Last Run"))
   gg <- gg + ggtitle(paste0("Simulated vs observed flow : ", siteId )) + theme_bw( base_size = 15)
   gg <- gg + geom_abline(intercept = 0, slope = 1) + coord_equal()+ xlim(0,maxval) + ylim(0,maxval)
   gg <- gg + xlab("Observed flow (cms)") + ylab ("Simulated flow (cms)")

   ggsave(filename=paste0(writePlotDir, "/", siteId, "_scatter.png"),
           plot=gg, units="in", width=8, height=8, dpi=300)

   rm(controlRun, lastRun, bestRun, obsStrDataPlot)


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



