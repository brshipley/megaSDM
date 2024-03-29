#' Construct ensemble models and project habitat suitability to current, past, and future climates
#'
#' This function conducts ensemble modelling on all replicates of the maxent model (see \code{MaxEntModel}
#' for information on the different replicates the function provides), by calculating the
#' median habitat suitability for each pixel across all replicates. Next, the function generates
#' binary presence/absence maps by applying a given threshold to the data. These processes are
#' repeated for each scenario/time period combination provided, to predict past or future
#' species distributions.
#'
#' @param input the full path name of the directory holding all model runs and outputs.
#' Outputs should be divided into sub-directories based on the species, named with the species
#' name (e.g. "./Results/Canis_lupus/"). The maxent.jar file should also be in this
#' directory (gets copied over if \code{MaxEntModel} is run).
#' @param time_periods a vector of the years in which the projection will occur. The first
#' element of this vector should be the original year (the year in which the model was
#' generated). If no precise years are available (e.g., using data from the Last Glacial Maximum),
#' order from current to least current (farthest into the future/past) and give character
#' strings for the years (e.g., "LGM").
#' @param scenarios a vector of character strings detailing the different climate models
#' used in the forecasted/hindcasted species distribution models. If only one distinct
#' climate scenario is used, still provide a name for that scenario. If no projection is
#' needed, set to NA (default).
#' @param study_dir The directory where all of the current/modern study area environmental rasters
#' are (should be the the same time period that the model was trained on).
#' @param predict_dirs A list of vectors: Each vector should include the directories of
#' the environmental rasters for all __time periods__ in a single __climate scenario__.
#' For example:
#' \code{list(Scenario1 = c(Time1, Time2), Scenario2 = c(Time1, Time2))}.
#'
#' If there is only one climate scenario, still create the list before the vector:
#' \code{list(Scenario1 = c(Time1,Time2))}.
#' @param ThreshMethod Which threshold to use for binary mapping; see MAXENT manual for details.
#'  Default is "Maximum.test.sensitivity.plus.specificity".
#' @param aucval (numeric or vector) Minimum AUC value necessary for each run to be counted.
#' AUC values estimate the predictive ability of a model, usually ranging from 0.5 (random)
#' to 1 (perfect). Can be a single number (same AUC threshold applied to all species) or
#' a vector with the same length as \code{input}. If set to \code{NA} (default), all replicates
#' will be used for ensembling and projection.
#' @param output The directory name where the outputs of the ensemble modelling and prediction
#' will be placed.
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @export
#' @return Returns binary (presence/absence) and ensembled distribution maps for all given species,
#' time periods, and climate scenarios provided. Each projected replicate of the maxent model is
#' placed in a newly-created folder (provided by the \code{projections} argument) within the directory given by the \code{input}
#' argument. The ensembled and binary maps are placed in the directory given by the \code{output} argument.

MaxEntProj <- function(input, time_periods, scenarios = NA, study_dir, predict_dirs,
                        output, ThreshMethod = "Maximum.test.sensitivity.plus.specificity",
                        aucval = NA, ncores = 1) {

  if(!dir.exists(output)) {
    dir.create(output)
  }

  studylayers <- list.files(study_dir, pattern = ".grd$", full.names = TRUE)
  #Ensure that all study area rasters have the same projection and extent
  projstudy <- rep(NA, len = length(studylayers))
  extstudy <- rep(NA, len = length(studylayers))
  for (i in 1:length(studylayers)) {
    projstudy[i] <- as.character(terra::crs(terra::rast(studylayers[[i]])))
    extstudy[i] <- as.character(terra::ext(terra::rast(studylayers[[i]])))
  }

  if (length(unique(projstudy)) > 1) {
    stop("Not all of the study area environmental rasters are in the same projection")
  } else if (length(unique(extstudy)) > 1) {
    stop("Not all of the study area environmental rasters have the same extent")
  }

  #Make sure that the study layers have a CRS that is not NA
  studystack <- terra::rast(studylayers)
  if (is.na(terra::crs(studystack))) {
    stop("study area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  } else {
    desiredCRS <- terra::crs(studystack)
  }
  
  #Format the matrix from which to parallelize the species
  ListSpp <- list.dirs(input, full.names = FALSE, recursive = FALSE)

  if (length(which(stats::complete.cases(aucval))) > 0) {
    #Removes species with all replicates less than the AUC threshold given by aucval
    print(paste0("    Removing species with Test AUC Values < AUC threshold from subsequent analyses"))
    AUCRetain <- c()
    if (length(aucval) == 1) {
      aucval <- rep(aucval, length(ListSpp))
    }
    for(i in 1:length(ListSpp)) {
      results <- utils::read.csv(paste0(input, "/", ListSpp[i], "/maxentResults.csv"))
      auc <- results$Test.AUC
      if (!is.na(aucval[i])) {
        aucAbove <- which(auc > aucval[i])
      } else {
        aucAbove <- 1:length(auc)
      }
      if (length(aucAbove) > 0){
        AUCRetain <- c(AUCRetain, i)
      } else {
        message(paste0("Removing: ", ListSpp[i]))
      }
    }
    ListSpp <- ListSpp[c(AUCRetain)]
  } else {
    ListSpp <- ListSpp
    results <- utils::read.csv(paste0(input, "/", ListSpp[1], "/maxentResults.csv"))
  }

  if (nrow(results) > 1) {
    nrep <- nrow(results) - 1
  } else {
    nrep <- 1
  }
  
  if (length(ListSpp) == 0) {
    stop(paste0("No species have AUC values higher than ", aucval))
  }
  if (length(ListSpp) < ncores) {
    ncores <- length(ListSpp)
  }
  ListSpp <- matrix(data = ListSpp, ncol = ncores)
  #Gets the date of the current year, the number of years and scenarios, and the total number of
  #distinct projections
  currentYear <- time_periods[1]
  numYear <- length(time_periods)
  if (!is.na(scenarios)[1]) {
    numScenario <- length(scenarios)
  } else {
    numScenario <- 0
  }
  nproj <- 1 + (numScenario * (numYear - 1))
  #Renames threshold to add in "logistic.threshold" to the mix
  ThreshMethodLower <- tolower(ThreshMethod)
  if(!grepl(".logistic.threshold", ThreshMethodLower)) {
    ThreshMethod <- paste0(ThreshMethod, ".Logistic.threshold")
  } else {
    ThreshMethod <- paste0(substr(ThreshMethod, 1, nchar(ThreshMethod) - 19), ".Logistic.threshold")
  }

  #Functions----------------------------
  #Creates binary rasters out of ensembled MaxEnt outputs
  threshold <- function(path, rasters, replicates, Scenario, time) {
    rasters <- matrix(rasters, nrow = replicates)
    rasterNames <-  c()

    ensemble.stack <- c()
    curmodel <- file.path(input, path)
    
    #Find which AUC value corresponds to the species of interest
    aucvalSpecies <- aucval[which(ListSpp == path)]
    
    #If species have been removed from the species list (e.g., not enough high-quality replicates),	
    #only take the first auc value provided	
    aucvalSpecies <- aucvalSpecies[1]	
    
          
    #Reads in the results file for each of the runs
    results <- utils::read.csv(paste0(curmodel, "/maxentResults.csv"))
    for (j in 1:replicates) {
      #Reads the threshold value for that results file
      if (!is.na(as.numeric(ThreshMethod))) {
        thresh <- as.numeric(ThreshMethod)
      } else {
        thresh <- results[j, ThreshMethod]
      }

      if(is.null(thresh)) {
        stop("The threshold for creating binary maps is not found in 'maxentResults.csv'. Revise in config")
      }

      if (!is.na(aucvalSpecies)) {
        tryCatch({
          #Finds the AUC value of the run
          auc <- results$Test.AUC[j]
          #If the AUC value is below the AUC threshold, the run is not used
          if (auc >= aucvalSpecies) {
            temp <- terra::rast(rasters[j])
            temp[temp >= thresh] <- 1
            temp[temp < thresh] <- 0
            rasterNames <- c(rasterNames, names(temp))
            #Creates a stack of output runs to be ensembled
            ensemble.stack <- c(ensemble.stack, temp)
          }
        }, error = function(err) {
          print(paste("MY_ERROR: ", path, " ", err," j= ", j,
                      " Scenario= ", Scenario, " time= ", time))
        })
      } else {
        temp <- terra::rast(rasters[j])
        temp[temp >= thresh] <- 1
        temp[temp < thresh] <- 0
        rasterNames <- c(rasterNames, names(temp))
        #Creates a stack of output runs to be ensembled
        ensemble.stack <- c(ensemble.stack, temp)
      }
    }
    #combine list of rasters into a single object
    ensemble.stack <- terra::rast(ensemble.stack)
    
    if (length(ensemble.stack) > 0) {
      #Takes the mean of the binary rasters
      #If the mean is >= 0.5 (more than half of the replicates show presence), set to 1
      ensemble.calc <- (terra::mean(ensemble.stack))
      rasterNames <- c(rasterNames, names(ensemble.calc))
      ensemble.calc[ensemble.calc >= 0.5] <- 1
      ensemble.calc[ensemble.calc < 0.5] <- 0
      terra::crs(ensemble.calc) <- desiredCRS
      #Writes out binary rasters

      if (!dir.exists(file.path(output, path, "map_pdfs"))) {
        dir.create(file.path(output, path, "map_pdfs"))
      }

      if (time == currentYear) {
        terra::crs(ensemble.calc) <- desiredCRS
        terra::writeRaster(ensemble.calc,
                    file.path(output, path, paste0(time, "_binary.grd")),
                    overwrite = TRUE)
        grDevices::pdf(file = file.path(output, path, "map_pdfs", paste0(time, "_binary.pdf")))
        terra::plot(ensemble.calc, main = paste0(path, "_", time, "_binary"))
        grDevices::dev.off()
      } else {
        if (!dir.exists(file.path(output, path, Scenario))) {
          dir.create(file.path(output, path, Scenario))
        }
        terra::crs(ensemble.calc) <- desiredCRS
        terra::writeRaster(ensemble.calc,
                    file.path(output, path, Scenario, paste0(time, "_", Scenario, "_binary.grd")),
                    overwrite = TRUE)
        grDevices::pdf(file = file.path(output, path, "map_pdfs", paste0(time, "_", Scenario, "_binary.pdf")))
        terra::plot(ensemble.calc, main = paste0(path, "_", time, "_", Scenario, "_binary"))
        grDevices::dev.off()
      }
      return(ensemble.calc)
    } else {
      message(paste0("No replicates of ", path, " had a test AUC value above ", aucval))
    }
    rm(temp, ensemble.stack)
    gc()
  }

  medianensemble <- function(path, rasters, replicates, Scenario, decade) {
    rasters <- matrix(rasters, nrow = replicates)
    #Lists all of the folders within "outputs"
    message("starting median ensemble")
    curmodel <- paste0(input, "/", path)
    ensemble.stack <- c()
    #Reads in the results file for each of the runs
    results <- utils::read.csv(paste0(curmodel, "/maxentResults.csv"))
    
    #Find which AUC value corresponds to the species of interest
    aucvalSpecies <- aucval[which(ListSpp == path)]
    
    #If species have been removed from the species list (e.g., not enough high-quality replicates),	
    #only take the first auc value provided	
    aucvalSpecies <- aucvalSpecies[1]	
    
    for (j in 1:replicates) {
      if (!is.na(aucvalSpecies)) {
        tryCatch({
          #Finds the AUC value of the run
          auc <- results$Test.AUC[j]
          #If the AUC value is below the AUC threshold, the run is not used
          if (auc >= aucvalSpecies) {
            temp <- terra::rast(rasters[j])
            #Creates a stack of output runs to be ensembled
            ensemble.stack <- c(ensemble.stack, temp)
          }
        }, error = function(err) {
          print(paste("MY_ERROR: ", path, " ", err, "j=", j,
                      " Scenario= ", Scenario, " Decade= ", decade))
        })
      } else {
        temp <- terra::rast(rasters[j])
        #Creates a stack of output runs to be ensembled
        ensemble.stack <- c(ensemble.stack, temp)
      }
    }
    
    #Combine list of rasters into single object
    ensemble.stack <- terra::rast(ensemble.stack)
    
    #Calculates the median of the stacked rasters
    if (terra::nlyr(ensemble.stack) > 1) {
      ensemble.calc <- terra::app(ensemble.stack, stats::median, na.rm = TRUE)
    } else {
      ensemble.calc <- ensemble.stack
    }
    
    if (decade == currentYear) {
      terra::crs(ensemble.calc) <- desiredCRS
      terra::writeRaster(ensemble.calc,
                  file.path(output, path, paste0(decade, "_ensembled.grd")),
                  overwrite = TRUE)
      grDevices::pdf(file = file.path(output, path, "map_pdfs", paste0(decade, "_ensembled.pdf")))
      terra::plot(ensemble.calc, main = paste0(path, "_", decade))
      grDevices::dev.off()
    } else {
      terra::crs(ensemble.calc) <- desiredCRS
      terra::writeRaster(ensemble.calc,
                  file.path(output, path, Scenario, paste0(decade, "_", Scenario, "_ensembled.grd")),
                  overwrite = TRUE)
      grDevices::pdf(file = file.path(output, path, "map_pdfs", paste0(decade, "_", Scenario, "_ensembled.pdf")))
      terra::plot(ensemble.calc, main = paste0(path, "_", decade, "_", Scenario))
      grDevices::dev.off()
    }
    rm(ensemble.stack)
    gc()
    return(ensemble.calc)
  }

  getSize <- function(raster) {
    Freq <- data.frame(terra::freq(raster, digits = 0, value = 1))
    return(Freq$count)
  }

  t1nott2 <- function(t1, t2) {
    tnew <- t1
    tnew[tnew == 1] <- 2
    Combined <- tnew + t2
    Combined[Combined != 2] <- 0
    Combined[Combined == 2] <- 1
    return(Combined)
  }

  overlap <- function(t1, t2) {
    Combined <- t1 + t2
    Combined[Combined != 2] <- 0
    Combined[Combined == 2] <- 1
    return(Combined)
  }

  getCentroid <- function(CentRaster) {
    #convert raster to points and only take the presence points
    points <- terra::as.points(CentRaster, values = TRUE)
    points <- points[which(terra::values(points) == 1),]
    
    #Get the coordinates of the points
    points <- data.frame(terra::geom(points))[, c("x", "y")]
    
    #average latitude (y)
    Clat <- mean(points[, 2], na.rm = TRUE)

    #average longitude (x)
    Clong <- mean(points[, 1], na.rm = TRUE)

    #returns the longitude & latitude of the Centroid
    return(c(Clong, Clat))
  }

  getStats <- function(cur, j, decade, binary, stats, modern.size, modern.binary) {
    #Adds all of the variables into the stats data frame
    stats$Projection[j] <- cur
    stats$NumberCells[j] <- getSize(binary)
    stats$CellChange[j] <- getSize(binary) - modern.size
    t1nott2_raster <- t1nott2(modern.binary, binary)
    t2nott1_raster <- t1nott2(binary, modern.binary)
    stats$T1notT2[j] <- getSize(t1nott2_raster)
    stats$T2notT1[j] <- getSize(t2nott1_raster)
    overlap_raster <- overlap(modern.binary, binary)
    stats$Overlap[j] <- getSize(overlap_raster)
    binary.centroid <- getCentroid(binary)
    stats$CentroidX[j] <- binary.centroid[1]
    stats$CentroidY[j] <- binary.centroid[2]
    return(stats)
  }

  run <- function(CurSpp) {
    spp.name <- CurSpp
    if(!dir.exists(file.path(output, spp.name))) {
      dir.create(file.path(output, spp.name))
    }

    #Creates a stats table to be filled in & passed to next scripts
    stats <- as.data.frame(cbind(Projection = rep(0, times = nproj),
                                 NumberCells = rep(0, times = nproj),
                                 CellChange = rep(0, times = nproj),
                                 T1notT2 = rep(0, times = nproj),
                                 T2notT1 = rep(0, nproj),
                                 Overlap = rep(0, times = nproj),
                                 CentroidX = rep(0, times = nproj),
                                 CentroidY = rep(0, times = nproj)))

    #Gets modern data stats first to use for comparison
    dir.create(file.path(output, spp.name, "projections"))
    dir.create(file.path(output, spp.name, "projections", time_periods[1]))
    for (g in 1:nrep) {
      #If there is more than 1 replicate, the lambdas file is named differently
      if (nrep == 1) {
        LambdaFile <- file.path(input, spp.name, paste0(spp.name, ".lambdas"))
        nrep_new <- 1
      } else {
        LambdaFile <- file.path(input, spp.name, paste0(spp.name, "_", g-1, ".lambdas"))
        if (!file.exists(LambdaFile)) {
          message(paste0("Model replicate ", g, " of species ", spp.name, " not found!"))
          message(paste0("nrep may be greater than number of unique occurrences."))
          message(paste0("Using fewer replicates for this species."))
          nrep_new <- g-1
          break()
        } else {
          nrep_new <- g
        }
      }
      #THIS IS THE PROJECTIONS COMMAND>>
      system(paste0("java -mx900m -cp ", file.path(input, "maxent.jar"), " density.Project ",
                    # location of lambdas file
                    file.path(LambdaFile), " ",
                    #location of folder with map rasters
                    file.path(study_dir), " ",
                    # where to output files (Run#)
                    file.path(output, spp.name, "projections", time_periods[1], paste0("RUN_", g-1)),
                    " noaskoverwrite nowarnings -a"))
    }

    modern.rasters <- c()
    r <- list.files(path = file.path(output, spp.name, "projections", time_periods[1]),
                    pattern = paste0(".asc$"),
                    full.names = TRUE)

    for (k in 1:length(r)) {
      if (!grepl("clamping", r[k])) {
        modern.rasters <- c(modern.rasters, r[k])
      }
    }
    modern.rasters <- gtools::mixedsort(modern.rasters)

    #Runs "threshold" function, which creates and ensembles binary maps
    modern.binary <- threshold(spp.name, modern.rasters, nrep_new, "", currentYear)
    gc()

    if (length(modern.binary) > 0) {
      #Runs medianensemble function, which ensembles maps using a median function
      modern.median <- medianensemble(spp.name, modern.rasters, nrep_new, "", currentYear)
      gc()

      statsrow <- 1
      #Fills in the stats row for the modern species
      stats$Projection[statsrow] <- currentYear
      modern.size <- getSize(modern.binary)
      stats$NumberCells[statsrow] <- modern.size
      modern.centroid <- getCentroid(modern.binary)
      stats$CentroidX[statsrow] <- modern.centroid[1]
      stats$CentroidY[statsrow] <- modern.centroid[2]

      #Generates projections for hindcasted/forecasted climate layers
      if (numScenario > 0) {
        for (ScenIndex in 1:numScenario) {
          dir.create(file.path(output, spp.name, "projections", scenarios[ScenIndex]))
          for (YearIndex in 2:numYear) {
            dir.create(file.path(output, spp.name, "projections", scenarios[ScenIndex], time_periods[YearIndex]))
          }
        }

        for (ScenIndex in 1:numScenario) {
          focusScen <- scenarios[ScenIndex]
          for (YearIndex in 2:numYear) {
            statsrow <- statsrow + 1
            focusDate <- time_periods[YearIndex]
            rasterLocation <- predict_dirs[[ScenIndex]][YearIndex - 1]

            for (g in 1:nrep_new) {
              #If there is more than 1 replicate, the lambdas file is named differently
              if (nrep == 1) {
                LambdaFile <- file.path(input, spp.name, paste0(spp.name, ".lambdas"))
              } else {
                LambdaFile <- file.path(input, spp.name, paste0(spp.name, "_", g-1, ".lambdas"))
              }
              #THIS IS THE PROJECTIONS COMMAND>>
              system(paste0("java -mx900m -cp ", file.path(input, "maxent.jar"), " density.Project ",
                            # location of lambdas file
                            LambdaFile, " ",
                            #location of folder with map rasters
                            rasterLocation, " ",
                            # where to output files (Run#)
                            output, "/", spp.name, "/projections/", focusScen, "/",
                            focusDate, "/RUN_", g-1,
                            " noaskoverwrite nowarnings -a"))
            }

            #Makes a list of the projections files to use for threshold/ensemble functions
            cur.rasters <- list.files(path = file.path(output, spp.name, "projections",
                                                    focusScen, focusDate),
                                      pattern = paste0('\\.asc$'),
                                      full.names = TRUE)
            cur.proj <- c()
            for (k in 1:length(cur.rasters)) {
              if (!grepl("clamping", cur.rasters[k])) {
                cur.proj <- c(cur.proj, cur.rasters[k])
              }
            }
            length(cur.proj)
            #Thresholds & gets medians of the projected rasters
            cur.binary <- threshold(spp.name, cur.proj, nrep_new, focusScen, focusDate)
            gc()
            cur.median <- medianensemble(spp.name, cur.proj, nrep_new, focusScen, focusDate)
            gc()
            #Fills in the stats table for the projected rasters
            stats <- getStats(paste0(focusScen, "_", focusDate), statsrow, focusDate, cur.binary,
                              stats, modern.size, modern.binary)
            gc()
          }
        }
      }
      #Writes the stats table to be used later
      utils::write.csv(stats, file = file.path(output, spp.name, "Results.csv"))
    } else {
      message(paste0("Removing ", spp.name, " from further analysis"))
    }

  }
  
  if (ncores == 1) {
    ListSpp <- as.vector(ListSpp)
    out <- sapply(ListSpp, function(x) run(x))
  } else {
    clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)
    
    parallel::clusterExport(clus, varlist = c("run", "threshold", "medianensemble", "getSize",
                                              "t1nott2", "overlap", "getCentroid","getStats", "nrep",
                                              "currentYear", "numYear", "numScenario", "nproj",
                                              "input", "time_periods", "scenarios", "predict_dirs",
                                              "study_dir", "ThreshMethod", "aucval", "output", "ncores",
                                              "ListSpp", "desiredCRS"), envir = environment())
    
    parallel::clusterEvalQ(clus, library(gtools))
    parallel::clusterEvalQ(clus, library(terra))
    for (i in 1:nrow(ListSpp)) {
      out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
      gc()
    }
    parallel::stopCluster(clus)
  }
}
