#' Constrain Modelled Species Distributions by Dispersal Ability
#'
#' This function incorporates the ability of a species to disperse over time into
#' projected habitat suitability models and presence/absence maps. The probability
#' of dispersal per year as a function of distance is modelled using an exponential
#' distribution, and summed together to create a probability of dispersal for the
#' intervals between each provided time step. Dispersal-constrained binary (presence
#' and absence) maps are generated, as well as continuous maps of "invadable suitability"
#' (see https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05450).
#'
#' NOTE: dispersal rate analyses are only informative for predicting species
#' distributions into the future (forecasting) rather than predicting past
#' distributions (hindcasting), as range contractions and extirpations are not
#' limited by dispersal rate.
#'
#' NOTE: Running this function for data in a longlat projection will take somewhat
#' longer than using data in equal area projections. Results are the same, however.
#'
#' @param result_dir the directory where the ensembled and binary maps are placed.
#' Each species should have its own sub-directory, and the forecasted/hindcasted
#' binary maps should be placed into directories like so: Species/Scenario/Time.
#' If \code{MaxEntProj} was used to make these maps, this is probably the same
#' as the \code{output} argument in that function.
#' @param dispersaldata either a dataframe or the complete path name of a
#' .csv file with two columns:
#'
#'   Column 1: species name (same as the name used for modelling).
#'
#'   Column 2: average dispersal rate of species in kilometers/year.
#' @param time_periods a vector of the years in which the projection will occur. The first
#' element should be the original year (the year in which the model was generated).
#' @param scenarios a vector of character strings detailing the different climate models
#' used in the forecasted/hindcasted species distribution models.
#' @param contiguous TRUE/FALSE: when constraining by dispersal rate, should only the contiguous
#' areas around each occurrence point be used for the first time step? setting this argument to
#' \code{TRUE} will mitigate the effects of overprediction in areas where a species has
#' not been observed. Default is \code{FALSE}.
#' @param point_data If contiguous = \code{TRUE}, a file path to the directory holding all
#' occurrence data used for the model generation. If previous steps of megaSDM have been run
#' (e.g., \code{MaxEntProj}, \code{MaxEntModel}, \code{OccurrenceManagement}), this will likely
#' be equal to \code{occ_output} in the vignette. If not, the occurrence data should have
#' coordinates in "x" and "y" columns.
#' @param hindcast TRUE/FALSE: is this a hindcasted model? If TRUE, dispersal rate calculations
#' will start at the first time period, not the current one.
#' @param startpoint if hindcast is \code{TRUE} then startpoint is a vector of length 2
#' describing a scenario/time combination, with the first argument as the scenario name
#' and the second as the time period desired as a starting point for the dispersal
#' simulations.
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @export
#' @return returns rasters and .pdf files of the projected species ranges for future time
#' periods, given different climate scenarios, when the ability of the species to disperse
#' is taken into account. Dispersal-constrained presence/absence maps and "invadable suitability"
#' maps are provided as outputs.

dispersalRate <- function(result_dir, dispersaldata, time_periods,
                          scenarios, contiguous = FALSE,
                          point_data = NA, hindcast = FALSE, 
                          startpoint = c(NA, NA), ncores = 1) {

  #Gets list of species from the directories given by "result_dir"
  spp.list <- list.dirs(result_dir, full.names = FALSE, recursive = FALSE)

  if (length(spp.list) == 0) {
    stop(paste0("No projected models found in 'result_dir': Ensure that 'result_dir' provides a path to the proper location"))
  }
  #Calculates number of scenarios, future time periods
  numScenario <- length(scenarios)
  numYear <- length(time_periods)

  #Reads in dispersal data
  if (class(dispersaldata) == "character") {
    dispersal <- utils::read.csv(dispersaldata, stringsAsFactors = FALSE)
    dispersal[, 1] <- gsub("_", " ", dispersal[, 1])
  } else {
    dispersal <- dispersaldata
    dispersal[, 1] <- gsub("_", " ", dispersal[, 1])
  }


  ListSpp <- c()
  #If no dispersal rate data found for a given species, prints a message
  for (w in 1:length(spp.list)) {
    FocSpec <- gsub("_", " ", spp.list[w])
    DispSpec <- grep(paste0("^", FocSpec, "$"), dispersal[, 1])
    if (length(DispSpec) == 0) {
      message(paste0("No dispersal rate values found for ", FocSpec, ": skipping dispersal rate analysis"))
    } else {
      ListSpp <- c(ListSpp, spp.list[w])
    }
  }
  if (length(ListSpp) < ncores) {
    ncores <- length(ListSpp)
  }
  ListSpp <- matrix(ListSpp, ncol = ncores)

  #Functions-------------------

  getSize <- function(raster) {
    Freq <- data.frame(terra::freq(raster, digits = 0, value = 1))
    return(Freq$count)
  }

  getCentroid <- function(raster) {
    #convert raster to points and only take the presence points
    points <- terra::as.points(raster, values = TRUE)
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


  #Creates distance rasters from original projection (studyarea environmental rasters)
  DistanceRaster <- function(spp, Time, Scen, CurrentBinary, TimeMap) {
    #Loads and reclassifies the binary maps
    CurrentPresence <- CurrentBinary
    CurrentPresence[which(terra::values(CurrentPresence) == 0)] <- NA
    #Trims the time map as an extent template for faster calculations
    TimeMap2 <- TimeMap
    TimeMap2[which(terra::values(TimeMap2) == 0)] <- NA
    TimeMap2 <- terra::trim(TimeMap2)
    TimeMap2[which(is.na(terra::values(TimeMap2)))] <- 0
    #Trims "CurrentPresence" to the extent of the time maps
    CurrentPresence <- terra::crop(CurrentPresence, terra::ext(TimeMap2))

    #Calculates the distances from each pixel to the nearest presence
    CurrDist <- terra::distance(CurrentPresence)
    #Extends the raster back out to full study area extent
    CurrDist <- terra::extend(CurrDist, terra::ext(CurrentBinary))
    CurrDist[which(is.na(terra::values(CurrDist)))] <- max(terra::values(CurrDist), na.rm = TRUE) + 1
    DistFinal <- terra::mask(CurrDist, CurrentBinary)

    #if the data have units that are not degrees or meters, convert distance values to meters.
    if(terra::linearUnits(DistFinal) != 0) {
      unitmult <- terra::linearUnits(DistFinal)
      DistFinal <- DistFinal * unitmult
    }

    #Converts distance (now in meters) to kilometers (for dispersal rate)
    DistFinal <- (DistFinal / 1000)
    #writes distance raster
    # terra::writeRaster(DistFinal,
    #             filename = paste0(result_dir, "/", spp, "/", "distance_", Time, "_", Scen, ".grd"),
    #             overwrite = TRUE)
    rm(CurrentPresence, CurrDist)
    gc()
    return(DistFinal)
  }

  #Creates dispersal probability raster from distance raster
  DispersalProbRaster <- function(rate, DistRaster, Elapsed) {
    #Calculates lambda of exponential distribution
    Lambda <- 1 / rate

    #When exponential distributions are added, they convolute to a gamma distribution
    GammaProbFun <- function(x) {
      1 - stats::pgamma(x, shape = Elapsed, rate = Lambda)
    }

    #Relates distance raster to dispersal probability
    DistProb <- terra::app(DistRaster, fun = function(x){GammaProbFun(x)})
    return(DistProb)
  }

  #Gets presence pixels in raster 1 but not raster 2
  t1nott2 <- function(t1, t2) {
    return(terra::mask(t1, t2, inverse = FALSE, maskvalue = 1, updatevalue = 0))
  }

  #Calculates overlap between raster 1 and raster 2 presences
  overlap <- function(t1, t2) {
    return(terra::mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
  }

  #Conducts the actual dispersal rate analyses
  FinalDispersal <- function(CurSpp) {
    #Gets species name and relevant dispersal rate
    speciesName = gsub("_", " ", CurSpp)
    if (length(grep(paste0(speciesName), dispersal[, 1])) > 0) {
      #Finds species-specific dispersal rate
      dispRateColumn <- which(unlist(lapply(dispersal, is.numeric)))
      dispersal_rate <- dispersal[grep(paste0(speciesName, "\\s*$"), dispersal[, 1]), dispRateColumn]
      if (length(dispersal_rate) > 1) {
        stop(paste0("More than one dispersal rate found for ", speciesName))
      }
      
      #Checking to see if the time periods pass through the year the model is projected on
      #or if they are historical time periods
      TrueCurr <- time_periods[1]
      if(hindcast == TRUE) {
        time_periods <- sort(time_periods)
      }
      
      if (!is.na(dispersal_rate)) {
        CurrentTime <- time_periods[1]

        #If desired, removes patches of the binary map at the first time step that do not border
        #the occurrence points used for the model. This is to avoid extrapolation errors.
        if (contiguous == TRUE) {

          #Read in the binary and ensembled rasters for the species at the time period the model was trained on
          CurrentBinary <- terra::rast(file.path(result_dir, CurSpp, paste0(CurrentTime, "_binary.grd")))
          CurrentEnsemble <- terra::rast(file.path(result_dir, CurSpp, paste0(CurrentTime, "_ensembled.grd")))

          CurrentPatch <- terra::patches(CurrentBinary, zeroAsNA = TRUE)

          #Get a list of occurrences and convert to SpatVector
          occurrence <- read.csv(file.path(point_data, paste0(CurSpp, ".csv")))
          occurrence <- terra::vect(occurrence, geom = c("x", "y"))
          terra::crs(occurrence) <- terra::crs(CurrentPatch)

          #Extract patch ID numbers that have points in them
          occ_extract <- terra::extract(CurrentPatch, occurrence, ID = FALSE)
          truepatches <- unique(occ_extract$patches)
          truepatches <- truepatches[which(!is.na(truepatches))]

          #Convert unoccupied patches to 0, occupied patches to 1
          CurrentPatch <- terra::match(CurrentPatch, truepatches)
          CurrentPatch[!is.na(CurrentPatch)] <- 1
          CurrentPatch[is.na(CurrentPatch)] <- 0
          CurrentPatch <- terra::mask(CurrentPatch, CurrentBinary)

          #Clip ensembled raster to boundary of occupied patches
          Ensemble2 <- terra::mask(CurrentEnsemble, CurrentPatch, maskvalue = 0)
          Ensemble2[is.na(Ensemble2)] <- 0
          Ensemble2 <- terra::mask(Ensemble2, CurrentBinary)

          #Write out rasters
          terra::writeRaster(CurrentBinary, file.path(result_dir, CurSpp,
                                                      paste0(CurrentTime, "_binary_original.grd")),
                             overwrite = TRUE)
          terra::writeRaster(CurrentEnsemble, file.path(result_dir, CurSpp,
                                                        paste0(CurrentTime, "_ensembled_original.grd")),
                             overwrite = TRUE)

          terra::writeRaster(CurrentPatch, file.path(result_dir, CurSpp,
                                                     paste0(CurrentTime, "_binary.grd")), overwrite = TRUE)
          terra::writeRaster(Ensemble2, file.path(result_dir, CurSpp,
                                                     paste0(CurrentTime, "_ensembled.grd")), overwrite = TRUE)
        }
        
        if (hindcast == FALSE) {
          CurrentBinary <- terra::rast(paste0(result_dir, "/", CurSpp, "/", CurrentTime, "_binary.grd"))
        } else {
          CurrentBinary <- terra::rast(paste0(result_dir, "/", CurSpp, "/", 
                                              startpoint[1], "/", startpoint[2], "_",
                                              startpoint[1],"_binary.grd"))
        }
        
        if(is.na(terra::linearUnits(CurrentBinary))) {
          stop("The spatial units of the data cannot be found!
           Choose a different coordinate system for all data or add units into the existing one.")
        }

        #Creates variables for stats
        Projection <- c("Current")
        NumberCells <- getSize(CurrentBinary)
        CellChange <- c(0)
        T1notT2 <- c(0)
        T2notT1 <- c(0)
        Overlap <- c(0)
        CentroidX <- c(getCentroid(CurrentBinary)[1])
        CentroidY <- c(getCentroid(CurrentBinary)[2])

        #Creates dispersal rasters and PDFs for each Scenario + time
        for (s in 1:length(scenarios)) {
          CurScen <- scenarios[s]
          curdir <- file.path(result_dir, CurSpp, CurScen)
          DispersalNames <- c()
          TimeMap <- terra::rast(file.path(result_dir, CurSpp, "TimeMapRasters", paste0("binary", CurScen, ".grd")))
          for (y in 2:length(time_periods)) {
            CurYear <- time_periods[y]

            #Calculates distance from current distribution
            if (y == 2) {
              DistanceRastersExist <- list.files(path = file.path(result_dir, CurSpp),
                                                 pattern = paste0("distance_", time_periods[1], "_Current.grd$"))
              if (length(DistanceRastersExist) == 0) {
                SppDistance <- DistanceRaster(CurSpp, CurrentTime, "Current", CurrentBinary, TimeMap)
                OriginalDistance <- SppDistance
              } else {
                OriginalDistance <- terra::rast(file.path(result_dir, CurSpp, DistanceRastersExist))
                SppDistance <- OriginalDistance
              }
            } else {
              FocusTime <- time_periods[y - 1]

              #Reads in the current distribution (from the previous time step)
              CurrentDistribution <- Binary_Dispersal

              #Sizes down the raster for faster distance measuring
              TimeMap2 <- TimeMap
              TimeMap2[which(terra::values(TimeMap2) == 0)] <- NA
              TimeMap2 <- terra::trim(TimeMap2)
              TimeMap2[which(is.na(terra::values(TimeMap2)))] <- 0

              #Highlights the places where species lived in the previous time step
              CurrentDistribution <- terra::crop(CurrentDistribution, terra::ext(TimeMap2))
              CurrentDistribution[which(terra::values(CurrentDistribution) == 0)] <- NA

              #Creates new distance raster
              CurrentDistance <- terra::distance(CurrentDistribution) / 1000
              SppDistance <- terra::extend(CurrentDistance, terra::ext(Binary_Dispersal))
              SppDistance[which(is.na(terra::values(SppDistance)))] <- max(terra::values(SppDistance), na.rm = TRUE) + 1

              if ((terra::ext(SppDistance) != terra::ext(CurrentBinary)) | (terra::ncell(SppDistance) != terra::ncell(CurrentBinary))) {
                message("Raster extents are not consistent: only the intersection of the rasters will be analysed")
                SppDistance <- terra::intersect(SppDistance, CurrentBinary)
                SppDistance <- terra::resample(SppDistance, CurrentBinary, method = "bilinear")
              }
              SppDistance <- terra::mask(SppDistance, CurrentBinary)
            }

            #Calculates the dispersal probability for the given time step
            TimeDiff <- abs(CurYear - time_periods[y - 1])
            SppDispProb <- DispersalProbRaster(dispersal_rate, SppDistance, TimeDiff)
            
            if(CurYear == TrueCurr) {
              RasterList <- list.files(file.path(result_dir, CurSpp), pattern = paste0(".grd$"))
              EnsembleNum <- grep(paste0(CurYear, "_ensembled.grd"), RasterList)
              EnsembleSD <- terra::rast(file.path(result_dir, CurSpp,
                                                  RasterList[EnsembleNum]))
            } else {
              RasterList <- list.files(path = curdir, pattern = paste0(".grd$"))
              EnsembleNum <- grep(paste0(CurYear, "_", CurScen, "_ensembled.grd"), RasterList)
              EnsembleSD <- terra::rast(file.path(curdir, RasterList[EnsembleNum]))
            }
            

            #Creates an ensembled raster that incorporates dispersal rate
            #Calculates the ensembled dispersal probability * habitat suitability to make "invadable suitability"
            if ((terra::ext(SppDispProb) != terra::ext(EnsembleSD)) | (terra::ncell(SppDispProb) != terra::ncell(EnsembleSD))) {
              message("Raster extents are not consistent: only the intersection of the rasters will be analysed")
              SppDispProb <- terra::intersect(SppDispProb, EnsembleSD)
              SppDispProb <- terra::resample(SppDispProb, EnsembleSD, method = "bilinear")
            }

            Ensemble_Dispersal <- SppDispProb * EnsembleSD
            if (!is.na(terra::crs(EnsembleSD))) {
              terra::crs(Ensemble_Dispersal) <- terra::crs(EnsembleSD)
            }

            #Writes the ensembled dispersal rate raster
            terra::writeRaster(Ensemble_Dispersal,
                        filename = file.path(curdir, paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate.grd")),
                        overwrite = TRUE)

            DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate"))

            #Creates a binary raster from the ensembled raster
            if(CurYear == TrueCurr) {
              BinaryNum <- grep(paste0(CurYear, "_binary.grd"), RasterList)
              BinarySD <- terra::rast(file.path(result_dir, CurSpp,
                                                  RasterList[EnsembleNum]))
            } else {
              BinaryNum <- grep(paste0(CurYear, "_", CurScen, "_binary.grd"), RasterList)
              BinarySD <- terra::rast(file.path(curdir, RasterList[BinaryNum]))
            }
            
            
            Binary_Dispersal <- (SppDispProb * BinarySD)
            Binary_Dispersal[Binary_Dispersal >= 0.5] <- 1
            Binary_Dispersal[Binary_Dispersal < 0.5] <- 0

            if (!is.na(terra::crs(EnsembleSD))) {
              terra::crs(Binary_Dispersal) <- terra::crs(EnsembleSD)
            }

            #Writes the binary raster
            terra::writeRaster(Binary_Dispersal,
                        filename = file.path(curdir, paste0(CurYear, "_", CurScen, "_binary_dispersalRate.grd")),
                        overwrite = TRUE)
            DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_binary_dispersalRate"))

            #Fills out the stats table
            Projection <- c(Projection, paste0(CurScen, "_", CurYear))
            FocusNCells <- getSize(Binary_Dispersal)
            NumberCells <- c(NumberCells, FocusNCells)
            CellChange <- c(CellChange, NumberCells[length(NumberCells)] - NumberCells[1])
            T1notT2_rast <- t1nott2(CurrentBinary, Binary_Dispersal)
            T2notT1_rast <- t1nott2(Binary_Dispersal, CurrentBinary)
            T1notT2 <- c(T1notT2, getSize(T1notT2_rast))
            T2notT1 <- c(T2notT1, getSize(T2notT1_rast))
            Overlap_rast <- overlap(Binary_Dispersal, CurrentBinary)
            Overlap <- c(Overlap, getSize(Overlap_rast))
            CentroidX <- c(CentroidX, getCentroid(Binary_Dispersal)[1])
            CentroidY <- c(CentroidY, getCentroid(Binary_Dispersal)[2])
          }

          #Writes PDFs
          DispersalRasters <- list.files(path = curdir, pattern = paste0("dispersalRate.grd$"), full.names = TRUE)
          if(hindcast == TRUE) {
            currentrastind <- grep(paste0("/\\", CurrentTime), DispersalRasters)
            if(length(currentrastind) > 0) {
              DispersalRasters <- DispersalRasters[-currentrastind]
            }
          }
          DispersalRasters <- gtools::mixedsort(DispersalRasters)
          DispersalNames <- gtools::mixedsort(DispersalNames)

          dir.create(file.path(result_dir, CurSpp, "map_pdfs"))
          for (d in 1:length(DispersalRasters)) {
            if (grepl("binary", DispersalRasters[d])) {
              title <- DispersalNames[d]
              grDevices::pdf(file = file.path(result_dir, CurSpp, "map_pdfs", paste0(CurSpp, "_", DispersalNames[d], ".pdf")))
              terra::plot(terra::rast(DispersalRasters[d]), legend = FALSE, main = title)
              graphics::legend("bottomright", legend = c("Absence", "Presence"), fill = c("white", "forestgreen"))
              grDevices::dev.off()
            } else {
              title <- DispersalNames[d]
              grDevices::pdf(file = file.path(result_dir, CurSpp, "map_pdfs", paste0(CurSpp, "_", DispersalNames[d], ".pdf")))
              terra::plot(terra::rast(DispersalRasters[d]), main = title)
              grDevices::dev.off()
            }
          }
        }

        #Fills out stats table
        stats <- data.frame(Projection = Projection,
                            NumberCells = NumberCells,
                            CellChange = CellChange,
                            T1notT2 = T1notT2,
                            T2notT1 = T2notT1,
                            Overlap = Overlap,
                            CentroidX = CentroidX,
                            CentroidY = CentroidY)

        utils::write.csv(stats, file = file.path(result_dir, CurSpp, "Results_Dispersal.csv"))
        rm(CurrentBinary)
        gc()
      } else {
        message(paste0("No dispersal rate values found for ", speciesName, ": skipping dispersal rate analysis"))
      }
    } else {
      message(paste0("No dispersal rate data found for ", speciesName, ": skipping dispersal rate analysis"))
    }
  }

  if (ncores == 1) {
    ListSpp <- as.vector(ListSpp)
    out <- sapply(ListSpp, function(x) FinalDispersal(x))

  } else {
    clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)

    parallel::clusterExport(clus, varlist = c("FinalDispersal", "getCentroid", "DistanceRaster",
                                              "DispersalProbRaster", "t1nott2", "overlap", "getSize",
                                              "numScenario", "numYear", "dispersal",
                                              "result_dir", "time_periods", "scenarios",
                                              "contiguous", "point_data", "ncores", "ListSpp"), envir = environment())

    parallel::clusterEvalQ(clus, library(gtools))
    parallel::clusterEvalQ(clus, library(terra))

    for(i in 1:ncol(ListSpp)) {
      for(j in 1:nrow(ListSpp)) {
        FinalDispersal(ListSpp[j,i])
      }
    }

    for (i in 1:nrow(ListSpp)) {
      out <- parallel::parLapply(clus, ListSpp[i, ], function(x) FinalDispersal(x))
      gc()
    }
    parallel::stopCluster(clus)
  }
}
