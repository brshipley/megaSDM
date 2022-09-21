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
#' @param result_dir the directory where the ensembled and binary maps are placed.
#' Each species should have its own sub-directory, and the forecasted/hindcasted
#' binary maps should be placed into directories like so: Species/Scenario/Time.
#' If \code{projectSuit} was used to make these maps, this is probably the same
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
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @export
#' @return returns rasters and .pdf files of the projected species ranges for future time
#' periods, given different climate scenarios, when the ability of the species to disperse
#' is taken into account. Dispersal-constrained presence/absence maps and "invadable suitability"
#' maps are provided as outputs.

dispersalRate <- function(result_dir, dispersaldata, time_periods,
                          scenarios, ncores = 1) {

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

  #Calculates the range centroid (average latitude, longitude)
  getCentroid <- function(raster) {
    # A matrix with three columns: x, y, and v (value)
    points <- raster::rasterToPoints(raster, fun = function(x){x == 1}, spatial = FALSE)

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
    CurrentBinary <- CurrentBinary
    CurrentPresence <- raster::reclassify(CurrentBinary, c(0, 0, NA), include.lowest = TRUE)
    #Trims the time map as an extent template for faster calculations
    TimeMap2 <- TimeMap
    TimeMap2[which(raster::values(TimeMap2) == 0)] <- NA
    TimeMap2 <- raster::trim(TimeMap2)
    TimeMap2[is.na(raster::values(TimeMap2))] <- 0
    #Trims "CurrentPresence" to the extent of the time maps
    CurrentPresence <- raster::crop(CurrentPresence, raster::extent(TimeMap2))
    #Calculates the distances from each pixel to the nearest presence
    CurrDist <- raster::distance(CurrentPresence, doEdge = TRUE)
    #Extends the raster back out to full study area extent
    CurrDist <- raster::extend(CurrDist, raster::extent(CurrentBinary), value = max(raster::values(CurrDist), na.rm = TRUE) + 1)
    CurrDist[is.na(raster::values(CurrDist))] <- max(raster::values(CurrDist), na.rm = TRUE) + 1
    DistFinal <- raster::mask(CurrDist, CurrentBinary)
    #Converts distance (in meters) to kilometers (for dispersal rate)
    DistFinal <- DistFinal / 1000
    #writes distance raster
    # writeRaster(DistFinal,
    #             filename = paste0(result_dir, "/", spp, "/", "distance_", Time, "_", Scen, ".bil"),
    #             overwrite = TRUE,
    #             format = format,
    #             prj = TRUE)
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
    DistProb <- raster::calc(DistRaster, fun = function(x){GammaProbFun(x)})
    return(DistProb)
  }

  #Gets presence pixels in raster 1 but not raster 2
  t1nott2 <- function(t1, t2) {
    return(raster::mask(t1, t2, inverse = FALSE, maskvalue = 1, updatevalue = 0))
  }

  #Calculates overlap between raster 1 and raster 2 presences
  overlap <- function(t1, t2) {
    return(raster::mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
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
      if (!is.na(dispersal_rate)) {
        CurrentTime <- time_periods[1]
        CurrentBinary <- raster::raster(paste0(result_dir, "/", CurSpp, "/", CurrentTime, "_binary.bil"))
        #Creates variables for stats
        Projection <- c("Current")
        NumberCells <- c(raster::cellStats(CurrentBinary, stat = sum))
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
          TimeMap <- raster::raster(file.path(result_dir, CurSpp, "TimeMapRasters", paste0("binary", CurScen, ".bil")))
          for (y in 2:length(time_periods)) {
            CurYear <- time_periods[y]

            #Calculates distance from current distribution
            if (y == 2) {
              DistanceRastersExist <- list.files(path = file.path(result_dir, CurSpp),
                                                 pattern = paste0("distance_", time_periods[1], "_Current.bil$"))
              if (length(DistanceRastersExist) == 0) {
                SppDistance <- DistanceRaster(CurSpp, CurrentTime, "Current", CurrentBinary, TimeMap)
                OriginalDistance <- SppDistance
              } else {
                OriginalDistance <- raster::raster(file.path(result_dir, CurSpp, DistanceRastersExist))
                SppDistance <- OriginalDistance
              }
            } else {
              FocusTime <- time_periods[y - 1]

              #Reads in the current distribution (from the previous time step)
              CurrentDistribution <- Binary_Dispersal

              #Sizes down the raster for faster distance measuring
              TimeMap2 <- TimeMap
              TimeMap2[which(raster::values(TimeMap2) == 0)] <- NA
              TimeMap2 <- raster::trim(TimeMap2)
              TimeMap2[is.na(raster::values(TimeMap2))] <- 0

              #Highlights the places where species lived in the previous time step
              CurrentDistribution <- raster::crop(CurrentDistribution, raster::extent(TimeMap2))
              CurrentDistribution[which(raster::values(CurrentDistribution) == 0)] <- NA

              #Creates new distance raster
              CurrentDistance <- raster::distance(CurrentDistribution) / 1000
              SppDistance <- raster::extend(CurrentDistance, raster::extent(Binary_Dispersal),
                                            value = max(raster::values(SppDistance), na.rm = TRUE) + 1)

              if ((raster::extent(SppDistance) != raster::extent(CurrentBinary)) | (raster::ncell(SppDistance) != raster::ncell(CurrentBinary))) {
                message("Raster extents are not consistent: only the intersection of the rasters will be analysed")
                SppDistance <- raster::intersect(SppDistance, CurrentBinary)
                SppDistance <- raster::resample(SppDistance, CurrentBinary, method = "bilinear")
              }
              SppDistance <- raster::mask(SppDistance, CurrentBinary)
            }

            #Calculates the dispersal probability for the given time step
            TimeDiff <- abs(CurYear - time_periods[y - 1])
            SppDispProb <- DispersalProbRaster(dispersal_rate, SppDistance, TimeDiff)
            RasterList <- list.files(path = curdir, pattern = paste0(".bil$"))

            #Creates an ensembled raster that incorporates dispersal rate
            #Calculates the ensembled dispersal probability * habitat suitability to make "invadable suitability"

            EnsembleNum <- grep(paste0(CurYear, "_", CurScen, "_ensembled.bil"), RasterList)
            EnsembleSD <- raster::raster(file.path(curdir, RasterList[EnsembleNum]))
            if ((raster::extent(SppDispProb) != raster::extent(EnsembleSD)) | (raster::ncell(SppDispProb) != raster::ncell(EnsembleSD))) {
              message("Raster extents are not consistent: only the intersection of the rasters will be analysed")
              SppDispProb <- raster::intersect(SppDispProb, EnsembleSD)
              SppDispProb <- raster::resample(SppDispProb, EnsembleSD, method = "bilinear")
            }

            Ensemble_Dispersal <- SppDispProb * EnsembleSD
            if (!is.na(raster::crs(EnsembleSD))) {
              raster::crs(Ensemble_Dispersal) <- raster::crs(EnsembleSD)
            }

            #Writes the ensembled dispersal rate raster
            raster::writeRaster(Ensemble_Dispersal,
                        filename = file.path(curdir, paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate.bil")),
                        overwrite = TRUE,
                        format = "EHdr",
                        prj = TRUE)

            DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate"))

            #Creates a binary raster from the ensembled raster
            BinaryNum <- grep(paste0(CurYear, "_", CurScen, "_binary.bil"), RasterList)
            BinarySD <- raster::raster(file.path(curdir, RasterList[BinaryNum]))
            Binary_Dispersal <- (SppDispProb * BinarySD)
            Binary_Dispersal[Binary_Dispersal >= 0.5] <- 1
            Binary_Dispersal[Binary_Dispersal < 0.5] <- 0

            if (!is.na(raster::crs(EnsembleSD))) {
              raster::crs(Binary_Dispersal) <- raster::crs(EnsembleSD)
            }

            #Writes the binary raster
            raster::writeRaster(Binary_Dispersal,
                        filename = file.path(curdir, paste0(CurYear, "_", CurScen, "_binary_dispersalRate.bil")),
                        overwrite = TRUE,
                        format = "EHdr",
                        prj = TRUE)
            DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_binary_dispersalRate"))

            #Fills out the stats table
            Projection <- c(Projection, paste0(CurScen, "_", CurYear))
            NumberCells <- c(NumberCells, raster::cellStats(Binary_Dispersal, stat = sum))
            CellChange <- c(CellChange, NumberCells[length(NumberCells)] - NumberCells[1])
            T1notT2 <- c(T1notT2, raster::cellStats(t1nott2(CurrentBinary, Binary_Dispersal), stat = sum))
            T2notT1 <- c(T2notT1, raster::cellStats(t1nott2(Binary_Dispersal, CurrentBinary), stat = sum))
            Overlap <- c(Overlap, raster::cellStats(overlap(Binary_Dispersal, CurrentBinary), stat = sum))
            CentroidX <- c(CentroidX, getCentroid(Binary_Dispersal)[1])
            CentroidY <- c(CentroidY, getCentroid(Binary_Dispersal)[2])
          }

          #Writes PDFs
          DispersalRasters <- list.files(path = curdir, pattern = paste0("dispersalRate.bil$"), full.names = TRUE)
          DispersalRasters <- gtools::mixedsort(DispersalRasters)
          DispersalNames <- gtools::mixedsort(DispersalNames)

          dir.create(file.path(result_dir, CurSpp, "map_pdfs"))
          for (d in 1:length(DispersalRasters)) {
            if (grepl("binary", DispersalRasters[d])) {
              title <- DispersalNames[d]
              grDevices::pdf(file = file.path(result_dir, CurSpp, "map_pdfs", paste0(CurSpp, "_", DispersalNames[d], ".pdf")))
              raster::plot(raster::raster(DispersalRasters[d], native = TRUE), legend = FALSE, main = title)
              graphics::legend("bottomright", legend = c("Absence", "Presence"), fill = c("white", "forestgreen"))
              grDevices::dev.off()
            } else {
              title <- DispersalNames[d]
              grDevices::pdf(file = file.path(result_dir, CurSpp, "map_pdfs", paste0(CurSpp, "_", DispersalNames[d], ".pdf")))
              raster::plot(raster::raster(DispersalRasters[d], native = TRUE), main = title)
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
                                              "DispersalProbRaster", "t1nott2", "overlap",
                                              "numScenario", "numYear", "dispersal",
                                              "result_dir", "time_periods", "scenarios",
                                              "ncores", "ListSpp"), envir = environment())
    
    parallel::clusterEvalQ(clus, library(gtools))
    parallel::clusterEvalQ(clus, library(raster))
    
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
