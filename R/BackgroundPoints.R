#' Generate background points for species distribution modelling
#'
#' This function generates a set of species-specific background points using a variety
#' of methods. These points can be generally across a given training area, or if
#' environmental data are provided, environmental subsampling (sensu Varela et al. 2014)
#' can be conducted. If a list of buffers around the occurrence points of each species
#' are provided, \code{BackgroundPoints} can conduct spatially-constrained sampling
#' within the buffer (\code{spatial_weights} = 1) or a hybrid method (\code{spatial_weights} between 0 and 1).
#'
#' @param spplist a vector of species names (used for the generation of output file names).
#' @param envdata a SpatRaster or list of raster files corresponding to the
#' area the model will be trained on. If environmental subsampling is desired (method = \code{"Varela"}),
#' all environmental layers used for the modelling should be included.
#' @param output a file folder where the output background point files will be placed.
#' @param nbg either a single integer or a vector of the same length as \code{spplist},
#' corresponding to the number of background points generated for each set of background points.
#' Should be in the same order as the species list.
#' @param spatial_weights a number between 0 and 1 determining the percentage of background points
#' sampled within the given buffer. A 1 indicates that background points are exclusively sampled
#' from within the buffers, whereas a 0 indicates random sampling throughout the training area.
#' @param buffers (optional) a list of shapefiles (SpatVec) or a list of file paths
#' for the buffers. Required if \code{spatial_weights} is not 0. Should be in the same order
#' as the species list. If only one species is needed, can be a single SpatialPolygon* object.
#' @param method either "Varela" or "random", where "Varela" incorporates environmental subsampling
#' (see Varela et al. 2014) and "random" simply randomly samples background points.
#' Note that if method = "Varela", exact numbers of background points cannot be guaranteed.
#' @param PCA (optional). The number of PC axes to use when binning climate data for environmental
#' subsampling. Can be a number from 1 to the number of environmental variables. Default selects the
#' number of axes that account for 95% of the environmental variation. If set to \code{NA},
#' PCA will not be run (for use when there are categorical variables, for instance).
#' @param ncores the number of computer cores to parallelize the background point generation on. Default is 1.
#'     Using one fewer core than the computer has is usually optimal.
#' @export
#' @return writes .csv files of background points for each species in \code{spplist} to a
#' directory provided by the \code{output} argument. These files have the coordinates of the points and the
#' environmental values at each point. In addition, if \code{method = "Varela"}, a .csv file is written
#' out with the number of background points generated for each species, as the exact number of
#' points may vary.


BackgroundPoints <- function(spplist, envdata, output,
                             nbg = 5000,
                             spatial_weights = 0,
                             buffers,
                             method = "Varela",
                             PCA = "Y",
                             ncores = 1) {

  if(!dir.exists(output)) {
    dir.create(output)
  }

  #Ensures that buffers are provided if spatial_weights are greater than 0
  if (spatial_weights > 0) {
    if (!methods::hasArg(buffers)) {
      stop("Spatially-constrained background points wanted, but no buffers given")
    } else if (length(buffers) != length(spplist)) {
      stop("The number of buffers does not match the number of species")
    }
  } else {
    buffers <- NA
  }

  #Ensures that nbg is the same length as spplist.
  if ((length(nbg) == 1)) {
    nbg <- rep(nbg, length(spplist))
  } else if (length(nbg) != length(spplist)) {
    stop("length(nbg) is not equal to length(spplist): Ensure that each element of spplist
         is associated with exactly one element of nbg")
  }

  #If the input training layers are not in SpatRaster form, ensure that they have the same projection/extent
  if (class(envdata) != "SpatRaster") {

    #Ensure that all training area rasters have the same projection and extent
    projtrain <- rep(NA, len = length(envdata))
    exttrain <- rep(NA, len = length(envdata))
    for (i in 1:length(envdata)) {
      projtrain[i] <- as.character(terra::crs(terra::rast(envdata[[i]])))
      exttrain[i] <- as.character(terra::crs(terra::rast(envdata[[i]])))
    }

    if (length(unique(projtrain)) > 1) {
      stop("Not all of the training area environmental rasters are in the same projection")
    } else if (length(unique(exttrain)) > 1) {
      stop("Not all of the training area environmental rasters have the same extent")
    }
    envstack <- terra::rast(envdata)
  } else {
    envstack <- envdata
  }

  #Make sure that the training layers have a CRS that is not NA
  
  if (is.na(terra::crs(envstack))) {
    stop("training area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }

  #Reads in the Varela sampling function and provides a folder for stats
  if (method == "Varela") {
    #creates a BG_Stats directory
    if (!dir.exists(paste0(output, "/BG_Stats"))) {
      dir.create(paste0(output, "/BG_Stats"))
    }

    StatsLoc <- paste0(output, "/BG_Stats")

    #Function for Varela sampling
    VarelaSample <- function (occurrences, env, no_bins, PCA, PCAxes) {

      occurrences <- terra::vect(occurrences,
                                 geom = c("Longitude", "Latitude"),
                                 crs = terra::crs(env))

      EnvOccur <- terra::extract(env, occurrences)
      EnvOccur <- EnvOccur[, 2:ncol(EnvOccur)]
      EnvCoords <- data.frame(terra::geom(occurrences))
      EnvOccur <- data.frame(x = EnvCoords$x, y = EnvCoords$y, EnvOccur)
      ClimOccur <- EnvOccur
      ClimOccur <- ClimOccur[stats::complete.cases(ClimOccur), ]
      if (PCA == "Y") {
        PCAEnv <- stats::prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
        PCAImp <- summary(PCAEnv)$importance
        #Determine the number of PC axes to use for subsampling
        if (is.numeric(PCAxes)) {
          NumberAxes <- PCAxes
        } else {
          NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
        }

        #Add PCA values to the unsubsampled data frame
        EnvOccur <- data.frame(cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes]))
      }

      #make a landing spot for the bin membership vectors

      #cycle through all of the environmental variables (columns 3 to end)
      nsamples <- c()
      for (j in 1:length(no_bins)) {
        out_ptz <- EnvOccur[, 1:2]
        for(i in 3:length(names(EnvOccur))) {
          #make a data frame that is this variable with no NA values
          k <- EnvOccur[!is.na(EnvOccur[, i]), i]
          #calculate the observed range of this variable
          rg <- range(k)
          #figure out the resolution from the number of bins
          res <- (rg[2] - rg[1]) / no_bins[j]
          #rescale the axis by the range and bin size, so the value is just a
          #number from 1 to no_bins for its bin membership
          d <- (EnvOccur[, i] - rg[1]) / res
          #d is now a vector of values ranging from 0 to no_bins
          f <- ceiling(d)
          #f is a vector of bin membership
          f[f == 0] <- 1 #move the zeros into the 1 bin
          #correct the name of the vector, so it will carry over to the output
          names(f) <- names(EnvOccur)[i]
          #add the bin membership vector to the output df for this section
          out_ptz <- cbind(out_ptz, f)
          #get the names correct
          names(out_ptz)[length(names(out_ptz))] <- names(EnvOccur)[i]
        }

        #subsample the bin membership df to come up with the filled bins
        sub_ptz <- dplyr::distinct(out_ptz[, -1:-2])
        #count the number of filled bins
        no_grps <- nrow(sub_ptz)
        #add a column with the group membership number; this number is arbitrary
        sub_ptz$grp <- c(1:no_grps)

        #join the out_ptz with the subsample to capture the group membership info
        #note: join() will automatically match the variable names from these two dfs
        out_ptz <- suppressMessages(dplyr::left_join(out_ptz, sub_ptz))
        #out_ptz now has a group membership  for each input point

        #select a random point for each group -- this is an additional improvement on the
        #Varela et al. function, because it does not pick the same points each time.

        #make a landing spot for the data
        final_out <- data.frame(x = numeric(), y = numeric())

        #cycle through each group
        for(i in 1:no_grps) {

          #subset to the members of the ith group, keep only the Latitude and Longitude
          grp_mbrs <- out_ptz[out_ptz$grp == i, c(1, 2)]

          #pick one of these group members to output
          grp_out <- grp_mbrs[sample(1:nrow(grp_mbrs), 1), ]
          #bind this sampled row to the output df
          final_out <- rbind(final_out, grp_out)
        }

        #return the subsampled points as a df of Latitude and Longitude values
        final_out <- data.frame(x = final_out[, 1], y = final_out[, 2])
        final_out <- merge(final_out, ClimOccur, by = c("x", "y"), all.x = TRUE)
        nsamples <- c(nsamples, nrow(final_out))
      }
      if (length(no_bins) == 1) {
        return(final_out)
      } else {
        return(data.frame(NumberofSamples = nsamples, NumberOfClimateBins = no_bins))
      }
    }
    nbins <- 25
  } else {
    #We need these variables for the parallelization, but they won't be used
    nbins <- NA
    StatsLoc <- NA
    VarelaSample <- NA
  }

  if (length(spplist) < ncores) {
    ncores <- length(spplist)
  }

  ListSpp <- matrix(data = spplist, ncol = ncores)

  run <- function(CurSpp) {
    s <- grep(CurSpp, spplist)
    
    envstack <- terra::rast(envstack)
    
    #If the number of background points wanted is more than the number of cells, prints warning, sets nbg-->ncells
    if (nbg[s] > terra::ncell(envstack)) {
      nbg[s] <- terra::ncell(envstack)
      message(paste0("Warning: the number of background points wanted for ", spplist[s],
                     " is more than the number of unique cells in the environmental raster"))
      message(paste0("    Changing the number of background points to ", terra::ncell(envstack)))

    }

    nbgBuff <- round(nbg[s] * spatial_weights)
    nbgFull <- round(nbg[s] - nbgBuff)

    #Sampling from the full training area--------------------
    if (spatial_weights < 1){
      if (method == "Varela") {

        if (is.numeric(PCA)) {
          PCAxes <- as.numeric(PCA)
          if (PCAxes > terra::nlyr(envstack)) {
            PCAxes <- terra::nlyr(envstack)
            message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
          }
          PCA <- "Y"
        } else if (!is.na(PCA)) {
          PCAxes <- "Y"
        } else {
          PCA <- "N"
          PCAxes <- ""
        }
        
        if (nbgFull * 10 > terra::ncell(envstack[[1]])) {
          RastFullBG <- terra::ncell(envstack[[1]])
        } else {
          RastFullBG <- nbgFull * 10
        }

        RandomTrain <- terra::spatSample(envstack, RastFullBG, na.rm = TRUE, xy = TRUE)
        names(RandomTrain)[1:2] <- c("Longitude", "Latitude")
        nbinscount <- nbins

        NUniqueClim <- nrow(unique(RandomTrain[, 3:ncol(RandomTrain)]))

        if (NUniqueClim < nbgFull) {
          nbgFull <- NUniqueClim
          message(paste0("Warning! The number of unique climates within the buffer is "))
          message(paste0("less than the number of background points wanted"))
          message(paste0("for species ", spplist[s]))

          message(paste0("Only ", nbgFull, " background points within the buffer can be generated"))
          nbinscount <- 98
        }

        FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envstack, nbinscount, PCA, PCAxes)

        if (nrow(FullPointsEnv) >= nbgFull) {

          while(nrow(FullPointsEnv) > nbgFull && nbinscount > 2) {
            nbinscount <- nbinscount - 1
            FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envstack, nbinscount, PCA, PCAxes)
          }

          FullPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envstack, (nbinscount + 1), PCA, PCAxes)

          diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
          diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)

          if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
          }

        } else {

          while(nrow(FullPointsEnv) < nbgFull && nbinscount < max(nbins, 99)) {
            nbinscount <- nbinscount + 1
            FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envstack, nbinscount, PCA, PCAxes)
          }

          FullPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envstack, (nbinscount - 1), PCA, PCAxes)

          diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
          diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)

          if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
          }

        }

        FullPointsEnv <- FullPointsEnv[stats::complete.cases(FullPointsEnv), ]

      } else if (method == "random") {
        FullPointsEnv <- terra::spatSample(x = envstack, size = nbgFull, na.rm = TRUE, xy = TRUE)
        names(FullPointsEnv)[1:2] <- c("Longitude", "Latitude")
      }

      if (spatial_weights == 0) {
        background <- data.frame("Species" = rep(spplist[s], nrow(FullPointsEnv)), FullPointsEnv)

        utils::write.csv(background, paste0(output, "/", gsub(" ", "_", spplist[s]), "_background.csv"), row.names = FALSE)
        if (method == "Varela") {
          utils::write.csv(data.frame("Species" = spplist[s],
                               "FullTrainPoints" = nrow(background),
                               "BuffPoints" = 0,
                               "FullPoints" = nrow(background)),
                    paste0(StatsLoc, "/", spplist[s], "_stats.csv"), row.names=FALSE)
        }
        return()
      }
    }

    #Sampling from the buffer area------------------

    if (length(spplist) > 1) {
      bufferSHP <- buffers[[s]]
    } else {
      bufferSHP <- buffers
    }

    if (!grepl("SpatVector", class(bufferSHP))) {
      bufferSHP <- terra::vect(bufferSHP)
    }

    if (as.character(terra::crs(bufferSHP)) != as.character(terra::crs(envstack))) {
      bufferSHP <- terra::project(bufferSHP, terra::crs(envstack))
    }

    if (method == "Varela") {

      if (is.numeric(PCA)) {
        PCAxes <- as.numeric(PCA)
        if (PCAxes > terra::nlyr(envstack)) {
          PCAxes <- terra::nlyr(envstack)
          message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
        }
      } else if (!is.na(PCA)) {
        PCAxes <- "Y"
      } else {
        PCAxes <- ""
      }

      if (nbgBuff * 10 > terra::ncell(envstack[[1]])) {
        RastBuffBG <- terra::ncell(envstack[[1]])
      } else {
        RastBuffBG <- nbgBuff * 10
      }

      #randomly samples nbg*5 number of points from the shapefile
      RandomBuff <- terra::spatSample(bufferSHP, RastBuffBG, "random")

      RandomBuff <- terra::geom(RandomBuff)
      RandomBuff <- data.frame(RandomBuff[, c("x", "y")])
      names(RandomBuff) <- c("Longitude", "Latitude")
      nbinscount <- nbins

      #Prints a warning if the nubmer of unique climates is less than the buffer background points wanted
      BuffClim <- terra::extract(envstack, RandomBuff, ID = FALSE)
      NUniqueClim <- nrow(unique(BuffClim))

      if (NUniqueClim < nbgBuff) {
        nbgBuff <- NUniqueClim
        message("Warning! The number of unique climates within the buffer is less than the number of background points wanted")
        message(paste0("Only ", nbgBuff, " background points within the buffer can be generated"))
        nbinscount <- 98
      }

      BuffPointsEnv <- VarelaSample(RandomBuff, envstack, nbinscount, PCA, PCAxes)

      if (nrow(BuffPointsEnv) >= nbgBuff) {

        while(nrow(BuffPointsEnv) > nbgBuff && nbinscount > 2) {
          nbinscount <- nbinscount - 1
          BuffPointsEnv <- VarelaSample(RandomBuff[, 1:2], envstack, nbinscount, PCA, PCAxes)
        }

        BuffPointsEnv2 <- VarelaSample(RandomBuff[, 1:2], envstack, (nbinscount + 1), PCA, PCAxes)

        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)

        if (diff2 < diff1) {
          BuffPointsEnv <- BuffPointsEnv2
        }

      } else {

        while(nrow(BuffPointsEnv) < nbgBuff && nbinscount < max(nbins, 99)) {
          nbinscount <- nbinscount + 1
          BuffPointsEnv <- VarelaSample(RandomBuff[, 1:2], envstack, nbinscount, PCA, PCAxes)
        }

        BuffPointsEnv2 <- VarelaSample(RandomBuff[, 1:2], envstack, (nbinscount - 1), PCA, PCAxes)

        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)

        if (diff2 < diff1) {
          BuffPointsEnv <- BuffPointsEnv2
        }

      }

      BuffPointsEnv <- BuffPointsEnv2[stats::complete.cases(BuffPointsEnv2), ]

    } else if (method == "random") {
      BuffPoints <- terra::spatSample(bufferSHP, nbgBuff * 10, "random")
      BuffPointsCoords <- terra::geom(BuffPoints)
      BuffPointsCoords <- data.frame(BuffPointsCoords[, c("x", "y")])
      names(BuffPointsCoords) <- c("Longitude", "Latitude")
      
      BuffPointsEnv <- terra::extract(envstack, BuffPointsCoords, ID = FALSE)
      BuffPointsEnv <- cbind(BuffPointsCoords, BuffPointsEnv)
      BuffPointsEnv <- stats::na.omit(BuffPointsEnv)
      BuffPointsEnv <- BuffPointsEnv[sample(c(1:nrow(BuffPointsEnv)), nbgBuff, replace = FALSE), ]
    }
    
    if (spatial_weights == 0) {
      background <- FullPointsEnv
      BuffPointsEnv <- FullPointsEnv[c(), ]
    } else if (spatial_weights == 1) {
      background <- BuffPointsEnv
      FullPointsEnv <- BuffPointsEnv[c(), ]
    } else {
      background <- rbind(FullPointsEnv, BuffPointsEnv)
    }
    
    background <- data.frame("Species" = rep(spplist[s], nrow(background)), background)

    if (method == "Varela") {
      utils::write.csv(data.frame("Species" = spplist[s],
                           "FullTrainPoints" = nrow(FullPointsEnv),
                           "BuffPoints" = nrow(BuffPointsEnv),
                           "FullPoints" = nrow(background)),
                paste0(StatsLoc, "/", gsub(" ", "_", spplist[s]), "_stats.csv"), row.names = FALSE)
    }

    utils::write.csv(background, paste0(output, "/", gsub(" ", "_", spplist[s]), "_background.csv"), row.names = FALSE)

  }
  
  #Convert envstack to a "packed" raster so that it can be parallelized
  envstack <- terra::wrap(envstack)
  
  if (ncores == 1) {
    ListSpp <- as.vector(ListSpp)
    out <- sapply(ListSpp, function(x) run(x))
  } else {
    clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)
    
    parallel::clusterExport(clus, varlist = c("VarelaSample", "run", "spplist", "envstack",
                                              "output", "nbg", "spatial_weights",
                                              "buffers", "method", "PCA", "ListSpp",
                                              "nbins", "StatsLoc"), envir = environment())
    
    parallel::clusterEvalQ(clus, library(dplyr))
    parallel::clusterEvalQ(clus, library(terra))
    
    for (i in 1:nrow(ListSpp)) {
      out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
    }
    
    parallel::stopCluster(clus)
  }
  
}
