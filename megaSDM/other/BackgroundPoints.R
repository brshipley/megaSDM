#' Generate Background Points for Species Distribution Modelling
#'
#' This function generates a set of species-specific background points using a variety
#' of methods. These points can be generally across a given training area, or if
#' environmental data are provided, environmental subsampling (sensu Varela et al. 2014)
#' can be conducted. If a list of buffers around the occurrence points of each species
#' are provided, \code{BackgroundPoints} can conduct spatially-constrained sampling
#' within the buffer (\code{spatial_weights} = 1) or a hybrid method (\code{spatial_weights} < 1 and > 0).
#'
#' @param spplist a vector of species names (used for the generation of output file names)
#' @param envdata a RasterStack or list of raster files corresponding to the
#' area the model will be trained on. If environmental subsampling is desired (method = \code{"Varela"}),
#' all environmental layers used for the modelling should be included
#' @param output a file folder where the output background point files will be placed.
#' @param nbg either a single integer or a vector of the same length as \code{spplist},
#' corresponding to the number of background points generated for each set of background points.
#' Should be in the same order as the species list.
#' @param spatial_weights a number between 0 and 1 determining the percentage of background points
#' sampled within the given buffer. A 1 indicates that background points are exclusively sampled
#' from within the buffers, whereas a 0 indicates random sampling throughout the training area.
#' @param buffers (optional) a list of shapefiles (SpatialPolygons*) or a list of file paths
#' for the buffers. Required if \code{spatial_weights} is not 0. Should be in the same order
#'  as the species list.If only one species is needed, can be a single SpatialPolygon* object.
#' @param method either "Varela" or "random", where "Varela" incorporates environmental subsampling
#' (see Varela et al. 2014) and "random" simply randomly samples background points.
#' Note that if method = "Varela", exact numbers of background points cannot be guaranteed.
#' @param PCA (optional). The number of PC axes to use when binning climate data for environmental
#' subsampling. Can be a number from 1 to the number of environmental variables. Default selects the
#' number of axes that account for 95% of the environmental variation. If set to \code{NA},
#' PCA will not be run (for use when there are categorical variables, for instance).
#' @export
#' @return writes .csv files of background points for each species in \code{spplist} to a
#' directory provided by \code{output}. These files have the coordinates of the points and the
#' environmental values at each point. In addition, if method = "Varela", a .csv file is written
#' out with the number of background points generated for each species, as the exact number of
#' points may vary.


BackgroundPoints <- function(spplist, envdata, output,
                             nbg = 5000,
                             spatial_weights = 0,
                             buffers,
                             method = "Varela",
                             PCA = "Y") {
  library(raster)
  library(sp)
  library(dplyr)

  #Ensures that buffers are provided if spatial_weights are greater than 0
  if (spatial_weights > 0) {
    if (!hasArg(buffers)) {
      stop("Spatially-constrained background points wanted, but no buffers given")
    } else if (length(buffers) != length(spplist)) {
      stop("The number of buffers does not match the number of species")
    }
  }

  #Ensures that nbg is the same length as spplist.
  if ((length(nbg) == 1)) {
    nbg <- rep(nbg, length(spplist))
  } else if (length(nbg) != length(spplist)) {
    stop("length(nbg) is not equal to length(spplist): Ensure that each element of spplist
         is associated with exactly one element of nbg")
  }

  if (method == "Varela") {
    #creates a BG_Stats directory
    if (!dir.exists(paste0(output, "/BG_Stats"))) {
      dir.create(paste0(output, "/BG_Stats"))
      StatsLoc <- paste0(output, "/BG_Stats")
    }

    #Function for Varela sampling
    VarelaSample <- function (occurrences, env, no_bins, PCA, PCAxes) {

      occurrences <- sp::SpatialPoints(occurrences, proj4string = crs(env))

      EnvOccur <- raster::extract(env, occurrences)
      EnvOccur <- data.frame(x = occurrences@coords[, 1], y = occurrences@coords[, 2], EnvOccur)
      ClimOccur <- EnvOccur
      ClimOccur <- ClimOccur[complete.cases(ClimOccur), ]
      if (PCA == "Y") {
        PCAEnv <- prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
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
    nbins <- 20
  }

  for (s in 1:length(spplist)) {

    #If the number of background points wanted is more than the number of cells, prints warning, sets nbg-->ncells
    if (nbg[s] > ncell(envdata)) {
      nbg[s] <- ncell(envdata)
      message(paste0("Warning: the number of background points wanted for ", spplist[s],
                     " is more than the number of unique cells in the environmental raster"))
      message(paste0("    Changing the number of background points to ", ncell(envdata)))

    }

    nbgBuff <- round(nbg[s] * spatial_weights)
    nbgFull <- round(nbg[s] - nbgBuff)

    #Sampling fron the full training area--------------------
    if (spatial_weights < 1){
      if (method == "Varela") {

        if (is.numeric(PCA)) {
          PCAxes <- as.numeric(PCA)
          if (PCAxes > nlayers(envdata)) {
            PCAxes <- nlayers(envdata)
            message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
          }
        } else if (!is.na(PCA)) {
          PCAxes <- "Y"
        } else {
          PCA <- "N"
          PCAxes <- ""
        }

        if (nbgFull * 10 > raster::ncell(envdata[[1]])) {
          RastFullBG <- raster::ncell(envdata[[1]])
        } else {
          RastFullBG <- nbgFull * 10
        }

        RandomTrain <- raster::sampleRandom(envdata, RastFullBG, na.rm = TRUE, xy = TRUE)
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

        FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCA, PCAxes)

        if (nrow(FullPointsEnv) >= nbgFull) {

          while(nrow(FullPointsEnv) > nbgFull && nbinscount > 2) {
            nbinscount <- nbinscount - 1
            FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCA, PCAxes)
          }

          FullPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, (nbinscount + 1), PCA, PCAxes)

          diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
          diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)

          if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
          }

        } else {

          while(nrow(FullPointsEnv) < nbgFull && nbinscount < max(nbins, 99)) {
            nbinscount <- nbinscount + 1
            FullPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCA, PCAxes)
          }

          FullPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, (nbinscount - 1), PCA, PCAxes)

          diff1 <- abs(nrow(FullPointsEnv) - nbgFull)
          diff2 <- abs(nrow(FullPointsEnv2) - nbgFull)

          if (diff2 < diff1) {
            FullPointsEnv <- FullPointsEnv2
          }

        }

        FullPointsEnv <- FullPointsEnv[complete.cases(FullPointsEnv), ]

      } else if (method == "random") {
        FullPointsEnv <- sampleRandom(envdata, nbgFull, na.rm = TRUE, xy = TRUE)
      }

      if (spatial_weights == 0) {
        background <- FullPointsEnv
        write.csv(background, paste0(output, "/", sub(" ", "_", spplist[s]), ".csv"), row.names = FALSE)
        if (method == "Varela") {
          write.csv(data.frame("Species" = spplist[s],
                               "FullTrainPoints" = nrow(background),
                               "BuffPoints" = 0,
                               "FullPoints" = nrow(background)),
                    paste0(StatsLoc, "/", spplist[s], "_stats.csv"), row.names=FALSE)
        }

        next()
      }
    }
    #Sampling from the buffer area------------------

    if (length(spplist) > 1) {
      bufferSHP <- buffers[[1]]
    } else {
      bufferSHP <- buffers
    }

    if (!grepl("Spatial", class(bufferSHP))) {
      bufferSHP <- raster::shapefile(bufferSHP)
    }

    if (as.character(raster::crs(bufferSHP)) != as.character(raster::crs(envdata))) {
      bufferSHP <- sp::spTransform(bufferSHP, CRSobj = raster::crs(envdata))
    }

    if (method == "Varela") {

      if (is.numeric(PCA)) {
        PCAxes <- as.numeric(PCA)
        if (PCAxes > nlayers(envdata)) {
          PCAxes <- nlayers(envdata)
          message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
        }
      } else if (!is.na(PCA)) {
        PCAxes <- "Y"
      } else {
        PCAxes <- ""
      }

      if (nbgBuff * 10 > raster::ncell(envdata[[1]])) {
        RastBuffBG <- raster::ncell(envdata[[1]])
      } else {
        RastBuffBG <- nbgBuff * 10
      }

      #randomly samples nbg*5 number of points from the shapefile
      RandomBuff <- sp::spsample(bufferSHP, RastBuffBG, "random", iter = 15)

      RandomBuff <- RandomBuff@coords

      nbinscount <- nbins

      #Prints a warning if the nubmer of unique climates is less than the buffer background points wanted
      BuffClim <- raster::extract(envdata, RandomBuff)
      NUniqueClim <- nrow(unique(BuffClim))

      if (NUniqueClim < nbgBuff) {
        nbgBuff <- NUniqueClim
        message("Warning! The number of unique climates within the buffer is less than the number of background points wanted")
        message(paste0("Only ", nbgBuff, " background points within the buffer can be generated"))
        nbinscount <- 98
      }



      BuffPointsEnv <- VarelaSample(RandomBuff, envdata, nbinscount, PCA, PCAxes)

      if (nrow(BuffPointsEnv) >= nbgBuff) {

        while(nrow(BuffPointsEnv) > nbgBuff && nbinscount > 2) {
          nbinscount <- nbinscount - 1
          BuffPointsEnv <- VarelaSample(RandomBuff[, 1:2], envdata, nbinscount, PCA, PCAxes)
        }

        BuffPointsEnv2 <- VarelaSample(RandomBuff[, 1:2], envdata, (nbinscount + 1), PCA, PCAxes)

        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)

        if (diff2 < diff1) {
          BuffPointsEnv <- BuffPointsEnv2
        }

      } else {

        while(nrow(BuffPointsEnv) < nbgBuff && nbinscount < max(nbins, 99)) {
          nbinscount <- nbinscount + 1
          BuffPointsEnv <- VarelaSample(RandomBuff[, 1:2], envdata, nbinscount, PCA, PCAxes)
        }

        BuffPointsEnv2 <- VarelaSample(RandomBuff[, 1:2], envdata, (nbinscount - 1), PCA, PCAxes)

        diff1 <- abs(nrow(BuffPointsEnv) - nbgBuff)
        diff2 <- abs(nrow(BuffPointsEnv2) - nbgBuff)

        if (diff2 < diff1) {
          BuffPointsEnv <- BuffPointsEnv2
        }

      }

      BuffPointsEnv <- BuffPointsEnv2[complete.cases(BuffPointsEnv2), ]

    } else if (method == "random") {
      BuffPoints <- sp::spsample(bufferSHP, nbgBuff, "random", iter = 15)
      BuffPointsCoords <- BuffPoints@coords
      BuffPointsEnv <- raster::extract(envdata, BuffPointsCoords)
      BuffPointsEnv <- cbind(BuffPointsCoords, BuffPointsEnv)
    }

    background <- rbind(FullPointsEnv, BuffPointsEnv)

    if (method == "Varela") {
      write.csv(data.frame("Species" = spplist[s],
                           "FullTrainPoints" = nrow(FullPointsEnv),
                           "BuffPoints" = nrow(BuffPointsEnv),
                           "FullPoints" = nrow(background)),
                paste0(StatsLoc, "/", spplist[s], "_stats.csv"), row.names = FALSE)
    }

    write.csv(background, paste0(output, "/", sub(" ", "_", spplist[s]), ".csv"), row.names = FALSE)
  }
}
