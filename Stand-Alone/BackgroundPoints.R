#Background Points

#nbg: number of pseudo-absence (background) points wanted
  #note: for method="Varela", this number will not be exact
#envdata: RasterLayer or RasterStack object of environmental data
#buffer (optional): shapefile polygon
#method: How should background points be generated:
  #Varela
  #random
#nbins (if method = Varela): provides a starting point for the Varela method to optimize for the number of points wanted (defaults to 10)
#PCA (if method = Varela): should climate data be binned by Principal Component Axes? 
  #If given an integer (n), will use n PC Axes for Varela subsampling
  #If given "Y", will use the number of axes required to explain 95% of the variance
  #If not given, the actual climate layers will be used for Varela subsampling


BackgroundPoints <- function(nbg, envdata, buffer, method, nbins, PCA) {
  require(raster)
  require(sp)
  
  size <- 10
  
  if (method == "random") {
    
    if (hasArg(buffer)) {
      BufferPoints <- spsample(buffer, nbg / 2, "random")
      BufferCoords <- BufferPoints@coords
      BufferEnv <- extract(envdata, BufferCoords)
      BufferPointsEnv <- cbind(BufferCoords, BufferEnv)
    } else {
      message("No buffer mask found, skipping spatially constrained background point selection")
    }
    
    if (nbg > ncell(envdata)) {
      RastPoints <- ncell(envdata)
    } else {
      RastPoints <- nbg
    }
    
    if (hasArg(buffer)) {
      FullAreaPointsEnv <- sampleRandom(envdata, RastPoints / 2, na.rm = TRUE, xy = TRUE)
    } else {
      FullAreaPointsEnv <- sampleRandom(envdata, RastPoints, na.rm = TRUE, xy = TRUE)
    }
    
    if (hasArg(buffer)) {
      background <- rbind(BufferPointsEnv, FullAreaPointsEnv)
    } else {
      background <- FullAreaPointsEnv
    }
    return(background)
    
  } else if (method == "Varela") {
    VarelaSample <- function (occurrences, env, no_bins, PCA, PCAxes) {
      require(dplyr)
      occurrences <- SpatialPoints(occurrences, proj4string = crs(env))
      occurrences <- spTransform(occurrences, crs(env))
      
      EnvOccur <- extract(env, occurrences)
      EnvOccur <- data.frame(x = occurrences@coords[, 1], y = occurrences@coords[, 2], EnvOccur)
      ClimOccur <- EnvOccur
      ClimOccur <- ClimOccur[complete.cases(ClimOccur), ]
      if (PCA == "Y") {
        PCAEnv <- prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
        PCAImp <- summary(PCAEnv)$importance
        #Determine the number of PC axes to use for subsampling
        if (!is.na(PCAxes)) {
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
        sub_ptz <- distinct(out_ptz[, -1:-2])
        #count the number of filled bins
        no_grps <- nrow(sub_ptz)
        #add a column with the group membership number; this number is arbitrary
        sub_ptz$grp <- c(1:no_grps)
        
        #join the out_ptz with the subsample to capture the group membership info
        #note: join() will automatically match the variable names from these two dfs
        out_ptz <- left_join(out_ptz, sub_ptz)
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
        nsamples=c(nsamples, nrow(final_out))
      }
      if (length(no_bins) == 1) {
        return(final_out)
      } else {
        plot(nsamples ~ no_bins, pch = 16, xlab = "Number of Climate Bins", ylab = "Number of Subsampled Points")
        return(data.frame(NumberofSamples = nsamples, NumberOfClimateBins = no_bins))
      }
    }
    
    if (nbg > ncell(envdata)) {
      nbg <- ncell(envdata)
    } else {
      nbg <- nbg
    }
    
    if (hasArg(nbins)){
      nbins <- nbins
    } else {
      nbins <- 10
    }
    
    if (hasArg(PCA)){
      PCAStep <- "Y"
      if (PCA != "Y") {
        PCAxes <- as.numeric(PCA)
        if (PCAxes > nlayers(envdata)) {
          PCAxes <- nlayers(envdata)
          message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
        }
      } else {
        PCAxes <- NA
      }
    } else {
      PCAStep <- "N"
      PCAxes <- ""
    }
    
    if (hasArg(buffer)) {
      BufferPoints <- spsample(buffer, nbg * (size / 2), "random", iter = 15)
      BufferCoords <- BufferPoints@coords
      BufferPointsEnv <- VarelaSample(BufferCoords, envdata, nbins, PCAStep, PCAxes)
      nbinscount <- nbins
      if (nrow(BufferPointsEnv) >= nbg / 2) {
        while(nrow(BufferPointsEnv) > nbg / 2 && nbinscount > 2) {
          nbinscount <- nbinscount - 1
          BufferPointsEnv <- VarelaSample(BufferCoords, envdata, nbinscount, PCAStep, PCAxes)
        }
        BufferPointsEnv2 <- VarelaSample(BufferCoords, envdata, nbinscount + 1, PCAStep, PCAxes)
        if (abs(nrow(BufferPointsEnv[complete.cases(BufferPointsEnv), ]) - nbg / 2) > abs(nrow(BufferPointsEnv2[complete.cases(BufferPointsEnv2), ]) - nbg / 2)) {
          BufferPointsEnv <- BufferPointsEnv2
        } else {
          BufferPointsEnv <- BufferPointsEnv
        }
      } else {
        while(nrow(BufferPointsEnv) < nbg / 2 && nbinscount < max(nbins, 99)) {
          nbinscount <- nbinscount + 1
          BufferPointsEnv <- VarelaSample(BufferCoords, envdata, nbinscount, PCAStep, PCAxes)
        }
        BufferPointsEnv2 <- VarelaSample(BufferCoords, envdata, nbinscount - 1, PCAStep, PCAxes)
        if (abs(nrow(BufferPointsEnv[complete.cases(BufferPointsEnv), ]) - nbg / 2) > abs(nrow(BufferPointsEnv2[complete.cases(BufferPointsEnv2), ]) - nbg / 2)) {
          BufferPointsEnv <- BufferPointsEnv2
        } else {
          BufferPointsEnv <- BufferPointsEnv
        }
      }
      BufferPointsEnv <- BufferPointsEnv[complete.cases(BufferPointsEnv), ]
      message(paste0("Number of Points in Buffer: ", nrow(BufferPointsEnv)))
    } else {
      message("No buffer mask found, skipping spatially constrained background point selection")
    }
    
    if (nbg * 10 > ncell(envdata[[1]])) {
      RastTrainBG <- ncell(envdata[[1]])
    } else {
      RastTrainBG <- nbg * 10
    }
    
    RandomTrain <- sampleRandom(envdata, RastTrainBG, na.rm = TRUE, xy = TRUE)
    
    if (hasArg(buffer)) {
      FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
      if (nrow(FullAreaPointsEnv) >= nbg / 2) {
        while(nrow(FullAreaPointsEnv) > nbg / 2 && nbinscount > 2) {
          nbinscount <- nbinscount - 1
          FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
        }
        FullAreaPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount + 1, PCAStep, PCAxes)
        if (abs(nrow(FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]) - nbg / 2) > abs(nrow(FullAreaPointsEnv2[complete.cases(FullAreaPointsEnv2), ]) - nbg / 2)) {
          FullAreaPointsEnv <- FullAreaPointsEnv2
        } else {
          FullAreaPointsEnv <- FullAreaPointsEnv
        }
      } else {
        while(nrow(FullAreaPointsEnv) < nbg / 2 && nbinscount < max(nbins, 99)) {
          nbinscount <- nbinscount + 1
         FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
        }
        FullAreaPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount - 1, PCAStep, PCAxes)
        if (abs(nrow(FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]) - nbg / 2) > abs(nrow(FullAreaPointsEnv2[complete.cases(FullAreaPointsEnv2), ]) - nbg / 2)) {
          FullAreaPointsEnv <- FullAreaPointsEnv2
        } else {
          FullAreaPointsEnv <- FullAreaPointsEnv
        }
      }
      FullAreaPointsEnv <- FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]
      message(paste0("Number of Points in Full Area: ", nrow(FullAreaPointsEnv)))
    } else {
      nbinscount <- nbins
      FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
      if (nrow(FullAreaPointsEnv) >= nbg) {
        while(nrow(FullAreaPointsEnv) > nbg && nbinscount > 2) {
          nbinscount <- nbinscount - 1
          FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
        }
        FullAreaPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount + 1, PCAStep, PCAxes)
        if (abs(nrow(FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]) - nbg) > abs(nrow(FullAreaPointsEnv2[complete.cases(FullAreaPointsEnv2), ]) - nbg)) {
          FullAreaPointsEnv <- FullAreaPointsEnv2
        } else {
          FullAreaPointsEnv <- FullAreaPointsEnv
        }
      } else {
        while(nrow(FullAreaPointsEnv) < nbg && nbinscount < max(nbins, 99)) {
          nbinscount <- nbinscount + 1
          FullAreaPointsEnv <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount, PCAStep, PCAxes)
        }
        FullAreaPointsEnv2 <- VarelaSample(RandomTrain[, 1:2], envdata, nbinscount - 1, PCAStep, PCAxes)
        if (abs(nrow(FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]) - nbg) > abs(nrow(FullAreaPointsEnv2[complete.cases(FullAreaPointsEnv2), ]) - nbg)) {
          FullAreaPointsEnv <- FullAreaPointsEnv2
        } else {
          FullAreaPointsEnv <- FullAreaPointsEnv
        }
      }
      FullAreaPointsEnv <- FullAreaPointsEnv[complete.cases(FullAreaPointsEnv), ]
      message(paste0("Number of Points: ", nrow(FullAreaPointsEnv)))
    }

    if (hasArg(buffer)) {
      background <- rbind(BufferPointsEnv, FullAreaPointsEnv)
    } else {
      background <- rbind(FullAreaPointsEnv)
    }
    return(background)
  }
}
