#Varela Subsampling
#occurrences: a data frame or table with columns [long,lat] in decimal degree form
#env: A RasterLayer or RasterStack object with environmental characteristics
#nbins: number of bins to divide each environmental layer into for stratified subsampling
  #If nbins is set to a range: c() or seq(), will provde a graph of the number of subsampled points for a range of bin numbers
#PCA: should climate data be binned by Principal Component Axes? 
  #"N" = no
  #"Y" = yes
#PCAxes (if PCA = "Y"): A number detailing the quantity of PC Axes climate should be binned into 
  #if not given, the number of axes required to explain 95% of the variance will be calculated


VarelaSample <- function (occurrences, env, nbins, PCA, PCAxes) {
  require(dplyr)
  occurrences <- SpatialPoints(occurrences, proj4string = crs(env))
  #  occurrences <- SpatialPoints(occurrences, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  occurrences <- spTransform(occurrences, crs(env))
  
  EnvOccur <- extract(env, occurrences)
  # EnvOccur <- data.frame(x = occurrences@coords[, 1], y = occurrences@coords[, 2], EnvOccur)
  ClimOccur <- EnvOccur
  ClimOccur <- ClimOccur[complete.cases(ClimOccur), ]
  if (PCA == "Y") {
    PCAEnv <- prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
    PCAImp <- summary(PCAEnv)$importance
    #Determine the number of PC axes to use for subsampling
    if (hasArg(PCAxes)) {
      NumberAxes <- PCAxes
    } else {
      NumberAxes <- max(2, min(which(PCAImp[3, ] > 0.95)))
    }
    
    if (NumberAxes > nlayers(train)) {
      NumberAxes <- nlayers(train)
      message(paste0("Number of PC Axes will be changed to number of environmental layers: ", NumberAxes))
    }
    
    #Add PCA values to the unsubsampled data frame
    EnvOccur <- data.frame(cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes]))
  }
  
  #make a landing spot for the bin membership vectors
  
  #cycle through all of the environmental variables (columns 3 to end)
  nsamples <- c()
  for (j in 1:length(nbins)) {
    out_ptz <- EnvOccur[, 1:2]
    for(i in 3:length(names(EnvOccur))) {
      #make a data frame that is this variable with no NA values
      k <- EnvOccur[!is.na(EnvOccur[, i]), i]
      #calculate the observed range of this variable
      rg <- range(k)
      #figure out the resolution from the number of bins
      res <- (rg[2] - rg[1]) / nbins[j]
      #rescale the axis by the range and bin size, so the value is just a 
      #number from 1 to nbins for its bin membership
      d <- (EnvOccur[, i] - rg[1]) / res
      #d is now a vector of values ranging from 0 to nbins
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
    nsamples <- c(nsamples, nrow(final_out))
  }
  if (length(nbins) == 1) {
    return(final_out)
  } else {
    plot(nsamples ~ nbins, pch = 16, xlab = "Number of Climate Bins", ylab = "Number of Subsampled Points")
    return(data.frame(NumberofSamples = nsamples, NumberOfClimateBins = nbins))
  }
}