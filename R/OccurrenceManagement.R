#' Manage and environmentally filter occurrence points
#'
#' This function takes a set of occurrence points (whether downloaded from
#' GBIF using \code{OccurrenceCollection} or provided by the user),
#' standardizes the column headings for effective use in species distribution
#' modelling, and, if requested, extracts the values of each environmental
#' variable used in the modelling for each point and environmentally subsamples
#' the data using the method developed by Varela et al. (2014).
#'
#' @param occlist a list of .csv file names, each containing the occurrence
#' points of a given taxon. The files should be named the same as the taxon
#' of interest (e.g.,: ".../Canis_lupus.csv"). The occurrence points should
#' be given in latlong coordinates (-90 to 90, -180 to 180).
#' @param output the directory where the managed occurrence points will be placed.
#' @param envextract (logical (\code{TRUE} or \code{FALSE})) should the
#' environments at each occurrence point be extracted?
#' @param envsample (logical (\code{TRUE} or \code{FALSE})) should
#' environmental (Varela) subsampling be conducted on the occurrence points?
#' @param envdata a RasterStack or list of raster files
#' corresponding to the area the model will be trained on. All environmental
#' variables should be provided.
#' @param nbins (integer) the number of bins each climate combination
#' will be placed in for Varela subsampling of occurrence points. See
#' Varela et al. (2014) and Castellanos et al. (2018) for discussion on this
#' method. Default is 25 bins.
#' @param PCA (optional). The number of PC axes to use when binning climate
#' data for environmental subsampling. Can be a number from 1 to the number
#' of environmental variables. Default selects the number of axes that account
#' for 95% of the environmental variation. If set to \code{NA}, PCA will not be
#' run (for use when there are categorical variables, for instance).
#' @export
#' @return a set of .csv files of occurrence points in the directory indicated by
#' the \code{output} argument with columns relating to the species name, the coordinates
#' (in the projection of the training area rasters), and the value of each
#' environmental layer at each point.


OccurrenceManagement <- function(occlist,
                                 output,
                                 envextract = FALSE,
                                 envsample = FALSE,
                                 envdata,
                                 nbins = 25,
                                 PCA = "Y") {


  if(!dir.exists(output)) {
    dir.create(output)
  }

  #If the input training layers are not in rasterstack form, ensure that they have the same projection/extent
  if (class(envdata) != "RasterStack") {

    #Ensure that all training area rasters have the same projection and extent
    projtrain <- rep(NA, len = length(envdata))
    exttrain <- rep(NA, len = length(envdata))
    for (i in 1:length(envdata)) {
      projtrain[i] <- as.character(raster::crs(raster::raster(envdata[[i]])))
      exttrain[i] <- as.character(raster::extent(raster::raster(envdata[[i]])))
    }

    if (length(unique(projtrain)) > 1) {
      stop("Not all of the training area environmental rasters are in the same projection")
    } else if (length(unique(exttrain)) > 1) {
      stop("Not all of the training area environmental rasters have the same extent")
    }
  }

  #Make sure that the training layers have a CRS that is not NA
  envstack <- raster::stack(envdata)

  #Provide names for the raster layers in the study area raster stack
  if (class(envdata) != "RasterStack") {
    EnvNames <- rep(NA, length = length(envdata))
    for(i in 1:length(EnvNames)) {
      focname <- unlist(strsplit(envdata[i], "/"))
      focname <- focname[length(focname)]
      focname <- unlist(strsplit(focname, "\\."))[1]
      EnvNames[i] <- focname
    }

    names(envstack) <- EnvNames
  }

  if (is.na(raster::crs(envstack))) {
    stop("training area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }

  if (envsample == "TRUE") {

    if (is.numeric(PCA)) {
      PCAxes <- as.numeric(PCA)
      if (PCAxes > raster::nlayers(envdata)) {
        PCAxes <- raster::nlayers(envdata)
        message(paste0("Setting Number of PC Axes to Number of Layers: ", PCAxes))
      }
    } else if (!is.na(PCA)) {
      PCAxes <- "Y"
    } else {
      PCA <- "N"
      PCAxes <- ""
    }

    VarelaSample <- function (occurrences, env, no_bins, PCA, PCAxes) {

      occurrences <- sp::SpatialPoints(occurrences, proj4string = raster::crs(env))

      EnvOccur <- raster::extract(env, occurrences)
      EnvOccur <- data.frame(x = occurrences@coords[, 1], y = occurrences@coords[, 2], EnvOccur)
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
  }

  for (s in 1:length(occlist)) {
    #Reads the occurrence file
    SpeciesOcc <- utils::read.csv(occlist[s])

    #Determines the species name from the name of the occurrence file
    SpeciesSplit <- unlist(strsplit(occlist[s], "/", fixed = TRUE))
    SpeciesName <- substr(SpeciesSplit[length(SpeciesSplit)], 1,
                          nchar(SpeciesSplit[length(SpeciesSplit)]) - 4)

    #Matches the species name to the correct column and renames it "Species"
    SpecCol <- grep(paste0("^", SpeciesName, "$|", "^", gsub("_", " ", SpeciesName), "$"),
                    SpeciesOcc[1,])

    #Adds a new column with species names if it doesn't already exist
    if (length(SpecCol) == 1) {
      names(SpeciesOcc)[SpecCol] <- "Species"
    } else {
      SpeciesOcc$Species <- rep(SpeciesName, nrow(SpeciesOcc))
    }

    #Matches the lat/long columns and renames them to standardize
    names(SpeciesOcc)[c(grep("lon", tolower(names(SpeciesOcc))),
                        grep("^x$", tolower(names(SpeciesOcc))))] <- "Longitude"
    names(SpeciesOcc)[c(grep("lat", tolower(names(SpeciesOcc))),
                        grep("^y$", tolower(names(SpeciesOcc))))] <- "Latitude"

    #Takes only the species name and the coordinates for reprojection
    SpeciesCoords <- SpeciesOcc[, c("Species", "Longitude", "Latitude")]

    if ((max(SpeciesCoords$Longitude) > 180) | (max(SpeciesCoords$Latitude) > 90) |
       (min(SpeciesCoords$Longitude) < -180) | (min(SpeciesCoords$Latitude) < -90)) {
      stop("The species occurrence points are not given in latitude/longitude form. Reproject species occurrences to a latlon projection")
    }

    SpeciesCoordsSP <- sp::SpatialPointsDataFrame(coords = SpeciesCoords[, c("Longitude", "Latitude")],
                                           data = SpeciesCoords,
                                           proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
    SpeciesCoordsSP <- sp::spTransform(SpeciesCoordsSP, CRSobj = raster::crs(envstack))

    SpeciesCoords <- data.frame("Species" = SpeciesCoordsSP$Species, SpeciesCoordsSP@coords)

    #If required, extracts environmental data from rasters. Otherwise, adds environmental data back
    #to the projected occurrence points
    if (envextract == "FALSE") {
      for (l in 1:raster::nlayers(envstack)) {
        env_var <- names(envstack)[l]
        EnvCol <- grep(paste0("^", env_var, "$"), names(SpeciesOcc))
        if (length(EnvCol) == 0) {
          stop(paste0("Environmental layer ", env_var, " not found in occurrence file for ",
                      SpeciesName, ". Ensure that all environmental variables are provided in the occurrence file"))
        } else {
          if (l == 1) {
            SpeciesEnv <- SpeciesCoords
          } else {
            SpeciesEnv <- cbind(SpeciesEnv, SpeciesOcc[EnvCol])
          }
        }
      }
    } else {
      SpeciesEnv <- raster::extract(envstack, SpeciesCoordsSP)
      SpeciesEnv <- cbind(SpeciesCoords, SpeciesEnv)
      SpeciesEnv <- stats::na.omit(SpeciesEnv)
      if (nrow(SpeciesEnv) == 0) {
        stop("Environmental extraction failed:
             ensure that the points and the raster have overlapping extents")
      }
    }

    if (envsample == "TRUE") {
      SpeciesFinal <- VarelaSample(SpeciesEnv[, c("Longitude", "Latitude")], envdata, nbins, PCA, PCAxes)
      SpeciesFinal <- data.frame("Species" = rep(SpeciesName, nrow(SpeciesFinal)), SpeciesFinal)
    } else {
      SpeciesFinal <- SpeciesEnv
      names(SpeciesFinal)[1:3] <- c("Species", "x", "y")
      SpeciesFinal$Species <- gsub(" ", "_", SpeciesFinal$Species)
    }

    utils::write.csv(SpeciesFinal, paste0(output, "/", SpeciesName, ".csv"), row.names = FALSE)
  }
}
