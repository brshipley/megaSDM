#' Create buffers for spatially-constrained background point generation
#'
#' This function takes a list of occurrence point files and generates buffer shapefiles
#' around each set of points. These buffers will be used if spatially-constrained background
#' points are required. The radius of the buffer can be defined as a single value for all
#' species or as a distinct value for each species. If no radius values are given, the
#' distances between the occurrence points themselves inform the buffer radius (see \code{buff_distance}).
#'
#' @param occlist a list of .csv files containing the occurrence points for all species
#' intending to be modelled. Either the complete occurrence data or the subsampled data
#' may be used although results may vary depending on which data set is used. The files
#' should be named the same as the taxon of interest (e.g.,: ".../Canis_lupus.csv"), and
#' should be in the same coordinate reference system as the projected environmental data
#' (see \code{megaSDM::OccurrenceManagement}).
#' @param envdata a SpatRaster or list of raster files corresponding to the training
#' region (the region the SDM model(s) will be trained on). This provides the CRS
#' arguments to project the created buffer shapefiles.
#' @param output the full path to the directory where the buffer shapefiles will be written
#' out to. It is recommended that this is a unique directory from the one background
#' points will be written to.
#' @param buff_distance how wide should the buffers be around the occurrence points? If a single
#' number is given, buffers for all species will have a radius of that number in map units.
#' Buffers for longlat projection should be given in meters. If you are unsure what spatial units the
#' data have, run \code{terra::linearUnits()} on \code{envdata}, which will give 0 if units are in degrees,
#' or a number corresponding to the number of meters within one map unit (e.g., map units of km will have 1000).
#' If a vector of the same length as \code{occlist} is given, the buffer radius for each species will
#' correspond to the matching value in the vector. Finally, if set to \code{NA} (default),
#' 2*the 95% quantile of the minimum distance between each point is taken as the radius.
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @export
#' @return shapefiles (.shp) that are buffers around the occurrence points of each species
#' are written out to the directory indicated by the \code{output} argument. Also prints a .csv file for each
#' species containing the width of the buffer in the projection units.
#'

BackgroundBuffers <- function(occlist,
                              envdata,
                              output,
                              buff_distance = NA,
                              ncores = 1) {

  if(!dir.exists(output)) {
    dir.create(output)
  }

  ListSpp <- c()
  for (i in 1:length(occlist)) {
    CurSpp <- occlist[i]
    SpeciesSplit <- unlist(strsplit(CurSpp, "/", fixed = TRUE))
    SpeciesName <- substr(SpeciesSplit[length(SpeciesSplit)], 1,
                          nchar(SpeciesSplit[length(SpeciesSplit)]) - 4)
    ListSpp <- c(ListSpp, SpeciesName)
  }
  if (length(ListSpp) < ncores) {
    ncores <- length(ListSpp)
  }
  ListSpp <- matrix(ListSpp, ncol = ncores)

  #If the input training layers are not in SpatRaster form, ensure that they have the same projection/extent
  if (class(envdata) != "SpatRaster") {

    #Ensure that all training area rasters have the same projection and extent
    projtrain <- rep(NA, len = length(envdata))
    exttrain <- rep(NA, len = length(envdata))
    for (i in 1:length(envdata)) {
      projtrain[i] <- as.character(terra::crs(terra::rast(envdata[[i]])))
      exttrain[i] <- as.character(terra::ext(terra::rast(envdata[[i]])))
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

  #Convert envstack to a "packed" raster so that it can be parallelized
  envstack <- terra::wrap(envstack)

  run <- function(CurSpp) {
    envstack <- terra::rast(envstack)
    species <- occlist[grep(paste0(CurSpp, ".csv"), occlist)]
    species <- utils::read.csv(species)

    #Converts species occurrence points into SpatialPoints with desiredCRS
    coordinates <- species[, c("x", "y")]
    Coordinates2 <- terra::vect(coordinates, geom = c("x", "y"), crs = terra::crs(envstack))

    if(is.na(terra::linearUnits(Coordinates2))) {
      stop("The spatial units of the data cannot be found!
           Choose a different coordinate system for all data or add units into the existing one.")
    }

    #Calculates a buffer width based on quantiled minimum distance among all occurrence points
    Locations <- as.data.frame(terra::geom(Coordinates2)[, c("x", "y")])

    if (nrow(Locations) > 5000) {
      message(paste0("Warning: the number of occurrences for ", CurSpp, " is > 5000"))
      message(paste0("    Generating backgrond buffers may take some time."))
    }

    #If different sizes of buffers are wanted for the different species, pull the one that matches the current species
    if (length(buff_distance) == 1) {
      cur_dist <- buff_distance
    } else {
      cur_dist <- buff_distance[grep(paste0(CurSpp, ".csv"), occlist)]
    }

    if (is.na(cur_dist)) {
      #Uses the 95% quantile of the minimum distance between each point
      Distance <- terra::distance(Coordinates2)
      Distance <- as.matrix(Distance)
      mindist <- c()
      for (q in 1:ncol(Distance)) {
        DistanceZero <- Distance[which(Distance[, q] > 0), q]
        mindist <- c(mindist, min(DistanceZero))
      }
      BufferWidth <- 2 * stats::quantile(mindist, 0.95)
      rm(Distance, mindist)
      gc()
    } else {
      BufferWidth <- cur_dist
    }

    #if the data have units that are not degrees or meters, convert buffer widths to meters.
    if(terra::linearUnits(Coordinates2) != 0) {
      unitmult <- terra::linearUnits(Coordinates2)
      BufferWidth <- BufferWidth / unitmult
    }

    #Makes a buffer (width dependent) around the first occurrence point
    combinedPolygon <- terra::buffer(Coordinates2, BufferWidth)
    combinedPolygon <- terra::aggregate(combinedPolygon)

    #Reprojects the buffer into the desired CRS
    combinedPolygon <- terra::project(combinedPolygon, terra::crs(envstack))

    #Write shapefile out of the buffer
    terra::writeVector(combinedPolygon, filename = file.path(output, paste0(CurSpp, ".shp")), filetype = "ESRI Shapefile", overwrite = TRUE)
    rm(combinedPolygon)
    gc()

    #Fill out background point stats table with information on the Buffer Width (m)
    BGPStats <- BufferWidth
    dir.create(file.path(output, CurSpp))
    utils::write.csv(BGPStats, file = file.path(output, CurSpp, "BackgroundPoints_stats.csv"))
  }

  if (ncores == 1) {
    ListSpp <- as.vector(ListSpp)
    out <- sapply(ListSpp, function(x) run(x))
  } else {
    clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)

    parallel::clusterExport(clus, varlist = c("occlist", "run", "envstack",
                                              "output", "buff_distance", "ListSpp"), envir = environment())
    parallel::clusterEvalQ(clus, library(terra))

    for (i in 1:nrow(ListSpp)) {
      out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
    }

    parallel::stopCluster(clus)
  }

}
