#' Project/clip training and study environmental layers
#'
#' This function manages environmental layers for the training and study areas of an SDM analysis.
#' The training area is the area where the model will be trained on (i.e., where the occurrence and
#' background points for the model generation are located). The study area is the region of interest
#' (i.e., where the parameters of the model will be applied and habitat suitability will be predicted
#' for both the current time period and future/past time periods). This function takes a list or stack
#' of rasters, projects and clips them to the desired projection, resolution, and extent. Next, if a
#' smaller study area is required than will be trained on, it clips the study area from the training
#' rasters and outputs a list containing both sets of rasters.
#'
#' @param input_TA RasterStack or a list of raster files for the training area (full directory path).
#' @param input_SA (optional) RasterStack or a list of raster files for the study area. Provide if
#' the study area cannot be clipped from the trainingarea. Defaults to NA (same list of raster files).
#' @param desiredCRS The coordinate system to project training/test data into given in PROJ4 notation.
#' Defaults to NA (no projection).
#' @param resolution The desired resolution of the raster data in the units of the __target projection__.
#'
#' NOTE: if CRS is NA, but resolution has a value, the rasters will be resampled to the given resolution.
#' @param clipTrain Extent object or vector of desired training extent in form c(xmin, xmax, ymin, ymax).
#' The extent should be given in latlong coordinates.
#' @param clipStudy Extent object or vector of desired study extent in form c(xmin, xmax, ymin, ymax).
#' The extent should be given in latlong coordinates. Alternatively, if \code{clipStudy = "train"}, the
#' training area extent will be used.
#' @param output If the rasters are to be written to the computer, the full path of the directory where
#' they will be written out to. If set to \code{NA} (the default), the rasters will not be written out
#' and will be returned as the value of this function.
#' @export
#' @return Returns a list of two rasterStacks: training "(\code{$training})" and study area "(\code{$study})"
#' environmental layers. if \code{output != NA}, the rasters will also be written out as ".bil" files.

TrainStudyEnv <- function(input_TA, input_SA, desiredCRS = NA,
                          resolution = NA, clipTrain = NA,
                          clipStudy = NA, output = NA) {

  #If the input training layers are not in rasterstack form, ensure that they have the same projection/extent
  if (class(input_TA) != "RasterStack") {

    #Ensure that all training area rasters have the same projection and extent
    projtrain <- rep(NA, len = length(input_TA))
    exttrain <- rep(NA, len = length(input_TA))
    for (i in 1:length(input_TA)) {
      projtrain[i] <- as.character(raster::crs(raster::raster(input_TA[[i]])))
      exttrain[i] <- as.character(raster::extent(raster::raster(input_TA[[i]])))
    }

    if (length(unique(projtrain)) > 1) {
      stop("Not all of the training area environmental rasters are in the same projection")
    } else if (length(unique(exttrain)) > 1) {
      stop("Not all of the training area environmental rasters have the same extent")
    }
  }

  #Make sure that the training layers have a CRS that is not NA
  envstack <- raster::stack(input_TA)
  if (is.na(raster::crs(envstack))) {
    stop("training area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }


  if (hasArg(input_SA)) {
    #If the input training layers are not in rasterstack form, ensure that they have the same projection/extent
    if (class(input_SA) != "RasterStack") {

      #Ensure that all training area rasters have the same projection and extent
      projstudy <- rep(NA, len = length(input_SA))
      extstudy <- rep(NA, len = length(input_SA))
      for (i in 1:length(input_SA)) {
        projstudy[i] <- as.character(raster::crs(raster::raster(input_SA[[i]])))
        extstudy[i] <- as.character(raster::extent(raster::raster(input_SA[[i]])))
      }

      if (length(unique(projstudy)) > 1) {
        stop("Not all of the training area environmental rasters are in the same projection")
      } else if (length(unique(extstudy)) > 1) {
        stop("Not all of the training area environmental rasters have the same extent")
      }
    }

    #Make sure that the training layers have a CRS that is not NA
    studystack <- raster::stack(input_SA)
    if (is.na(raster::crs(studystack))) {
      stop("training area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
    }
  }

  if (!is.na(desiredCRS)) {
    if (is.na(resolution)) {
      envstack <- raster::projectRaster(envstack, crs = desiredCRS)
    } else {
      envstack <- raster::projectRaster(envstack, crs = desiredCRS, res = resolution, method = "bilinear")
    }

  } else if (!is.na(resolution)) {
    envstack <- raster::projectRaster(envstack, crs = raster::crs(envstack), res = resolution, method = "bilinear")
  }

  if (class(clipTrain) != "logical") {
    clipTrain <- raster::extent(clipTrain)
    ExtentTA <- sp::SpatialPoints(clipTrain, proj4string = raster::crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
    ExtentProjTA <- sp::spTransform(ExtentTA, raster::crs(envstack))

    NewExtentTA <-raster::extent(ExtentProjTA@coords)

    #Weird projections cause the reprojected layers to have issues (dealt with here)
    if (NewExtentTA[1] > NewExtentTA[2] | NewExtentTA[3] > NewExtentTA[4]) {
      stop("The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
    }

    if (NewExtentTA > raster::extent(envstack)) {
      message("Warning: the desired training extent is larger than the original raster")
      message("Only the intersection of the two extents will be used")
    }

    envstack <- raster::crop(envstack, raster::extent(max(NewExtentTA[1], raster::extent(envstack)[1]),
                                              min(NewExtentTA[2], raster::extent(envstack)[2]),
                                              max(NewExtentTA[3], raster::extent(envstack)[3]),
                                              min(NewExtentTA[4], raster::extent(envstack)[4])))
    print("Training area complete!")
  }

  if (!exists("studystack")) {
    studystack <- envstack
  }
  if (class(clipStudy) != "logical") {
    clipStudy <- raster::extent(clipStudy)
    ExtentSA <- sp::SpatialPoints(clipStudy, proj4string = raster::crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
    ExtentProjSA <- sp::spTransform(ExtentSA, raster::crs(envstack))

    NewExtentSA <- raster::extent(ExtentProjSA@coords)

    if (NewExtentSA > raster::extent(studystack)) {
      message("Warning: the desired study extent is larger than the study rasters")
      message("Only the intersection of the two extents will be used")
    }

    #Weird projections cause the reprojected layers to have issues (dealt with here)
    if (NewExtentSA[1] > NewExtentSA[2] | NewExtentSA[3] > NewExtentSA[4]) {
      stop("The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
    }

    studystack <- raster::crop(envstack, raster::extent(max(NewExtentSA[1], raster::extent(studystack)[1]),
                                              min(NewExtentSA[2], raster::extent(studystack)[2]),
                                              max(NewExtentSA[3], raster::extent(studystack)[3]),
                                              min(NewExtentSA[4], raster::extent(studystack)[4])))
    print("Study area complete!")
  }

  if (!hasArg(input_SA) & (class(clipStudy) == "logical")) {
    studystack <- envstack
  }
  envstack <- raster::trim(envstack)
  studystack <- raster::trim(studystack)
  EnvRasters <- list("training" = raster::stack(envstack), "study" = raster::stack(studystack))

  if (!is.na(output)) {
    if (!dir.exists(paste0(output, "/trainingarea"))) {
      dir.create(paste0(output, "/trainingarea"))
    }
    if (!dir.exists(paste0(output, "/studyarea"))) {
      dir.create(paste0(output, "/studyarea"))
    }
    raster::writeRaster(EnvRasters$training, paste0(output, "/trainingarea/", names(EnvRasters$training), ".bil"),
                        bylayer = TRUE, format = "EHdr", overwrite = TRUE)
    raster::writeRaster(EnvRasters$study, paste0(output, "/studyarea/", names(EnvRasters$study), ".bil"),
                        bylayer = TRUE, format = "EHdr", overwrite = TRUE)
  }

  return(EnvRasters)
}
