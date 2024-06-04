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
#' @param input_TA SpatRaster or a list of raster files for the training area (full directory path).
#' @param input_SA (optional) SpatRaster or a list of raster files for the study area. Provide if
#' the study area cannot be clipped from the trainingarea. Defaults to NA (same list of raster files).
#' @param desiredCRS The coordinate system to project training/test data into given in PROJ4 notation.
#' Defaults to NA (no projection).
#' @param resolution The desired resolution of the raster data in the units of the __target projection__.
#' NOTE: if CRS is NA, but resolution has a value, the rasters will be resampled to the given resolution.
#' NOTE: different x and y resolutions (i.e., rectangular pixels) are acceptable for megaSDM in general but MaxEntProj cannot be run.
#' @param clipTrain Extent object or vector of desired training extent in form c(xmin, xmax, ymin, ymax).
#' The extent should be given in latlong coordinates.
#' @param clipStudy Extent object or vector of desired study extent in form c(xmin, xmax, ymin, ymax).
#' The extent should be given in latlong coordinates. Alternatively, if \code{clipStudy = "train"}, the
#' training area extent will be used.
#' @param output If the rasters are to be written to the computer, the full path of the directory where
#' they will be written out to. If set to \code{NA} (the default), the rasters will not be written out
#' and will be returned as the value of this function.
#' @param maxentproj TRUE/FALSE: Will the MaxEntProj step be run on these data? If so, rectangular pixels
#' without a defined resolution will be resampled to square pixels using the longer of the two sides. In
#' addition, the NA value of the raster will be set to the maximum raster value + 0.01 instead of the
#' nun-numeric NaN when written out.
#' @export
#' @return Returns a list of two SpatRasters: training "(\code{$training})" and study area "(\code{$study})"
#' environmental layers. if \code{output != NA}, the rasters will also be written out as ".grd" files.

TrainStudyEnv <- function(input_TA, input_SA, desiredCRS = NA,
                          resolution = NA, clipTrain = NA,
                          clipStudy = NA, output = NA, maxentproj = TRUE) {

  #If the input training layers are not in SpatRaster form, ensure that they have the same projection/extent
  if (class(input_TA) != "SpatRaster") {

    #Ensure that all training area rasters have the same projection and extent
    projtrain <- rep(NA, len = length(input_TA))
    exttrain <- rep(NA, len = length(input_TA))
    for (i in 1:length(input_TA)) {
      projtrain[i] <- as.character(terra::crs(terra::rast(input_TA[[i]])))
      exttrain[i] <- as.character(terra::ext(terra::rast(input_TA[[i]])))
    }

    if (length(unique(projtrain)) > 1) {
      stop("Not all of the training area environmental rasters are in the same projection")
    } else if (length(unique(exttrain)) > 1) {
      stop("Not all of the training area environmental rasters have the same extent")
    }
    envstack <- terra::rast(input_TA)
  } else {
    envstack <- input_TA
  }

  #Make sure that the training layers have a CRS that is not NA
  if (is.na(terra::crs(envstack))) {
    stop("training area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }

  #If study area rasters are provided, repeat for study area
  if (methods::hasArg(input_SA)) {
    #If the input study area layers are not in SpatRaster form, ensure that they have the same projection/extent
    if (class(input_SA) != "SpatRaster") {

      #Ensure that all training area rasters have the same projection and extent
      projstudy <- rep(NA, len = length(input_SA))
      extstudy <- rep(NA, len = length(input_SA))
      for (i in 1:length(input_SA)) {
        projstudy[i] <- as.character(terra::crs(terra::rast(input_SA[[i]])))
        extstudy[i] <- as.character(terra::ext(terra::rast(input_SA[[i]])))
      }

      if (length(unique(projstudy)) > 1) {
        stop("Not all of the study area environmental rasters are in the same projection")
      } else if (length(unique(extstudy)) > 1) {
        stop("Not all of the study area environmental rasters have the same extent")
      }
    }

    #Make sure that the study area layers have a CRS that is not NA
    studystack <- terra::rast(input_SA)
    if (is.na(terra::crs(studystack))) {
      stop("study area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
    } else if (as.character(terra::crs(studystack)) != as.character(terra::crs(envstack))) {
      studystack <- terra::project(studystack, terra::crs(envstack))
      message("Warning: study area has different coordinate projection than training area: reprojecting to training area CRS")
    }

  }

  #If pixels are rectangular and maxentprojection needs to be run, resample to coarser resolution
  #Also, force resample to avoid tiny resolution errors (see GitHub Issue #2)
  if (maxentproj) {
    #Set desired resolution to the longest axis of the pixels
    #If projection to a new CRS is required, first project, then get resolution

    if (is.na(max(resolution)) & is.na(desiredCRS)) {
      resolution <- max(terra::res(envstack))
    } else if (is.na(max(resolution))) {
      TestProject <- terra::project(envstack[[1]], desiredCRS)
      resolution <- max(terra::res(TestProject))
    }
    if(length(resolution) > 1) {
      resolution <- max(resolution)
      message("The desired resolution will lead to rectangular pixels: resampling to square pixels for use in MaxEnt projection")
    }
    if(!identical(terra::res(envstack)[1], terra::res(envstack)[2])) {
      resolution <- max(terra::res(envstack))
      message("The original training raster pixels are rectangular: resampling to square pixels for use in MaxEnt projection")
    }
    if(exists("studystack")) {
      if(!identical(terra::res(studystack)[1], terra::res(studystack)[2])) {
        resolution <- max(c(terra::res(envstack), terra::res(studystack)))
        message("The original study area raster pixels are rectangular: resampling to square pixels for use in MaxEnt projection")
      }
    }
  }

  #Clip rasters to desired extent
  #Having the clipping step first makes the process faster
  if (class(clipTrain) != "logical") {
    clipTrain <- terra::ext(clipTrain)
    ExtentTA <- terra::as.points(clipTrain)
    terra::crs(ExtentTA) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    ExtentProjTA <- terra::project(ExtentTA, terra::crs(envstack))

    NewExtentTA <- terra::ext(ExtentProjTA)

    #Weird projections cause the reprojected layers to have issues (dealt with here)
    if (NewExtentTA[1] > NewExtentTA[2] | NewExtentTA[3] > NewExtentTA[4]) {
      stop("The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
    }

    if (NewExtentTA > terra::ext(envstack)) {
      message("Warning: the desired training extent is larger than the original raster")
      message("Only the intersection of the two extents will be used")
    }

    envstack <- terra::crop(envstack, terra::ext(max(NewExtentTA[1], terra::ext(envstack)[1]),
                                                      min(NewExtentTA[2], terra::ext(envstack)[2]),
                                                      max(NewExtentTA[3], terra::ext(envstack)[3]),
                                                      min(NewExtentTA[4], terra::ext(envstack)[4])))
    print("Training area clipped!")
  }

  if (!exists("studystack")) {
    studystack <- envstack
  }

  if (class(clipStudy) != "logical") {
    clipStudy <- terra::ext(clipStudy)
    ExtentSA <- terra::as.points(clipStudy)
    terra::crs(ExtentSA) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    ExtentProjSA <- terra::project(ExtentSA, terra::crs(envstack))

    NewExtentSA <- terra::ext(ExtentProjSA)

    if (NewExtentSA > terra::ext(studystack)) {
      message("Warning: the desired study extent is larger than the study rasters")
      message("Only the intersection of the two extents will be used")
    }

    #Weird projections cause the reprojected layers to have issues (dealt with here)
    if (NewExtentSA[1] > NewExtentSA[2] | NewExtentSA[3] > NewExtentSA[4]) {
      stop("The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
    }

    studystack <- terra::crop(envstack, terra::ext(max(NewExtentSA[1], terra::ext(studystack)[1]),
                                                        min(NewExtentSA[2], terra::ext(studystack)[2]),
                                                        max(NewExtentSA[3], terra::ext(studystack)[3]),
                                                        min(NewExtentSA[4], terra::ext(studystack)[4])))
    print("Study area clipped!")
  }

  #Reproject and/or resample rasters to desired crs and resolution
  #TODO allow for "method" argument of resampling to be vared based on layer
  if (!is.na(desiredCRS)) {
    if (is.na(resolution)) {
      envstack <- terra::project(envstack, desiredCRS)
      gc(full = TRUE)
      if (exists("studystack")) {
        studystack <- terra::project(studystack, desiredCRS)
        gc(full = TRUE)
      }

    } else {
      envstack <- terra::project(envstack, desiredCRS)
      envstack_res <- terra::rast(extent = terra::ext(envstack), resolution = resolution, crs = terra::crs(envstack))
      envstack <- terra::resample(envstack, envstack_res, method = "bilinear")
      gc(full = TRUE)
      if (exists("studystack")) {
        studystack <- terra::project(studystack, desiredCRS)
        studystack_res <- terra::rast(extent = terra::ext(studystack), resolution = resolution, crs = terra::crs(studystack))
        studystack <- terra::resample(studystack, studystack_res, method = "bilinear")
        gc(full = TRUE)
      }
    }

  } else if (!is.na(resolution)) {
    if (exists("studystack")) {
      studystack_res <- terra::rast(extent = terra::ext(studystack), resolution = resolution, crs = terra::crs(studystack))
      studystack <- terra::resample(studystack, studystack_res, method = "bilinear")
      gc(full = TRUE)
    }
    envstack_res <- terra::rast(extent = terra::ext(envstack), resolution = resolution, crs = terra::crs(envstack))
    envstack <- terra::resample(envstack, envstack_res, method = "bilinear")
    gc(full = TRUE)
  }

  envstack <- terra::trim(envstack)
  studystack <- terra::trim(studystack)
  EnvRasters <- list("training" = envstack, "study" = studystack)

  if (!is.na(output)) {
    if (!dir.exists(paste0(output, "/trainingarea"))) {
      dir.create(paste0(output, "/trainingarea"), recursive = TRUE)
    }
    if (!dir.exists(paste0(output, "/studyarea"))) {
      dir.create(paste0(output, "/studyarea"), recursive = TRUE)
    }
    #If maxent projection is necessary: change NA value of each training raster to maximum value + 0.01
    if (maxentproj) {
      for(e in 1:terra::nlyr(EnvRasters$training)) {
        FocusRast <- EnvRasters$training[[e]]
        MaxValue <- max(terra::values(FocusRast), na.rm = TRUE)
        FocusRast[which(is.na(terra::values(FocusRast)))] <- as.numeric(MaxValue + 0.01)
        terra::writeRaster(FocusRast,
                           filename = paste0(output, "/trainingarea/", names(FocusRast), ".grd"),
                           overwrite = TRUE, NAflag = as.numeric(MaxValue + 0.01))
      }

      for(e in 1:terra::nlyr(EnvRasters$study)) {
        FocusRast <- EnvRasters$study[[e]]
        MaxValue <- max(terra::values(FocusRast), na.rm = TRUE)
        FocusRast[which(is.na(terra::values(FocusRast)))] <- as.numeric(MaxValue + 0.01)
        terra::writeRaster(FocusRast,
                           filename = paste0(output, "/studyarea/", names(FocusRast), ".grd"),
                           overwrite = TRUE, NAflag = as.numeric(MaxValue + 0.01))
      }

    } else {

      terra::writeRaster(EnvRasters$training, paste0(output, "/trainingarea/", names(EnvRasters$training), ".grd"),
                         overwrite = TRUE)
      terra::writeRaster(EnvRasters$study, paste0(output, "/studyarea/", names(EnvRasters$study), ".grd"),
                         overwrite = TRUE)
    }

  }
  return(EnvRasters)
}
