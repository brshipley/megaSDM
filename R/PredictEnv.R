#' Project, clip, and store forecasted/hindcasted environmental rasters for SDM prediction
#'
#' This function takes lists of RasterStacks that correspond to future or past time periods
#' of a single climate model (e.g., RCP4.5, CCSM3), ensures that the environmental variables
#' are the same as those that the model will be developed on, and projects, clips, and resamples
#' these layers to the characteristics of a given study region. If more than one climate scenario
#' is required, run this function multiple times for each climate scenario.
#'
#' @param studylayers A single RasterStack or list of raster files that constitute all environmental
#' variables and parameters (e.g., extent, resolution) used for projecting the modelled relationship
#' (see \code{$study} in \code{megaSDM::TrainStudyEnv}).
#' @param futurelayers A list of rasterstacks or vectors of file-names corresponding to the environmental
#' variables at the different time periods the model will be forecasted/hindcasted to.
#' @param time_periods a vector of the time periods the models will be projected to, with the first element
#' as the year the model will be trained on (usually the current data). The projected time periods should be
#' given in the same order as the list given in \code{futurelayers}. If no precise years are available
#' (e.g., using data from the Last Glacial Maximum), order the \code{futurelayers} from current to least
#' current (farthest into the future/past) and give character strings for the years (e.g., "LGM"). If
#' running dispersal analyses, \code{time_periods} must be numeric (e.g., -21000 instead of "LGM").
#' @param output If the rasters are to be written to the computer, the full path of the directory where
#' they will be written out to. If there are multiple climate scenarios wanted in this SDM analysis,
#' give the output directory the name of the climate scenario (e.g., ".../output/RCP4.5). If set to \code{NA}
#' (the default), the rasters will not be written out and will be returned as the value of this function.
#' @param scenario_name (optional) If the rasters are to be written to the disk, a character string with the
#' name of the climate model/scenario. A sub-directory will be created within \code{output} and files will
#' be placed in there. If only one distinct climate scenario is needed, still give it a name for reference
#' in other functions within \code{megaSDM}.
#'
#' @export
#' @return Returns a list of environnmental layer RasterStacks with each RasterStack corresponding to a time period.

PredictEnv <- function(studylayers, futurelayers,
                       time_periods, output = NA,
                       scenario_name = NA) {

  if (class(studylayers) != "RasterStack") {

    #Ensure that all study area rasters have the same projection and extent
    projstudy <- rep(NA, len = length(studylayers))
    extstudy <- rep(NA, len = length(studylayers))
    for (i in 1:length(studylayers)) {
      projstudy[i] <- as.character(raster::crs(raster::raster(studylayers[[i]])))
      extstudy[i] <- as.character(raster::extent(raster::raster(studylayers[[i]])))
    }

    if (length(unique(projstudy)) > 1) {
      stop("Not all of the current study area environmental rasters are in the same projection")
    } else if (length(unique(extstudy)) > 1) {
      stop("Not all of the current study area environmental rasters have the same extent")
    }
  }

  studystack <- raster::stack(studylayers)
  
  #Provide names for the raster layers in the study area raster stack
  if (class(studylayers) != "RasterStack") {
    SANames <- rep(NA, length = length(studylayers))
    for(i in 1:length(SANames)) {
      focname <- unlist(strsplit(studylayers[i], "/"))
      focname <- focname[length(focname)]
      focname <- unlist(strsplit(focname, "\\."))[1]
      SANames[i] <- focname
    }
    
    names(studystack) <- SANames
  }
  
  if (is.na(raster::crs(studystack))) {
    stop("study area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }

  time_periods <- time_periods[2:length(time_periods)]
  if(length(time_periods) != length(futurelayers)) {
    stop("The number of time_periods given does not match the number of environmental layer sets")
  }

  for (j in 1:length(time_periods)) {
    focuslayers <- futurelayers[[j]]
    if (class(focuslayers) != "RasterStack") {
      #Ensure that all future/past rasters have the same projection and extent
      projstudy <- rep(NA, len = length(focuslayers))
      extstudy <- rep(NA, len = length(focuslayers))
      for (i in 1:length(focuslayers)) {
        projstudy[i] <- as.character(raster::crs(raster::raster(focuslayers[[i]])))
        extstudy[i] <- as.character(raster::extent(raster::raster(focuslayers[[i]])))
      }

      if (length(unique(projstudy)) > 1) {
        stop("Not all of the current study area environmental rasters are in the same projection")
      } else if (length(unique(extstudy)) > 1) {
        stop("Not all of the current study area environmental rasters have the same extent")
      }
    }

    focusstack <- raster::stack(focuslayers)
    
    #Provide names for the raster layers in the study area raster stack
    if (class(focuslayers) != "RasterStack") {
      FocusNames <- rep(NA, length = length(focuslayers))
      for(i in 1:length(FocusNames)) {
        focname <- unlist(strsplit(focuslayers[i], "/"))
        focname <- focname[length(focname)]
        focname <- unlist(strsplit(focname, "\\."))[1]
        FocusNames[i] <- focname
      }
      
      names(focusstack) <- FocusNames
    }
    
    
    if (is.na(raster::crs(focusstack))) {
      stop("study area raster crs = NA: Ensure all raster layers have a defined coordinate projection")
    }

    if (!setequal(names(focusstack), names(studystack))) {
      message("Warning: the environmental layer names do not match between current and future/past raster data")
    }

    if (raster::res(focusstack)[1] > raster::res(studystack)[1]) {
      message("Warning: the future/past raster data have coarser resolution than the current raster data")
      message(paste0("Changing the resolution of the current raster data to ", raster::res(focusstack)[1], "is recommended"))
    }

    if (as.character(raster::crs(focusstack)) != as.character(raster::crs(studystack))) {
      focusstack <- raster::projectRaster(focusstack,
                                           crs = raster::crs(studystack),
                                           res = raster::res(studystack),
                                           method = "bilinear")
    }

    if (raster::extent(focusstack) != raster::extent(studystack)) {
      focusstack <- raster::crop(focusstack, raster::extent(studystack))
    }

    if (!setequal(raster::res(focusstack), raster::res(studystack))) {
      focusstack <- raster::resample(focusstack, studystack, method = "bilinear")
    }

    if (j == 1) {
      PredictEnv <- focusstack
    } else {
      PredictEnv <- c(PredictEnv, focusstack)
    }

  }
  names(PredictEnv) <- time_periods

  if (!is.na(output)) {
    for (i in 1:length(PredictEnv)) {
      if (!is.na(scenario_name)) {
        if (!dir.exists(paste0(output, "/", scenario_name))) {
          dir.create(paste0(output, "/", scenario_name))
        }
      }
      time <- names(PredictEnv)[i]
      if (!dir.exists(paste0(output, "/", scenario_name ,"/", time))) {
        dir.create(paste0(output, "/", scenario_name, "/", time))
      }
      raster::writeRaster(PredictEnv[[i]], paste0(output, "/", scenario_name, "/", time, "/", names(PredictEnv[[i]]), ".grd"),
                          bylayer = TRUE, format = "raster", overwrite = TRUE)
    }
  }
  return(PredictEnv)
}
