% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrainStudyEnv.R
\name{TrainStudyEnv}
\alias{TrainStudyEnv}
\title{Project/clip training and study environmental layers}
\usage{
TrainStudyEnv(
  input_TA,
  input_SA,
  desiredCRS = NA,
  resolution = NA,
  clipTrain = NA,
  clipStudy = NA,
  output = NA,
  maxentproj = TRUE
)
}
\arguments{
\item{input_TA}{SpatRaster or a list of raster files for the training area (full directory path).}

\item{input_SA}{(optional) SpatRaster or a list of raster files for the study area. Provide if
the study area cannot be clipped from the trainingarea. Defaults to NA (same list of raster files).}

\item{desiredCRS}{The coordinate system to project training/test data into given in PROJ4 notation.
Defaults to NA (no projection).}

\item{resolution}{The desired resolution of the raster data in the units of the \strong{target projection}.
NOTE: if CRS is NA, but resolution has a value, the rasters will be resampled to the given resolution.
NOTE: different x and y resolutions (i.e., rectangular pixels) are acceptable for megaSDM in general but MaxEntProj cannot be run.}

\item{clipTrain}{Extent object or vector of desired training extent in form c(xmin, xmax, ymin, ymax).
The extent should be given in latlong coordinates.}

\item{clipStudy}{Extent object or vector of desired study extent in form c(xmin, xmax, ymin, ymax).
The extent should be given in latlong coordinates. Alternatively, if \code{clipStudy = "train"}, the
training area extent will be used.}

\item{output}{If the rasters are to be written to the computer, the full path of the directory where
they will be written out to. If set to \code{NA} (the default), the rasters will not be written out
and will be returned as the value of this function.}

\item{maxentproj}{TRUE/FALSE: Will the MaxEntProj step be run on these data? If so, rectangular pixels
without a defined resolution will be resampled to square pixels using the longer of the two sides. In
addition, the NA value of the raster will be set to the maximum raster value + 0.01 instead of the
nun-numeric NaN when written out.}
}
\value{
Returns a list of two SpatRasters: training "(\code{$training})" and study area "(\code{$study})"
environmental layers. if \code{output != NA}, the rasters will also be written out as ".grd" files.
}
\description{
This function manages environmental layers for the training and study areas of an SDM analysis.
The training area is the area where the model will be trained on (i.e., where the occurrence and
background points for the model generation are located). The study area is the region of interest
(i.e., where the parameters of the model will be applied and habitat suitability will be predicted
for both the current time period and future/past time periods). This function takes a list or stack
of rasters, projects and clips them to the desired projection, resolution, and extent. Next, if a
smaller study area is required than will be trained on, it clips the study area from the training
rasters and outputs a list containing both sets of rasters.
}
