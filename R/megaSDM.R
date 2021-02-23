#' megaSDM: Integrating Dispersal and Time-step Analyses into Species Distribution Models
#' 
#' @description
#' This package can efficiently create and project species distribution models using the
#' MaxEnt framework and parallel processing. It can find and download occurrence data 
#' for a list of species on GBIF (Global Biodiversity Information Facility), 
#' environmentally subsample the occurrences to mitigate spatial bias, generate background
#' (pseudo-absence) points, train the model and project it to different times (incorporating
#' dispersal rate of each species and intermediate range fluctuations), and create species
#' richness maps for each time period and taxon. 
#' 
#' ## Contributors
#' Benjamin Shipley, Renee Bach, Younje Do, Heather Strathearn, Jenny McGuire, Bistra Dilkina
#' 
#' ## Package Dependencies
#' dplyr (1.0.4), gtools (3.8.2), plotfunctions (1.4), raster (3.4), rgbif (3.5.2)
#' rgdal (1.5.23), rgeos (0.5.5), sampSurf (0.7.5), sp (1.4.5)
#' 
#' ## Vignette
#' We have provided an example vignette ("vignettes/megaSDM") using occurrence data from GBIF
#' for 6 North American mammal species. Environmental data are from the WorldClim database
#' (Hijmans et al. 2005; https://www.worldclim.org/), and data on dispersal rate (in km/year)
#' were collected by HS. The vignette displays the entire functionality of megaSDM in a cohesive
#' workflow, from data collection to the analysis and presentation of results.
