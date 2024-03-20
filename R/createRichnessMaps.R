#' Create regular and dispersal-constrained richness maps from stacked SDMs
#'
#' This function stacks binary (presence/absence) species distribution maps to
#' create richness maps for a list of species. If higher taxa are provided
#' (i.e., mammals), separate richness maps for each higher taxon will be created
#' in addition to the full species richness maps. Given hindcasted/forecasted
#' binary maps, future/past species richness will also be calculated. Rasters
#' and .pdf files of each are generated, and the legend is standardized across
#' all .pdf maps for effective comparisons. Finally, provided distribution maps
#' that are constrained by dispersal rate, this function compares between the
#' dispersal-constrained and regular richness maps.
#'
#' @param result_dir the directory where the ensembled and binary maps are placed.
#' Each species should have its own sub-directory, and the forecasted/hindcasted
#' binary maps should be placed into directories like so: Species/Scenario/Time.
#' If \code{projectSuit} was used to make these maps, this is probably the same
#' as the \code{output} argument in that function.
#' @param time_periods a vector of the years in which the projection will occur with the first element
#' as the original year (the year in which the model was generated). If no precise years
#' are available for future/past time periods (e.g., using data from the Last Glacial Maximum),
#' order the remaining time periods from most current to least current (farthest into the future/past)
#' and give character strings for the years (e.g., "LGM").
#' @param scenarios a vector of character strings detailing the different climate models
#' used in the forecasted/hindcasted species distribution models. If no projection is
#' needed, set to NA (default).
#' @param dispersal (logical \code{TRUE} or \code{FALSE}) Are dispersal rate analyses
#' needed? If set to \code{TRUE}, this function constructs dispersal-constrained richness
#' maps for all time periods/scenarios given, and compares those maps to the non-dispersal
#' constrained maps. If these analysis are not needed, or if the \code{megaSDM::dispersalRate}
#' function has yet to be run, this should be set to \code{FALSE} (default).
#' @param taxonlist (optional) data.frame or .csv file path name.
#' If taxon-specific species richness maps are required, this should correspond to
#' a list of all species modelled (Column 1) and a selected higher taxon
#' (Column 2; e.g., Order, Class) to split the richness maps on. If taxon-specific
#'  richness maps are not required, \code{taxonlist} should be set to \code{FALSE} (default).
#' @export
#' @return writes .pdf files and rasters (.grd) of species richness for all given
#' taxa, time periods, and scenarios and dispersal constraints to the directory
#' "\code{result_dir}/RichnessMaps".

createRichnessMaps <- function(result_dir, time_periods, scenarios = NA,
                               dispersal = FALSE, taxonlist = FALSE) {

  if (is.na(scenarios[1])) {
    numScenario <- 0
    numYear <- 0
  } else {
    numScenario <- length(scenarios)
    numYear <- length(time_periods)
  }

  spp.list <- list.dirs(result_dir, full.names = FALSE, recursive = FALSE)
  spp.list <- spp.list[grep("_", spp.list)]
  if (length(spp.list) == 0) {
    stop(paste0("No projected models found in 'result_dir': Ensure that 'result_dir' provides a path to the proper location"))
  }
  spp.list <- data.frame(Species = spp.list)
  if (taxonlist == TRUE) {
    if (class(taxonlist) == "character") {
      taxonlist <- utils::read.csv(taxonlist, stringsAsFactors = FALSE)
      taxonlist[, 1] <- gsub("_", " ", taxonlist[, 1])
    } else {
      taxonlist[, 1] <- gsub("_", " ", taxonlist[, 1])
    }

    ListSpp <- merge(spp.list, taxonlist, by.x = "Species", by.y = c(colnames(taxonlist)[1]))

    if((nrow(ListSpp) != nrow(spp.list)) | (nrow(ListSpp) != nrow(taxonlist))) {
      message("Warning! The list of species modelled is not the same as the list of species given in 'taxonlist'")
      message("Ensure that the taxonomy in consistent and that 'taxonlist' contains all modelled species")
    }

  } else {
    spp.list$Taxon <- rep("All_Taxa", nrow(spp.list))
    ListSpp <- spp.list
  }

  #Finds and removes species with AUC Values less than desired threshold
  DeleteSP <- c()
  for (sp in 1:nrow(ListSpp)) {
    focusspp <- gsub(" ", "_", ListSpp[sp, 1])
    curfocus <- list.files(path = file.path(result_dir, focusspp), pattern = paste0("binary.grd$"))
    #If no binary maps were generated, AUC < threshold
    if (!length(curfocus) > 0) {
      message(paste0(focusspp, " was not modelled (likely due to low AUC values)"))
      DeleteSP <- c(DeleteSP, sp)
    }
  }

  #Deletes species from the taxonlist if AUC values are too small
  if (length(DeleteSP) > 0) {
    ListSpp <- ListSpp[-DeleteSP, ]
  } else {
    ListSpp <- ListSpp
  }

  #Gets the different broad taxa to do species richness on
  taxa <- c(levels(ListSpp[, 2]), "All_Taxa")
  ntaxa <- length(unique(taxa))
  if (ntaxa <= 2) {
    taxa <- taxa[1]
  }
  ntaxa <- length(taxa)

  #Create species richness directory
  if(!dir.exists(file.path(result_dir, "RichnessMaps"))) {
    dir.create(file.path(result_dir, "RichnessMaps"))
  }

  RichnessMaps <- file.path(result_dir, "RichnessMaps")

  if (numScenario == 0) {
    scenarios <- "Current"
  }

  #Functions------------------------------------

  #Uploads and stacks binary SDM rasters of particular time, taxon, and scenario
  FutureRichnessMaps <- function(taxon, y, s, disp) {
    futstack <- c()
    FocusYear <- time_periods[y]
    FocusScenario <- scenarios[s]

    if (disp == TRUE) {
      filepath <- "_dispersalRate.grd"
    } else {
      filepath <- ".grd"
    }

    for (sp in 1:length(taxonlist2)) {
      focusspp <- taxonlist2[sp]
      path <- paste0(result_dir, "/", focusspp)
      r <- list.files(path = path, pattern = paste0(FocusYear, "_", FocusScenario, "_binary", filepath, "$"), recursive = TRUE)
      futraster <- terra::rast(paste0(path, "/", r))
      futstack <- c(futraster, futstack)
    }


    #Calculates sum (total species richness) for given parameters, makes rasters and PDFs
    FutureRichness <- terra::app(futstack, fun = sum)
    terra::writeRaster(FutureRichness,
                filename = file.path(RichnessMaps, paste0("Richness_", taxon, "_", FocusYear, "_", FocusScenario, filepath)),
                overwrite = TRUE)
  }

  DispersalDiff <- function(taxon, y, s) {
    FocusYear <- time_periods[y]
    FocusScenario <- scenarios[s]
    #Loads Non-dispersal constrained richness maps
    NonDispersalSR <- terra::rast(paste0(RichnessMaps, "/Richness_", taxon, "_", FocusYear, "_", FocusScenario, ".grd"))

    #Loads dispersal-constrained richness maps for comparison
    DispersalSR <- terra::rast(paste0(RichnessMaps, "/Richness_", taxon, "_", FocusYear, "_", FocusScenario, "_dispersalRate.grd"))

    #Calculates difference between dispersal and non-dispersal species richness maps
    DifferenceSR <- DispersalSR - NonDispersalSR
    terra::writeRaster(DifferenceSR,
                filename = file.path(RichnessMaps, paste0("Richness Difference_", taxon, "_", FocusYear, "_", FocusScenario, ".grd")),
                overwrite = TRUE)

    #Sets graphical parameters for dispersal-difference PDFs
    DiffVec <- terra::unique(DifferenceSR)
    DiffVec <- DiffVec[,1]
    if (min(DiffVec) < 0 & max(DiffVec) > 0) {
      color <- grDevices::colorRampPalette(c("red", "grey91", "blue"))(length(DiffVec))
    } else if (min(DiffVec) < 0 & max(DiffVec) <= 0) {
      color <- grDevices::colorRampPalette(c("red", "gold1", "grey91"))(length(DiffVec))
    } else {
      color <- grDevices::colorRampPalette(c("grey91", "deepskyblue", "blue"))(length(DiffVec))
    }

    #Creates dispersal-difference PDFs
    grDevices::pdf(file = file.path(result_dir, "RichnessMaps", paste0(taxon, "_", FocusYear, "_", FocusScenario, "_", "dispersalDifference.pdf")))
    terra::plot(DifferenceSR, col = color, xlab = "", ylab = "", legend = FALSE, main = paste0(taxon, " ", FocusScenario, " ", FocusYear, " Dispersal Difference"))
    plotfunctions::gradientLegend(c(min(DiffVec):max(DiffVec)), color = color, pos = 0.125, side = 4, n.seg = 2, dec = 0, fit.margin = TRUE, inside = TRUE)
    grDevices::dev.off()
  }

  MakePDFMaps <- function(taxon) {
    #Upload all raster files, calculate overall maximum richness (regardless of time, scenario)
    RasterList <- list.files(path = RichnessMaps, pattern = paste0("\\.grd$"), full.names = TRUE)
    FocusRasters <- RasterList[grep(paste0("Richness_", taxon), RasterList)]
    FocusStack <- terra::rast(FocusRasters)
    MaximumRich <- max(terra::global(FocusStack, fun = max, na.rm = TRUE))

    for (e in 1:terra::nlyr(FocusStack)) {
      Raster1 <- FocusStack[[e]]
      MaxRaster1 <- terra::global(Raster1, fun = max, na.rm = TRUE)
      MaxRaster1 <- as.numeric(MaxRaster1)
      ValueFreq <- terra::freq(Raster1, bylayer = FALSE)
      #If the maximum species richness occurs in too few pixels to display in pdf, only plots the next highest species richness
      while ((ValueFreq[which(ValueFreq[, 1] == MaxRaster1), 2] < floor(0.00005 * terra::ncell(Raster1))) && (MaxRaster1 > 1)) {
        Raster1[Raster1 == MaxRaster1] <- MaxRaster1 - 1
        message(paste0("The maximum species richness in the study area, ", MaxRaster1, ", covers too few cells to display on PDF:"))
        message("Examine the created raster for an accurate count of pixels.")
        MaxRaster1 <- MaxRaster1 - 1
      }

      #Gets color scheme for legend
      breakpoints <- seq(from = 0, to = MaximumRich, length = MaximumRich)
      color <- c("lightgrey", grDevices::colorRampPalette(c("springgreen", "dodgerblue", "darkblue"))(length(breakpoints)))

      #Creates pdf of the species richness raster
      grDevices::pdf(file = file.path(paste0(substr(FocusRasters[e], 1, (nchar(FocusRasters[e]) - 4)), ".pdf")))
      terra::plot(legend = FALSE, Raster1, col = color[1:(MaxRaster1 + 1)], xlab = "", ylab = "", main = names(FocusStack[[e]]))
      plotfunctions::gradientLegend(c(0:MaximumRich), color = color, pos = 0.125, side = 4, n.seg = 2, dec = 0, fit.margin = TRUE, inside = TRUE)
      grDevices::dev.off()
    }
    rm(FocusStack)
    gc()
  }

  #Makes richness maps by taxon
  for (t in 1:ntaxa) {
    #Subsets out taxa of interest
    taxon <- taxa[t]
    if (taxon == "All_Taxa") {
      taxonlist2 <- gsub(" ", "_", ListSpp[, 1])
    } else {
      taxonlist2 <- gsub(" ", "_", ListSpp[which(ListSpp[, 2] == taxon), 1])
    }
    print(paste0("    Creating current species richness maps for ", taxon, ", with ", length(taxonlist2), " total species"))

    #Uploads and stacks current binary rasters of all species within the taxa
    #Removes species with AUC values < threshold
    DeleteSP <- c()
    for (sp in 1:length(taxonlist2)) {
      focusspp <- taxonlist2[sp]
      if (dir.exists(paste0(result_dir, "/", focusspp))) {
        curfocus <- list.files(path = file.path(result_dir, focusspp), pattern = paste0("binary.grd$"), full.names = TRUE)
        if (length(curfocus) > 0){
          if(!exists("curstack")) {
            curstack <- terra::rast(curfocus)
          } else {
            curraster <- terra::rast(curfocus)
            curstack <- c(curstack, curraster)
          }
        } else {
          message(paste0(focusspp, " will be removed (was modelled but had no replicates of a high enough AUC value)"))
          DeleteSP <- c(DeleteSP, sp)
        }
      }
    }

    #Deletes species with AUC values < threshold
    if (length(DeleteSP) > 0) {
      taxonlist2 <- taxonlist2[-DeleteSP]
    } else {
      taxonlist2 <- taxonlist2
    }

    #Calculates the sum of the rasters (species richness), and writes rasters
    CurrentRichness <- terra::app(curstack, fun = sum)
    terra::crs(CurrentRichness) <- terra::crs(curstack)
    terra::writeRaster(CurrentRichness,
                filename = file.path(RichnessMaps, paste0("Richness_", taxon, "_", time_periods[1], ".grd")),
                overwrite = TRUE)

    #Creates species richness maps for other time frames and scenarios
    if (numYear > 1) {
      for (y in 2:numYear) {
        for (s in 1:numScenario) {
          FutureRichnessMaps(taxon, y, s, FALSE)
        }
      }
      print(paste0("    Creating species richness maps for forecasted/hindcasted scenarios"))
    }

    #Calculates species richness with dispersal constraints
    if (dispersal == TRUE) {
      #Highlights species with AUC values < threshold or no dispersal data, adds others to taxonlist
      for (sp in 1:length(taxonlist2)) {
        focusspp <- taxonlist2[sp]
        curfocus <- list.files(path = file.path(result_dir, focusspp, scenarios[1]),
                               pattern = paste0("binary_dispersalRate.grd$"), full.names = TRUE)
        if (!length(curfocus) > 0) {
          message(paste0(focusspp, " will be removed (no dispersal rate rasters were found)"))
          DeleteSP <- c(DeleteSP, sp)
        }
      }

      #Deletes species with AUC values < threshold
      if (length(DeleteSP) > 0) {
        taxonlist2 <- taxonlist2[-DeleteSP]
      } else {
        taxonlist2 <- taxonlist2
      }

      if (length(taxonlist2) > 0) {
        #Creates species richness maps for time frames and scenarios
        for (y in 2:numYear) {
          for (s in 1:numScenario) {
            FutureRichnessMaps(taxon, y, s, TRUE)
            DispersalDiff(taxon, y, s)
          }
        }
        print(paste0("    Creating species richness maps for forecasted/hindcasted scenarios (including dispersal rate)"))
      }
    }
  }
  print("    Making PDF Maps...")
  for(t in 1:ntaxa) {
    taxon <- taxa[t]
    MakePDFMaps(taxon)
  }
}
