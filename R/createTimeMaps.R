#' Create maps describing species range shifts across many time periods
#'
#' Creates maps describing species range shifts across multiple time periods. These
#' maps detail the step-wise expansions and contractions of the species distribution
#' through those time-steps, allowing for the visualization of both unidirectional
#' range shifts (e.g., a range expansion across all time-steps) and more complex
#' dynamics (e.g., a range expansion followed by a range contraction).
#'
#' @param result_dir the directory where the ensembled and binary maps are placed.
#' Each species should have its own sub-directory, and the forecasted/hindcasted
#' binary maps should be placed into directories like so: Species/Scenario/Time.
#' If \code{projectSuit} was used to make these maps, this is probably the same
#' as the \code{output} argument in that function.
#' @param time_periods a vector of the years in which the projection will occur. The first
#' element should be the original year (the year in which the model was generated). If no
#' precise years are available (e.g., using data from the Last Glacial Maximum), order data
#' from current to least current (furthest into the future/past) and give character strings
#' for the years (e.g., "LGM").
#' @param scenarios a vector of character strings detailing the different climate models
#' used in the forecasted/hindcasted species distribution models.
#' @param dispersal (logical \code{TRUE} or \code{FALSE}). Are the binary maps given
#' dispersal-constrained maps or non-dispersal-constrained maps? If dispersal rate
#' analyses are not needed, or if the \code{megaSDM::dispersalRate} function has yet
#' to be run, this should be set to \code{FALSE} (the default).
#'
#' NOTE: if \code{dispersal = TRUE}, then \code{time_periods} have to be numeric.
#' @param dispersaldata either a dataframe or the complete path name of a .csv file with two columns:
#'
#'   Column 1: species name (same as the name used for modelling).
#'
#'   Column 2: average dispersal rate of species in kilometers/year.
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @export
#' @return writes both rasters and pdfs of the maps showing the expansion and contraction
#' of each species range across all time-steps for each climate scenario.

createTimeMaps <- function(result_dir, time_periods, scenarios,
                           dispersal, dispersaldata, ncores) {

  spp.list <- list.dirs(result_dir, full.names = FALSE, recursive = FALSE)
  spp.list <- spp.list[grep("_", spp.list)]
  if (length(spp.list) == 0) {
    stop(paste0("No projected models found in 'result_dir': Ensure that 'result_dir' provides a path to the proper location"))
  }
  #Generates the species list for parallelization
  if (dispersal == "TRUE") {

    #Reads in dispersal data
    if (class(dispersaldata) == "character") {
      dispersaldata <- utils::read.csv(dispersaldata, stringsAsFactors = FALSE)
      dispersaldata[, 1] <- gsub("_", " ", dispersaldata[, 1])
    } else {
      dispersaldata[, 1] <- gsub("_", " ", dispersaldata[, 1])
    }

    ListSpp <- c()
    for (i in 1:length(spp.list)) {
      curspecies <- spp.list[i]
      if (file.exists(paste0(result_dir, "/", curspecies, "/Results_Dispersal.csv"))) {
        ListSpp <- c(ListSpp, spp.list[i])
      }
    }
    for (w in 1:length(spp.list)) {
      FocSpec <- gsub("_", " ",spp.list[w])
      DispSpec <- grep(paste0("^", FocSpec, "$"), dispersaldata[, 1])
      if (length(DispSpec) == 0) {
        message(paste0("No dispersal rate values found for ", FocSpec, ": skipping dispersal rate analysis"))
      }
    }
  } else {
    ListSpp <- spp.list
    dispersaldata <- NA
  }
  if (length(ListSpp) < ncores) {
    ncores <- length(ListSpp)
  }
  ListSpp <- matrix(ListSpp, ncol = ncores)

  #Lists Scenarios and time_periods
  numScenario <- length(scenarios)
  numYear <- length(time_periods)
  #Creates a vector of years
  largestNum <- "1"
  for (year in 1:numYear) {
    largestNum <- paste0(largestNum, "0")
  }

  #Checking to see if the time periods pass through the year the model is projected on
  largestNum <- as.numeric(largestNum)
  if (is.numeric(time_periods)) {
    if (time_periods[length(time_periods)] > time_periods[1]) {
      if (!identical(sort(time_periods), time_periods)) {
        timeSort <- sort(time_periods)
      } else {
        timeSort <- time_periods
      }
    } else if (time_periods[length(time_periods)] < time_periods[1]) {
      if (!identical(sort(time_periods, decreasing = TRUE), time_periods)) {
        timeSort <- sort(time_periods)
      } else {
        timeSort <- time_periods
      }
    }
  }

  #Functions---------------------

  #Overlaps two rasters
  overlap <- function(t1, t2) {
    return(raster::mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
  }

  #Removes some temporary raster files for disk management
  removeTempRasterFiles <- function(rasterNames) {
    for (i in 1:length(rasterNames)) {
      removingName <- paste0(substr((rasterNames[i]), 0, (nchar(rasterNames[i]) - 4)), ".gri")
      file.remove(removingName)
      file.remove(rasterNames[i])
    }
    gc()
  }

  #Determines if there is a consistent trend (expanding or contracting) among the time periods
  #If there is more than one change from 0 to 1 or vice versa, "FALSE"
  consecutiveCheck <- function(val) {
    change <- 0
    valString<- unlist(strsplit(as.character(val), ""))
    firstVal <- as.numeric(valString[1])
    #Counts the number of times a "0" switches to a "1" or vice-versa
    for (i in 2:length(valString)) {
      if (as.numeric(valString[i]) != firstVal) {
        firstVal <- as.numeric(valString[i])
        change <- change + 1
      }
    }
    if (change > 1) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  #Sets the values of the binary (0, 1) raster to the values of the given state (e.g.: 101)
  setbinary <- function(r, val) {
    return (raster::calc(r, fun = function(x) {x * val}))
  }

  #Converts binary numbers to decimal numbers
  BinToDec <- function(x) {
    sum(2 ^ (which(rev((unlist(strsplit(as.character(x), "")) == 1))) - 1))
  }

  run <- function(CurSpp) {
    #Creates new directories for the Time Maps
    if (!dir.exists(file.path(result_dir, CurSpp, "TimeMapRasters"))) {
      dir.create(file.path(result_dir, CurSpp, "TimeMapRasters"))
    }

    if (!dir.exists(file.path(result_dir, CurSpp, "TimeMaps"))) {
      dir.create(file.path(result_dir, CurSpp, "TimeMaps"))
    }

    #Lists current binary raster files
    modernData <- list.files(path = file.path(result_dir, CurSpp), pattern = paste0("binary.bil$"), full.names = TRUE)
    rasterCRS <- raster::crs(raster::raster(modernData[1]))
    directories <- list.dirs(path = file.path(result_dir, CurSpp))
    correctDirectories <- c()

    #Gets only the forecasted/hindcasted file folders
    for (w in 1:length(scenarios)) {
      cordir <- grep(paste0("^", result_dir, "/", CurSpp, "/", scenarios[w], "$"), directories, perl = TRUE)
      if(length(cordir) == 1) {
        correctDirectories <- c(correctDirectories, directories[cordir])
      } else {
        stop("Directory with name ", scenarios[w], " not found! Revise binary map file management to Species/Scenario/Time")
      }
    }

    if (dispersal == TRUE) {
      filepattern <- "binary_dispersalRate.bil$"
    } else {
      filepattern <- "binary.bil$"
    }

    #Creates a matrix with the correct file paths to the current and future rasters
    collength <- length(list.files(path = file.path(correctDirectories[1]),
                                   pattern = paste0(filepattern),
                                   full.names = TRUE))

    results <- matrix(data = NA,
                      nrow = collength + 1,
                      ncol = length(correctDirectories),
                      byrow = FALSE,
                      dimnames = NULL)

    sortmodern <- which(timeSort == time_periods[1])

    for (i in 1:ncol(results)) {
      results[sortmodern, i] <- modernData
    }

    for(i in 1:length(correctDirectories)) {
      files <- list.files(path = file.path(correctDirectories[i]),
                          pattern = paste0(filepattern),
                          full.names = TRUE)

      for (j in 2:nrow(results)) {
        sortfut <- which(timeSort == time_periods[j])
        futfile <- files[grep(time_periods[j], files)]
        results[sortfut, i] <- futfile
      }
    }

    #Stacks the list of projected raster files and gives the stack the name %Scenarios[i]%
    for(i in 1:ncol(results)) {
      if (raster::extent(raster::raster(results[1, 1])) == raster::extent(raster::raster(results[2, 1]))) {
        assign(scenarios[i], raster::stack(results[, i]))
      } else {
        #if the extents don't match, crop the rasters to ensure that they do
        assign(scenarios[i], NULL)
        for(j in 1:nrow(results)) {
          minX <- max(raster::extent(raster::raster(results[1, 1]))[1], raster::extent(raster::raster(results[2, 1]))[1])
          maxX <- min(raster::extent(raster::raster(results[1, 1]))[2], raster::extent(raster::raster(results[2, 1]))[2])
          minY <- max(raster::extent(raster::raster(results[1, 1]))[3], raster::extent(raster::raster(results[2, 1]))[3])
          maxY <- min(raster::extent(raster::raster(results[1, 1]))[4], raster::extent(raster::raster(results[2, 1]))[4])
          Extent2 <- rbind(c(minX, maxX), c(minY, maxY))
          if (j == 1) {
            assign(scenarios[i], raster::crop(raster::raster(results[j, i]), Extent2))
          } else {
            assign(scenarios[i], raster::stack(get(scenarios[i]), raster::crop(raster::raster(results[j, i]), Extent2)))
          }
        }
      }
    }

    #Makes a list of possible states of presence/absence through the multiple years
    l <- rep(list(0:1), numYear)
    possible <- expand.grid(l)
    possible <- possible[-1, ]

    ScenariosCalc <- c()
    allRasterNames <- list()
    for (i in 1:numScenario) {
      ScenariosCalc[[i]] <- list()
    }

    for(i in 1:numScenario) {
      #Defines variables
      rasterNames <- c()
      consRasterNames <- c()
      allRasters <- get(scenarios[i])
      breakpoints <- c(0)
      allRasterNames <- c()
      FinalPrintRast <- c()
      #For each possible state
      for(j in 1:nrow(possible)) {
        name <- ""
        binary <- possible[j, ]
        #Overlaps each binary raster (with a state of "1") for the given scenario and state
        for (col in 1:ncol(binary)) {
          name <- paste0(name, binary[col])
          if (binary[col] == 1) {
            if (!exists("computedRaster")) {
              computedRaster <- allRasters[[col]]
            } else {
              computedRaster <- overlap(computedRaster, allRasters[[col]])
              allRasterNames <- c(allRasterNames, raster::filename(computedRaster))
            }
          }
        }

        #Removes the temporary raster files
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        allRasterNames <- c()

        #creates a mask of the overlapped raster (removes areas that are state "0" for other times)
        for (k in 1:ncol(binary)) {
          if (binary[k] == 0) {
            computedRaster <- raster::mask(computedRaster, allRasters[[k]], inverse = FALSE, maskvalue = 1, updatevalue = 0)
            allRasterNames <- c(allRasterNames, raster::filename(computedRaster))
          }
        }

        #Removes the temporary raster files
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        allRasterNames <- c()

        #Creates a composite raster with the binary values of computedRaster
        PrintedRaster <- computedRaster
        PrintedRaster[PrintedRaster@data@values == 1] <- as.numeric(name)
        if (j == 1) {
          FinalPrintRast <- PrintedRaster
        } else {
          PrintStack <- raster::stack(PrintedRaster, FinalPrintRast)
          FinalPrintRast <- raster::calc(PrintStack, fun = max)
        }

        raster::crs(FinalPrintRast) <- rasterCRS

        #Seperates the rasters showing a consistent trend (consecutiveCheck) and those that don't
        if (consecutiveCheck(name)) {
          breakpoints <- c(breakpoints, as.numeric(name) - 0.1)
          consRasterNames <- c(consRasterNames, paste0("raster_", name))
          assign(paste0("raster_", name), setbinary(computedRaster, as.numeric(name)))
          removeTempRasterFiles(raster::filename(computedRaster))
          rm(computedRaster)
          gc()
        } else {
          #If the first and last time periods are a "0", the species distribution expanded and then contracted ("expand-contract")
          if ((binary[, 1] == 0) && (binary[, ncol(binary)] == 0)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 2)))
            breakpoints <- c(breakpoints, largestNum * 2 - 0.1)
            removeTempRasterFiles(raster::filename(computedRaster))
            rm(computedRaster)
            gc()
            #If the first and last time periods are a "1", the species distribution contracted and then expanded ("contract-expand")
          } else if ((binary[, 1] == 1) && (binary[, ncol(binary)] == 1)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 4)))
            breakpoints <- c(breakpoints, largestNum * 4 - 0.1)
            removeTempRasterFiles(raster::filename(computedRaster))
            rm(computedRaster)
            gc()
            #If the first and last time periods are different and there is no consistent trend ("mixed")
          } else {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 3)))
            breakpoints <- c(breakpoints, largestNum * 3 - 0.1)
            removeTempRasterFiles(raster::filename(computedRaster))
            rm(computedRaster)
            gc()
          }
        }
      }

      #write the output raster
      if (dispersal == TRUE){
        raster::writeRaster(FinalPrintRast,
                    filename = file.path(result_dir, CurSpp, "TimeMapRasters", paste0("binary", scenarios[i], "_dispersal.bil")),
                    overwrite = TRUE,
                    format = "EHdr",
                    prj = TRUE)
      } else {
        raster::writeRaster(FinalPrintRast,
                            filename = file.path(result_dir, CurSpp, "TimeMapRasters", paste0("binary", scenarios[i], ".bil")),
                            overwrite = TRUE,
                            format = "EHdr",
                            prj = TRUE)
      }

      #Creates breakpoints (used to delineate colors in the PDFs)
      breakpoints <- c(breakpoints, max(breakpoints) + 1)
      breakpoints <- unique(breakpoints)
      #Makes names for the legend
      consRasterNames <- sort(consRasterNames)
      rasterNames <- gtools::mixedsort(rasterNames)
      stackedRaster <- get(consRasterNames[1])
      for (index in 2:length(consRasterNames)) {
        stackedRaster <- raster::stack(stackedRaster, get(consRasterNames[index]))
      }
      if (length(rasterNames) > 0) {
        for (index in 1:length(rasterNames)) {
          stackedRaster <- raster::stack(stackedRaster, get(rasterNames[index]))
        }
      }

      #Gets colors for the graph
      color <- c()
      color <- c(color, "lightgrey")
      dbtolb <- grDevices::colorRampPalette(c("darkblue", "dodgerblue"))(floor(length(consRasterNames) / 2))
      redtoyellow <- grDevices::colorRampPalette(c("red", "yellow"))(ceiling(length(consRasterNames) / 2))

      #The amount of "mixed" depends on the number of time steps
      if (length(time_periods) == 3) {
        pinktopurple <- grDevices::colorRampPalette(c("magenta3", "pink"))(2)
      } else {
        pinktopurple <- grDevices::colorRampPalette(c("magenta3", "pink"))(3)
      }
      if (length(time_periods) == 2) {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
      } else {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
        color <- c(color, pinktopurple)
      }

      #Combines all of the rasters into a single one
      ScenariosCalc <- c(ScenariosCalc, raster::calc(stackedRaster, fun = max))
      r <- ScenariosCalc[[length(ScenariosCalc)]]

      #Gets the full legend captions
      bin <- c()
      for (index in 1:length(consRasterNames)) {
        bin <- c(bin, unlist(strsplit(consRasterNames[index], "_"))[2])
      }
      if (length(rasterNames) > 0) {
        for (index in 1:length(rasterNames)) {
          bin <- c(bin, unlist(strsplit(rasterNames[index], "_"))[2])
        }
      }

      BinaryPoss <- rep(c(), nrow(possible))
      for(p in 1:nrow(possible)) {
        q <- 1
        BinaryPoss[p] <- possible[p, q]
        for (q in 2:ncol(possible)) {
          BinaryPoss[p] <- paste0(BinaryPoss[p], (possible[p, q]))
        }
      }

      #Makes "expand-contract", "contract-expand", and "mixed" as legend captions
      for(binIndex in 1:length(bin)) {
        if (as.numeric(bin[binIndex]) > largestNum) {
          for (p in 1:length(BinaryPoss)) {
            if (grepl(BinaryPoss[p], bin[binIndex]) == TRUE) {
              if ((possible[p, 1] == 0) && (possible[p, ncol(possible)] == 0)) {
                bin[binIndex] <- paste0("expand-contract")
              } else if ((possible[p, 1] == 1) && (possible[p, ncol(possible)] == 1)) {
                bin[binIndex] <- paste0("contract-expand")
              } else {
                bin[binIndex] <- paste0("mixed")
              }
            }
          }
        }
      }

      bin <- unique(bin)
      for(deleteIndex in 1:length(names(stackedRaster))) {
        removeTempRasterFiles(raster::filename(stackedRaster[[deleteIndex]]))
      }
      rm(stackedRaster)
      rm(list = consRasterNames)
      rm(list = rasterNames)
      gc()

      if (dispersal == TRUE) {
        #creates pdf of the Time Maps
        grDevices::pdf(file = file.path(result_dir, CurSpp, "TimeMaps", paste0(scenarios[i], "_dispersalRate_TimeMap.pdf")))
        raster::plot(legend = FALSE,
             breaks = breakpoints,
             r,
             col = color,
             xlab = "",
             ylab = "",
             main = paste0(CurSpp, " ", scenarios[i], " dispersalRate"))
        graphics::legend("bottomright", legend = rev(bin), fill = rev(color), cex = 0.6)
        grDevices::dev.off()
      } else {
        #creates pdf of the Time Maps
        grDevices::pdf(file = file.path(result_dir, CurSpp, "TimeMaps", paste0(scenarios[i], "_TimeMap.pdf")))
        raster::plot(legend = FALSE,
             breaks = breakpoints,
             r,
             col = color,
             xlab = "",
             ylab = "",
             main = paste0(CurSpp, " ", scenarios[i]))
        graphics::legend("bottomright", legend = rev(bin), fill = rev(color), cex = 0.6)
        grDevices::dev.off()
      }

      removeTempRasterFiles(raster::filename(r))
      rm(r)
      gc()
    }
    raster::removeTmpFiles(h = 0)
    rm(ScenariosCalc)
    gc()
  }
  if (ncores == 1) {
    ListSpp <- as.vector(ListSpp)
    out <- sapply(ListSpp, function(x) run(x))
  } else {
    clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)
    
    parallel::clusterExport(clus, varlist = c("run", "overlap", "removeTempRasterFiles",
                                              "consecutiveCheck", "setbinary", "BinToDec",
                                              "numScenario", "timeSort", "numYear", "largestNum",
                                              "result_dir", "time_periods", "scenarios",
                                              "dispersal", "dispersaldata", "ncores",
                                              "ListSpp"), envir = environment())
    
    parallel::clusterEvalQ(clus, library(gtools))
    parallel::clusterEvalQ(clus, library(raster))
    
    for (i in 1:nrow(ListSpp)) {
      out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
      gc()
    }
    parallel::stopCluster(clus)
  }
}
