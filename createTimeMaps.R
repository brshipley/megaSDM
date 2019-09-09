#Initializations-------------------------------------------------
numYear <- df[, "numYear"]
numScenario <- df[, "numScenario"]
result_dir <- df[, "result_dir"]
test <- df[, "test"]
rastertype <- df[, "rastertype"]
dispersalStep <- df[, "dispersalStep"]
ncores <- df[, "ncores"]
dispersalRan <- df[, "dispersalRan"]
if (dispersalStep == "Y") {
  dispersal <- read.csv(list.files(df[, "dispersalRate_dir"], full.names = TRUE)[1])
}
outfile <- "statsout.txt"
ScenarioNames <- as.character(df[grep("^Scenario", colnames(df))])

setwd(result_dir)
spp.list <- substr(spp_batch, 1, nchar(spp_batch) - 4)
nspp <- length(spp.list)

library(parallel)
library(raster)
library(gtools)

#Functions------------------------------------------------------
#Sets the values of the binary (0, 1) raster to the values of the given state (e.g.: 101) 
setbinary <- function(r, val) {
  return (calc(r, fun = function(x) {x * val}))
}

#This function determines if there is a consistent trend (waxing or waning) among the time periods
#If there is more than one change from 0 to 1 or vice versa, "FALSE"
consecutiveCheck <- function(val) {
  change <- 0
  valString<- unlist(strsplit(as.character(val), ""))
  firstVal <- as.numeric(valString[1])
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

#This function removes some raster files for disk management
removeTempRasterFiles <- function(rasterNames) {
  for (i in 1:length(rasterNames)) {
    removingName <- paste0(substr((rasterNames[i]), 0, (nchar(rasterNames[i]) - 4)), ".gri")
    file.remove(removingName)
    file.remove(rasterNames[i])
  }
  gc()
}

#This function calculates the number of ones (presence) for each raster 
numberOfOnes <- function(val) {
  return(length(which((unlist(strsplit(as.character(val), "")) == 1))))
}

#This function overlaps two rasters 
overlap <- function(t1, t2) {
  return(mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
}

#This function converts binary numbers to decimal numbers
BinToDec <- function(x) {
  sum(2 ^ (which(rev((unlist(strsplit(as.character(x), "")) == 1))) - 1))
}

#This function actually creates the Time Maps (without considering dispersal)
timeMapsNoDispersal <- function(spp) {
  tryCatch({
    setwd(result_dir)
    setwd(spp)
    dir.create("TimeMapRasters")
    dir.create("TimeMaps")
    setwd(result_dir)
    
    #Lists current binary raster files 
    modernData <- list.files(path = spp, pattern = paste0("binary", rastertype, "$"), full.names = TRUE)
    directories <- list.dirs(path = spp)
    correctDirectories <- c()
    
    #Gets only the projection raster files 
    for (w in 1:length(ScenarioNames)) {
      for(i in 1:length(directories)) {
        if(grepl(paste0("^", spp, "/", ScenarioNames[w], "$"), directories[i], perl = TRUE)) {
          correctDirectories <- c(correctDirectories, directories[i])
        }
      }
    }
    collength <- length(list.files(path = paste0(result_dir, "/", correctDirectories[1]), 
                                   pattern = paste0("binary", rastertype, "$"), 
                                   full.names = TRUE))
    results <- matrix(data = NA, 
                      nrow = collength + 1, 
                      ncol = length(correctDirectories), 
                      byrow = FALSE,
                      dimnames = NULL)
    
    for (i in 1:ncol(results)) {
      results[1, i] <- modernData
    }
    for(i in 1:length(correctDirectories)) {
      files <- list.files(path = paste0(result_dir, "/", correctDirectories[i]), 
                          pattern = paste0("binary", rastertype, "$"), 
                          full.names = TRUE)
      for (j in 2:nrow(results)) {
        results[j, i] <- files[j - 1]
      }
    }
    
    #Stacks the list of projected raster files and gives the stack the name %Scenarios[i]%
    for(i in 1:ncol(results)) {
      if (extent(raster(results[1, 1])) == extent(raster(results[2, 1]))) {
        assign(Scenarios[i], stack(results[, i]))
      } else {
        assign(Scenarios[i], NULL)
        for(j in 1:nrow(results)) {
          minX <- max(extent(raster(results[1, 1]))[1], extent(raster(results[2, 1]))[1])
          maxX <- min(extent(raster(results[1, 1]))[2], extent(raster(results[2, 1]))[2])
          minY <- max(extent(raster(results[1, 1]))[3], extent(raster(results[2, 1]))[3])
          maxY <- min(extent(raster(results[1, 1]))[4], extent(raster(results[2, 1]))[4])
          Extent2 <- rbind(c(minX, maxX), c(minY, maxY))
          assign(Scenarios[i], stack(crop(raster(results[j, i]), Extent2), get(Scenarios[i])))
        }
      }
    }
    
    #Makes a list of possible states of presence/absence through the multiple years 
    l <- rep(list(0:1), numYear)
    possible <- expand.grid(l)
    possible <- possible[-1, ]
    
    allConsRasterNames <- c()
    ScenariosCalc <- c()
    
    for(i in 1:numScenario) {
      #Defines variables
      rasterNames <- c()
      consRasterNames <- c()
      allRasters <- get(Scenarios[i])
      breakpoints <- c(0)
      allRasterNames <- c()
      FinalPrintRast <- c()
      
      #For each possible state:
      for(j in 1:nrow(possible)) {
        name <- ""
        binary <- possible[j, ]
        #Overlaps each binary raster (with a state of "1") for the given scenario and state, 
        for (col in 1:ncol(binary)) {
          name <- paste0(name, binary[col])
          if (binary[col] == 1) {
            if (!exists("computedRaster")) {
              computedRaster <- allRasters[[col]]
            } else {
              computedRaster <- overlap(computedRaster, allRasters[[col]])
              allRasterNames <- c(allRasterNames, filename(computedRaster))
            }
          }
        }
        
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        allRasterNames <- c()
        
        #creates a mask of the overlapped raster (removes areas that are state "0" for other times)
        for (k in 1:ncol(binary)) {
          if (binary[k] == 0) {
            computedRaster <- mask(computedRaster, allRasters[[k]], inverse = FALSE, maskvalue = 1, updatevalue = 0)
            allRasterNames <- c(allRasterNames, filename(computedRaster))
          }
        }
        
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        
        allRasterNames <- c()
        print("HERE3")
        print(name)
        print(class(name))
        print(as.character(name))
        
        #Creates a composite raster with the binary values of computedRaster
        PrintedRaster <- computedRaster
        PrintedRaster[PrintedRaster@data@values == 1] <- as.numeric(name)
        if (j == 1) {
          FinalPrintRast <- PrintedRaster
        } else {
          PrintStack <- stack(PrintedRaster, FinalPrintRast)
          FinalPrintRast <- calc(PrintStack, fun = max)
        }
        
        #Seperates the rasters showing a consistent trend (consecutiveCheck) and those that don't
        if (consecutiveCheck(name)) {
          breakpoints <- c(breakpoints, as.numeric(name) - 0.1)
          consRasterNames <- c(consRasterNames, paste0("raster_", name))
          assign(paste0("raster_", name), setbinary(computedRaster, as.numeric(name)))
          removeTempRasterFiles(filename(computedRaster))
          rm(computedRaster)
          gc()
        } else {
          #If the first and last time periods are a "0", the species distribution expanded and then contracted ("wax-wane")
          if ((binary[, 1] == 0) && (binary[, ncol(binary)] == 0)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 2)))
            breakpoints <- c(breakpoints, largestNum * 2 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          #If the first and last time periods are a "1", the species distribution contracted and then expanded ("wane-wax")
          } else if ((binary[, 1] == 1) && (binary[, ncol(binary)] == 1)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 4)))
            breakpoints <- c(breakpoints, largestNum * 4 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          #If the first and last time periods are different and there is no consistent trend, ("mixed")
          } else {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 3)))
            breakpoints <- c(breakpoints, largestNum * 3 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          }
        }
      }
      
      #Write the output raster
      writeRaster(FinalPrintRast, 
                  filename = paste0(result_dir, "/", spp,"/TimeMapRasters", "/binary", Scenarios[i], rastertype),
                  overwrite = TRUE, 
                  format = format,
                  prj = TRUE)
      #create breakpoints (used to delineate colors in the time maps)
      breakpoints <- c(breakpoints, max(breakpoints) + 1)
      breakpoints <- unique(breakpoints)
      #Make names for the legend
      consRasterNames <- sort(consRasterNames)
      rasterNames <- mixedsort(rasterNames)
      stackedRaster <- get(consRasterNames[1])
      for (index in 2:length(consRasterNames)) {
        stackedRaster <- stack(stackedRaster, get(consRasterNames[index]))
      }
      if (length(rasterNames) > 0) {
        for (index in 1:length(rasterNames)) {
          stackedRaster <- stack(stackedRaster, get(rasterNames[index]))
        }
      }
      # Getting colors for the graph
      color <- c()
      color <- c(color, "lightgrey")
      dbtolb <- colorRampPalette(c("darkblue", "dodgerblue"))(floor(length(consRasterNames) / 2))
      redtoyellow <- colorRampPalette(c("red", "yellow"))(ceiling(length(consRasterNames) / 2))
      
      if (length(years) == 3) {
        pinktopurple <- colorRampPalette(c("magenta3", "pink"))(2)
      } else {
        pinktopurple <- colorRampPalette(c("magenta3", "pink"))(3)
      }
      
      if (length(years) == 2) {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
      } else {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
        color <- c(color, pinktopurple)
      }

      ScenariosCalc <- c(ScenariosCalc, calc(stackedRaster, fun = max))
      r <- ScenariosCalc[[length(ScenariosCalc)]]
      r <- as.factor(r)
      
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
      
      for(binIndex in 1:length(bin)) {
        if (as.numeric(bin[binIndex]) > largestNum) {
          for (p in 1:length(BinaryPoss)) {
            if (grepl(BinaryPoss[p], bin[binIndex]) == TRUE) {
              if ((possible[p, 1] == 0) && (possible[p, ncol(possible)] == 0)) {
                bin[binIndex] <- paste0("wax-wane")
              } else if ((possible[p, 1] == 1) && (possible[p, ncol(possible)] == 1)) {
                bin[binIndex] <- paste0("wane-wax")
              } else {
                bin[binIndex] <- paste0("mixed")
              }
            }
          }
        }
      }
      
      bin <- unique(bin)
      
      for(deleteIndex in 1:length(names(stackedRaster))) {
        removeTempRasterFiles(filename(stackedRaster[[deleteIndex]]))
      }
      rm(stackedRaster)
      rm(FinalPrintRast)
      rm(list = consRasterNames)
      rm(list = rasterNames)
      gc()
      
      #creates raster and pdf of the Time Maps
      pdf(file = paste0(result_dir, "/", spp, "/TimeMaps/", Scenarios[i],"_TimeMap.pdf"))
      plot(legend = FALSE, breaks = breakpoints, r, col = color, xlab = "", ylab = "", main = paste0(spp, " ", Scenarios[i]))
      legend("bottomright", legend = rev(bin), fill = rev(color), cex = 0.6)
      dev.off()
      removeTempRasterFiles(filename(r))
      rm(r)
      gc()
    }
    removeTmpFiles(h = 0)
    rm(ScenariosCalc)
    gc()
  }, error = function(err) {
    print(paste("MY_ERROR: ", spp, " ", err))
    return(paste0("error: ", err))
  })
}

timeMapsDispersal <- function(spp) {
  tryCatch({
    setwd(result_dir)
    setwd(spp)
    dir.create("TimeMapRasters")
    setwd(result_dir)
    #Finds and lists current binary rasters
    modernData <- list.files(path = spp, pattern = paste0("binary", rastertype, "$"), full.names = TRUE)
    modernData <- paste0(result_dir, "/", modernData)
    directories <- list.dirs(path = spp)
    correctDirectories <- c()
    #Locates the projection files
    for (w in 1:length(ScenarioNames)) {
      for(i in 1:length(directories)) {
        if(grepl(paste0("^", spp, "/", ScenarioNames[w], "$"), directories[i], perl = TRUE)) {
          correctDirectories <- c(correctDirectories, directories[i])
        }
      }
    }
    collength <- length(list.files(path = paste0(result_dir, "/", correctDirectories[1]), 
                                   pattern = paste0("binary_dispersalRate", rastertype, "$"), 
                                   full.names = TRUE))
    results <- matrix(data = NA, 
                      nrow = collength + 1, 
                      ncol = length(correctDirectories), 
                      byrow = FALSE,
                      dimnames = NULL)
    for (i in 1:ncol(results)) {
      results[1, i] <- modernData
    }
    for(i in 1:length(correctDirectories)) {
      files <- list.files(path = paste0(result_dir, "/", correctDirectories[i]), 
                          pattern = paste0("binary_dispersalRate", rastertype, "$"), 
                          full.names = TRUE)
      for (j in 2:nrow(results)) {
        results[j, i] <- files[j - 1]
      }
    }
    #Stacks the list of projected raster files and gives the stack the name %Scenarios[i]%
    for(i in 1:ncol(results)) {
      if (extent(raster(results[1, 1])) == extent(raster(results[2, 1]))) {
        assign(Scenarios[i], stack(results[, i]))
      } else {
        assign(Scenarios[i], NULL)
        for(j in 1:nrow(results)) {
          minX <- max(extent(raster(results[1, 1]))[1], extent(raster(results[2, 1]))[1])
          maxX <- min(extent(raster(results[1, 1]))[2], extent(raster(results[2, 1]))[2])
          minY <- max(extent(raster(results[1, 1]))[3], extent(raster(results[2, 1]))[3])
          maxY <- min(extent(raster(results[1, 1]))[4], extent(raster(results[2, 1]))[4])
          Extent2 <- rbind(c(minX, maxX), c(minY, maxY))
          assign(Scenarios[i], stack(crop(raster(results[j, i]), Extent2), get(Scenarios[i])))
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
      allRasters <- get(Scenarios[i])
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
              allRasterNames <- c(allRasterNames, filename(computedRaster))
            }
          }
        }
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        allRasterNames <- c()
        #creates a mask of the overlapped raster
        for (k in 1:ncol(binary)) {
          if (binary[k] == 0) {
            computedRaster <- mask(computedRaster, allRasters[[k]], inverse = FALSE, maskvalue = 1, updatevalue = 0)
            allRasterNames <- c(allRasterNames, filename(computedRaster))
          }
        }
        allRasterNames <- allRasterNames[-length(allRasterNames)]
        if (!is.null(allRasterNames) && length(allRasterNames) != 0) {
          removeTempRasterFiles(allRasterNames)
        }
        allRasterNames <- c()
        
        PrintedRaster <- computedRaster
        PrintedRaster[PrintedRaster@data@values == 1] <- as.numeric(name)	
        if (j == 1) {		
          FinalPrintRast <- PrintedRaster
        } else {		
          PrintStack <- stack(PrintedRaster, FinalPrintRast)		
          FinalPrintRast <- calc(PrintStack, fun = max)
        }		
        
        #Seperates the rasters showing a consistent trend (consecutiveCheck) and those that don't
        if (consecutiveCheck(name)) {
          breakpoints <- c(breakpoints, as.numeric(name) - 0.1)
          consRasterNames <- c(consRasterNames, paste0("raster_", name))
          assign(paste0("raster_", name), setbinary(computedRaster, as.numeric(name)))
          removeTempRasterFiles(filename(computedRaster))
          rm(computedRaster)
          gc()
        } else {
          #If the first and last time periods are a "0", the species distribution expanded and then contracted ("wax-wane")
          if ((binary[, 1] == 0) && (binary[, ncol(binary)] == 0)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 2)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 2)))
            breakpoints <- c(breakpoints, largestNum * 2 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          #If the first and last time periods are a "1", the species distribution contracted and then expanded ("wane-wax")
          } else if ((binary[, 1] == 1) && (binary[, ncol(binary)] == 1)) {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 4)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 4)))
            breakpoints <- c(breakpoints, largestNum * 4 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          #If the first and last time periods are different and there is no consistent trend ("mixed")
          } else {
            rasterNames <- c(rasterNames, paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))))
            assign(paste0("raster_", as.character(as.numeric(name) + (largestNum * 3)), "_", BinToDec(as.numeric(name))), setbinary(computedRaster, (largestNum * 3)))
            breakpoints <- c(breakpoints, largestNum * 3 - 0.1)
            removeTempRasterFiles(filename(computedRaster))
            rm(computedRaster)
            gc()
          }
        }
      }
      
      #write raster of time maps
      writeRaster(FinalPrintRast,  
                  filename = paste0(result_dir, "/", spp,"/TimeMapRasters", "/binary", Scenarios[i], "_dispersal", rastertype),
                  overwrite = TRUE, 
                  format = format,
                  prj = TRUE)
      
      #create breakpoints (used to delineate colors in the time maps)
      breakpoints <- c(breakpoints, max(breakpoints) + 1)
      breakpoints <- unique(breakpoints)
      #Make names for the legend
      consRasterNames <- sort(consRasterNames)
      rasterNames <- mixedsort(rasterNames)
      stackedRaster <- get(consRasterNames[1])
      for (index in 2:length(consRasterNames)) {
        stackedRaster <- stack(stackedRaster, get(consRasterNames[index]))
      }
      if (length(rasterNames) > 0) {
        for (index in 1:length(rasterNames)) {
          stackedRaster <- stack(stackedRaster, get(rasterNames[index]))
        }
      }
      # Getting colors for the graph
      color <- c()
      color <- c(color, "lightgrey")
      dbtolb <- colorRampPalette(c("darkblue", "dodgerblue"))(floor(length(consRasterNames) / 2))
      redtoyellow <- colorRampPalette(c("red", "yellow"))(ceiling(length(consRasterNames) / 2))
      
      #The amount of "mixed" depends on the number of time steps
      
      if (length(years) == 3) {
        pinktopurple <- colorRampPalette(c("magenta3", "pink"))(2)
      } else {
        pinktopurple <- colorRampPalette(c("magenta3", "pink"))(3)
      }
      if (length(years) == 2) {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
      } else {
        color <- c(color, dbtolb)
        color <- c(color, redtoyellow)
        color <- c(color, pinktopurple)
      }
      
      ScenariosCalc <- c(ScenariosCalc, calc(stackedRaster, fun = max))
      r <- ScenariosCalc[[length(ScenariosCalc)]]
      r <- as.factor(r)
      
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
      
      for(binIndex in 1:length(bin)) {
        if (as.numeric(bin[binIndex]) > largestNum) {
          for (p in 1:length(BinaryPoss)) {
            if (grepl(BinaryPoss[p], bin[binIndex]) == TRUE) {
              if ((possible[p, 1] == 0) && (possible[p, ncol(possible)] == 0)) {
                bin[binIndex] <- paste0("wax-wane")
              } else if ((possible[p, 1] == 1) && (possible[p, ncol(possible)] == 1)) {
                bin[binIndex] <- paste0("wane-wax")
              } else {
                bin[binIndex] <- paste0("mixed")
              }
            }
          }
        }
      }
      
      bin <- unique(bin)
      for(deleteIndex in 1:length(names(stackedRaster))) {
        removeTempRasterFiles(filename(stackedRaster[[deleteIndex]]))
      }
      rm(stackedRaster)
      rm(list = consRasterNames)
      rm(list = rasterNames)
      gc()
      #creates pdf of the Time Maps		
      pdf(file = paste0(result_dir, "/",  spp, "/TimeMaps/", Scenarios[i], "_dispersalRate_TimeMap.pdf"))
      plot(legend = FALSE, 
           breaks = breakpoints, 
           r, 
           col = color, 
           xlab = "", 
           ylab = "",
           main = paste0(spp, " ", Scenarios[i], " dispersalRate"))
      legend("bottomright", legend = rev(bin), fill = rev(color), cex = 0.6)
      dev.off()
      removeTempRasterFiles(filename(r))
      rm(r)
      gc()
    }
    removeTmpFiles(h = 0)
    rm(ScenariosCalc)
    gc()
  }, error = function(err) {
    print(paste("MY_ERROR: ", spp, " ", err))
    return(paste0("error: ", err))
  })
}

#Run-------------------------------
#Creates a vector of years
years <- c()
largestNum <- "1"
for (year in 1:numYear) {
  largestNum <- paste0(largestNum, "0")
  if(year != 1) {
    currentYear <- paste0("Year.", year - 1)
  } else {
    currentYear <- "Year"
  }
  workingYear <- as.numeric(df[currentYear])
  years[year] <- workingYear
}
largestNum <- as.numeric(largestNum)

#Creates a vector of scenarios
Scenarios <- c()
for (ScenIndex in 1:numScenario) {
  if(ScenIndex != 1) {
    cScenario <- paste0("Scenario.", ScenIndex - 1)
  } else {
    cScenario <- "Scenario"
  }
  Scenario <- df[, cScenario]
  Scenarios[ScenIndex] <- Scenario[[1]]
}

clus <- makeCluster(ncores, outfile = outfile)
clusterExport(clus, varlist = c("timeMapsNoDispersal", "timeMapsDispersal",
                                "test", "nsubsamp", "nrep", "reptype","result_dir",
                                "ncores", "consecutiveCheck", "setbinary",
                                "removeTempRasterFiles", "numberOfOnes", "numYear", 
                                "years", "numScenario", "Scenarios", "ScenarioNames", "BinToDec", 
                                "overlap", "largestNum", "rastertype", "format"))
clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(gtools))

Sys.time()
if (dispersalRan == "N") {
  out <- parLapply(clus, spp.list, function(x) timeMapsNoDispersal(x))
}

if ((dispersalStep == "Y") && (dispersalRan == "Y")) {
  for (w in 1:length(spp.list)) {
    FocSpec <- gsub("_", " ", spp.list[w])
    DispSpec <- grep(paste0("^", FocSpec, "$"), dispersal[, 1])
    if (length(DispSpec) == 0) {
      message(paste0("No dispersal rate values found for ", FocSpec, ": skipping dispersal rate analysis"))
    }
  }
  out <- parLapply(clus, spp.list, function(x) timeMapsDispersal(x))
}

Sys.time()
#frees computer cores
stopCluster(clus)