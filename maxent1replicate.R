####maxent1replicate.R####
##This is the main maxent modeling script (for 1 replicate)
##It does both the initial modeling and the projections

#Initializations------------------------
#Loads the necessary packages
library(gtools)
library(parallel)

#Loads the necessary variables from "df"
test <- df[, "test"]
proj <- df[, "proj_studyarea"]
samples <- "samples"
numScenario <- as.numeric(df[, "numScenario"])
numYear <- as.numeric(df[, "numYear"])

if (numScenario > 0) {
  predictenv <- df[, "proj_predictenv"]
} else {
  predictenv <- c()
}

#Raster Type --> Method dictionary
rastertype <- df[, "rastertype"]

if (rastertype == ".asc") {
  format <- "ascii"
} else if (rastertype == ".bil") {
  format <- "EHdr"
} else if (rastertype == ".grd") {
  format <- "raster"
} else if (rastertype == ".tif") {
  format <- "GTiff"
} else if (rastertype == ".img") {
  format <- "HFA"
} else if (rastertype == ".rst") {
  format <- "IDRISI"
} else {
  message("Error: Raster type unknown")
}

#number of cores to use in parallel
ncores <- as.numeric(df[, "ncores"]) 
#name of output file for parallel
outfile <- "modelout.txt"
#number of times to repeat ENTIRE PROCESS
nsubsamp <- as.numeric(df[, "nsubsamp"]) 
#number of replications WITHIN MAXENT (crossvalidation/subsampling/etc)
nrep <- as.numeric(df[, "nrep"])
#type of replicates (Subsample, etc.; see MAXENT manual)
reptype <- df[, "reptype"]
#percentage of occurrences to hold back for cross-validation/accuracy testing
test_percent <- df[, "test_percent"]
#whether or not to use hinge features when modelling
hinge <- tolower(df[, "hinge"])
#regularization parameter (penalty for more complex models)
regularization <- as.numeric(df[, "regularization"])
#prefix string for categorical data
Categorical <- df[, "Categorical"]
#Which threshold to use
ThreshMethod <- df[, "threshold"]
#result directory
result_dir <- df[, "result_dir"]
#AUC cutoff value. We will only use models with AUC values greater than this value.
aucval <- df[, "aucval"]
#number of projections used
nproj <- ((numScenario) * (numYear - 1) + 1)

DeleteIndex <- c()

UrbanAnalysis  <- df[, "UrbanAnalysis"]
ProtectedAnalysis <- df[, "ProtectedAnalysis"]

if (UrbanAnalysis == "Y") {
  proj_urbanized_dir <- df[, "proj_urbanized_dir"]
  setwd(proj_urbanized_dir)
  urblist <- list.files(path = getwd(), pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
} else {
  proj_urbanized_dir <- c()
  urblist <- c()
}

if (ProtectedAnalysis == "Y") {
  proj_protected_dir <- df[, "proj_protected_dir"]
} else {
  proj_protected_dir <- c()
}

#Should all secondary MaxEnt outputs be generated?
AllOutputs <- df[, "AllOutputs"]
if (AllOutputs == "Y") {
  alloutputs <- "true"
} else {
  alloutputs <- "false"
}

#Creates results directory
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

#Functions------------------------------
#Calculates total number of cells in the raster
getSize <- function(raster) {
  return(freq(raster, digits = 0, value = 1, useNA = 'no', progress = ''))
}

#Calculates the number of pixels covered by one raster but not the other
t1nott2 <- function(t1, t2) {
  tnew <- t1
  tnew[tnew == 1] <- 2
  Combined <- tnew + t2
  Combined[Combined != 2] <- 0
  Combined[Combined == 2] <- 1
  return(Combined)
}

#Calculates the number of pixels covered by both rasters
overlap <- function(t1, t2) {
  Combined <- t1 + t2
  Combined[Combined != 2] <- 0
  Combined[Combined == 2] <- 1
  return(Combined)
}

#Calculates the number of pixels covered by a protected area
getProtected <- function(raster, decade, prot_files) {
  correctIndex <- 0
  if (length(prot_files) > 1) {
    correctIndex <- grep(decade, prot_files)
  } else {
    correctIndex <- 1
  }
  protected <- shapefile(prot_files[correctIndex])
  m <- mask(raster, protected)
  return(freq(m, digits = 0, value = 1, useNA = 'no', progress = ''))
}

#Calculates the number of pixels covered by an urban area
getUrbanized <- function(focusraster, decade, urban_files, urbanized) {
  #get correct decade
  correctIndex <- 0
  if (length(urban_files) > 1) {
    correctIndex <- grep(decade, urban_files)
  } else {
    correctIndex <- 1
  }
  cur <- urbanized[[correctIndex]]
  if (length(unique(cur)) > 2) {
    UrbThresh <- median(cur, na.rm = TRUE)
  } else {
    UrbThresh <- cellStats(cur, stat = max)
  }
  cur2 <- cur >= UrbThresh
  cur2 <- crop(cur2, extent(focusraster))
  focusraster <- crop(focusraster, extent(cur2))
  m <- mask(focusraster, cur2, inverse = TRUE, maskvalue = 1, updateValue = 0)
  return (freq(m, digits = 0, value = 1, useNA = 'no', progress = ''))
}

#Calculates the centroid of the binary distribution
getCentroid <- function(CentRaster) {
  #A matrix with three columns: x, y, and v (value)
  points <- rasterToPoints(CentRaster, fun = function(x){x == 1}, spatial = FALSE)
  
  #average latitude (y)
  Clat <- mean(points[, 2], na.rm = TRUE)
  
  #average longitude (x)
  Clong <- mean(points[, 1], na.rm = TRUE)
  
  #returns the longitude & latitude of the Centroid
  return(c(Clong, Clat))
}

#Creates binary rasters out of ensembled MaxEnt outputs
threshold <- function(path, rasters, replicates, Scenario, decade) {
  rasters <- matrix(rasters, nrow = nrep)
  setwd(test)
  rasterNames <-  c()
  #Lists all of the folders within "outputs"
  runs <- list.dirs(path = paste0("outputs/", path), full.names = TRUE, 
                    recursive = FALSE)
  nruns <- length(runs)
  ensemble.stack <- c()
  for (i in 1:nruns) {
    curmodel <- runs[i]
    #Reads in the results file for each of the runs
    results <- read.csv(paste0(curmodel, "/maxentResults.csv"))
    for (j in 1:replicates) {
      #Reads the threshold value for that results file
      if (!is.na(as.numeric(ThreshMethod))) {
        thresh <- as.numeric(ThreshMethod)
      } else {
        thresh <- results[j, ThreshMethod]
      }
      if(is.null(thresh)) {
        stop("The threshold for creating binary maps is not found in 'maxentResults.csv'. Revise in config")
      }
      tryCatch({
        #Finds the AUC value of the run
        auc <- results$Test.AUC[j]
        #If the AUC value is below the AUC threshold, the run is not used
        if (auc >= aucval) {
          temp <- raster(rasters[j, i]) 
          temp[temp >= thresh] <- 1
          temp[temp < thresh] <- 0
          rasterNames <- c(rasterNames, filename(temp))
          #Creates a stack of output runs to be ensembled
          ensemble.stack <- stack(c(ensemble.stack, temp))
        }
      }, error = function(err) {
        print(paste("MY_ERROR: ", path, " ", err, " i= ", i, " j= ", j, 
                    " Scenario= ", Scenario, " Decade= ", decade))
      })
    }
  }
  
  if (length(ensemble.stack) > 0) {
    #Takes the mean of the binary rasters
    #If the mean is >= 0.5 (more than half of the replicates show presence), set to 1
    ensemble.calc <- (mean(ensemble.stack))
    rasterNames <- c(rasterNames, filename(ensemble.calc))
    ensemble.calc[ensemble.calc >= 0.5] <- 1
    ensemble.calc[ensemble.calc < 0.5] <- 0
    
    #Writes out binary rasters
    setwd(result_dir)
    if(decade == currentYear) {
      writeRaster(ensemble.calc, 
                  filename = paste0(path, "/", decade, "_binary", rastertype),
                  overwrite = TRUE,
                  format = format,
                  prj = TRUE)
    } else {
      if (!dir.exists(paste0(path, "/", Scenario))) {
        dir.create(paste0(path, "/", Scenario))
      }
      writeRaster(ensemble.calc, 
                  filename = paste0(path, "/", Scenario,"/", decade, "_", Scenario, "_binary", rastertype),
                  overwrite = TRUE, 
                  format = format,
                  prj = TRUE)
      
    }
    setwd(test)
    return(ensemble.calc)
  } else {
    message(paste0("No replicates of ", path, " had a test AUC value above ", aucval))
  } 
  rm(temp, ensemble.stack)
  gc()
}

#Ensembles multiple runs on non-binary presence outputs from MaxEnt
medianensemble <- function(path, rasters, replicates, Scenario, decade) {
  rasters <- matrix(rasters, nrow = nrep)
  setwd(test)
  #Lists all of the folders within "outputs"
  message("starting median ensemble")
  runs <- list.dirs(path = paste0("outputs/", path), 
                    full.names = TRUE, 
                    recursive = FALSE)
  nruns <- length(runs)
  ensemble.stack <- c()
  for (i in 1:nruns) {
    curmodel <- runs[i]
    #Reads in the results file for each of the runs
    results <- read.csv(paste0(curmodel, "/maxentResults.csv"))
    for (j in 1:replicates) {
      tryCatch({
        #Finds the AUC value of the run
        auc <- results$Training.AUC[j]
        #If the AUC value is below the AUC threshold, the run is not used
        if (auc >= aucval) {
          temp <- raster(rasters[j, i])
          #Creates a stack of output runs to be ensembled
          ensemble.stack <- stack(c(ensemble.stack, temp))
        }
      }, error = function(err) {
        print(paste("MY_ERROR: ", path, " ", err, " i= ", i, "j=", j, 
                    " Scenario= ", Scenario, " Decade= ", decade))
      })
    }
  }
  #Calculates the median of the stacked rasters
  if (nruns > 1) {
    ensemble.calc <- calc(ensemble.stack, median, na.rm = TRUE)
  } else {
    ensemble.calc <- ensemble.stack[[1]]
  }
  setwd(result_dir)
  if (decade == currentYear) {
    writeRaster(ensemble.calc, 
                filename = paste0(path, "/", decade, "_ensembled", rastertype), 
                overwrite = TRUE, 
                format = format,
                prj = TRUE)
  } else {
    writeRaster(ensemble.calc, 
                filename = paste0(path, "/", Scenario, "/", decade, "_", Scenario, "_ensembled", rastertype), 
                overwrite = TRUE, 
                format = format,
                prj = TRUE)
  }
  rm(ensemble.stack)
  gc()
  setwd(test)
  return(ensemble.calc)
}

#Creates the stats file
getStats <- function(cur, j, decade, binary, stats, modern.size, modern.binary) {
  #Adds all of the variables into the stats data frame 
  stats$Projection[j] <- cur
  stats$NumberCells[j] <- getSize(binary)
  stats$CellChange[j] <- getSize(binary) - modern.size
  t1nott2_raster <- t1nott2(modern.binary, binary)
  t2nott1_raster <- t1nott2(binary, modern.binary)
  stats$T1notT2[j] <- getSize(t1nott2_raster)
  stats$T2notT1[j] <- getSize(t2nott1_raster)
  overlap_raster <- overlap(modern.binary, binary)
  stats$Overlap[j] <- getSize(overlap_raster)
  binary.centroid <- getCentroid(binary)
  stats$CentroidX[j] <- binary.centroid[1]
  stats$CentroidY[j] <- binary.centroid[2]
  return(stats)
}

#Error handling
checkError <- function(result, spp.name, failed_runs, run) {
  if (grepl("error", result)) {
    cur <- paste0(spp.name, ", failed run: ", run)
    failed_runs <- c(failed_runs, cur)
    return(failed_runs)
  } else {
    return(failed_runs)
  }
}

#Projects algorithm onto study area/future climate layers
projectSpeciesHabitat <- function(CurSpp) {
  setwd(test)
  #Identifies the species/set of species being run
  spp.name <- CurSpp
  spp.name <- substr(spp.name, 9, nchar(spp.name)) #9 corresponds to "outputs/"
  print(spp.name)
  setwd(result_dir)
  dir.create(spp.name)
  setwd(test)
  #Creates a stats table to be filled in & passed to next scripts
  stats <- as.data.frame(cbind(Projection = rep(0, times = nproj), 
                               NumberCells = rep(0, times = nproj), 
                               CellChange = rep(0, times = nproj), 
                               T1notT2 = rep(0, times = nproj), 
                               T2notT1 = rep(0, nproj), 
                               Overlap = rep(0, times = nproj), 
                               CentroidX = rep(0, times = nproj), 
                               CentroidY = rep(0, times = nproj)))
  
  if (UrbanAnalysis == "Y") {
    stats$Urbanized <- rep(0, times = nproj)
  }
  
  if (ProtectedAnalysis == "Y") {
    stats$Protected <- rep(0, times = nproj)
  }
  
  #Gets modern data stats first to use for comparison  
  setwd(test)
  modern <- list.dirs(path = paste0("outputs/", spp.name), recursive = FALSE)
  modern.rasters <- c()
  #Creates and re-orders a list of modern rasters
  for (j in 1:length(modern)) {
    r <- list.files(path = modern[j], 
                    pattern = paste0((strsplit(proj, "/"))[[1]][length((strsplit(proj, "/"))[[1]])], '.asc$'), 
                    full.names = TRUE)
    modern.rasters <- c(modern.rasters, (r[1]))
  }
  modern.rasters <- mixedsort(modern.rasters)
  #Runs "threshold" function, which creates and ensembles binary maps
  modern.binary <- threshold(spp.name, modern.rasters, nrep, "", currentYear)
  gc()
  
  if (length(modern.binary) > 0) {
    #Runs medianensemble function, which ensembles maps using a median function
    modern.median <- medianensemble(spp.name, modern.rasters, nrep, "", currentYear)
    gc()
    
    #Fills in the stats row for the modern species
    stats$Projection[1] <- currentYear
    modern.size <- getSize(modern.binary)
    stats$NumberCells[1] <- modern.size
    modern.centroid <- getCentroid(modern.binary)
    stats$CentroidX[1] <- modern.centroid[1]
    stats$CentroidY[1] <- modern.centroid[2]
    if (ProtectedAnalysis == "Y") {
      setwd(proj_protected_dir)
      ProtectList <- paste0("/", list.files(path = ".", pattern = paste0("\\", ".shp", "$"), full.names = FALSE))
      stats$Protected[1] <- getProtected(modern.binary, currentYear, ProtectList)
    }
    if (UrbanAnalysis == "Y") {
      setwd(proj_urbanized_dir)
      urblist <- list.files(path = getwd(), pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
      urbanized <- stack(urblist)
      stats$Urbanized[1] <- getUrbanized(modern.binary, currentYear, urblist, urbanized)
    }
    #Generates projections for hindcasted/forecasted climate layers
    if (numScenario > 0) {
      #Lists the species to be projected
      spp.name <- substr(CurSpp, 9, nchar(CurSpp))
      print(spp.name)
      setwd(test)
      dir.create(paste0("projections"))
      dir.create(paste0("projections/", spp.name))
      models <- list.dirs(path = paste0("outputs/", spp.name), full.names = TRUE, recursive = FALSE)
      for (fileIndex in 1:length(futlistfull)) {
        #Creates folders for each scenario/future date combo
        futscenario <- strsplit(futlist[fileIndex], "/")[[1]][2]
        futdate <- strsplit(futlist[fileIndex], "/")[[1]][3]
        dir.create(paste0(test, "/projections/", spp.name, "/", futscenario))
        dir.create(paste0(test, "/projections/", spp.name, "/", futscenario, "/",futdate))
      }
      for (dataIndex in 1:length(futlistfull)) {
          for (i in 1:nsubsamp) {
            #THIS IS THE PROJECTIONS COMMAND>>
            setwd(test)
            system(paste0("java -mx900m -cp maxent.jar density.Project ", 
                          # location of lambdas file
                          models[i], "/", spp.name, ".lambdas ", 
                          #location of folder with map rasters
                          gsub(" ", "\\\\ ", futlistfull[dataIndex]),
                          # where to outputfiles (Run#)
                          " projections/", spp.name, futlist[dataIndex], "/", 
                          (strsplit(models[i], "/")[[1]][3]),
                          " noaskoverwrite nowarnings -a"))
          }
          setwd(test)
          #Parses out the scenario & date values
          futscenario <- strsplit(futlist[dataIndex], "/")[[1]][2]
          futdate <- strsplit(futlist[dataIndex], "/")[[1]][3]
          
          #Makes a list of the projections files to use for threshold/ensemble functions
          cur.rasters <- list.files(path = paste0(test, "/projections/", spp.name, "/", futscenario, "/", futdate), 
                                    pattern = paste0('\\.asc$'), 
                                    full.names = TRUE)
          cur.proj <- c()
          for (k in 1:length(cur.rasters)) {
            if (!grepl("clamping", cur.rasters[k])) {
              cur.proj <- c(cur.proj, cur.rasters[k])
            }
          }
          length(cur.proj)
          #Thresholds & gets medians of the projected rasters
          cur.binary <- threshold(spp.name, cur.proj, nrep, futscenario, futdate)
          gc()
          cur.median <- medianensemble(spp.name, cur.proj, nrep, futscenario, futdate)
          gc()
          #Fills in the stats table for the projected rasters
          stats <- getStats(paste0(futscenario, "_", futdate), dataIndex + 1, futdate, cur.binary, 
                            stats, modern.size, modern.binary)
          gc()
          
          #Conducts urban analysis (if requested)
          if (UrbanAnalysis  == "Y") {
            setwd(proj_urbanized_dir)
            urblist <- list.files(pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
            urbanized <- stack(urblist)
            stats$Urbanized[dataIndex + 1] <- getUrbanized(cur.binary, futdate, urblist, urbanized)
          }
          
          #Conducts protected analysis (if requested)
          if (ProtectedAnalysis == "Y") {
            setwd(proj_protected_dir)
            ProtectList <- list.files(pattern = ("\\.shp$"), full.names = TRUE)
            if (length(ProtectList) > 1) {
              stats$Protected[dataIndex + 1] <- getProtected(cur.binary, futdate, ProtectList)
            } else {
              stats$Protected[dataIndex + 1] <- getProtected(cur.binary, currentYear, ProtectList)
            }
          }
      }
    }
    setwd(paste0(result_dir, "/", spp.name))
    #Writes the stats table to be used later
    write.csv(stats, file = "Results.csv")
  } else {
    message(paste0("Removing ", spp.name, " from further analysis"))
  }
  rm(modern.binary, modern.rasters)
  gc()
}

#Trains and calculates the species distribution model for training data
maxent <- function(CurSpp) {
  spp.name <- substr(CurSpp, 9, nchar(CurSpp)) #9 corresponds to "samples/"
  print(spp.name)
  
  #Creates directories for habitat suitability models
  starttime <- Sys.time()
  print(paste0("Start time for ", spp.name, ": ", starttime))
  setwd(test)
  dir.create(paste0("outputs"))
  dir.create(paste0("outputs/", spp.name))
  all_runs <- list.files(path = CurSpp, full.names = TRUE, pattern = "\\.csv$")
  all_runs <- mixedsort(all_runs)
  failed_runs <- c() 
  #Models the species habitat suitability
  for (i in 1:nsubsamp) {
    run <- i
    dir.create(paste0("outputs/", spp.name, "/RUN_", run))
    SppRun <- all_runs[i]
    #THIS IS THE MODEL CREATION COMMAND>>
    model.out <- tryCatch({
      #can turn off a lot of output writing for the final experiment 
      #(jackknifing, write plot pngs) press help in maxent for details
      system(paste0("java -mx900m -jar maxent.jar -e ", "backgrounds/", spp.name, 
                    "/", spp.name, "_background_", run, ".csv -s ", SppRun, 
                    " -J -o outputs/", spp.name, "/RUN_", run, 
                    " noaskoverwrite logistic threshold -X ", test_percent, " replicates=", nrep,
                    " betamultiplier=", regularization, " writeclampgrid=", alloutputs, 
                    " writemess=", alloutputs, " nowarnings writeplotdata=", alloutputs ," -j ", 
                    gsub(" ", "\\\\ ", proj), " -a ", reptype, " hinge=", hinge, " togglelayertype=", Categorical))
    }, error = function(err) {
      print(paste("MY_ERROR: ", spp.name, " ", err))
      return(paste0("error: ", err))
    })
    failed_runs <- checkError(model.out, spp.name, failed_runs, 
                              paste0("model creation number ", i))
    if (!is.null(failed_runs)) {
      message(failed_runs)
    }
  }
  endtime <- Sys.time()
  print(paste0("End time for ", spp.name, ": ", endtime))
  return(failed_runs)
}

#Run------------------------------------------------------------------------
#NOTE: folders need to be in correct format: folder/scenario/time period
#Creates a list of locations for all scenarios/time periods
if (numScenario > 0) {
  predictenvdir <- list.dirs(path = predictenv, recursive = TRUE)
} else {
  predictenvdir <- c()
}

futlistfull <- c()
futlist <- c()
Scenarios <- c()

years <- c()
years <- rep(NA, length = numYear)
years <- as.numeric(config[grep("^Year", row.names(config)), ])
currentYear <- years[1]

#Formats forecasted/hindcasted climate data directories
if (numScenario > 0) {
  setwd(predictenv)
  #Lists all directories and subdirectories in "predictenv"
  directory <- paste0("/", list.dirs(path = ".", full.names = FALSE, recursive = TRUE))
  nfolders <- length(directory)
  #Selects only the addresses of the subdirectory locations
  fut <- directory[grep("/([^:]+)/", directory)]
  
  Scenarios <- rep(NA, length = numScenario)
  Scenarios <- config[grep("^Scenario", row.names(config)), ]
  
  for (s in 1:length(Scenarios)) {
    FutTest1 <- fut[grep(Scenarios[s], fut)]
    if (length(FutTest1) == 0) {
      stop(paste0("Error: folder for cliamte scenario '", Scenarios[s], "' not found within the 'proj_predictenv' folders"))
    }
    
    for (y in 2:length(years)) {
      FutTest2 <- grep(years[y], FutTest1)
      if (length(FutTest2) == 0) {
        stop(paste0("Error: folder for time period '", years[y], "' not found within ", Scenarios[s]))
      }
    }
  }
  
  if (dir.exists(paste0(proj, "/maxent.cache"))) {
    CacheFile <- list.dirs(proj)
    unlink(CacheFile[grep("maxent.cache", CacheFile)], recursive = TRUE)
  }
  
  #creates full directory addresses, including "predictenv" location
  for (i in 1:length(fut)) {
    for (ScenIndex in 1:numScenario) {
      if (grepl(Scenarios[ScenIndex], fut[i], ignore.case = TRUE)) {
        for(yearIndex in 1:numYear) {
          if (grepl(years[yearIndex], fut[i])) {
            futlist <- c(futlist, fut[i])
            futlistfull <- c(futlistfull, paste0(predictenv, fut[i]))
          }
        }
      }
    }
  }
}
futlist <- unique(futlist)
futlistfull <- unique(futlistfull)
setwd(test)

#Generates the species list for parallelization
setwd(test)
ListSpp <- c()
speciesWorked <- spp_batch

for (i in 1:length(speciesWorked)) {
  ListSpp[i] <- paste0(samples, "/", substr(speciesWorked[i], 1, (nchar(speciesWorked[i]) - 4)))
}
ListSpp <- unique(ListSpp)
print("    Will evaluate species:")
print(paste0("        ", spp_batch))

#Parallelization
clus <- makeCluster(ncores, outfile = outfile, setup_timeout = 0.5)
clusterExport(clus, varlist = c("projectSpeciesHabitat", "checkError", "maxent", "proj", "years",
                                "futlistfull", "test", "predictenv", "predictenvdir", "samples",   
                                "nsubsamp", "nrep", "reptype", "regularization", "hinge", "result_dir", "numScenario",
                                "getCentroid", "getProtected", "getSize", "DeleteIndex",
                                "getStats", "getUrbanized", "medianensemble", "ListSpp", "alloutputs", 
                                "overlap", "t1nott2", "threshold", "proj_protected_dir",
                                "proj_urbanized_dir", "test", "aucval", "nproj", "ncores", "ThreshMethod",
                                "urblist", "UrbanAnalysis", "ProtectedAnalysis", "currentYear", 
                                "rastertype", "test_percent", "format", "futlist",  "Categorical"))
clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(gtools))

#Prints where the output files can be found
print("    Beginning maxent files in directory:")
print(paste0("        ", test, "/outputs/%Species_Name%"))

#Runs maxent
out <- parLapply(clus, ListSpp, function(x) maxent(x))
print("    Completing maxent files.")
Sys.time()

#Ensures that all species were modelled correctly
setwd(test)
ListSpp <- c()
speciesWorked <- spp_batch

for (i in 1:length(speciesWorked)) {
  ListSpp[i] <- paste0("outputs/", substr(speciesWorked[i], 1, (nchar(speciesWorked[i]) - 4)))
}
ListSpp <- unique(ListSpp)

#Removes species with low AUC values from subsequent analysis
AUCRetain <- c()
print(paste0("    Removing species with Test AUC Values < ", df[,"aucval"], " from subsequent analyses:"))
for (l in 1:length(ListSpp)) {
  runs <- list.dirs(path = paste0(ListSpp[l]), full.names = TRUE, 
                    recursive = FALSE)
  runs <- runs[grep("RUN", runs)]
  nruns <- length(runs)
  HighAUC <- c()
  for (i in 1:nruns) {
    curmodel <- runs[i]
    #Reads in the results file for each of the runs
    results <- read.csv(paste0(curmodel, "/maxentResults.csv"))
    #Reads the threshold value for that results file
    auc <- results$Test.AUC
    aucAbove <- which(auc > df[, "aucval"])
    HighAUC <- c(HighAUC, aucAbove)
  }
  if (length(HighAUC) > 0){
    AUCRetain <- c(AUCRetain, l)
  } else {
    message(paste0("Removing: ", ListSpp[l]))
  }
}

#Adds in the subsample index to the results file
for (l in 1:length(ListSpp)) {
  runs <- list.dirs(path = paste0(ListSpp[l]), full.names = TRUE, 
                    recursive = FALSE)
  runs <- runs[grep("RUN", runs)]
  nruns <- length(runs)
  HighAUC <- c()
  for (i in 1:nruns) {
    curmodel <- runs[i]
    #Reads in the results file for each of the runs
    results <- read.csv(paste0(curmodel, "/maxentResults.csv"))
    results$Subsample <- rep(i, (nrep + 1))
    write.csv(results, paste0(curmodel, "/maxentResults.csv"), row.names = FALSE)
  }
}


#Creates list of species to forecast/hindcast the maxent model
spp_batch <- spp_batch[c(AUCRetain)]
speciesWorked <- spp_batch
ListSpp <- rep(NA, length = length(speciesWorked))
for (i in 1:length(speciesWorked)) {
  ListSpp[i] <- paste0("outputs/", substr(speciesWorked[i], 1, nchar(speciesWorked[i]) - 4))
}
ListSpp <- unique(ListSpp)

#Parallelization
if (length(ListSpp) < ncores) {
  stopCluster(clus)
  clus <- makeCluster(length(ListSpp), outfile = outfile, setup_timeout = 0.5)
  clusterExport(clus, varlist = c("projectSpeciesHabitat", "checkError", "maxent", "proj",
                                  "futlistfull", "test", "predictenv", "samples", "predictenvdir", "years", 
                                  "nsubsamp", "nrep", "reptype", "result_dir", "numScenario",
                                  "getCentroid", "getProtected", "getSize", "DeleteIndex",
                                  "getStats", "getUrbanized", "medianensemble", "hinge", "alloutputs", 
                                  "overlap", "t1nott2", "threshold", "proj_protected_dir",
                                  "proj_urbanized_dir", "test", "aucval", "nproj", "ncores", "ListSpp",
                                  "UrbanAnalysis", "ProtectedAnalysis", "currentYear", "rastertype", 
                                  "urblist", "test_percent", "format", "futlist", "ThreshMethod"))
  clusterEvalQ(clus, library(raster))
  clusterEvalQ(clus, library(gtools))
  ListSpp <- ListSpp
}

nruns <- length(runs)
ensemble.stack <- c()

#Runs the "projectSpeciesHabitat" function
if (length(spp_batch) > 0) {
  print("    Beginning projection files in directory:")
  print(paste0("        ", test, "/projections/%Species_Name%"))
  Sys.time()
  out <- parLapply(clus, ListSpp, function(x) projectSpeciesHabitat(x))
}

print("    Completing project files.")
Sys.time()
#Frees computer cores
stopCluster(clus)