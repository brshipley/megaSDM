#Script for megaSDM_run
#Initializations-------------------------------------
#sets the number of DLLs able to be loaded to 500
Sys.setenv("R_MAX_NUM_DLLS" = 500)

#set the working directory by navigating to and clicking on this file (or any other script file in the same directory)
setwd(dirname(file.choose()))

#Reads and formats the "config.txt" file

source("format.R")

library(raster)
library(parallel)
library(rgdal)

#Variable Creation------------------------------------
#creates variables for folder names

scripts <- df[, "scripts"]
test <- df[, "test"]
occurrences <- df[, "occurrences"]
result_dir <- df[, "result_dir"]

#creates variables for Y/N steps 
ClipEnvDataStep <- df[, "ClipEnvDataStep"]
CoordinateProjectionStep <- df[, "CoordinateProjectionStep"]
gbifstep <- df[, "gbifstep"]
OccurEnvFeaturesStep <- df[, "OccurEnvFeaturesStep"]
subsampleVarelaStep <- df[, "subsampleVarelaStep"]
speciesBufferStep <- df[, "speciesBufferStep"]
backgroundPointsStep <- df[, "backgroundPointsStep"]
dispersalStep <- df[, "dispersalStep"]
UrbanAnalysis <- df[, "UrbanAnalysis"]
ProtectedAnalysis <- df[, "ProtectedAnalysis"]
RichnessStep <- df[, "RichnessStep"]

#creates variables for data directories
buff_dir <- df[, "buff_dir"]
samples <- "samples"
train <- df[, "trainingarea"]
projtrain <- df[, "proj_trainingarea"]
if (projtrain == "NA") {
  projtrain <- train
}
sa <- df[, "studyarea"]
projsa <- df[, "proj_studyarea"]
if (projsa == "NA") {
  projsa <- sa
}
if (df[, "numScenario"] > 0) {
  projpredictenv <- df[, "proj_predictenv"]
  predictenv <- df[, "predictenv"]
}
if (ProtectedAnalysis == "Y") {
  protected_dir <- df[, "protected_dir"]
  proj_protected_dir <- df[, "proj_protected_dir"]
}
if (UrbanAnalysis == "Y") {
  urbanized_dir <- df[, "urbanized_dir"]
  proj_urbanized_dir <- df[, "proj_urbanized_dir"]
}
rastertype <- df[, "rastertype"]

#Raster Management------------------------------------
#This only needs to be run once, as long as study area and coordinate reference system remain constant
#Runs envlayerprojection.R if necessary: IF ALL OF THE DATA ARE NOT IN THE SAME PROJECTION OR NOT IN THE DESIRED PROJECTION,THIS SHOULD BE "Y"
#When examining another study area with the same training area and species, 
  #replace climate files in the "studyarea" directory and re-run projection (if needed) and clip

{
#If the units of meters are not specified in the desired projection, reformat desired projection to have units of meters  
projcond <- length(grep("\\+units *= *m", df[, "desiredCRS"])) == 0
unitscond <- length(grep("\\+units", df[, "desiredCRS"])) == 0

if (projcond && unitscond) {
  df[, "desiredCRS"] <- paste(df[, "desiredCRS"], " +units=m")
  CoordinateProjectionStep <- "Y"
} else if (projcond){
  df[, "desiredCRS"] <- sub("\\+units.*\\+","", df[, "desiredCRS"])
  df[, "desiredCRS"] <- paste(df[, "desiredCRS"], " +units=m")
  CoordinateProjectionStep <- "Y"
}

trainingEnv <- list.files(path = trainingarea, full.names = TRUE, pattern = paste0("\\", rastertype, "$"))

#Ensure that all training environmental rasters are in the same projection
ProjEnv <- rep(NA, len = length(trainingEnv))

for (i in 1:length(trainingEnv)) {
  ProjEnv[i] <- as.character(crs(raster(trainingEnv[[i]])))
}

ProjUnique <- unique(ProjEnv)
if (length(ProjUnique) > 1) {
  stop(message ("Error! Not all of the training environmental rasters are in the same projection!"))
}

#Ensure that the raster-type is correct
setwd(trainingarea)
trainstack <- stack(list.files(trainingarea, pattern = paste0(rastertype, "$")))
if (length (trainstack) == 0) {
  stop(message(paste0("No ", rastertype, " files found in Training Area Directory: add files or change rastertype")))
}

#Ensure that the rasters have some sort of projection system already
if (is.na(crs(trainstack))) {
  stop(message("crs = NA: Define a Coordinate Reference System for all raster layers"))
}

#Forces all input rasters to be in ".bil" form
  #This avoids issues that arise in MaxEnt's software with ".asc" files
if (rastertype != ".bil") {
  CoordinateProjectionStep <- "Y"
}

if (CoordinateProjectionStep == "Y") {
  setwd(scripts)
  source("envlayerprojection.R")
  print("Layers projection complete!")
} else {
  setwd(train)
  
  #If the directories for projected training area rasters are empty, copy over the original training area rasters 
  if (length(dir(projtrain)) == 0) {
    for (i in 1:length(dir(train))) {
      file.copy(from = dir(train)[i], to = paste0(projtrain, "/", dir(train)[i]))
    }
  }
  
  #Define the desired coordinate reference system as the CRS of the copied (or already projected) rasters
  df[, "desiredCRS"] <- as.character(crs(raster(dir(projtrain)[1])))
  print(paste0("desired CRS is:", df[, "desiredCRS"], ""))
  
  #If the directories for projected study area rasters are empty, copy over the original study area rasters 
  setwd(sa)
  if (length(dir(projsa)) == 0) {
    for (i in 1:length(dir(sa))) {
      file.copy(from = dir(sa)[i], to = paste0(projsa, "/", dir(sa)[i]))
    }
  }
  
  #If the directories for projected future (or past) rasters are empty, copy over the original rasters
  if (df[, "numScenario"] > 0) {
    setwd(predictenv)
    if (length(dir(projpredictenv)) == 0) {
      for (i in 1:length(dir(predictenv))) {
        file.copy(from = dir(predictenv)[i], to = paste0(projpredictenv, "/", dir(predictenv)[i]))
      }
    }
  }
  
  #If the directory for projected shapefiles of protected areas is empty, copy over the original shapefiles
  if (ProtectedAnalysis == "Y") {
    setwd(protected_dir)
    if (length(dir(proj_protected_dir)) == 0) {
      for (i in 1:length(dir(protected_dir))) {
        file.copy(from = dir(protected_dir)[i], to = paste0(proj_protected_dir, "/", dir(protected_dir)[i]))
      }
    }
  }
  
  #If the directory for projected urbanization rasters is empty, copy over the original rasters
  if (UrbanAnalysis == "Y") {  
    setwd(urbanized_dir)
    if (length(dir(proj_urbanized_dir)) == 0) {
      for (i in 1:length(dir(urbanized_dir))) {
        file.copy(from = dir(urbanized_dir)[i], to = paste0(proj_urbanized_dir, "/", dir(urbanized_dir)[i]))
      }
    }
  }
  
}

gc()

#clips the rasters to the study site if necessary

if (CoordinateProjectionStep == "Y") {
  ClipEnvDataStep <- "Y"
}


if (ClipEnvDataStep == "Y") {
  setwd(scripts)
  source("clip.R")
  print("Clip complete!")
} else {
  setwd(projsa)
  
  #Defines a standard raster to which the other rasters will be resampled
  env_files <- list.files(path = ".", pattern = paste0("\\.bil$"), full.names = TRUE)
  resample_raster <- raster(env_files[1])
  
  #Are there future/past environmental layers? If so, resample them to standard resolution
  if (df[, "numScenario"] > 0) {
    setwd(projpredictenv)
    predictenvdir <- list.dirs(path = projpredictenv, full.names = T)
    
    for (i in 1:length(predictenvdir)) {
      correctDir <- list.dirs(path = predictenvdir[i], full.names = T)
      
      if (length(correctDir) == 1) {
        
        #Makes a list of all future environmental layer files
        setwd(correctDir)
        future <- list.files(path = ".", pattern = paste0("\\.bil$"), full.names = TRUE)
        
        #Resamples to the environmental layer raster
        future_res <- c()
        for (i in 1:length(future)) {
          future_raster <- raster(future[[i]])
          if ((future_raster@extent != resample_raster@extent) | (future_raster@ncols != resample_raster@ncols)){
            future_res <- c(future_res, resample(future_raster, resample_raster, method = "ngb"))
          }
        }
        
        #Writes the rasters
        if (length(future_res) > 0) {
          future_res <- stack(future_res)
          setwd(projpredictenv)
          for (k in 1:nlayers(future_res)){
              writeRaster(future_res[[k]], filename = paste0((names(future_res)[[k]]), ".bil"), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
          }
        }
      }
    }
  }
  
  #Resample urbanization rasters (if provided)
  if (UrbanAnalysis == "Y") {
    setwd(proj_urbanized_dir)
    urb_files <- list.files(path = ".", pattern = paste0("\\.bil$"), full.names = TRUE)
    urban_res <- c()
    
    for (i in 1:length(urb_files)) {
      urban_raster <- raster(urb_files[[i]])
      if ((urban_raster@extent != resample_raster@extent) | (urban_raster@ncols != resample_raster@ncols)){
        urban_res <- c(urban_res, resample(urban_raster, resample_raster, method = "ngb"))
      }
    }
    
    #Writes the rasters
    if (length(urban_res) > 0){
      urban_res <- stack(urban_res)
      setwd(proj_urbanized_dir)
      for (j in 1:nlayers(urban_res)){
        if (length(grep(paste(years,collapse="|"), names(urban_res))) == nlayers(urban_res)){
          writeRaster(urban_res[[j]], filename = paste0("SA_urb", years[j], ".bil", sep=""), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
        } else {
          writeRaster(urban_res[[j]], filename = paste0("SA_urb", j, ".bil", sep=""), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
        }
      }
    }
  }
  gc()
}
}

#Occurrence Data Management---------------------------- 
#runs rgbif.R if necessary
if (gbifstep == "Y") {
  setwd(scripts)
  source("rgbif.R")
} else {
  #Copies maxent.jar from the data directory to the "test" folder
  setwd(test)
  file.copy(paste0(occurrences, "/maxent.jar"), test)
  
  #Creates a list of occurrence CSV files found, copies them into "test" folder 
  setwd(occurrences)
  ListSpp <- list.files(pattern = '\\.csv', full.names = TRUE)
  nspp <- length(ListSpp)
  for(j in 1:nspp) {
    file.copy(paste0(occurrences, "/", ListSpp[j]), test)
  }
  
  setwd(test)
  dir.create("species")
  ListSpp <- list.files(pattern = '\\.csv', full.names = TRUE)
  ListSpp <- ListSpp[1:length(ListSpp)]
  
  #Re-formats the occurrence layers for use in subsequent steps
  setwd(test)
  CreateCounts <- c()
  ListofSpeciesFiles <- c()
  
  for (G in 1:length(ListSpp)) {
    FocusSpecies <- ListSpp[G]
    SpeciesName <- gsub("_", " ", substr(FocusSpecies, 3, (nchar(FocusSpecies) - 4)))
    CurSpp <- read.csv(FocusSpecies)
    CurSpp <- data.frame(lapply(CurSpp, as.character), stringsAsFactors = FALSE)
    
    #If scientific names of the occurrences already have underscores
    if (length(grep(SpeciesName, CurSpp[1, ])) > 0) {
      names(CurSpp)[grep(SpeciesName, CurSpp[1, ])[1]] <- "Species"
    } else {
      SpeciesName <- substr(FocusSpecies, 3, (nchar(FocusSpecies) - 4))
      names(CurSpp)[grep(SpeciesName, CurSpp[1, ])[1]] <- "Species"
    }
    
    if (tolower(names(CurSpp)[1]) == "x" && length(grep("^y$", tolower(names(CurSpp)))) == 0) {
      CurSpp <- CurSpp[, 2:ncol(CurSpp)]
    }
    
    #Locate the coordinates of the provided occurrence data
    names(CurSpp)[c(grep("lon", tolower(names(CurSpp))), grep("^x$", tolower(names(CurSpp))))] <- "Longitude"
    names(CurSpp)[c(grep("lat",tolower(names(CurSpp))), grep("^y$", tolower(names(CurSpp))))] <- "Latitude"
    CurSpp <- CurSpp[, c("Species", "Longitude", "Latitude")]
    s <- as.character(CurSpp$Species[1])
    s <- gsub(" ", "_" , s)
    write.csv(CurSpp, file = paste0("species/", s, ".csv"), row.names = FALSE)
    CreateCounts <- c(CreateCounts, (nrow(unique(CurSpp))))
    ListofSpeciesFiles <- c(ListofSpeciesFiles, SpeciesName)
  }
  
  #Writes the species counts file if it hasn't already been written
  CountsDF <- data.frame(Species = ListofSpeciesFiles, Occurrences = CreateCounts)
  write.csv(CountsDF, file = df[, "counts"], row.names = FALSE)
}

#Lists occurrence csv files
all1 <- list.files(path = paste0(test, "/species/"), pattern = "\\.csv")
spplist <- read.csv(df[, "spplist"], stringsAsFactors = FALSE)

if(length(all1) == 0) {
  stop(paste0("No species .csv files found in occurrences directory -- ", occurrences))
}

#Collect species occurrence counts
counts <- read.csv(df[, "counts"])
counts <- counts[!is.na(counts$Species), ]
counts$Species <- as.character(counts$Species)
counts$Species <- gsub(" ", "_", counts$Species)
counts <- counts[order(counts$Occurrences, decreasing = FALSE), ]

#Filters out unwanted species and list which species to evaluate
spp_total <- c()
spp_counts <- c()
setwd(occurrences)

#Remove species with < 3 occurrence points, print warning if species has < 20 occurrence points
print("The following species will be evaluated:")
for (i in 1:nrow(spplist)) {
  cur <- spplist[i,2]
  row <- (grep(paste0(gsub(" ", "_", cur), ".csv$"), all1))
  rowCount <- (grep(paste0(gsub(" ", "_", cur), "$"), counts$Species))
  if (length(row) > 0) {
    if (counts$Occurrences[rowCount] <= 3) {
      message(paste0("Warning! species ", cur, " has too few occurrences to model (", counts$Occurrences[rowCount],  ") and will be removed from the analysis"))
    } else if (counts$Occurrences[rowCount] <= 20) {
      message(paste0("Warning! ", cur, " has only ", counts$Occurrences[rowCount], " unique occurrences:"))
      message(paste0("       SDM Analysis may be unreliable"))
      spp_total <- c(spp_total,all1[row])
      spp_counts <- c(spp_counts, counts$Occurrences[rowCount])
      print(paste0("   ", cur))
    } else {
      spp_total <- c(spp_total,all1[row])
      spp_counts <- c(spp_counts, counts$Occurrences[rowCount])
      print(paste0("   ", cur))
    }
  }
}

spp_total <- spp_total[order(spp_counts, decreasing = FALSE)]
all_spp <- matrix(spp_total, nrow = as.numeric(df[, "ncores"]), byrow = FALSE)

#Modelling--------------------------------------------
if (backgroundPointsStep == "Y") {
  if (length(grep("x", tolower(nbg))) == 0) {
    setwd(scripts)
    print("Running backgroundPoints1.R")
    source("backgroundPoints1.R")
  }
  gc()
}

#Subsamples occurrence points and creates background points within buffers
for (speciesBatchIndex in 1:ncol(all_spp)) { 
  print(Sys.time())
  spp_batch <- all_spp[, speciesBatchIndex] #get the current chunk
  df[,"ncores"] <- nrow(all_spp)
  
  #If there are repeats due to ncores not being a factor of species number, revise the last number of cores
  if (speciesBatchIndex == ncol(all_spp)) {
    if (length(unique(spp_total)) < df[, "ncores"]) {
      df[, "ncores"] <- length(unique(spp_total))
    } else {
      Remainder <- (length(unique(spp_total)) %% nrow(all_spp))
      if(Remainder > 0) {
        spp_batch <- spp_batch[1:Remainder]
        df[, "ncores"] <- length(spp_batch)
      }
    }
  }
  
  print(paste0("Currently running species set ", speciesBatchIndex, " of ", ncol(all_spp), ":"))
  df[, "dispersalRan"] <- "N"
  for (z in 1:length(unique(spp_batch))) {
    print(paste0("   ", unique(spp_batch)[z]))
  }
  
  #Creates a results folder for each species
  for (j in 1:length(spp_batch)) {
    curspec <- substr(spp_batch[j], 1, nchar(spp_batch[j]) - 4)
    dir.create(paste0(result_dir, "/", curspec))
  }
  
  #Extracts the environmental features from the occurrence (and background if given) data
  if (gbifstep == "Y") {
    OccurEnvFeaturesStep <- "Y"
  }
  
  if (OccurEnvFeaturesStep == "Y") {
    setwd(scripts)
    print("Running getOccurEnvFeatures.R...")
    source("getOccurEnvFeatures.R")
  }
  print(Sys.time())
  gc()
  
  #Subsamples occurrence data using methods developed by Varela et al. 2014
  if (subsampleVarelaStep == "Y") {
    setwd(scripts)
    print("Running subsampleOccur.R...")
    source("subsampleOccur.R")
  } else {
    #Lists the species worked through in this batch
    setwd(test)
    spp.list <- c()
    speciesWorked <- spp_batch
    for (s in 1:length(speciesWorked)) {
      spp.list <- c(spp.list, list.files(path = samples, full.names = TRUE, pattern = speciesWorked[s]))
    }
    spp.list <- unique(spp.list)
    print("   Will evaluate species:")
    print(spp.list)
    
    #For each species, use all occurrence points
    for (cur in 1:length(spp.list)) {
      cur2 <- na.omit((read.csv(spp.list[cur])))
      #pull the species' name from this file (either in the csv itself or through the file name)
      if (length(grep("_", spp.list[cur])) > 0) {
        s <- as.character(cur2$Species[1])
        s <- gsub(" ", "_" , s)
      } else {
        s <- c(unlist(strsplit(as.character(cur2$Species[1]), " ")))[1]
      }
      curdir <- paste0(samples, "/", s, sep = "")
      if (!dir.exists(curdir)) {
        dir.create(curdir)
      }
      #Write out the CSV file
      for (r in 1:nsubsamp) {
        write.csv(cur2, file = paste0(curdir, "/OccurrenceSamplePoints_", r, ".csv", sep = ""), row.names = FALSE)
      }
    }
  }
  print(Sys.time()) 
  gc()
  
  #Creates buffer rasters for background sampling
  if (speciesBufferStep == "Y") {
    #Get a list of buffers already in the buffer directory
    BuffersList <- substr(spp_batch, 1, nchar(spp_batch) - 4)
    BuffIndex <- grep(paste(BuffersList, collapse = "|"), list.files(buff_dir, pattern = paste0(".shp$")))
    #If buffers have already been generated, skip this step
    if (length(BuffIndex) != length(spp_batch)) {
      setwd(scripts)
      print("Running createBackgroundBuffers.R...")
      source("createBackgroundBuffers.R") 
      print(Sys.time())
    }
  }
  
  if (speciesBufferStep == "Y") {
    backgroundPointsStep <- "Y"
  }
  gc()
  
  #Samples background points from buffer (if desired)
  BufferFiles <- list.files(buff_dir, pattern = paste0(rastertype, "$|.shp$"))
  if (((length(BufferFiles) > 0) | (length(grep("x", tolower(nbg))) > 0)) && (backgroundPointsStep == "Y")) {
    #If buffer files exist or if the number of desired background points is species-dependent, run this step
    setwd(scripts)
    print("Running backgroundPoints2.R...")
    source("backgroundPoints2.R")
    print(Sys.time())
  } else if (backgroundPointsStep == "Y") {
    #Use background points generated in the first step
    for (s in 1:length(spp_batch)) {
      currentspec <- substr(spp_batch[s], 1, nchar(spp_batch[s]) - 4)
      dir.create(paste0(test, "/backgrounds/", currentspec))
      setwd(paste0(test, "/backgrounds/", currentspec))
      for (i in 1:nsubsamp) {
        Full_BGPoints <- read.csv(paste0(test,"/backgrounds/Train_Background_", i, ".csv"))
        write.csv(Full_BGPoints, file = paste0(getwd(), "/", currentspec, "_background_", i, ".csv"), row.names = FALSE)
      }
    }
  } else {
    #Use background points provided by user 
    setwd(paste0(test, "/backgrounds/"))
    BackgroundFiles <- list.files(paste0(getwd()), pattern=paste0(".csv$"), full.names = TRUE, recursive = TRUE)
    for (s in 1:length(spp_batch)) {
      currentspec <- substr(spp_batch[s], 1, nchar(spp_batch[s]) - 4)
      SpecIndex <- grep(currentspec, BackgroundFiles)
      if (!dir.exists(paste0(test, "/backgrounds/", currentspec))) {
        dir.create(paste0(test, "/backgrounds/", currentspec))
      }  
      FocusBGFiles <- BackgroundFiles[SpecIndex]
      if (length(FocusBGFiles) != nsubsamp) {
        message("Warning! the number of background files does not match the requested number of subsamples")
        FocusBGFiles <- rep(FocusBGFiles, len = nsubsamp)
        message("Resizing files to number of subsamples")
      }
      setwd(paste0(test, "/backgrounds/", currentspec))
      for (i in 1:nsubsamp) {
        Full_BGPoints <- read.csv(paste0(FocusBGFiles[i]))
        write.csv(Full_BGPoints, file = paste0(getwd(), "/", currentspec,"_background_", i, ".csv"), row.names = FALSE)
      }
    }
    gc()
  }
}

#Runs MaxEnt, creates SDMs and calculates statistics
for (speciesBatchIndex in 2:ncol(all_spp)) {
  print(Sys.time())
  spp_batch <- all_spp[, speciesBatchIndex]
  df[, "ncores"] <- nrow(all_spp)
  
  #If there are repeats due to ncores not being a factor of species number, revise the last number of cores
  if (speciesBatchIndex == ncol(all_spp)){
    if (length(unique(spp_total)) < df[, "ncores"]) {
      df[, "ncores"] <- length(unique(spp_total))
    } else {
      Remainder <- (length(unique(spp_total)) %% nrow(all_spp))
      if(Remainder > 0) {
        spp_batch <- spp_batch[1:Remainder]
        df[, "ncores"] <- length(spp_batch)
      }
    }
  }
  
  #Message
  print(paste0("Currently running species set ", speciesBatchIndex, " of ", ncol(all_spp), ":"))
  df[, "dispersalRan"] <- "N"
  for (z in 1:length(unique(spp_batch))) {
    print(paste0("   ", unique(spp_batch)[z]))
  }
  
  spp_orig <- spp_batch
  
  for (j in 1:length(spp_batch)) {
    curspec <- substr(spp_batch[j], 1, nchar(spp_batch[j]) - 4)
    dir.create(paste0(result_dir, "/", curspec))
  }
  
  #Creates SDMs using MaxEnt with statistics and output map displays
  if (df[, "nrep"] == 1) {
    setwd(scripts)
    print("Running maxent1replicate.R...")
    source("maxent1replicate.R")
  } else {
    setwd(scripts)
    print("Running maxent2replicates.R...")
    source("maxent2replicates.R")
  }
  
  gc()
  
  copySpeciesFiles <- function(originalSpp) {
    setwd(test)
    newspp <- originalSpp
    for (b in 1:length(newspp)) {
      newspp[b] <- substr(newspp[b], 1, (nchar(newspp[b]) - 4))
    }
    
    #For each species batch analysed
    for (j in 1:length(newspp)) {
      spp.name <- newspp[j]
      #Create directories in %result_dir%
      print(paste0("Copying results for species ", gsub("_", " ", spp.name), " into:"))
      print(paste0(" ", result_dir, "/", spp.name))
      dir.create(paste0(result_dir, "/", spp.name, "/lambdas"))
      dir.create(paste0(result_dir, "/", spp.name, "/samples"))
      dir.create(paste0(result_dir, "/", spp.name, "/backgrounds"))
      dir.create(paste0(result_dir, "/", spp.name, "/logs"))
      dir.create(paste0(result_dir, "/", spp.name, "/outputs"))
      if (numScenario > 0) {
        dir.create(paste0(result_dir, "/", spp.name, "/projections"))
      }
      fol <- list.files(path = paste0(test, "/outputs/", spp.name), full.names = TRUE)
      allresults <- c()
      #Copy the lambda files and the MaxEnt logs into result_dir
      for (f in 1:length(fol)) {
        directory = paste0("outputs/", spp.name, "/RUN_", f)
        if (alloutputs == "Y") {
          lambdas = list.files(path = directory, pattern = "\\.lambdas", full.names = TRUE)
          logs = list.files(path = directory, pattern = "\\.log", full.names = TRUE)
          if (length(lambdas) > 0) {
            for (k in 1:length(lambdas)) {
              file.rename(from = lambdas[k], to = paste0(result_dir, "/", spp.name, "/lambdas/RUN_", f, ".lambdas"))
              file.rename(from = logs[k], to = paste0(result_dir, "/", spp.name, "/logs/RUN_", f, ".log"))
            }
          }
        }
        res <- paste0(directory, "/maxentResults.csv")
        allresults <- c(res, allresults)
      }
      
      #Create maxentResults.csv in result_dir
      r <- c()
      for (f in 1:length(fol)) {
        r <- rbind(r, read.csv(allresults[f]))
      }
      write.csv(r, paste0(result_dir, "/", spp.name, "/maxentResults.csv"))
      
      #Copy samples, backgrounds, specific outputs, and speciic projection files (if desired)
      if (AllOutputs == "Y") {
        samples.csv = list.files(path = paste0(test, "/samples/", spp.name), full.names=TRUE)
        for (f in 1:length(samples.csv)) {
          file.copy(from = samples.csv[f], 
                    to = paste0(result_dir, "/", spp.name, "/samples"),
                    overwrite = TRUE,
                    recursive = TRUE)
        }
        
        bg.csv = list.files(path=paste0(test, "/backgrounds/", spp.name), full.names = TRUE)
        for (f in 1:length(bg.csv)) {
          file.copy(from = bg.csv[f], 
                    to = paste0(result_dir, "/", spp.name, "/backgrounds"),
                    overwrite = TRUE,
                    recursive = TRUE)
        }
      }
      
      if (AllOutputs == "Y"){
        file.copy(from = paste0(test, "/outputs/", spp.name), 
                  to = paste0(result_dir, "/", spp.name, "/outputs"),
                  overwrite = TRUE, recursive = TRUE)
        if (dir.exists(paste0(test, "/projections/", spp.name))) {
          if (numScenario > 0) {
            file.copy(from = paste0(test, "/projections/", spp.name),
                      to = paste0(result_dir, "/", spp.name, "/projections"),
                      overwrite = TRUE, recursive = TRUE)
          }
        }
      }
      
      #Delete the copied (or not) folders from the test directory
      unlink(paste0(test, "/outputs/", spp.name),recursive = TRUE)
      
      if (backgroundPointsStep == "Y") {
        unlink(paste0(test, "/backgrounds/", spp.name),recursive = TRUE)
      }
      if (subsampleVarelaStep == "Y") {
        unlink(paste0(test, "/samples/", spp.name),recursive = TRUE)
      }
      if (numScenario > 0) {
        unlink(paste0(test, "/projections/", spp.name), recursive = TRUE)
      }
    }
    
    if (speciesBatchIndex == ncol(all_spp)){
      #Reading in the Taxon-Species list
      specieslist <- read.csv(df[,"spplist"])
      speciesfolders <- list.dirs(result_dir, recursive = FALSE)
      speciesfolders <- speciesfolders[grep("_", speciesfolders)]
      
      #Making a list of all species that have folders
      sppfold <- c()
      for (i in 1:length(speciesfolders)){
        if (length(list.files(speciesfolders[i]) > 0)) {
          sppfold <- c(sppfold,speciesfolders[i])
        }
      }
      
      sppfoldlist <- c()
      for (i in 1:length(sppfold)) {
        foldsplit <- unlist(strsplit(sppfold[i], "/"))
        sppfoldlist <- c(sppfoldlist, foldsplit[length(foldsplit)])
        sppfoldlist[i] <- gsub("_", " ", sppfoldlist[i])
      }
      SppFold <- data.frame(Species = sppfoldlist)
      
      #Merging the two lists
      taxonlist <- merge(specieslist, SppFold, by.x = c(colnames(specieslist)[2]), by.y = "Species", all.y = TRUE)
      colnames(taxonlist)[1] <- "Species"
      taxonlist$Species <- as.character(taxonlist$Species)
      
      #Finding and removing species with AUC Values less than desired threshold
      DeleteSP <- c()
      for (sp in 1:nrow(taxonlist)) {
        focusspp <- gsub(" ", "_", taxonlist[sp, 1])
        setwd(paste0(result_dir, "/", focusspp))
        curfocus <- list.files(pattern = paste0("binary", rastertype, "$"))
        if (!length(curfocus) > 0){
          message(paste0(focusspp, " will be removed (no replicates with an AUC > ", aucval, ")"))
          DeleteSP <- c(DeleteSP, sp)
        }
      }
      
      if (length(DeleteSP) > 0) {
        AnalysedSppList <- taxonlist[-DeleteSP,]
      } else {
        AnalysedSppList <- taxonlist
      }
      
      write.csv(AnalysedSppList, paste0(result_dir, "/AnalysedSpecies.csv"), row.names = FALSE)
      print(paste0(" Writing CSV file of fully analysed species to: ", result_dir,"/AnalysedSpecies.csv"))
    }
  }
  
  if (length(unique(spp_batch)) == 0) {
    copySpeciesFiles(spp_orig)
    next
  } else if (length(unique(spp_batch)) < df[, "ncores"]) {
    df[, "ncores"] <- length(unique(spp_batch))
  }
  
  print(Sys.time())
    
  setwd(scripts)
  print("Running additionalStats.R...")
  source("additionalStats.R")

  print(Sys.time())
  
  if (numScenario > 0) {
    setwd(scripts)
    print("Running createTimeMaps.R...")
    source("createTimeMaps.R")
    print(Sys.time())
  }

  gc()
  
  #Incorporates dispersal rate into projected species distributions
  if (numScenario == 0) {
    dispersalStep == "N"
  }
  
  if (dispersalStep == "Y") {
    setwd(scripts)
    print("Running dispersalRate.R...")
    source("dispersalRate.R")

    setwd(scripts)
    print("Running additionalStats.R...")
    source("additionalStats.R")

    setwd(scripts)
    print("Running createTimeMaps.R...")
    source("createTimeMaps.R")
    print(Sys.time())
  }
  
  gc()
  
  #Creates PDFs of output maps and graphs 
  setwd(scripts)
  print("Running resultsToPDF.R...")
  source("resultsToPDF.R")
  
  alloutputs <- df[, "AllOutputs"]
  
  #File Management: moves all files in "test" to "results" directories, clears test directory
  copySpeciesFiles(unique(spp_orig))
  
  print(Sys.time())
  rm(spp_batch)
  gc(full = TRUE)
}

#Reset number of cores
df[,"ncores"] <- nrow(all_spp)

#creates Rasters and PDF Maps of Species Richness
if (RichnessStep == "Y") {
  setwd(scripts)
  source("createRichnessMaps.R")
}

#Finalizations----------------------------------------
unlink(test, recursive = TRUE)
print(paste0("All runs complete: ", Sys.time()))