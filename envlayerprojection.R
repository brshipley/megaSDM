####envlayerprojection.R####
##ensures that all data are in the required projections

#Initializations---------------------------
#Loads the necessary packages
library(raster)
library(rgdal)

#Loads the necessary variables from "df"
test <- df[, "test"]
rastertype <- df[, "rastertype"]
original <- df[, "trainingarea"]
proj <- df[, "proj_trainingarea"]
originalsa <- df[, "studyarea"]
projsa <- df[, "proj_studyarea"]
desiredCRS <- df[, "desiredCRS"]
numScenario <- df[, "numScenario"]
resolution <- df[, "resolution"]
buff_dir <- df[, "buff_dir"]

UrbanAnalysis <- df[, "UrbanAnalysis"]
ProtectedAnalysis <- df[, "ProtectedAnalysis"]

if (numScenario > 0) {
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

setwd(proj)
#Raster type -> format dictionary
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

#Current Data---------------------------
#Creates a stack of training-area environmental raster files
trainingEnv <- list.files(path = original, pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
trainingEnv <- stack(trainingEnv)

#If the training area doesn't need to be clipped, makes the extent of the training stack the desired extent
if (TrainingAreaClip == "N") {
  TAExtent <- extent(trainingEnv)
  df[, "TrainClipLongitude"] <- paste0(ceiling(TAExtent[1]), ",", floor(TAExtent[2]))
  df[, "TrainClipLatitude"] <- paste0(ceiling(TAExtent[3]), ",", floor(TAExtent[4]))
}

#Projects each training-area environmental layer file to desired CRS
setwd(proj)
for (i in 1:nlayers(trainingEnv)) {
  trainbio <- projectRaster(trainingEnv[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
  writeRaster(trainbio, paste0(names(trainingEnv[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
}
print("Projecting Rasters (Training Area)")
studyEnv = list.files(path = originalsa, full.names = TRUE, pattern = paste0("\\", rastertype, "$"))

if (length(studyEnv) >= 1) {
  
  #Ensure that all study area rasters have the same projection
  ProjEnv <- rep(NA, len = length(studyEnv))
  for (i in 1:length(studyEnv)) {
    ProjEnv[i] <- as.character(crs(raster(studyEnv[[i]])))
  }
  
  ProjUnique <- unique(ProjEnv)
  if (length(ProjUnique) > 1) {
    stop("Not all of the study area environmental rasters are in the same projection")
  }
  
  studyEnv2 <- stack(studyEnv) 
  setwd(projsa)
  if (is.na(crs(studyEnv2))) {
    stop("study raster crs = NA: Ensure all raster layers have a defined coordinate projection")
  }
  #Converts each study area environmental layer file to the correct CRS
  for (i in 1:nlayers(studyEnv2)) {
    secbio <- projectRaster(studyEnv2[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(secbio, paste0(names(studyEnv2[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
  }
  print("Projecting Rasters (Study Area)")
} else {
  ClipEnvDataStep <- "Y"
}

#Future Data---------------------------
#Converts each forecasted/hindcasted scenario into the correct CRS
if (numScenario > 0) {
  #Finds the file faths for the scenarios
  setwd(predictenv)
  predictenvdir <- list.dirs(path = predictenv, full.names = TRUE)
  
  for (i in 1:length(predictenvdir)) {
    correctDir <- list.dirs(path = predictenvdir[i], full.names = TRUE)
    #If the directory is an "end" directory with no subfolders
    if (length(correctDir) == 1) {
      #Stacks climate data
      setwd(correctDir)
      futureenv <- list.files(full.names = TRUE, pattern = paste0("\\", rastertype, "$"))
      if (length(futureenv) == 0) {
        stop(paste0("Forecasted/hindcasted climate rasters not found or not in ", rastertype, " format!", "\n", 
                    "  Ensure that all predicted climate rasters are in the correct locations"))
      } else if (length(futureenv) != nlayers(trainingEnv)) {
        stop(paste0("Number of forecasted/hindcasted climate rasters does not equal number of current climate variables", "\n",
                    "  Ensure that every predicted climate raster is in the correct location"))
      }
      
      #Ensures that all environmental rasters are in the same projection
      ProjEnv <- rep(NA, len = length(futureenv))
      for (i in 1:length(futureenv)) {
        ProjEnv[i] <- as.character(crs(raster(futureenv[[i]])))
      }
      
      ProjUnique <- unique(ProjEnv)
      if (length(ProjUnique) > 1) {
        stop("Not all of the forecasted/hindcasted environmental rasters are in the same projection")
      }
      
      futureenv <- stack(futureenv)
      #If prediction layers don't have crs, prints error
      if (is.na(crs(futureenv))) {
        stop(" forecasted/hindcasted raster crs = NA: Ensure all raster layers have a defined coordinate projection")
      }
      setwd(projpredictenv)
      Directories <- unlist(strsplit(correctDir, "/"))
      newDir <- Directories[length(Directories) - 1]
      if (!dir.exists(paste0(projpredictenv, "/", newDir))) {
        dir.create(newDir)
      }
      setwd(newDir)
      if (!dir.exists(paste0(getwd(), "/", Directories[length(Directories)]))) {
        dir.create(Directories[length(Directories)])
      }
      setwd(Directories[length(Directories)])
      FutFilesRemove <- list.files()
      unlink(FutFilesRemove)
      #Projects and writes forecasted/hindcasted data
      for (j in 1:nlayers(futureenv)) {
        secbio <- projectRaster(futureenv[[j]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
        writeRaster(secbio, paste0(names(futureenv[[j]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
      }
    }
  }
  print("Projecting Rasters (Future Env)")
}

#Buffer Rasters---------------------------
#Projects buffer rasters (if provided)
BuffRastList <- list.files(buff_dir, pattern = paste0(rastertype, "$"), full.names = TRUE)
if (length(BuffRastList) > 0) {
  for(i in 1:length(BuffRastList)) {
    currast <- raster(BuffRastList[i])
    currast <- projectRaster(currast, crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(currast, paste0(names(currast), rastertype), format = format, overwrite = TRUE, prj = TRUE)
    
  }
  print(paste0("Projecting Buffer Rasters"))
  rm(currast)
}

#Urban + Protected Analysis
#Projects protected area shapefile(s) into desired CRS
if (df$ProtectedAnalysis == "Y") {
  setwd(protected_dir)
  ProtectedAreas = list.files(path = protected_dir, full.names = TRUE, pattern = paste0("\\", ".shp", "$"))
  protected <- shapefile(ProtectedAreas)
  projectedProtected <- spTransform(protected, desiredCRS)
  setwd(proj_protected_dir)
  writeOGR(projectedProtected, dsn = proj_protected_dir, layer = "protected_areas_sa", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  print("Projecting Protected Areas")
}

#Projects urban intensity raster(s) into desired CRS
if (df$UrbanAnalysis == "Y") {  
  setwd(urbanized_dir)
  urbanized <- list.files(path = urbanized_dir, pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
  urbanized <- stack(urbanized)
  
  setwd(proj_urbanized_dir)
  for (i in 1:nlayers(urbanized)) {
    urbanizedRaster <- projectRaster(urbanized[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(urbanizedRaster, paste0(names(urbanized[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
  }
  print("Projecting Urban Layers")
}