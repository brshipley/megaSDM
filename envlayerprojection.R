#Initializations---------------------------
library(rgdal)
library(raster)

#Variables from config.txt
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
trainingEnv <- list.files(path = original, full.names = TRUE, pattern = paste0("\\", rastertype, "$"))
trainingEnv <- stack(trainingEnv)

#If the training area doesn't need to be clipped, make the extent of the training stack the desired extent
if (TrainingAreaClip == "N") {
  TAExtent <- extent(trainingEnv)
  df[, "TrainClipLongitude"] <- paste0(ceiling(TAExtent[1]), ",", floor(TAExtent[2]))
  df[, "TrainClipLatitude"] <- paste0(ceiling(TAExtent[3]), ",", floor(TAExtent[4]))
}

# Projects each training-area environmental layer file to desired CRS
setwd(proj)
for (i in 1:nlayers(trainingEnv)) {
  trainbio <- projectRaster(trainingEnv[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
  writeRaster(trainbio, paste0(names(trainingEnv[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
  print("Writing File (Training Area)")
}

studyEnv = list.files(path = originalsa, full.names = TRUE, pattern = paste0("\\", rastertype, "$"))

if (length(studyEnv) >= 1) {
  studyEnv2 <- stack(studyEnv)  
  # Converts each study area environmental layer file to the correct CRS
  setwd(projsa)
  
  for (i in 1:nlayers(studyEnv2)) {
    secbio <- projectRaster(studyEnv2[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(secbio, paste0(names(studyEnv2[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
    print("writing File (Study Area)")
  }
  
} else {
  ClipEnvDataStep <- "Y"
}

#Future Data---------------------------
#Converts each future scenario into the correct CRS
if (numScenario > 0) {
  setwd(predictenv)
  predictenvdir <- list.dirs(path = predictenv, full.names = TRUE)
  for (i in 1:length(predictenvdir)) {
    correctDir <- list.dirs(path = predictenvdir[i], full.names = TRUE)
    if (length(correctDir) == 1) {
      setwd(correctDir)
      futureenv <- list.files(correctDir, full.names = TRUE, pattern = paste0("\\", rastertype, "$"))
      futureenv <- stack(futureenv)
      setwd(projpredictenv)
      Directories <- unlist(strsplit(correctDir, "/"))
      newDir <- Directories[length(Directories) - 1]
      if (!dir.exists(newDir)) {
        dir.create(newDir)
      }
      setwd(newDir)
      dir.create(Directories[length(Directories)])
      setwd(Directories[length(Directories)])
      FutFilesRemove <- list.files()
      unlink(FutFilesRemove)
      for (j in 1:nlayers(futureenv)) {
        secbio <- projectRaster(futureenv[[j]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
        writeRaster(secbio, paste0(names(futureenv[[j]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
        print("Writing File (Future Environment)")
      }
    }
  }
}

#Buffer Rasters---------------------------
#Projects buffer rasters (if provided)
BuffRastList <- list.files(buff_dir , pattern = paste0(rastertype, "$"), full.names = TRUE)
if (length(BuffRastList) > 0) {
  for(i in 1:length(BuffRastList)) {
    currast <- raster(BuffRastList[i])
    currast <- projectRaster(currast, crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(currast, paste0(names(currast), rastertype), format = format, overwrite = TRUE, prj = TRUE)
    print(paste0("Writing Buffer Rasters"))
  }
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
  print("Writing File (Protected Areas)")
}

#Projects urban intensity raster(s) into desired CRS
if (df$UrbanAnalysis == "Y") {  
  setwd(urbanized_dir)
  urbanized <- list.files(pattern=paste0("\\", rastertype, "$"), full.names=TRUE)
  urbanized <- stack(urbanized)
  
  setwd(proj_urbanized_dir)
  for (i in 1:nlayers(urbanized)) {
    urbanizedRaster <- projectRaster(urbanized[[i]], crs = desiredCRS, res = c(resolution, resolution), method = "ngb")
    writeRaster(urbanizedRaster, paste0(names(urbanized[[i]]), ".bil"), format = "EHdr", overwrite = TRUE, prj = TRUE)
    print("Writing File (Urban Areas)")
  }
}