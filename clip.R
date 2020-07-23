####clip.R####
##Clips all spatial data to training/study area
#Initializations---------------------------
#Loads the necessary packages
library(raster)
library(rgdal)

#Loads the necessary variables from "df"
rastertype <- df[, "rastertype"]
trainingarea <- df[, "proj_trainingarea"]
minlat <- as.numeric(df[, "minlat"])
maxlat <- as.numeric(df[, "maxlat"])
minlong <- as.numeric(df[, "minlong"])
maxlong <- as.numeric(df[, "maxlong"])
studyarea <- df[, "proj_studyarea"]
UrbanAnalysis <- df[, "UrbanAnalysis"]
ProtectedAnalysis <- df[, "ProtectedAnalysis"]
numScenario <- df[, "numScenario"]
desiredCRS <- as.character(df[, "desiredCRS"])
defaultCRS <- df[, "defaultCRS"]
decimalLatitude <- df[, "decimalLatitude"]
decimalLongitude <- df[, "decimalLongitude"]

if (UrbanAnalysis == "Y") {
  proj_urbanized_dir <- df[, "proj_urbanized_dir"]
}
if (ProtectedAnalysis == "Y") {
  proj_protected_dir <- df[, "proj_protected_dir"]
}
if (numScenario > 0) {
  proj_predictenv_dir <- df[, "proj_predictenv"]
}

years <-  rep(NA, length = numYear)
years <-  as.numeric(config[grep("^Year", row.names(config)), ])


#Raster->format Dictionary
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

#Extent Checking----------------------------
#Converts Study Area extent boundaries to correct projection
Extent <- cbind(c(minlong, maxlong), c(minlat, maxlat))
colnames(Extent) <- c("long", "lat")
ExtentSP <- SpatialPoints(Extent, proj4string = CRS(defaultCRS))
ExtentProj <- spTransform(ExtentSP, desiredCRS)

#Makes data frame out of study area extent
NewExtent <-as.data.frame(ExtentProj@coords)
colnames(NewExtent) <- c("long", "lat")
#Weird projections cause the reprojected layers to have issues (dealt with here)
if (NewExtent$long[1] > NewExtent$long[2] | NewExtent$lat[1] > NewExtent$lat[2]) {
  NewExtent$long <- NewExtent$long[order(NewExtent$long, decreasing = FALSE)]
  NewExtent$lat <- NewExtent$lat[order(NewExtent$lat, decreasing = FALSE)]
  flipper <- "Y"
  message("Warning: The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
} else {
  flipper <- "N"
}

#Gets a list of files to clip from training area
setwd(trainingarea)
bioclim <- list.files(path = getwd(), pattern = paste0("\\.bil$"), full.names = FALSE)

#Determines if all of the extents are constant
ExtentClim <- matrix(data = NA, nrow = 4, ncol = length(bioclim))
for (i in 1:length(bioclim)) {
  ExtentClim[, i] <- extent(raster(bioclim[[i]]))[1:4]
}
ExtentUnique <- unique(ExtentClim, MARGIN = 2)

#Checks and deals with multiple raster extents
if (ncol(ExtentUnique) > 1) {
  message("Warning: the environmental rasters have different extents: megaSDM will use only the intersection of the rasters")
  df[, "TrainingAreaClip"] <- "Y"
  StackList <- c()
  
  #Stack all of the rasters with the same extents and create a list 
  for (j in 1:ncol(ExtentUnique)) {
    FocusStack <- c()
    for (p in 1:length(bioclim)) {
      if (identical(extent(raster(bioclim[[p]]))[1:4], ExtentUnique[, j])) {
        FocusStack <- c(FocusStack, p)
      }
    }
    StackList <- c(StackList, stack(bioclim[FocusStack]))
  }
  
  #Goes through each raster stack and intersects with the others
  #"pat" is template for resampling all of the rasters
  v <- 1
  pat <- StackList[[1]][[1]]
  while (v <= (length(StackList) - 1)) {
    SL_na <- c()
    if (v == 1) {
      #Resample all layers in the next rasterstack 
      for (d in 1:nlayers(StackList[[v + 1]])) {
        SL_na <- c(SL_na, resample(StackList[[v + 1]][[d]], pat, method = "ngb"))
      }
      
      #Intersect the two raster stacks
      StackList[[v + 1]] <- stack(SL_na)
      bioclim_int <- intersect(StackList[[v + 1]], StackList[[v]])
      print(paste0("Intersecting Rasters..."))
      bioclim <- stack(StackList[[v]], bioclim_int)
    } else {
      #Resample all layers in the next rasterstack
      for (d in 1:nlayers(StackList[[v + 1]])) {
        SL_na <- c(SL_na, resample(StackList[[v + 1]][[d]], pat, method = "ngb"))
      }
      
      #Intersect the two raster stacks
      StackList[[v + 1]] <- stack(SL_na)
      bioclim_int <- intersect(StackList[[v + 1]], bioclim)
      print(paste0("Intersecting Rasters..."))
      bioclim <- stack(bioclim, bioclim_int)
    }
    v <- v + 1
  }
} else {
  bioclim <- stack(bioclim) 
}

#Convert Training Area extent boundaries to correct projection
if (TrainingAreaClip == "N") {
  Extent_TA <- extent(bioclim[[1]])
  Extent_TA2 <- matrix(NA, nrow = 2, ncol = 2)
  Extent_TA2[, 1] <- Extent_TA[1:2]
  Extent_TA2[, 2] <- Extent_TA[3:4]
  Extent_TA <- Extent_TA2
} else {
  Extent_TA <- cbind(as.numeric(strsplit(TrainClipLongitude, ",")[[1]][1:2]),
                   as.numeric(strsplit(TrainClipLatitude, ",")[[1]][1:2]))
}

colnames(Extent_TA) <- c("long", "lat")

#If TrainingAreaClip is "N", then the CRS of the boudaries are in the correct projection
if (TrainingAreaClip == "N") {
  ExtentSP_TA <- SpatialPoints(Extent_TA, proj4string = CRS(desiredCRS))
} else {
  ExtentSP_TA <- SpatialPoints(Extent_TA, proj4string = CRS(defaultCRS))
}

#Convert the extent of the training area to the desiredCRS
ExtentProj_TA <- spTransform(ExtentSP_TA, crs(desiredCRS))
NewExtent_TA <- as.data.frame(ExtentProj_TA@coords)
colnames(NewExtent_TA) <- c("long", "lat")

#Weird projections cause the reprojected layers to have issues (dealt with here)
if ((NewExtent_TA$long[1] > NewExtent_TA$long[2]) | (NewExtent_TA$lat[1] > NewExtent_TA$lat[2])) {
  NewExtent_TA$long <- NewExtent_TA$long[order(NewExtent_TA$long, decreasing = FALSE)]
  NewExtent_TA$lat <- NewExtent_TA$lat[order(NewExtent_TA$lat, decreasing = FALSE)]
  flipperTA <- "Y"
  message("Warning: The desired CRS chosen may lead to incorrect raster clipping. Choose another CRS or clip the unprojected raster before projecting")
} else {
  flipperTA <- "N"
}

#Current Climate Data--------------------------------------
#Clip Training Area (if necessary)
if (TrainingAreaClip == "Y") {
  bioclim_t <- c()
  bioclim_a <- c()
  for (i in 1:nlayers(bioclim)) {
    if (flipperTA == "N") {  
      #Clips training layers
      bioclim_t <- c(bioclim_t, crop(bioclim[[i]], extent(max(NewExtent_TA$long[1], extent(bioclim)[1]), 
                                                            min(NewExtent_TA$long[2], extent(bioclim)[2]),
                                                            max(NewExtent_TA$lat[1], extent(bioclim)[3]),
                                                            min(NewExtent_TA$lat[2], extent(bioclim)[4]))))
      #Clips study layers
      bioclim_a <- c(bioclim_a, crop(bioclim[[i]], extent(max(NewExtent$long[1], extent(bioclim)[1]),
                                                            min(NewExtent$long[2], extent(bioclim)[2]),
                                                            max(NewExtent$lat[1], extent(bioclim)[3]),
                                                            min(NewExtent$lat[2], extent(bioclim)[4]))))
    } else {
      #Clips training layers
      bioclim_t <- c(bioclim_t, crop(bioclim[[i]], extent(min(NewExtent_TA$long[1], extent(bioclim)[1]), 
                                                          max(NewExtent_TA$long[2], extent(bioclim)[2]),
                                                          min(NewExtent_TA$lat[1], extent(bioclim)[3]),
                                                          max(NewExtent_TA$lat[2], extent(bioclim)[4]))))
      #Clips study layers
      bioclim_a <- c(bioclim_a, crop(bioclim[[i]], extent(min(NewExtent$long[1], extent(bioclim)[1]),
                                                          max(NewExtent$long[2], extent(bioclim)[2]),
                                                          min(NewExtent$lat[1], extent(bioclim)[3]),
                                                          max(NewExtent$lat[2], extent(bioclim)[4]))))
    }
  }
  
  #Writes and stores projected files in a temporary directory
  if (!dir.exists(paste0(DataDirectory, "/TEMP"))) {
    dir.create(paste0(DataDirectory, "/TEMP"))  
  }
  setwd(paste0(DataDirectory, "/TEMP"))
  for (k in 1:length(bioclim_t)) {
    writeRaster(bioclim_t[[k]], 
                filename = paste0(names(bioclim_t[[k]]), ".bil", sep = ""),
                bylayer = TRUE,
                overwrite = TRUE,
                format = "EHdr",
                prj = TRUE)
  }
  print("Clipping Rasters (Training Area)")
  rm(bioclim_t, bioclim)
  setwd(trainingarea)
} else {
  #Writes and stores training area layers in a temporary directory
  bioclim_t2 <- stack(bioclim)
  dir.create(paste0(DataDirectory, "/TEMP"))
  setwd(paste0(DataDirectory, "/TEMP")) 
  for (k in 1:nlayers(bioclim_t2)) {
    writeRaster(bioclim_t2[[k]], 
                filename = paste0(names(bioclim_t2[[k]]), ".bil"), 
                overwrite = TRUE, 
                bylayer = TRUE, 
                format = "EHdr",
                prj = TRUE)
  }
  print("Clipping Rasters (Training Area)")
  #Clips study area to desired extent
  bioclim_a <- c()
  for (k in 1:nlayers(bioclim_t2)) {
    bioclim_a <- c(bioclim_a, crop(bioclim_t2[[k]], extent(NewExtent$long[1], NewExtent$long[2], NewExtent$lat[1], NewExtent$lat[2])))
  }
  rm(bioclim_t2, bioclim)
}

#Resampling Study Area Environmental Rasters
#Uses the first raster as the resampling template
pat <- bioclim_a[[1]]
bioclim_na <- c(pat)

#Resamples clipped rasters to template
for (i in 2:length(bioclim_a)) {
  bioclim_na <- c(bioclim_na, resample(bioclim_a[[i]], pat, method = "ngb"))
}

bioclim_na2 <- stack(bioclim_na)
rm(bioclim_a)

#Deletes all previous files within the projected study area directory
setwd(studyarea)
FilesToClear <- list.files(studyarea) 
unlink(FilesToClear)

#Names and writes clipped rasters
for (k in 1:nlayers(bioclim_na2)) {
  writeRaster(bioclim_na2[[k]], filename = paste0(names(bioclim_na2)[k], ".bil"), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
}
print("Clipping Rasters (Study Area)")

#Clipping Future Environmental Layers----------------------------------
#Clips all of the future environmental rasters
if (numScenario > 0) {
  setwd(proj_predictenv_dir)
  
  #Separates out directories
  predictenvdir <- list.dirs(path = proj_predictenv_dir, full.names = TRUE)
  for (i in 1:length(predictenvdir)) {
    correctDir <- list.dirs(path = predictenvdir[i], full.names = TRUE)
    
    #If the directory is an endpoint, clips
    if (length(correctDir) == 1) {
      #Makes a list of all future environmental layer files
      setwd(correctDir)
      future <- list.files(path = getwd(), pattern=paste0("\\.bil$"), full.names = TRUE)
      if (length(future) == 0) {
        stop(paste0("Forecasted/hindcasted climate rasters not found or not in ", rastertype, " format!", "\n", 
                    "  Ensure that all predicted climate rasters are in the correct locations"))
      } else if (length(future) != length(bioclim_na)) {
        stop(paste0("Number of forecasted/hindcasted climate rasters does not equal number of current cliamte variables", "\n",
                    "  Ensure that every predicted climate raster is in the correct location"))
      }
      future_a <- c()
      
      #Ensures that all environmental rasters are in the same projection
      ProjEnv <- rep(NA, len = length(future))
      for (i in 1:length(future)) {
        ProjEnv[i] <- as.character(crs(raster(future[[i]])))
      }
      
      ProjUnique <- unique(ProjEnv)
      if (length(ProjUnique) > 1) {
        stop("Not all of the forecasted/hindcasted environmental rasters are in the same projection")
      }
      
      #Ensures that all rasters have a defined coordinate projection
      futurestack <- stack(future)
      if(is.na(crs(futurestack))) {
        stop("forecasted/hindcasted raster crs = NA: Ensure all rasters have a defined coordinate projection")
      }
      
      #Clipping each individual raster file
      for (j in 1:length(future)) {
        future_a <- c(future_a, crop(raster(future[j]), extent(pat)))
      }  
      
      #Resample to the environmental layer raster
      future_na <- c()
      for (j in 1:length(future)) {
        future_na <- c(future_na, resample(future_a[[j]], pat, method = "ngb"))
      }
      
      #Write the rasters
      future_na2 <- stack(future_na)
      setwd(correctDir)
      for (k in 1:nlayers(future_na2)){
        writeRaster(future_na2[[k]], filename = paste0((names(future_na2)[k]), ".bil"), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
        
      }
    }
  }
  print("Clipping Rasters (Future Env)")
}

#Urban and Protected Analysis-------------------------------
#Clip the urban and protected data
if (UrbanAnalysis == "Y") {
  
  #Making a list of the urban raster files
  setwd(proj_urbanized_dir)
  urban <- paste0("/", list.files(path = ".", pattern = paste0("\\.bil$"), full.names = FALSE))
  urban_a <- c()
  
  #Clipping each individual raster file
  for (i in 1:length(urban)) {
    urban_a <- c(urban_a, crop(raster(urban[i]), extent(pat)))
  }
  
  #resample to the environmental reference raster
  urban_na <- c()
  for (i in 1:length(urban)) {
    urban_na <- c(urban_na, resample(urban_a[[i]], pat, method="ngb"))
  }
  
  rm(urban_a)
  
  #overwriting the urbanized raster files
  urban_na <- stack(urban_na)
  setwd(proj_urbanized_dir)
  unlist(proj_urbanized_dir)
  for (j in 1:nlayers(urban_na)) {
    #If thehre are multiple urbanized raster (for different time periods), write out each one
    if (length(grep(paste(years, collapse = "|"), names(urban_na))) == nlayers(urban_na)) {
      writeRaster(urban_na[[j]], filename = paste0("SA_urb", years[j], rastertype, sep = ""), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
    } else {
      writeRaster(urban_na[[j]], filename = paste0("SA_urb", j, rastertype, sep = ""), overwrite = TRUE, bylayer = TRUE, format = "EHdr", prj = TRUE)
    }
  }
  print("Writing Urban Rasters")
}

if (ProtectedAnalysis == "Y") {
  #Lists the protected shapefiles (there usually should only be one)
  setwd(proj_protected_dir)
  protect <- paste0("/", list.files(path = ".", pattern = paste0("\\", ".shp", "$"), full.names = FALSE))
  protect_a <- c()
  
  #Clips the protected shapefiles
  for (i in 1:length(protect)) {
    protect_a <- c(protect_a, crop(shapefile(protect[i]), extent(pat)))
  }
  
  #Writes the shapefiles
  for (i in 1:length(protect_a)) {
    writeOGR(protect_a[[i]], dsn = proj_protected_dir, layer = paste0("protected_areas_sa_",i), overwrite_layer = TRUE, driver = "ESRI Shapefile")
  }
}

#Copies training area files from temporary directory to desired directory
rm(bioclim_na, bioclim_na2, pat)
unlink(trainingarea)
TrainingFiles <- list.files(paste0(DataDirectory, "/TEMP"))
for (f in 1:length(TrainingFiles)) {
  file.copy(from = paste0(DataDirectory, "/TEMP/", TrainingFiles[f]), to = paste0(trainingarea), overwrite = TRUE)
}

rm(TrainingFiles)
unlink(paste0(DataDirectory, "/TEMP"), recursive = TRUE)