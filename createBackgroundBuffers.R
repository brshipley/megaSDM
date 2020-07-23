####createBackgroundBuffers.R####
##Creates buffers to spatially weight the background points near the occurrences
#Initializations---------------------------
#Loads the necessary packages
library(parallel)
library(raster)
library(rgdal)
library(rgeos)
library(sampSurf)

#Loads the necessary variables from "df"
proj_trainingarea <- df[, "proj_trainingarea"]
buff_dir <- df[, "buff_dir"]
if (!dir.exists(buff_dir)) {
  dir.create(buff_dir)
}

result_dir <- df[, "result_dir"]
occurrences <- df[, "occurrences"]
test <- df[, "test"]
desiredCRS <- df[, "desiredCRS"]
defaultCRS <- df[, "defaultCRS"]
outfile <- "statsout.txt"
nsubsamp <- df[, "nsubsamp"]

setwd(test)
speciesWorked <- spp_batch

for (i in 1:length(speciesWorked)) {
  spp.list <- c(spp.list, list.files(path = samples, full.names = TRUE, pattern = speciesWorked[i]))
}

spp.list <- unique(spp.list)
print("   Building Buffers For:")
print(paste0("        ", paste0(spp_batch)))

#Functions------------------------
BuffFiles <- function(CurSpp) {
  #Species name management
  setwd(paste0(test, "/species"))
  species <- substr(CurSpp, 9, nchar(CurSpp))
  species <- read.csv(species)
  if (length(grep("_", CurSpp)) > 0) {
    name <- as.character(species$Species[1])
    name <- (gsub(" ", "_", name))
  } else {
    name <- c(unlist(strsplit(as.character(species$Species[1]), " ")))[1]
  }
  
  #Converts species occurrence points into SpatialPoints with desiredCRS
  coordinates <- species[, c("Longitude", "Latitude")]
  Coordinates2 <- SpatialPoints(coordinates, proj4string = crs(defaultCRS))
  Coordinates2 <- spTransform(Coordinates2, as.character(desiredCRS)) 
  
  #Calculates a buffer width based on quantiled minimum distance among all occurrence points
  Locations <- as.data.frame(Coordinates2@coords)    
  
  #Uses the first 1000 points (randomly sampled) to create buffers and distances
  if (nrow(Locations) > 1000) {
    BufferPointNumber <- 1000
  } else {
    BufferPointNumber <- nrow(Locations) - 1
  }
  
  SampleLocations <- Locations[sample(c(1:nrow(Locations)), size = BufferPointNumber, replace = FALSE), ]
  
  #Uses the 95% quantile of the minimum distance between each point
  Distance <- pointDistance(SampleLocations, lonlat = FALSE)
  mindist <- c()
  for (q in 1:ncol(Distance)) {
    DistanceZero <- Distance[which(Distance[, q] > 0), q]
    mindist <- c(mindist, min(DistanceZero))
  }
  BufferWidth <- 2 * quantile(mindist, 0.95) 
  rm(Distance, mindist)
  gc()
  
  #Makes a buffer (width dependent) around the first occurrence point
  combinedPolygon <- spCircle(BufferWidth, spUnits = crs(desiredCRS),
                              centerPoint = c(x = Locations[1, 1], y = Locations[1, 2]))$spCircle
  
  #Creates buffers aroudn the remaining points and combines them
  for (b in 2:nrow(Locations)) {
    circle <- spCircle(BufferWidth, spUnits = crs(desiredCRS), centerPoint = c(x = Locations[b, 1], y = Locations[b, 2]))
    binded <- bind(combinedPolygon, circle$spCircle)
    combinedPolygon <- gUnaryUnion(binded)
  }
  
  #Reprojects the buffer into the desired CRS
  proj4string(combinedPolygon) <- crs(desiredCRS)
  
  #Write shapefile out of the buffer
  setwd(buff_dir)
  sp_poly_df <- SpatialPolygonsDataFrame(combinedPolygon, data = data.frame(ID = 1), match.ID = FALSE)
  writeOGR(sp_poly_df, dsn = paste0(buff_dir, "/", name, ".shp"), layer = paste0(name), overwrite_layer = TRUE, driver = "ESRI Shapefile")
  rm(combinedPolygon, sp_poly_df)
  gc()
  
  #Fill out background point stats table with information on the Buffer Width (m)
  BGPStats <- rep(BufferWidth, nsubsamp)
  dir.create(paste0(result_dir, "/", name))
  write.csv(BGPStats,file = paste0(result_dir, "/", name, "/BackgroundPoints_stats.csv"))
}

#Run------------------------------
#Parallelization
clus <- makeCluster(ncores, outfile = outfile, setup_timeout = 0.5)
clusterExport(clus, varlist = c("buff_dir", "proj_trainingarea", "occurrences", "test", 
                              "desiredCRS", "defaultCRS", "spp.list", "BuffFiles",
                              "format", "nsubsamp","result_dir"))
clusterEvalQ(clus, library(maptools))
clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(rgdal))
clusterEvalQ(clus, library(rgeos))
clusterEvalQ(clus, library(sampSurf))

out<-parLapply(clus, spp.list, function(x) BuffFiles(x))
stopCluster(clus)