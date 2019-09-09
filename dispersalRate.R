#Initializations-----------------------
library(biomod2)
library(raster)
library(parallel)
library(gtools)

result_dir <- df[, "result_dir"]
dispersalRate_dir <- df[, "dispersalRate_dir"]

if (UrbanAnalysis =="Y" ) {
  proj_urbanized_dir <- df[, "proj_urbanized_dir"]
} else {
  proj_urbanized_dir <- c()
}

if (ProtectedAnalysis == "Y") {
  proj_protected_dir <- df[, "proj_protected_dir"]
} else {
  proj_protected_dir <- c()
}

rastertype <- df[, "rastertype"]
ncores <- as.numeric(df[, "ncores"])
outfile <- "statsout.txt"

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

#Gets vectors of the scenarios and years
numScenario <- df[, "numScenario"]
Scenarios <- df[, grep("^Scenario", colnames(df))]

numYear <- df[, "numYear"]
years <-  rep(NA, length = numYear)
years <-  as.numeric(config[grep("^Year", row.names(config)), ])

#Reads in the dispersal data
setwd(dispersalRate_dir)
dispersal <- read.csv(list.files(getwd(), pattern = ".csv"), stringsAsFactors = FALSE)

#Gets a list of the species analyzed on this go-around
ListSpp <- c()
speciesWorked <- spp_batch

for (i in 1:length(speciesWorked)) {
  ListSpp[i] <- paste0(substr(speciesWorked[i], 1, nchar(speciesWorked[i]) - 4))
}
ListSpp <- unique(ListSpp)
print("   Will evaluate species:")
print(ListSpp)

#Functions----------------------------
#Creates distance rasters from original projection (studyarea environmental rasters)
DistanceRaster <- function(spp, Time, Scen, CurrentBinary) {
  CurrentBinary <- CurrentBinary
  CurrentPresence <- reclassify(CurrentBinary, c(0, 0, NA),include.lowest = TRUE)
  CurrentPresence <- trim(CurrentPresence)
  CurrDist <- distance(CurrentPresence, doEdge = TRUE)
  CurrDist <- extend(CurrDist, extent(CurrentBinary))
  CurrDist[is.na(values(CurrDist))] <- max(values(CurrDist), na.rm = TRUE) + 1
  DistFinal <- mask(CurrDist, CurrentBinary)
  DistFinal <- DistFinal / 1000
  writeRaster(DistFinal, 
              filename = paste0(result_dir, "/", spp, "/", "distance_", Time, "_", Scen, rastertype), 
              overwrite = TRUE, 
              format = format, 
              prj = TRUE)
  rm(CurrentPresence, CurrDist)
  gc()
  return(DistFinal)
}

#Creates dispersal probability raster
DispersalProbRaster <- function(rate, DistRaster, Elapsed) {
  #Calculates lambda of exponential distribution
  Lambda <- 1 / rate
  
  #When exponential distributions are added, they convolute to a gamma distribution
  GammaProbFun <- function(x) {
    1 - pgamma(x, shape = Elapsed, rate = Lambda)
  }
  
  DistProb <- calc(DistRaster, fun = function(x){GammaProbFun(x)})
  print("Creating Dispersal Probability Raster")
  return(DistProb)
}

t1nott2 <- function(t1, t2) {
  return(mask(t1, t2, inverse = FALSE, maskvalue = 1, updatevalue = 0))
}

overlap <- function(t1, t2) {
  return(mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
}

getCentroid <- function(raster) {
  # A matrix with three columns: x, y, and v (value)
  points <- rasterToPoints(raster, fun = function(x){x == 1}, spatial = FALSE)
  
  #average latitude (y)
  Clat <- mean(points[, 2], na.rm = TRUE)
  
  #average longitude (x)
  Clong <- mean(points[, 1], na.rm = TRUE)
  
  #returns the longitude & latitude of the Centroid
  return(c(Clong, Clat))
}

getProtected <- function(raster, protected) {
  m <- mask(raster, protected)
  return(freq(m, digits = 0, value = 1, useNA = 'no', progress = ''))
}

#The amount of pixels covered by an urban area
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
  m <- mask(focusraster, cur2, inverse = TRUE, maskvalue = 1, updateValue = 0)
  return (freq(m, digits = 0, value = 1, useNA = 'no', progress = ''))
}

FinalDispersal <- function(spp) {
  #Gets species name and relevant dispersal rate
  speciesName = gsub("_", " ", spp)
  if (length(grep(paste0(speciesName), dispersal[, 1])) > 0) {
    dispRateColumn <- which(unlist(lapply(dispersal, is.numeric)))
    dispersalRate <- dispersal[grep(paste0(speciesName, "\\s*$"), dispersal[, 1]), dispRateColumn]
    if (!is.na(dispersalRate)) {
      CurrentTime <- years[1]
      CurrentBinary <- raster(paste0(result_dir, "/", spp, "/", CurrentTime, "_binary", rastertype))
      #For stats
      Projection <- c("Current")
      NumberCells <- c(cellStats(CurrentBinary, stat = sum))
      CellChange <- c(0)
      T1notT2 <- c(0)
      T2notT1 <- c(0)
      Overlap <- c(0)
      CentroidX <- c(getCentroid(CurrentBinary)[1])
      CentroidY <- c(getCentroid(CurrentBinary)[2])
      Protected <- c()
      Urbanized <- c()
      
      if (ProtectedAnalysis == "Y") {
        setwd(proj_protected_dir)
        prot_files <- list.files(path = ".", pattern = paste0("\\", ".shp", "$"), full.names = TRUE)
        protected <- shapefile(prot_files[1])
        Protected <- c(Protected, getProtected(CurrentBinary, protected))
      } 
      
      if (UrbanAnalysis == "Y") {
        setwd(proj_urbanized_dir)
        urblist <- list.files(pattern = paste0("\\.bil$"), full.names = TRUE)
        urbanized <- stack(urblist)
        Urbanized <- c(Urbanized, getUrbanized(CurrentBinary, CurrentTime, urblist, urbanized))
      }
      
      #Creates dispersal rasters and PDFs
      for (s in 1:length(Scenarios)) {
        CurScen <- Scenarios[s]
        curdir <- paste0(result_dir, "/", spp, "/", CurScen)
        setwd(curdir)
        DispersalNames <- c()
        TimeMap <- raster(paste0(result_dir, "/", spp, "/TimeMapRasters/binary", CurScen, rastertype))
        for (y in 2:length(years)) {
          CurYear <- years[y]
          
          #Calculates distance from current distribution
          if (y == 2) {
            DistanceRastersExist <- list.files(paste0(result_dir, "/", spp), 
                                               pattern = paste0("distance_", years[1], "_Current", rastertype, "$"))
            if (length(DistanceRastersExist) == 0) {
              SppDistance <- DistanceRaster(spp, CurrentTime, "Current", CurrentBinary)
              OriginalDistance <- SppDistance
            } else {
              OriginalDistance <- raster(paste0(result_dir, "/", spp, "/", DistanceRastersExist))
              SppDistance <- OriginalDistance
            }
          } else {
              FocusTime <- years[y - 1]
              #Reading in the current distribution
              CurrentDistribution <- Binary_Dispersal
              
              #Sizing down the raster for faster distance measuring
              TimeMap2 <- TimeMap
              TimeMap2[which(values(TimeMap2) == 0)] <- NA
              TimeMap2 <- trim(TimeMap2)
              TimeMap2[is.na(values(TimeMap2))] <- 0
              
              #Making the places where species lived in the previous time step
              CurrentDistribution <- crop(CurrentDistribution, extent(TimeMap2))
              CurrentDistribution[which(values(CurrentDistribution) == 0)] <- NA
              
              #Creating new distance raster
              CurrentDistance <- distance(CurrentDistribution) / 1000
              SppDistance <- extend(CurrentDistance, extent(Binary_Dispersal), value = max(values(SppDistance), na.rm = TRUE) + 1)
              SppDistance <- mask(SppDistance, CurrentBinary)
          }
          
          TimeDiff <- abs(CurYear - years[y - 1])
          SppDispProb <- DispersalProbRaster(dispersalRate, SppDistance, TimeDiff)
          setwd(curdir)
          RasterList <- list.files(curdir, pattern = paste0(rastertype, "$"))
          #Creates an ensembled raster that incorporates dispersal rate
          EnsembleNum <- grep(paste0(CurYear, "_", CurScen, "_ensembled", rastertype), RasterList)
          EnsembleSD <- raster(RasterList[EnsembleNum])
          if (extent(SppDispProb) != extent(EnsembleSD) | ncol(SppDispProb) != ncol(EnsembleSD)) {
            message("Raster extents are not consistent: only the intersection of the rasters will be analysed")
            SppDispProb <- intersect(SppDispProb, EnsembleSD)
            SppDispProb <- resample(SppDispProb, EnsembleSD, method = "bilinear")
          }
          Ensemble_Dispersal <- SppDispProb * EnsembleSD
          
          writeRaster(Ensemble_Dispersal,
                      filename = paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate", rastertype),
                      overwrite = TRUE, 
                      format = format,
                      prj = TRUE)
          
          DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_ensembled_dispersalRate"))
          
          #Creates a binary raster that incorporates dispersal rate
          BinaryNum <- grep(paste0(CurYear, "_", CurScen, "_binary", rastertype), RasterList)
          BinarySD <- raster(RasterList[BinaryNum])
          Binary_Dispersal <- BinaryTransformation(SppDispProb * BinarySD, 0.5)
          
          writeRaster(Binary_Dispersal,
                      filename = paste0(CurYear, "_", CurScen, "_binary_dispersalRate", rastertype),
                      overwrite = TRUE, 
                      format = format,
                      prj = TRUE)
          DispersalNames <- c(DispersalNames, paste0(CurYear, "_", CurScen, "_binary_dispersalRate"))
          
          #Fills out the stats table
          
          Projection <- c(Projection, paste0(CurScen, "_", CurYear))
          NumberCells <- c(NumberCells, cellStats(Binary_Dispersal, stat = sum))
          CellChange <- c(CellChange, NumberCells[length(NumberCells)] - NumberCells[1])
          T1notT2 <- c(T1notT2, cellStats(t1nott2(CurrentBinary, Binary_Dispersal), stat = sum))
          T2notT1 <- c(T2notT1, cellStats(t1nott2(Binary_Dispersal, CurrentBinary), stat = sum))
          Overlap <- c(Overlap, cellStats(overlap(Binary_Dispersal, CurrentBinary), stat = sum))
          CentroidX <- c(CentroidX, getCentroid(Binary_Dispersal)[1])
          CentroidY <- c(CentroidY, getCentroid(Binary_Dispersal)[2])
          
          if (UrbanAnalysis  == "Y") {
            setwd(proj_urbanized_dir)
            urblist <- list.files(pattern = paste0("\\", rastertype, "$"), full.names = TRUE)
            urbanized <- stack(urblist)
            Urbanized <- c(Urbanized, getUrbanized(Binary_Dispersal, CurYear, urblist, urbanized))
          }
          
          if (ProtectedAnalysis == "Y") {
            setwd(proj_protected_dir)
            ProtectList <- list.files(pattern = ("\\.shp$"), full.names = TRUE)
            if (!exists("protected")) {
              protected <- shapefile(ProtectList[1])
            }
            Protected <- c(Protected, getProtected(Binary_Dispersal, protected))
          }
        }
        #Writes PDFs
        DispersalRasters <- list.files(curdir, pattern = paste0("dispersalRate", rastertype, "$"), full.names = TRUE)
        DispersalRasters <- mixedsort(DispersalRasters)
        DispersalNames <- mixedsort(DispersalNames)
        setwd(paste0(result_dir, "/", spp))
        dir.create("map_pdfs")
        setwd("map_pdfs")
        for (d in 1:length(DispersalRasters)) {
          title <- DispersalNames[d]
          pdf(file = paste0(DispersalNames[d], ".pdf"))
          plot(raster(DispersalRasters[d], native = TRUE), main = title)
          dev.off()
        }
      }
      
      if ((UrbanAnalysis == "Y") | (ProtectedAnalysis == "Y")) {
        if (UrbanAnalysis == "N") {
          stats <- cbind(Projection, NumberCells, CellChange,                # Which columns do you want?
                      T1notT2, T2notT1, Overlap, CentroidX, CentroidY, Protected)
        } else if (ProtectedAnalysis == "N") {
          stats <- cbind(Projection, NumberCells, CellChange, 
                      T1notT2, T2notT1, Overlap, CentroidX, CentroidY, Urbanized)       # Which columns do you want?  
        } else {
          stats <- cbind(Projection, NumberCells, CellChange, 
                      T1notT2, T2notT1, Overlap, CentroidX, 
                      CentroidY, Urbanized, Protected)
        }
      } else {
        stats <- cbind(Projection, NumberCells, CellChange, T1notT2,
                    T2notT1, Overlap, CentroidX, CentroidY)
      }
      stats <- as.data.frame(stats)
      write.csv(stats, file = (paste0(result_dir, "/", spp, "/Results_Dispersal.csv")))
      rm(CurrentBinary)
      gc()
    } else {
      message(paste0("No dispersal rate values found for ", speciesName, ": skipping dispersal rate analysis"))
    }
  } else {
    message(paste0("No dispersal rate data found for ", speciesName, ": skipping dispersal rate analysis"))
  }
}

#Run-----------------------------------
for (w in 1:length(ListSpp)) {
  FocSpec <- gsub("_", " ", ListSpp[w])
  DispSpec <- grep(paste0("^", FocSpec, "$"), dispersal[, 1])
  if (length(DispSpec) == 0) {
    message(paste0("No dispersal rate values found for ", FocSpec, ": skipping dispersal rate analysis"))
  }
}

clus <- makeCluster(ncores, outfile = outfile)
clusterExport(clus, varlist = c("result_dir", "dispersalRate_dir", "format", "rastertype", "proj_urbanized_dir",
                              "proj_protected_dir", "numScenario", "Scenarios", "numYear", "years", "dispersal",
                              "ncores", "DistanceRaster", "DispersalProbRaster", "FinalDispersal", "ListSpp",
                              "getCentroid", "getProtected", "getUrbanized", "t1nott2", "overlap", "UrbanAnalysis",
                              "ProtectedAnalysis"))
clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(biomod2))
clusterEvalQ(clus, library(gtools))

out<-parLapply(clus, ListSpp, function(x) FinalDispersal(x))
stopCluster(clus)
df[, "dispersalRan"] <- "Y"