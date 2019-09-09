#Initializations-----------------------------
library(plotfunctions)
library(raster)
library(gtools)

result_dir <- df[, "result_dir"]
occ <- df[, "occurrences"]
spplist <- df[, "spplist"]
occurrences <- df[, "occurrences"]
rastertype <- df[, "rastertype"]
desiredCRS <- df[, "desiredCRS"]
numScenario <- df[, "numScenario"]
numYear <- df[, "numYear"]
if (numScenario > 0) {
  predictenv <- df[, "proj_predictenv"]
}

#Raster-Format Dictionary
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

#Reading in the Taxon-Species list
specieslist <- read.csv(spplist)
speciesfolders <- list.dirs(result_dir, recursive = FALSE)
speciesfolders <- speciesfolders[grep("_", speciesfolders)]

#Making a list of all species that have folders
sppfold <- c()
for (i in 1:length(speciesfolders)){
  if (length(list.files(speciesfolders[i]) > 0)) {
    sppfold <- c(sppfold, speciesfolders[i])
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
  if (!length(curfocus) > 0) {
    message(paste0(focusspp, " will be removed (no replicates with an AUC > ", aucval, ")"))
    DeleteSP <- c(DeleteSP, sp)
  }
}

if (length(DeleteSP) > 0) {
  InitialTaxonList <- taxonlist[-DeleteSP, ]
} else {
  InitialTaxonList <- taxonlist
}

#Gets the different broad taxa to do species richness on
taxa <- c(levels(InitialTaxonList[, 2]), "All_Taxa")
ntaxa <- length(taxa)
if (ntaxa == 2) {
  taxa <- taxa[1]
}
ntaxa <- length(taxa)

dir.create(paste0(result_dir, "/", "RichnessMaps"))
RichnessMaps <- paste0(result_dir, "/", "RichnessMaps")

#Gets a vector of Scenarios
if (numScenario > 0) {
  predictenvdir <- list.dirs(path = predictenv, recursive = TRUE)
  Scenarios <- c(as.character(df[, grep("^Scenario", colnames(df))]))
} else {
  Scenarios <- c("Current")
}
#Gets a vector of Years
years <-  as.numeric(df[, grep("^Year", colnames(df))])

#No Dispersal------------------------------------------------------------
FutureRichnessMaps <- function(taxon, y, s) {
  futstack <- c()
  #Uploads and stacks binary SDM rasters of particular time, taxon, and scenario
  for (sp in 1:length(taxonlist2)) {
    FocusYear <- years[y]
    FocusScenario <- Scenarios[s]
    focusspp <- taxonlist2[sp]
    path <- paste0(result_dir, "/", focusspp)
    r <- list.files(path = path, pattern = paste0(FocusYear, "_", FocusScenario, "_binary", rastertype, "$"), recursive = TRUE)
    futraster <- raster(paste0(path, "/", r), crs = desiredCRS)
    futstack <- stack(c(futraster, futstack))
  }
  
  #Calculates sum (total species richness) for given parameters, makes rasters and PDFs
  FutureRichness <- calc(futstack, fun = sum)
  print(paste0("Creating species richness maps for ", taxon, " for Scenario ", FocusScenario, " and time ", FocusYear))
  setwd(RichnessMaps)
  writeRaster(FutureRichness, 
              filename = paste0("Richness_", taxon, "_", FocusYear, "_", FocusScenario, rastertype),
              format = format,
              overwrite = TRUE,
              prj = TRUE)
}

#Makes richness maps by taxon
for (t in 1:ntaxa) {
  #Subsets out taxa of interest
  taxon <- taxa[t]
  if (taxon == "All_Taxa") {
    taxonlist2 <- gsub(" ", "_", InitialTaxonList[, 1])
  } else {
    taxonlist2 <- gsub(" ", "_", InitialTaxonList[which(InitialTaxonList[, 2] == taxon), 1])
  }
  print(paste0("Creating current species richness maps for ", taxon, ", with ", length(taxonlist2), " total species"))
  
  #Uploads and stacks binary rasters of all species within the taxa
  curstack <- c()
  DeleteSP <- c()
  for (sp in 1:length(taxonlist2)) {
    focusspp <- taxonlist2[sp]
    if (dir.exists(paste0(result_dir, "/", focusspp))) {  
      setwd(paste0(result_dir, "/", focusspp))
      curfocus <- list.files(pattern = paste0("binary", rastertype, "$"))
      if (length(curfocus) > 0){
        curraster <- raster(curfocus, crs = desiredCRS)
        curstack <- stack(c(curstack, curraster))
      } else {
        message(paste0(focusspp, " will be removed (no replicates with an AUC > ", aucval, ")"))
        DeleteSP <- c(DeleteSP, sp)
      }
    }
  }
  
  if (length(DeleteSP) > 0) {
    taxonlist2 <- taxonlist2[-DeleteSP]
  } else {
    taxonlist2 <- taxonlist2
  }
  
  #Calculates the sum of the rasters (and therefore species richness), and creates rasters and PDFs
  CurrentRichness <- calc(curstack, fun = sum)
  setwd(RichnessMaps)
  writeRaster(CurrentRichness, 
              filename = paste0("Richness_", taxon, "_", years[1], rastertype), 
              format = format, 
              overwrite = TRUE,
              prj = TRUE)
  
  #Creates species richness maps for other time frames and scenarios
  if (numYear > 1) {
    for (y in 2:numYear) {
      for (s in 1:numScenario) {
        FutureRichnessMaps(taxon, y, s)
      }
    }
  }
}

#Dispersal---------------------------------------------------------------
FutureRichnessMapsDisp <- function(taxon, y, s) {
  futstack <- c()
  #Uploads and stacks binary SDM rasters (dispersal rate included) of particular time, taxon, and scenario
  for (sp in 1:length(taxonlist2)) {
    FocusYear <- years[y]
    FocusScenario <- Scenarios[s]
    focusspp <- taxonlist2[sp]
    path <- paste0(result_dir, "/", focusspp)
    r <- list.files(path = path, 
                    pattern = paste0(FocusYear, "_", FocusScenario, "_binary_dispersalRate", rastertype, "$"),
                    recursive = TRUE)
    if (length(r) == 1) {
      futraster <- raster(paste0(path, "/", r), crs = desiredCRS)
      futstack <- stack(c(futraster, futstack))
    }
  }
  
  #Calculates sum (total species richness) for given parameters, makes rasters and PDFs
  FutureRichness <- calc(futstack, fun = sum)
  print(paste0("Creating species richness maps for ", taxon, " for Scenario ", FocusScenario, " and time ", FocusYear, " (including dispersal rate)"))
  setwd(RichnessMaps)
  writeRaster(FutureRichness, 
              filename = paste0("Richness_", taxon, "_", FocusYear, "_", FocusScenario, "_", "dispersalRate", rastertype),
              format = format,
              overwrite = TRUE,
              prj = TRUE)
  
  NonDispersalSR <- raster(paste0("Richness_", taxon, "_", FocusYear, "_", FocusScenario, rastertype))
  
  #Calculates difference between dispersal and non-dispersal species richness maps
  DifferenceSR <- FutureRichness - NonDispersalSR
  writeRaster(DifferenceSR, 
              filename = paste0("Richness Difference_", taxon, "_", FocusYear, "_", FocusScenario, rastertype),
              format = format,
              overwrite = TRUE,
              prj = TRUE)
  DiffVec <- unique(DifferenceSR)
  if (min(DiffVec) < 0 & max(DiffVec) > 0) {
    color <- colorRampPalette(c("red", "grey91", "blue"))(length(DiffVec))
  } else if (min(DiffVec) < 0 & max(DiffVec) <= 0) {
    color <- colorRampPalette(c("red", "gold1", "grey91"))(length(DiffVec))
  } else {
    color <- colorRampPalette(c("grey91", "deepskyblue", "blue"))(length(DiffVec))
  }
  
  pdf(file = paste0(result_dir, "/", "RichnessMaps/", taxon, "_", FocusYear, "_", FocusScenario, "_", "dispersalDifference.pdf"))
  plot(DifferenceSR, col = color, xlab = "", ylab = "", main = paste0(taxon, " ", FocusScenario, " ", FocusYear, " Dispersal Difference"))
  dev.off()
  rm(futstack)
}

if (dispersalStep == "Y" && length(taxonlist2) > 0) {
  for (t in 1:ntaxa) {
    #Subsets out taxa of interest
    taxon <- taxa[t]
    
    if (taxon == "All_Taxa") {
      taxonlist2 <- gsub(" ", "_", InitialTaxonList[, 1])
    } else {
      taxonlist2 <- gsub(" ", "_", InitialTaxonList[which(InitialTaxonList[, 2] == taxon), 1])
    }
    print(paste0("Creating current species richness maps for ", taxon, ", with ", length(taxonlist2), " total species"))
    
    #Uploads and stacks binary rasters of all species within the taxa
    curstack <- c()
    DeleteSP <- c()
    for (sp in 1:length(taxonlist2)) {
      focusspp <- taxonlist2[sp]
      setwd(paste0(result_dir, "/", focusspp, "/", Scenarios[1]))
      curfocus <- list.files(pattern = paste0("binary_dispersalRate", rastertype, "$"))
      if (!length(curfocus) > 0) {
        message(paste0(focusspp, " will be removed (no dispersal rate rasters were found)"))
        DeleteSP <- c(DeleteSP, sp)
      }
    }
    
    if (length(DeleteSP) > 0) {
      taxonlist2 <- taxonlist2[-DeleteSP]
    } else {
      taxonlist2 <- taxonlist2
    }
    
    if (length(taxonlist2) > 0) {
      #Creates species richness maps for time frames and scenarios
      for (y in 2:numYear) {
        for (s in 1:numScenario) {
          FutureRichnessMapsDisp(taxon, y, s)
        }
      }
    }
  }
}

#Making PDF Maps---------------------------------------------------------
MakePDFMaps <- function(taxon) {
  setwd(RichnessMaps)
  RasterList <- list.files(getwd(), pattern = paste0("\\", rastertype, "$"))
  FocusRasters <- RasterList[grep(paste0("Richness_", taxon), RasterList)] 
  FocusStack <- stack(FocusRasters)
  MaximumRich <- max(cellStats(FocusStack, stat = max, na.rm = TRUE))
  
  for (e in 1:nlayers(FocusStack)) {
    Raster1 <- FocusStack[[e]]
    MaxRaster1 <- cellStats(Raster1, stat = max, na.rm = TRUE)
    ValueFreq <- freq(Raster1)
    #If the maximum species richness occurs in too few pixels to display in pdf, only plots the next highest species richness
    while ((ValueFreq[which(ValueFreq[, 1] == MaxRaster1), 2] < floor(0.00005 * ncell(Raster1))) && (MaxRaster1 > 1)) {
      Raster1[Raster1 == MaxRaster1] <- MaxRaster1 - 1
      message(paste0("The maximum species richness in the study area, ", MaxRaster1, ", covers too few cells to display on PDF:"))
      message("Examine the created raster for an accurate count of pixels.")
      MaxRaster1 <- MaxRaster1 - 1
    }
    
    #Gets color scheme for legend
    breakpoints <- seq(from = 0, to = MaximumRich, length = MaximumRich)
    color <- c("lightgrey", colorRampPalette(c("springgreen", "dodgerblue", "darkblue"))(length(breakpoints)))
    
    #Creates pdf of the species richness raster
    pdf(file = paste0(result_dir, "/", "RichnessMaps/", substr(FocusRasters[e], 1, (nchar(FocusRasters[e]) - 4)), ".pdf"))
    plot(legend = FALSE, Raster1, col = color[1:(MaxRaster1 + 1)], xlab = "", ylab = "", main = names(FocusStack[[e]]))
    gradientLegend(c(0:MaximumRich), color = color, pos = 0.125, side = 4, n.seg = 2, fit.margin = TRUE, inside = TRUE)
    dev.off()
  }
  rm(FocusStack)
  gc()
}

print("Making PDF Maps...")
for(t in 1:ntaxa) {
  taxon <- taxa[t]
  MakePDFMaps(taxon)
}