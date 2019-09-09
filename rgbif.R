#Initializations-----------------------------------
library(rgbif)
#Loads in variables from the configuration data frame
spplist <- df[, "spplist"]
occurrences <- df[, "occurrences"]
sppcountsloc <- df[, "sppcountsloc"]
lat <- df[, "decimalLatitude"] 
long <- df[, "decimalLongitude"]
minlat <- as.numeric(df[, "minlat"])
maxlat <- as.numeric(df[, "maxlat"])
minlong <- as.numeric(df[, "minlong"])
maxlong <- as.numeric(df[, "maxlong"])
test <- df[, "test"]
result_dir <- df[, "result_dir"]
dispersalStep <- df[, "dispersalStep"]
if (dispersalStep == "Y") {
  dispersalRate_dir <- df[, "dispersalRate_dir"]
}
#Allows for connection to last longer than 60 seconds
options(timeout = 1000)
#Copies the maxent javascript to the test folder
setwd(test)
file.copy(paste0(DataDirectory, "/maxent.jar"), test)

#Functions------------------------------------------
#Converts headings and name of GBIF data to be used in subsequent scripts
convert <- function(CurSpp) {
  tryCatch({
    CurSpp2 <- (read.csv(CurSpp))
    CurSpp2 <- data.frame(lapply(CurSpp2, as.character), stringsAsFactors = FALSE)
    if (length(grep("_", CurSpp)) > 0) {
      s <- as.character(CurSpp2$species[1])
      s <- gsub(" ", "_" , s)
    } else {
      GenusLevels <- levels(as.factor(CurSpp2$species))
      SpeciesInGenus <- GenusLevels[grep(" ", GenusLevels)]
      s <- as.character(SpeciesInGenus)
      s <- gsub(" ", "_", s)
    }
    
    for (u in 1:length(s)) {
      sSpace <- gsub("_", " ", s)
      CurSppFinal <- CurSpp2[which(CurSpp2$species == sSpace[u]), ]
      names(CurSppFinal)[names(CurSppFinal) == "species"] <- "Species"
      names(CurSppFinal)[names(CurSppFinal) == "decimalLongitude"] <- "Longitude"
      names(CurSppFinal)[names(CurSppFinal) == "decimalLatitude"] <- "Latitude"
      CurSppFinal <- CurSppFinal[, c("Species", "Longitude", "Latitude")]
      write.csv(CurSppFinal, file = paste0("species/", s[u], ".csv"), row.names = FALSE)
    }
  }, error = function(err) {
    print(paste("MY_ERROR: ", spp.name, " ", err))
    return(paste0("error: ", err))
  })
}

#iterate through species
speciterate <- function(er, nspp) {
  FailedSpecies <- c()
  for (i in er:nspp) {
    tryCatch({
      #get species name
      s <- as.character(OurSpp$Scientific.Name[i]) 
      if (length(grep(" ", s)) > 0) {
        print(paste("Species ", i, " of ", nspp, ": ", s, sep = ""))
        print(paste0("   Beginning search: ", Sys.time()))
        #get taxon key
        data <- name_suggest(q = s, rank = "species", fields = c('key', 'canonicalName', 'rank', 'higherClassificationMap'))
        hierarchy <- data$hierarchy
        data <- data$data
        if (nrow(data) > 1) {
          for (h in 1:length(hierarchy)) {
            H1 <- hierarchy[[h]]
            CorrectTaxon <- grep(tolower(OurSpp$Taxon[i]), tolower(H1$name))
            if (length(CorrectTaxon) == 0) {
              data <- data[-h, ]
            }
          }
          if (nrow(data) == 0) {
            message("Taxon not found: check that taxon is consistent with GBIF Taxonomy")
            #Make a list of species failed 
            FailedSpecies <- c(FailedSpecies, OurSpp$Scientific.Name[i])
          }
        }
        OurSpp$OrigOccurrences[i] <- 0
      } else {
        print(paste("Genus ", i, " of ", nspp, ": ", s, sep = ""))
        print(paste0("   Beginning search: ", Sys.time()))
        #get taxon key
        data <- name_suggest(q = s, rank = "genus", fields = c('key', 'canonicalName', 'rank', 'higherClassificationMap')) ## data contains canonical names, keys, and rank
        hierarchy <- data$hierarchy
        data <- data$data
        for (h in 1:length(hierarchy)) {
          H1 <- hierarchy[[h]]
          CorrectTaxon <- grep(tolower(OurSpp$Taxon[i]), tolower(H1$name))
          if (length(CorrectTaxon) == 0) {
            data <- data[-h, ] 
          }
        }
        OurSpp$OrigOccurrences[i] <- 0
        OurSpp$Occurrences[i] <- 0
      }
      ## Begin as empty/null values. Will be updated as applicable in if statement that follows
      Occ <- data.frame()
      keyslist <- NULL
      nameslist <- NULL
      
      if (nrow(data) > 0) {
        for (k in 1:nrow(data)) {
          if((!is.na(data$canonicalName[k])) && (tolower(data$canonicalName[k]) == tolower(s))) {
            #perform search
            Occ1 <- occ_search(taxonKey = data$key[k], 
                               decimalLatitude = lat, 
                               decimalLongitude = long, 
                               return = 'data', 
                               hasCoordinate = TRUE, 
                               limit = 200000, 
                               fields = c('species', 'decimalLatitude', 'decimalLongitude', 'basisOfRecord', 'issues'))
            #if the search found anything
            ## Fix columns of Occ and Occ1 to match if need be before rbind
            if (!is.atomic(Occ1)) {
              #check if all columns have been returned, not to cause problems with rbind to Occ
              Occnames <- colnames(Occ)
              Occ1names <- colnames(Occ1)
              Occnames <- sort(Occnames)
              Occ1names <- sort(Occ1names)
              
              if (ncol(Occ) == 0) { 
              } else if (identical(Occnames, Occ1names)) {
              } else {
                boolean <- data.frame(x = rep(FALSE, ncol(Occ1)))
                
                for(q in 1:length(Occnames)) {
                  add <- TRUE
                  for (w in 1:length(Occ1names)) {
                    if (Occnames[q] == Occ1names[w]) {
                      add <-  FALSE
                      boolean[w,] <- TRUE
                    }  
                  }
                  
                  if(add) {
                    Occ1[, Occnames[q]] <- rep("NA", nrow(Occ1))
                  }
                }
                
                for (r in 1:nrow(boolean)) {
                  if (!boolean[r, ]) {
                    Occ[, Occ1names[r]] <- rep("NA", nrow(Occ))
                  }
                }
              }
              
              OurSpp$OrigOccurrences[i] <- OurSpp$OrigOccurrences[i] + nrow(Occ1)
              Occ <- rbind(Occ,Occ1)
              keyslist <- paste(keyslist, data$key[k])
              OurSpp$Keys[i] <- keyslist
              if (is.null(nameslist)) {
                nameslist <- paste(nameslist, Occ1$species[1], sep = "")
              } else {
                nameslist <- paste(nameslist, ", ", Occ1$species[1], sep = "")
              }
              
              OurSpp$Names[i] <- nameslist
            }
          } else {
            print(paste("While searching for", s, " ignored GBIF data for species ", data$canonicalName[k], " and key ", data$key[k], sep=""))
          }
        }
      }
      
      print(paste0("   Finishing search: ", Sys.time()))
      print(paste0("   Number original occurrences: ", nrow(Occ)))
      print(paste0("   Keys List: ", keyslist))
      
      if(nrow(Occ) == 0) {
        print(paste0("   Species failed, no search data found: ", s))
      } else if(!is.atomic(Occ)) {
        #don't want fossils
        Occ <- Occ[!Occ$basisOfRecord == "FOSSIL_SPECIMEN", ]
        #don't want geographical issues
        Occ <- Occ[!grepl("cdiv",  Occ$issues), ]
        Occ <- Occ[!grepl("cdout",  Occ$issues), ]
        Occ <- Occ[!grepl("cdrepf",  Occ$issues), ]
        Occ <- Occ[!grepl("cdreps",  Occ$issues), ]
        Occ <- Occ[!grepl("gdativ",  Occ$issues), ]
        Occ <- Occ[!grepl("preneglat",  Occ$issues), ]
        Occ <- Occ[!grepl("preneglon",  Occ$issues), ]
        Occ <- Occ[!grepl("preswcd",  Occ$issues), ]
        Occ <- Occ[!grepl("txmatnon",  Occ$issues), ]
        Occ <- Occ[!grepl("zerocd",  Occ$issues), ]
        #remove duplicates, log those removed
        Occ <- Occ[!duplicated(data.frame(Occ$decimalLatitude, Occ$decimalLongitude)), ]
        OurSpp$Occurrences[i] <- nrow(Occ)
        
        #clip to the study area to log counts, can change these lat/longs to anything to log counts here
        SA_Occ <- Occ[as.numeric(Occ$decimalLatitude) <= maxlat, ]
        SA_Occ <-SA_Occ[as.numeric(SA_Occ$decimalLatitude) >= minlat, ]
        SA_Occ <-SA_Occ[as.numeric(SA_Occ$decimalLongitude) <= maxlong, ]
        SA_Occ<- SA_Occ[as.numeric(SA_Occ$decimalLongitude) >= minlong, ]
        OurSpp$SA_Occur[i] <- nrow(SA_Occ) 
        print(paste0("   Number in the study area: ", nrow(SA_Occ)))
        
        OurSpp$SpeciesSearched[i] <- s
        # record details about species as found by GBIF
        
        sppkeys <- OurSpp$Keys[i]
        sppkeys <- unlist(strsplit(sppkeys, split = " "))
        sppkeys <- sppkeys[sppkeys != ""]
        specurl <- paste("http://api.gbif.org/v1/species/", trimws(sppkeys[1]), sep="")
        URLRead <- colnames(read.csv(specurl))
        gbifapidata <- gsub("\\.", "_", URLRead)
        if (length(grep("^species_", gbifapidata)) > 0) {
          OurSpp$Species[i] <- paste(as.character(substr(gbifapidata[grep("^species_", gbifapidata)], 9, nchar(gbifapidata[grep("^species_", gbifapidata)]))))
        } else {
          OurSpp$Species[i] <- as.character(substr(gbifapidata[grep("^genus_", gbifapidata)], 7, nchar(gbifapidata[grep("^genus_", gbifapidata)])))
        }
        if (length(grep("^class_", gbifapidata)) > 0) {
          OurSpp$Class[i] <- paste(as.character(substr(gbifapidata[grep("^class_", gbifapidata)], 7, nchar(gbifapidata[grep("^class_", gbifapidata)]) - 1)))
        }
        if (length(grep("^family_", gbifapidata)) > 0) {
          OurSpp$Family[i] <- paste(as.character(substr(gbifapidata[grep("^family_", gbifapidata)], 8, nchar(gbifapidata[grep("^family_", gbifapidata)]))))
        }
        if (length(grep("^genus_", gbifapidata)) > 0) {
          OurSpp$Genus[i] <- paste(as.character(substr(gbifapidata[grep("^genus_", gbifapidata)], 7, nchar(gbifapidata[grep("^genus_", gbifapidata)]))))
        }
        
        # Find status from IUCN
        if (length(grep(" ", s)) > 0) {
          spec <- gsub("_", "%20", tolower(OurSpp$Species[i]))
          url <- paste("http://apiv3.iucnredlist.org/api/v3/species/",spec,"?token=dbb7d8e67e97f2666db295c05334569d8cb16f42cc977eb9860814c99e7eb368", sep = "")
          iucndata <- colnames(read.csv(url))
          iucndata <- gsub("\\.", "_", iucndata)
          ## If page contains data
          if(length(grep("^category_", iucndata)) == 1) {
            ## Get data
            OurSpp$Status[i] <- substr(iucndata[grep("^category_", iucndata)], 10, nchar(iucndata[grep("^category_", iucndata)]))
          } else {
            OurSpp$Status[i] <- NA
          }
        }
        
        # write occurrences
        setwd(occurrences)
        # Use species name from our SppList file as occurrence file name for the csv
        write.csv(Occ, file = paste(gsub(" ", "_", s), ".csv", sep = "")) 
        print(paste0("   Finishing species: ", Sys.time()))
      } else {
        print(paste0("   Species failed, no search data found: ", s))
      }
      p <<- i
      OurSpp <<- OurSpp
      #Create csv file with list of species
      #ADD SPECIES WITH 0 OCCURRENCES TO FAILEDSPECIES
      if (OurSpp$Occurrences[i] == 0) {
        FailedSpecies <- c(FailedSpecies, OurSpp$Scientific.Name[i])
        message(paste0("Species Generated in Failed Species List: Check ", paste0(result_dir, "/", "FailedSpecies.csv"), " for failed species"))
      }
    }
    , error = function(e){
    })
  }
  write.csv(FailedSpecies, paste0(result_dir, "/", "FailedSpecies.csv"))
  rm(Occ1, Occ)
  gc()
}

#Run------------------------------------------------
#set up our excel file with species list
SppList <- read.csv(paste0(spplist), strip.white = TRUE, stringsAsFactors = FALSE)
OurSpp <- SppList
OurSpp <- data.frame(lapply(OurSpp, as.character),stringsAsFactors = FALSE)
OurSpp$Taxon <- OurSpp[, 1]
OurSpp$Scientific.Name <- OurSpp[, 2]
nspp <- nrow(OurSpp)
OurSpp <- as.data.frame(cbind(SpeciesSearched = rep(NA, times = nspp), 
                              OurSpp, 
                              OrigOccurrences = rep(0, times = nspp), 
                              Occurrences = rep(0, times = nspp), 
                              SA_Occur = rep(0, nspp), 
                              Keys = rep(NA, times = nspp), 
                              Names = rep(NA, times = nspp), 
                              Family = rep(NA, times = nspp), 
                              Genus = rep(NA, times = nspp), 
                              Status = rep(NA, times = nspp)))

setwd(occurrences)

# Renaming species
for (f in 1:nspp) {
  s <- as.character(OurSpp$Scientific.Name[f])
  OccP <- occ_search(scientificName = s, 
                     decimalLatitude = lat, 
                     decimalLongitude = long, 
                     return = 'data', 
                     hasCoordinate = TRUE, 
                     limit = 2, 
                     fields = c('species'))
  if (OccP[[1]][1] == "no data found, try a different search") {
    message("No Occurrence data found for ", s, ", check spelling or try searching a synonym")
    next
  } else if (is.na(OccP[2, 1])) {
    message("No Occurrence data found for ", s, ", check spelling or try searching a synonym")
    next
  }
  if (OccP[2, 1] != s) {
    if (length(grep(OccP[2, 1], OurSpp$Scientific.Name) > 0)) {
      message(paste0(s, " is the same species as ", OccP[2, 1], " (also listed in spplist), and will therefore be discarded from this analysis"))
    } else {
      message(paste0(s, " will be renamed to: ", OccP[2, 1]))
    }
    if (dispersalStep == "Y") {
      DispDataName <- list.files(dispersalRate_dir, pattern = ".csv")
      DispData <- read.csv(paste0(dispersalRate_dir, "/", DispDataName[1]), stringsAsFactors = FALSE)
      DispSpec <- grep(" ", DispData[1, ])
      Specloc <- grep(paste0(s, "$"), DispData[, DispSpec])
    }
    s <- OccP[2, 1]
    if (dispersalStep == "Y") {
      DispData[Specloc, DispSpec] <- s
    }
    OurSpp$Scientific.Name[f] <- s
    SppList[f, 2] <- s
  }
}

SppList <- unique(SppList)
OurSpp <- unique(OurSpp)
OurSpp <- OurSpp[which(SppList[, 2] != "NA"), ]
SppList <- SppList[which(SppList[, 2] != "NA"), ]
nspp <- nrow(OurSpp)

#Suppress URL warnings
options(warn = -1)

#Extract occurrences from GBIF
speciterate(1, nspp)
while (p < nspp) {
  speciterate(p, nspp)
}

#write to files so we can view
setwd(sppcountsloc)
for (h in 1:ncol(OurSpp)) {
  OurSpp[, h] <- unlist(as.vector(OurSpp[, h]))
}
write.csv(OurSpp, file = df[, "counts"])

setwd(occurrences)
ListSpp <- list.files(pattern = '\\.csv', full.names = TRUE)
nspp <- length(ListSpp)
for(j in 1:nspp) {
  file.copy(paste0(occurrences, "/", ListSpp[j]), test)
}
options(warn = 0)
setwd(test)
dir.create("species")
ListSpp <- list.files(pattern = '\\.csv', full.names = TRUE)
ListSpp <- ListSpp[1:length(ListSpp)]
out<- lapply(ListSpp, function(x) convert(x))
write.csv(SppList, file = paste0(df[, "spplist"]), row.names = FALSE)

if (!exists("DispData") && dispersalStep == "Y") {
  DispDataName <- list.files(dispersalRate_dir, pattern = ".csv")
  DispData <- read.csv(paste0(dispersalRate_dir, "/", DispDataName[1]), stringsAsFactors = FALSE)
}
if (dispersalStep == "Y") {
  write.csv(DispData, file = paste0(dispersalRate_dir, "/", DispDataName), row.names = FALSE)
}
gc()