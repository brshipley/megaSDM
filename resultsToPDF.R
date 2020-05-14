####resultsToPDF.R####
##Generates PDF maps for all raster data (habitat suitability maps/binary distribution maps)

#Initializations-------------------------
#Loads the necessary variables from "df"
test <- df[, "test"]
result_dir <- df[, "result_dir"]
numScenario <- df[, "numScenario"]
rastertype <- df[, "rastertype"]
if (numScenario > 0) {
  predictenv <- df[, "proj_predictenv"]
} else {
  predictenv <- c()
}

Scenarios <- c()
if (numScenario > 0) {
  if (length(predictenv) > 0) {
    for (ScenIndex in 1:numScenario) {
      if (ScenIndex != 1) {
        cScenario <- paste0("Scenario.", ScenIndex - 1)
      } else {
        cScenario <- "Scenario"
      }
      Scenario <- df[,cScenario]
      Scenarios[ScenIndex] <- Scenario[[1]]
    }
  }
}

#Gets list of species (all of the species)
setwd(test)
spp.list <- list.dirs(path = "samples", full.names = FALSE, recursive = FALSE)
nspp <- length(spp.list)

#Message
print("   Creating species pdf Scenario maps in:")
print(paste0("      ", result_dir, "/SPECIES_NAME/map_pdfs"))

#Run----------------------------------
for (i in 1:nspp) {
  #Creates PDF directory
  setwd(result_dir)
  setwd(spp.list[i])
  dir.create("map_pdfs")
  
  #Lists the paths of the current species distribution maps
  maps <- list.files(full.names = TRUE, pattern = paste0("\\", rastertype, "$"))
  
  if (length(maps) > 0) {
    if (numScenario > 0) {
      #Lists and appends the paths of the non-current species distribution maps
      for(ScenIndex in 1:numScenario) {
        files <- list.files(path = Scenarios[ScenIndex], full.names = TRUE, pattern = paste0("\\", rastertype, "$"))
        maps <- c(maps, files)
      }
    } else {
      maps <- maps
    }
    
    for (j in 1:length(maps)) {
      tryCatch({
        if (substr(maps[j], 1, 1) == ".") {
          #Creates pdf with same title as file name
          title <- substr(maps[j], 3, nchar(maps[j]) - 4)
          pdf(file = paste0("map_pdfs/", title, ".pdf"))
        } else {
          #Finds name of Scenario/time for pdf creation
          compScenario <- strsplit(maps[j], "/")[[1]][1]
          Scenariolength <- 0
          for (k in 1:numScenario) {
            if (compScenario == Scenarios[k]) {
              Scenariolength <- nchar(Scenarios[k])
            }
          }
          #Creates pdf with same title as file name
          title <- substr(maps[j], Scenariolength + 2, nchar(maps[j]) - 4)
          pdf(file = paste0("map_pdfs/", title, ".pdf"))
        }
        #plots the maps on the created pdf file
        plot(raster(maps[j], native = TRUE), main = title)
        dev.off()
        }, error = function(err) {
          print(paste0(err, maps[j]))
      })
    }
  }
}
