####additionalStats.R####
##constructs graphs and calculates statistics for SDM shifts

#Initializations-------------------------------
#Loads the necessary packages
library(graphics)
library(parallel)

#Loads the necessary variables from "df"
result_dir <- df[, "result_dir"]
UrbanAnalysis <- df[, "UrbanAnalysis"]
ProtectedAnalysis <- df[, "ProtectedAnalysis"]
numYear <- as.numeric(df[, "numYear"])
numScenario <- as.numeric(df[, "numScenario"])
dispersalStep <- df[, "dispersalStep"]
dispersalRan <- df[, "dispersalRan"]

if (dispersalStep == "Y") {
  dispersal <- read.csv(list.files(df[, "dispersalRate_dir"], full.names = TRUE)[1])
}

if (numScenario > 0) {
  predictenv <- df[, "proj_predictenv"]
} else {
  predictenv <- c()
}

outfile <- "outfile.txt"

#Creates a vector of years
years <- c()
years <-  as.numeric(df[, grep("^Year", colnames(df))])

#Creates a vector of scenarios
Scenarios <- c()

if (numScenario > 0) {
  #Lists all projected directories
  predictenvdir <- list.dirs(predictenv, recursive = TRUE)
  for (ScenIndex in 1:numScenario) {
    if(ScenIndex != 1) {
      cScenario <- paste0("Scenario.", ScenIndex - 1)
    } else {
      cScenario <- "Scenario"
    }
    Scenario <- df[, cScenario]
    Scenarios[ScenIndex] <- Scenario[[1]]
  }
}

#Functions-----------------------------------
#Creates a bar-graph of the number of cells at all time periods for all scenarios  
getCellsGraph <- function(spp, stats, dispersalApplied) {
  
  #If dispersal has been applies, prints out a new graph
  if (dispersalApplied == "Y") {
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/NumCells.pdf"))
  } else {
    pdf(file = paste0(result_dir, "/", spp, "/Additional Stats/NumCells.pdf"))
  }
  
  #Graphical parameters and barplot
  par(mar = c(8, 4, 4, 2), mfrow = c(1, 1), las = 2)
  ticks <- signif(seq(from = 0, to = max(stats$NumberCells), length.out = 10), digits = 3)
  par(yaxp = c(0, max(stats$NumberCells), 20))
  
  barplot(stats$NumberCells, 
          ylab = "Number of Cells", 
          axes = FALSE, 
          main = spp, 
          col = col13, 
          names.arg = stats$Projection, 
          cex.names = 0.7, 
          cex.lab = 0.7, 
          beside = TRUE)
  axis(2, ticks, cex.axis = 0.6)
  dev.off()
}

#Creates multiple bar-graphs showing number of cells for each time period (one for each scenario)
getDiffScenariosGraph <- function(spp, stats, dispersalApplied) {
  
  #If dispersal has been applied, prints out a new graph
  if (dispersalApplied == "Y") {
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/NumCells", numScenario, "Graphs.pdf"))
  } else {
    pdf(file = paste0(result_dir, "/", spp, "/Additional Stats/NumCells", numScenario, "Graphs.pdf"))
  }
  
  #graphical parameters and bar plot (for each climate scenario)
  old.par <- par(mfrow = c(2, ceiling((max(numScenario, 1) / 2))), las = 2, mar = c(8, 4, 4, 2))
  
  for (ScenIndex in 1:numScenario) {
    ticks <- signif(seq(from = 0, to = max(stats$NumberCells), length.out = 10), digits = 3)
    barplot(stats$NumberCells[c(1, (2 + ((ScenIndex - 1) * (numYear - 1))):((2 + ((ScenIndex - 1) * (numYear - 1))) + numYear - 2))], 
            ylab = "Number of Cells", 
            axes = FALSE, 
            main = Scenarios[ScenIndex], 
            names.arg = years, 
            cex.names = 0.7, 
            cex.axis = 0.5, 
            cex.lab = 0.7, 
            beside = TRUE, 
            col = colYear)
    axis(2, ticks, cex.axis = 0.6)
  }
  dev.off()
  par(old.par)
}

#Creates bar graphs showing change in distriubtion overlap with urban areas and protected areas
getUrbProtGraph <- function(spp, stats, dispersalApplied) {
  nstats <- nrow(stats)
  
  #Calculates percent of the species range protected
  if (ProtectedAnalysis == "Y") {
    PercentProtected <- c()
    for (n in 1:nstats) {
      ncells <- stats$NumberCells[n]
      npro <- stats$Protected[n]
      if (ncells > 0) {
        perpro <- (npro / ncells * 100)
      } else {
        perpro <- 0
      }
      PercentProtected <- c(PercentProtected, perpro)
    }
    PercentProtected <- unlist(rapply(list(PercentProtected), 
                                      f = function(x) ifelse(is.nan(x), 0, x), 
                                      how = "replace"))
  } else {
    PercentProtected <- rep(NA,length = nstats)
  }
  
  #Calculates percent of the species range urbanized
  if (UrbanAnalysis == "Y") {
    PercentUrbanized <- c()
    for (n in 1:nstats) {
      ncells <- stats$NumberCells[n]
      nurb <- stats$Urbanized[n]
      if (ncells > 0) {
        perurb <- (nurb / ncells * 100)
      } else {
        perurb <- 0
      }
      PercentUrbanized <- c(PercentUrbanized, perurb)
    }
    PercentUrbanized <- unlist(rapply(list(PercentUrbanized), 
                                      f = function(x) ifelse(is.nan(x), 0, x), 
                                      how = "replace"))
  } else {
    PercentUrbanized <- rep(NA,length = nstats)
  }
  
  #Combines protected and urbanized data into a matrix
  t <- matrix(c(PercentProtected, PercentUrbanized), ncol = 2)
  colnames(t) <- c("% Protected", "% Urbanized")
  rownames(t) <- stats$Projection

  #If dispersal has been applied, prints out a new graph
  if (dispersalApplied == "Y") {
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/UrbanizedProtected.pdf"))
  } else {
    pdf(file = paste0(result_dir, "/", spp, "/Additional Stats/UrbanizedProtected.pdf"))
  }
  
  #graphical parameters and barplot
  ticks <- ceiling(signif(seq(from = 0, to = max(t, na.rm = TRUE), length.out = 10), digits = 2))
  par(mar = c(8,4,4,2))
  
  UrbanizedProtectedPlot <- barplot(t,
                                    ylab = "Percentage in Range", 
                                    axes = FALSE,
                                    main = spp,
                                    cex.names = 0.7,
                                    cex.axis = 0.7, 
                                    cex.lab = 0.7, 
                                    beside = TRUE, 
                                    col = col13, 
                                    legend = rownames(t), 
                                    args.legend = list(x = "topleft", cex = 0.6))
  axis(2, ticks, cex.axis = 0.6)
  dev.off()
}

#Creates a bar-graph of the change in cells (percent of original) from original
getpercentDiffGraphs <- function(spp, stats, dispersalApplied) {
  nstats <- nrow(stats)
  PercentChange <- c()
  origcells <- stats$NumberCells[1]
  
  #Calculates percent change between time period range and original range
  for (n in 1:nstats) {
    nchange <- stats$CellChange[n]
    if (origcells > 0) {
      perchange <- (nchange / origcells * 100)
    } else {
      perchange <- 0
    }
    PercentChange <- c(PercentChange, perchange)
  }
  
  PercentChange <- matrix(PercentChange)
  rownames(PercentChange) <- stats$Projection
  
  #If dispersal has been applied, prints out a new graph
  if (dispersalApplied == "Y"){
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/CellChange.pdf"))
  } else {
    pdf(file = paste0(result_dir, "/", spp, "/Additional Stats/CellChange.pdf"))
  }
  
  #barplot
  CellChangePlot <- barplot(PercentChange[, 1], 
                            ylab = paste0("Percent Change from ", years[1]), 
                            main = spp, 
                            cex.names = 0.7, 
                            cex.axis = 0.7, 
                            cex.lab = 0.7, 
                            col = col13)
  dev.off()
}

#Creates bar-graph with average number of cells per time period (averages across scenarios)
getMinMaxGraphs <- function(spp, stats, dispersalApplied) {
  numlist <- c(stats$NumberCells[1])
  
  #Calculates average number of cells per time period across scenarios
  for (i in 2:numYear) {
    #Sums all of the cell numbers from each scenario within a time period
    for (j in 0:(max(numScenario, 1) - 1)) {
      if (j == 0) {
        num <- stats$NumberCells[j * (numYear - 1) + (i - 1) + 1]
      } else {
        num <- num + stats$NumberCells[j * (numYear - 1) + (i - 1) + 1]
      }
    }
    num <- num / max(numScenario, 1)
    numlist <- c(numlist, num)
  }
  
  #Calculates maximum range size per time period
  maxlist <- c(0)
  for (i in 2:numYear) {
    tempList <- c()
    for (j in 0:(max(numScenario, 1) - 1)) {
      tempList <- c(tempList, stats$NumberCells[j * (numYear - 1) + (i - 1) + 1])
    }
    maxlist <- c(maxlist, max(tempList))
  }
  
  #Calculates minimum range size per time period
  minlist <- c(0)
  for (i in 2:numYear) {
    tempList <- c()
    for (j in 0:(max(numScenario, 1) - 1)) {
      tempList <- c(tempList, stats$NumberCells[j * (numYear - 1) + (i - 1) + 1])
    }
    minlist <- c(minlist, min(tempList))
  }

  #If dispersal has been applied, prints out a new graph
  if (dispersalApplied == "Y") {
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/AvgNumberCells.pdf"))
  } else {
    pdf(file = paste0(result_dir, "/", spp, "/Additional Stats/AvgNumberCells.pdf"))
  }
  
  #graphical parameters and barplot
  ticks <- signif(seq(from = 0, to = max(numlist), length.out = 10), digits = 2)
  
  avg <- barplot(numlist, 
                 axes=FALSE, 
                 ylab = "Average # of Cells per Decade", 
                 main = spp, 
                 names.arg = years, 
                 cex.names = 0.9, 
                 cex.axis = 0.5, 
                 cex.lab = 0.9, 
                 beside=TRUE, 
                 col=colYear)
  axis(2, ticks, cex.axis = 0.6)
  segments(x0 = avg, y0 = minlist, x1 = avg, y1 = maxlist)
  dev.off()
}

#Creates a bar graph comparing area of suitable habitat between dispersal and regular analysis
DispersalCompare <- function(spp, stats, myvars) {
  setwd(result_dir)
  spp.name <- spp
  setwd(paste0(result_dir, "/", spp.name))
  
  #Reads in the results data frame (not dispersal-constrained)
  statsNoDisp <- read.csv("Results.csv") 
  statsNoDisp <- data.frame(lapply(statsNoDisp, as.character), stringsAsFactors = FALSE)
  statsNoDisp <- statsNoDisp[myvars]
  nstatsNoDisp <- nrow(statsNoDisp)
  statsNoDisp[, 2:ncol(statsNoDisp)] <- lapply(statsNoDisp[, 2:ncol(statsNoDisp)], as.numeric)
  
  for(scen in 1:numScenario) {
    DispComp <- matrix(rep(NA, 2 * (length(years) - 1)), nrow = 2, ncol = length(years) - 1)
    
    #Fill in "DispComp" matrix with area of suitable habitat (both dispersal and regular)
    for(y in 2:length(years)) {
      CurYear <- years[y]
      DispCells <- stats[grep(paste0(Scenarios[scen], "_", years[y], "$"), stats$Projection), 2]
      NoDispCells <- statsNoDisp[grep(paste0(Scenarios[scen], "_", years[y], "$"), statsNoDisp$Projection), 2]
      DispComp[1, (y - 1)] <- NoDispCells
      DispComp[2, (y - 1)] <- DispCells
    }
    
    #Name the output pdf
    rownames(DispComp) <- c("No Dispersal", "Dispersal")
    pdf(file = paste0(result_dir, "/", spp, "/Dispersal Applied Additional Stats/", 
                      Scenarios[scen], "_DispersalCompare.pdf"))
    
    #graphical parameters and barplot
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 8))
    barplot(DispComp,
            main = paste0("Dispersal vs. Non-Dispersal: ", Scenarios[scen]),
            ylab = "Number of Cells",
            names.arg = c(years[2:length(years)]),
            col = c("darkblue", "red"),
            legend = c("NoDisp", "Disp"),
            args.legend = list(x = "bottomright",bty = "n",inset = c(-0.25, 0)),
            beside = TRUE)
    dev.off()
  }
}

#Run-------------------------------------------
run <- function(CurSpp) {
  setwd(result_dir)
  spp.name <- substr(CurSpp, 1, nchar(CurSpp) - 4)
  setwd(paste0(result_dir, "/", spp.name))
  
  #Determine which variables are necessary given Urban and ProtectedAnalysis
  if ((UrbanAnalysis == "Y") | (ProtectedAnalysis == "Y")) {
    if (UrbanAnalysis == "N") {
      myvars <- c("Projection", "NumberCells", "CellChange",           # What columns do you want?
                "T1notT2", "T2notT1", "Overlap", "CentroidX", "CentroidY", "Protected")
    } else if (ProtectedAnalysis == "N") {
      myvars <- c("Projection", "NumberCells", "CellChange", 
                  "T1notT2", "T2notT1", "Overlap", "CentroidX", "CentroidY", "Urbanized")       # Which columns do you want?  
    } else {
      myvars <- c("Projection", "NumberCells", "CellChange",  
                  "T1notT2", "T2notT1", "Overlap", "CentroidX", 
                  "CentroidY", "Urbanized", "Protected")
    }
  } else {
    myvars <- c("Projection", "NumberCells", "CellChange", "T1notT2",
                "T2notT1", "Overlap", "CentroidX", "CentroidY")
  }
  
  #If dispersal has been applied, load the dispersal stats (if not, load regular stats)
  if ((dispersalStep == "Y") && (dispersalRan == "Y")) {
    dir.create("Dispersal Applied Additional Stats")
    stats <- read.csv("Results_Dispersal.csv")
    stats <- data.frame(lapply(stats, as.character), stringsAsFactors = FALSE)
    stats <- stats[myvars]
    nstats <- nrow(stats)
    stats[, 2:ncol(stats)] <- lapply(stats[, 2:ncol(stats)], as.numeric)
  } else {
    dir.create("Additional Stats")
    stats <- read.csv("Results.csv")
    stats <- data.frame(lapply(stats, as.character), stringsAsFactors = FALSE)
    stats <- stats[myvars]
    nstats <- nrow(stats)
    stats[, 2:ncol(stats)] <- lapply(stats[, 2:ncol(stats)], as.numeric)
  }
  
  #If urban analysis or protected analysis is required, construct graphs
  if ((UrbanAnalysis == "Y") | (ProtectedAnalysis == "Y")) {
    getUrbProtGraph(spp.name, stats, dispersalRan)
  }
  
  #construct graphs of area changes
  getCellsGraph(spp.name, stats, dispersalRan)
  getDiffScenariosGraph(spp.name, stats, dispersalRan)
  getpercentDiffGraphs(spp.name, stats, dispersalRan)
  
  if (numScenario > 0) {
    getMinMaxGraphs(spp.name, stats, dispersalRan)
  }
  
  if (dispersalRan == "Y") {
    DispersalCompare(spp.name, stats, myvars)
  }
}

#Graphical parameters for the barplots
colScenario <- colorRampPalette(c("blue", "darkred"))(max(numScenario, 1))  #COLOR SCHEME
colYear <- colorRampPalette(c("darkgreen", "brown"))(numYear)
col13 <- colorRampPalette(c("yellow", "darkgrey"))(max(numScenario, 1) * (numYear - 1) + 1)

setwd(paste0(test, "/samples"))
spp.list <- spp_batch
ListSpp <- c()

#Generates the species list for parallelization
if ((dispersalStep == "Y") && (dispersalRan == "Y")) {
  for (i in 1:length(spp.list)) {
    curspecies <- substr(spp.list[i], 1, nchar(spp.list[i]) - 4)
    if (file.exists(paste0(result_dir, "/", curspecies, "/Results_Dispersal.csv"))) {
      ListSpp <- c(ListSpp, spp.list[i])
    }
  }
  for (w in 1:length(spp.list)) {
    FocSpec <- gsub("_", " ", substr(spp.list[w], 1, nchar(spp.list[w]) - 4))
    DispSpec <- grep(paste0("^", FocSpec, "$"), dispersal[, 1])
    if (length(DispSpec) == 0) {
      message(paste0("No dispersal rate values found for ", FocSpec, ": skipping dispersal rate analysis"))
    }
  }
} else {
  ListSpp <- spp.list
}

#Parallelization
clus <- makeCluster(ncores, outfile = outfile)
clusterExport(clus, varlist = c("colScenario", "colYear", "col13", "result_dir",
                                "UrbanAnalysis", "ProtectedAnalysis", "numYear", "numScenario", "dispersalStep",
                                "dispersalRan", "DispersalCompare", "predictenv", "years","Scenarios",
                                "getCellsGraph", "getDiffScenariosGraph", "getUrbProtGraph",
                                "getpercentDiffGraphs", "getMinMaxGraphs", "run"))
clusterEvalQ(clus, library(graphics))



print("   Creating barplot pdfs in:")
print(paste0("      ", result_dir, "/SPECIES_NAME"))
out <- parLapply(clus, ListSpp, function(x) run(x))
stopCluster(clus)