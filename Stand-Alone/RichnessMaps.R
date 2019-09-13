#Species Richness Maps
  #speciesrasters: a list of RasterLayers or RasterStackobject
    #Binary distribution data: 0 for absence, 1 for presence
  #plotmap (optional): TRUE/FALSE should a species richness map be displayed in R
  #pdffile (optional): a character string of a filename for a PDF output (must include ".pdf")
  #title (optional): a character string of the PDF and/or displayed map title

RichnessMaps <- function(speciesrasters, plotmap, pdffile, title) {
  require(raster)
  if (is.list(speciesrasters) | is.vector(speciesrasters)) {
    SpeciesStack <- stack(speciesrasters)
  } else {
    SpeciesStack <- speciesrasters
  }
  
  SpeciesRichness <- calc(SpeciesStack, fun = sum)
  
  if((hasArg(plotmap)) & (plotmap == TRUE)) {
    if (hasArg(title)) {
      plot(SpeciesRichness, main = paste0(title))
    } else {
      plot(SpeciesRichness, main = "Species Richness")
    }
  }
  
  if(hasArg(pdffile)) {
    require(plotfunctions)
    maxrich <- cellStats(Raster, stat = max, na.rm = TRUE)
    ValueFreq <- freq(Raster)
    
    #If the maximum species richness occurs in too few pixels to display in pdf, only plots the next highest species richness
    if ((ValueFreq[which(ValueFreq[, 1] == maxrich), 2] < floor(0.00153 * ncell(Raster))) && (maxrich > 1)) {
      Raster[Raster == maxrich] <- maxrich - 1
      message(paste0("The maximum species richness in the study area, ", maxrich, ", covers too few cells to display on PDF:"))
      message("Examine the created raster for an accurate count of pixels.")
      maxrich <- maxrich - 1
    }
    
    #Gets color scheme for legend
    breakpoints <- seq(from = 0, to = maxrich, length = maxrich)
    color <- c("lightgrey", colorRampPalette(c("springgreen", "dodgerblue", "darkblue"))(length(breakpoints)))
    
    #Creates pdf of the species richness raster
    pdf(file = paste0(results, "/", "RichnessMaps/", taxon, "_", scenario, "_", year, ".pdf"))
    if (hasArg(title)) {
      plot(legend = FALSE, Raster, col = color, xlab = "", ylab = "", main = paste0(title))
    } else {
      plot(legend = FALSE, Raster, col = color, xlab = "", ylab = "", main = paste0("SpeciesRichness"))
    }
    gradientLegend(c(0:maxrich), color = color, pos = 0.125, side = 4, n.seg = 2, fit.margin = TRUE, inside = TRUE)
    dev.off()
  }
  
  return(SpeciesRichness)
}


