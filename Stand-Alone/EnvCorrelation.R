#Environmental Raster Correlation:
  #envlayers: a list of RasterLayers or a RasterStack object
  #threshold (optional): number between 0 and 1
    #if supplied, warning outputs will be generated for any two rasters with pearson correlation coefficients higher than %threshold% 

EnvCorrelation <- function(envlayers, threshold) {
  TrainingStack <- stack(envlayers)
  EnvCorr <- layerStats(TrainingStack, stat = "pearson", na.rm = TRUE)
  EnvCC <- EnvCorr$`pearson correlation coefficient` 
  
  if (hasArg(threshold)) {
    for (i in 1:nrow(EnvCC)) {
      EnvCC2 <- EnvCC[, -i]
      HighCor <- which(abs(EnvCC2[i, ]) > threshold)
      if (length(HighCor > 0)) {
        message(paste0("Warning! The environmental layer ", rownames(EnvCC2)[i], " is highly correlated with ", length(HighCor), " other layers:"))
        for (w in 1:length(HighCor)) {
          CorLayer <- names(HighCor[w])
          CorValue <- EnvCC2[i, HighCor[w]]
          message (paste0("    ", CorLayer, " (", round(CorValue, 3), ")"))
        }
      }
    }
  }
  return(EnvCC)
}