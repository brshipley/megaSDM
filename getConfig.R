#getConfig
#output: the destination location (where to save the downloaded file).
getConfig <- function(output) {
  ConfURL <- "https://raw.githubusercontent.com/brshipley/megaSDM/master/config.txt"
  download.file(ConfURL, paste0(output, "/config.txt"))
}