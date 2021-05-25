# read in CTD data, write to file

library(oce)
library(ggplot2)
library(purrr)
library(stringr)
library(lubridate)

# create vector of filenames
ctd_files <- dir(path = "./CTD/", pattern = ".cnv$")

# read CTD data to a list of ctd class objects
ctd <- lapply(paste0("CTD/", ctd_files), read.ctd)

ctd_dec <- lapply(ctd, ctdDecimate)

# function for converting ctd class objects to a data frame
# also extracting date+time, Brunt-vaisala, station ID from filenames
ctd_to_df <- function(ctd, filenames){
  ctd <- subset(ctd, sigmaT>10 & depth > 0)
  d <- as.data.frame(ctd@data)
  d$date_time <- ctd@metadata$startTime
  d$N2 <- swN2(ctd, derivs = "smoothing", df = 50)
  d$station_id <- str_match(filenames, "[.](.*)[.]")[,2]
  d$date <- as.Date(str_match(filenames, "^([^.]*)")[,1])

  return(d)
} 


# apply function to all objects and filenames
ctd_df <- purrr::map2(ctd, ctd_files, ctd_to_df)

ctd_df_dec <- purrr::map2(ctd_dec, ctd_files, ctd_to_df)

# convert list of data frames to 1 data frame
ctd_full <- map_dfr(ctd_df, rbind)
ctd_full_dec <- map_dfr(ctd_df_dec, rbind)

#readr::write_tsv(ctd_full_dec, "CTD_full_dec.tsv")
#readr::write_tsv(ctd_full, "CTD/CTD_full.tsv")
