# Save PRISM daily data to file (should be sourced from another script)
# 3 Dec 2025

library(httr2)
# library(dplyr)
# library(lubridate)
library(stringr)
library(tibble)
# library(terra)

# Create folder to hold temporary prism data
prism_folder <- "prismtmp/"
dir.create(prism_folder)

# Base URL for PRISM Web Service (see https://www.prism.oregonstate.edu/downloads/)
url_base <- "https://services.nacse.org/prism/data/get/us/4km"

# Select tmin or tmax
var <- tmp_variable

# Sequence of dates
current_month <- str_pad(month(today), width = 2, pad = "0")
current_day <- day(today)
# Pick a couple days ago to ensure PRISM data are available
recent_day <- str_pad(current_day - 2, width = 2, pad = "0")

dates <- seq.Date(as.Date(paste0(yr, "-01-01")),
                  as.Date(paste(yr, current_month, recent_day, sep = "-")), 
                  by = 1)

for (i in 1:length(dates)) {
  
  message(paste0("Downloading ", var, " data for ", dates[i]))
  
  date_nodash <- str_remove_all(dates[i], "-")
  
  # Create folder/zip names
  tmp_folder <- paste0("prism_", var, "_us_25m_", date_nodash)
  tmp_zip <- paste0(prism_folder, tmp_folder, ".zip")
  
  # Download weather data
  req <- request(url_base)
  
  req_tmp <- req %>%
    req_url_path_append(var, date_nodash)
  invisible(req_perform(req_tmp, path = tmp_zip))
  
  # Unzip folder
  suppressWarnings(
    utils::unzip(tmp_zip, exdir = paste0(prism_folder, tmp_folder))
  )
  
  # Remove zips
  invisible(file.remove(tmp_zip))
  
  # Read file in as SpatRaster
  tmp_file <- paste0(prism_folder, tmp_folder, "/", tmp_folder, ".tif")
  tmp <- terra::rast(tmp_file)
  names(tmp) <- paste0(var, "_", date_nodash)
  
  # If Jan 1, create new SpatRaster. Otherwise append.
  if (i == 1) {
    tmp_yr <- tmp
  } else {
    tmp_yr <- c(tmp_yr, tmp)
  }
  
} # i

# Convert degrees C to degrees F
tmp_yr <- tmp_yr * 9/5 + 32

# Save RDS file (slightly smaller than GeoTiff)
saveRDS(tmp_yr, paste0("PRISM/", var, "-", yr, ".rds"))

# Remove data files and folders
unlink(prism_folder, recursive = TRUE)
rm(tmp_yr)
rm(tmp)
