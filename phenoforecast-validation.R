# Validating PhenoForecasts: How well do PhenoForecast model predictions line up
# with Nature's Notebook observations?
# 3 Dec 2025

# This script is based in large part on previous work done by T. Crimmins and 
# M. Crimmins.

# Steps: 
# Download status observations for each species and phenophase
# Obtain AGDD value associated with observation location and date
# Use AGDD value to predict whether the species-phenophase is present or not
# Compare predictions with observations (confusion matrix)

library(rnpn)
library(lubridate)
library(stringr)
library(dplyr)
library(terra)
library(caret)

# Today
today <- Sys.Date()

# Year
yr <- year(today)
# yr <- 2025 # Alternatively, can select something other than current year

# Load csv with info on species, phenophase, gdd parameters
pf <- read.csv("phenoforecast-matrix.csv")
# Format a few things
pf <- pf %>%
  mutate(gdd_start = ydm(paste0(yr, "-", gdd_start))) %>%
  mutate(gdd_upper_threshold = ifelse(gdd_upper_threshold == "NULL",
                                      NA, gdd_upper_threshold),
         gdd_upper_threshold = as.numeric(gdd_upper_threshold))

# Files with temperature data
tmin_rast <- paste0("PRISM/tmin-", yr, ".rds")
tmax_rast <- paste0("PRISM/tmax-", yr, ".rds")

# Logical to identify whether temperature data should be overwritten
replace_prism <- FALSE

# First, we'll need daily minimum/maximum temperatures for the year of interest
# if they haven't already been downloaded or if they need to be overwritten. 
# We'll get data from PRISM at 4-km resolution.
if (!file.exists(tmin_rast) | replace_prism) {
  tmp_variable <- "tmin"
  source("prism-download.R")
}
if (!file.exists(tmax_rast) | replace_prism) {
  tmp_variable <- "tmax"
  source("prism-download.R")
}

# Load weather data
tmin <- readRDS(tmin_rast)
tmax <- readRDS(tmax_rast)

# TODO: 
# Add columns to pf dataframe that will hold results:
  # Number of NPN observations
  # Number of positive observations 
  # Confusion matrix elements

# Loop through each phenoforecast
# for (i in 1:nrow(pf)) {
i = 3

  # Download observations 
  obs <- npn_download_status_data(
    request_source = "erinz", 
    years = yr, 
    species_ids = pf$species_id[i]
  ) %>% data.frame()
  
  # Filter for phenophase of interest and start date for GDD accumulations and
  # remove any unknown status observations (status = -1)
  obs <- obs %>%
    mutate(observation_date = ymd(observation_date)) %>%
    filter(phenophase_description == pf$npn_phenophase[i]) %>%
    filter(observation_date >= pf$gdd_start[i]) %>%
    filter(phenophase_status >= 0)
  
  # If no observations (positive or negative) this year, skip
  if (nrow(obs) == 0) {
    message("There are 0 observations of ", pf$species[i], " in ", yr)
    next
  }

  # If there are no observations of live individuals this year, print message
  # but continue with calculations
  if (sum(obs$phenophase_status == 0)) {
    message("There are no observations of live ", pf$species[i], " ", 
            str_to_lower(pf$npn_phenophase[i]), " in ", yr)
  }
  
  # Summarize information by site
  locs <- obs %>%
    group_by(site_id, latitude, longitude) %>%
    summarize(last_date = max(observation_date), .groups = "keep")
  locsv <- vect(locs, geom = c("longitude", "latitude"), crs = "epsg:4326")
  
  # Calculate AGDD values (for entire raster I think....)
  # Extract AGDD values for each observation
  # Append AGDD values to observation dataframe
  # Make prediction based on pred_lower/upper_bounds
  # Calculate confusion matrix
  
# }

