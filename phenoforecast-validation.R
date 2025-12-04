# Validating PhenoForecasts: How well do PhenoForecast model predictions line up
# with Nature's Notebook observations?
# 4 Dec 2025

# This script is based in large part on previous work done by T. Crimmins and 
# M. Crimmins.

# Steps: 
# Download status observations for each species and phenophase
# Obtain AGDD value associated with observation location and date
# Use AGDD value to predict whether the species-phenophase is present or not
# Compare predictions with observations

library(rnpn)
library(lubridate)
library(stringr)
library(dplyr)
library(terra)
library(degday)

# Today
today <- Sys.Date()

# Year
yr <- year(today)
# yr <- 2025 # Alternatively, can select something other than current year

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

# Load csv with info on species, phenophase, gdd parameters
pf <- read.csv("phenoforecast-matrix.csv")

# Format GDD start date, gdd type, and pred upper bound
# Add columns to hold validation results
pf <- pf %>%
  mutate(gdd_start = ydm(paste0(yr, "-", gdd_start))) %>%
  mutate(dd_method = case_when(
    gdd_type == "simple average" ~ "simp_avg",
    gdd_type == "double-sine" ~ "dbl_sine",
    gdd_type == "single triangle" ~ "sng_tri",
    .default = NA
  )) %>%
  mutate(pred_upper_bound = ifelse(is.na(pred_upper_bound), 
                                      Inf, pred_upper_bound)) %>%
  mutate(n_obs = NA,
         n_positive = NA,
         n_correct = NA,
         accuracy = NA, # true pos + true neg / n
         n_pred1_obs0 = NA, # false positive
         n_pred0_obs1 = NA) # false negative

# Loop through each phenoforecast
for (i in 1:nrow(pf)) {
  
  cat("Evaluating data for", pf$species[i], " ", str_to_lower(pf$np_phenophase[i]))

  # Download observations 
  obs <- npn_download_status_data(
    request_source = "erinz", 
    years = yr, 
    species_ids = pf$species_id[i]
  ) %>% data.frame()
  
  # If no observations (positive or negative) this year, skip
  if (nrow(obs) == 0) {
    pf$n_obs[i] <- 0
    message("There are 0 observations of ", pf$species[i], " in ", yr)
    next
  }
  
  # Filter for phenophase of interest and start date for GDD accumulations and
  # remove any unknown status observations (status = -1)
  obs <- obs %>%
    mutate(obsdate = ymd(observation_date)) %>%
    filter(phenophase_description == pf$npn_phenophase[i]) %>%
    filter(obsdate >= pf$gdd_start[i]) %>%
    filter(phenophase_status >= 0)

  # If no observations for this phenophase (positive or negative) this year, skip
  if (nrow(obs) == 0) {
    pf$n_obs[i] <- 0
    message("There are 0 observations of ", pf$species[i], " ",
            str_to_lower(pf$npn_phenophase[i]), " in ", yr)
    next
  }

  # If there are no observations of live individuals this year, print message
  # but continue with calculations
  if (sum(obs$phenophase_status) == 0) {
    message("There are no observations of live ", pf$species[i], " ", 
            str_to_lower(pf$npn_phenophase[i]), " in ", yr)
  }
  
  # Summarize information by site
  locs <- obs %>%
    group_by(site_id, latitude, longitude) %>%
    summarize(last_date = max(obsdate), .groups = "keep")
  locsv <- vect(locs, geom = c("longitude", "latitude"), crs = "epsg:4326")
  locsv <- terra::project(locsv, crs(tmin))
  
  # Extract daily min temperatures for each site
  tmins_all <- terra::extract(tmin, locsv, ID = FALSE) %>%
    mutate(site_id = locsv$site_id) %>%
    tidyr::pivot_longer(cols = -site_id, 
                        names_to = "obsdate",
                        values_to = "tmin") %>%
    mutate(obsdate = str_remove_all(obsdate, "tmin_")) %>%
    mutate(obsdate = ymd(obsdate))
  
  # Extract daily max temperatures for each site
  tmaxs_all <- terra::extract(tmax, locsv, ID = FALSE) %>%
    mutate(site_id = locsv$site_id) %>%
    tidyr::pivot_longer(cols = -site_id, 
                        names_to = "obsdate",
                        values_to = "tmax") %>%
    mutate(obsdate = str_remove_all(obsdate, "tmax_")) %>%
    mutate(obsdate = ymd(obsdate))
  
  # Merge min and max temperatures
  tminmax <- left_join(tmins_all, tmaxs_all, by = c("site_id", "obsdate"))
  
  # Warning if any tmin > tmax 
  if (sum(tminmax$tmin > tminmax$tmax) > 0) {
    warning("One or more dates with tmin > tmax!")
  }
  
  # Loop through sites and calculate daily AGDD values
  for (j in 1:nrow(locs)) {
    minmax <- filter(tminmax, site_id == locs$site_id[j])
    agdd <- dd_calc(
      daily_min = minmax$tmin,
      daily_max = minmax$tmax,
      thresh_low = pf$gdd_base[i],
      thresh_up = pf$gdd_upper_threshold[i], # NA and NULL produce same results
      method = pf$dd_method[i],
      cumulative = TRUE
    )
    minmax$agdd <- agdd
    if (j == 1) {
      minmax_all <- minmax
    } else {
      minmax_all <- rbind(minmax_all, minmax)
    }
  }
  
  # Append AGDD values to observation dataframe
  obs <- obs %>%
    left_join(select(minmax_all, site_id, obsdate, agdd),
              by = c("site_id", "obsdate"))
  # Remove observations in last few days (that won't have AGDD values)
  last_temp_day <- ymd(str_remove(last(names(tmin)), "tmin_"))
  obs <- obs %>%
    filter(obsdate <= last_temp_day)

  # Predict whether species-phenophase should be present based on GDD model
  obs <- obs %>%
    mutate(pred = ifelse(agdd >= pf$pred_lower_bound[i] &
                           agdd <= pf$pred_upper_bound[i],
                         1, 0))
  
  # Identify when prediction was correct, false positives, false negatives
  obs <- obs %>%
    mutate(correct = ifelse(phenophase_status == pred, 1, 0),
           pred1_obs0 = ifelse(pred == 1 & phenophase_status == 0, 1, 0),
           pred0_obs1 = ifelse(pred == 0 & phenophase_status == 1, 1, 0))
  
  # Calculate accuracy rate, number of false positives/negatives and add
  # to phenoforecast dataframe
  pf$n_obs[i] <- nrow(obs)
  pf$n_positive[i] <- sum(obs$phenophase_status)
  pf$n_correct[i] <- sum(obs$correct)
  pf$n_pred1_obs0[i] <- sum(obs$pred1_obs0)
  pf$n_pred0_obs1[i] <- sum(obs$pred0_obs1)
}

pf <- pf %>%
  mutate(accuracy = round(n_correct/n_obs, 4))

pf
