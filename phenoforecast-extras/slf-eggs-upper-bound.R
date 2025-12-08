
library(rnpn)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)
library(terra)
library(degday)
library(ggplot2)

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
         n_pred0_obs1 = NA, # false negative
         n_pred0_obs0 = NA, # correct: negative
         n_pred1_obs1 = NA) # correct: positive

# Create table to hold results for SLF egg hatch
slf <- pf %>%
  filter(species == "spotted lanternfly" & npn_phenophase == "Nymphs")

# Download observations 
obs_orig <- npn_download_status_data(
  request_source = "erinz", 
  years = yr, 
  species_ids = slf$species_id[1]
) %>% data.frame()

# What do NN data look like?
table(obs_orig$phenophase_description, obs_orig$phenophase_status)

# Filter for phenophase of interest and start date for GDD accumulations and
# remove any unknown status observations (status = -1)
obs <- obs_orig %>%
  mutate(obsdate = ymd(observation_date)) %>%
  filter(phenophase_description == slf$npn_phenophase[1]) %>%
  filter(obsdate >= slf$gdd_start[1]) %>%
  filter(phenophase_status >= 0)

# If there are no observations of live individuals this year, print message
# but continue with calculations
if (sum(obs$phenophase_status) == 0) {
  message("There are no observations of live ", slf$species[1], " ", 
          str_to_lower(slf$npn_phenophase[1]), " in ", yr)
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
  gdd <- dd_calc(
    daily_min = minmax$tmin,
    daily_max = minmax$tmax,
    thresh_low = slf$gdd_base[1],
    thresh_up = slf$gdd_upper_threshold[1], # NA and NULL produce same results
    method = slf$dd_method[1],
    cumulative = FALSE
  )
  minmax$gdd <- gdd
  minmax <- minmax %>%
    mutate(agdd_working = ifelse(obsdate < slf$gdd_start[1], 0, gdd)) %>%
    mutate(agdd = cumsum(agdd_working)) %>%
    select(-agdd_working)
  
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
obs <- obs %>%
  filter(!is.na(agdd))

# Predict whether species-phenophase should be present based on GDD model
obs <- obs %>%
  mutate(pred = ifelse(agdd >= slf$pred_lower_bound[1] &
                         agdd <= slf$pred_upper_bound[1],
                       1, 0))

# Identify when prediction was correct, false positives, false negatives
obs <- obs %>%
  mutate(correct = ifelse(phenophase_status == pred, 1, 0),
         pred1_obs0 = ifelse(pred == 1 & phenophase_status == 0, 1, 0),
         pred0_obs1 = ifelse(pred == 0 & phenophase_status == 1, 1, 0),
         pred0_obs0 = ifelse(pred == 0 & phenophase_status == 0, 1, 0),
         pred1_obs1 = ifelse(pred == 1 & phenophase_status == 1, 1, 0))

select(obs, site_id, observation_date, phenophase_status, agdd, pred, correct)

# One site had a bunch of records in Sep-Nov when AGDD values are very high 
# (>2800) with no positive observations, so they were all classified as 
# incorrect (pred = 1, obs = 0). If there were an upper bound on this SLF stage 
# class, then these would be classified as correct (pred = 0, obs = 0) and 
# accuracy would go from ~70% to ~93%.

