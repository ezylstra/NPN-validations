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
           pred0_obs1 = ifelse(pred == 0 & phenophase_status == 1, 1, 0),
           pred0_obs0 = ifelse(pred == 0 & phenophase_status == 0, 1, 0),
           pred1_obs1 = ifelse(pred == 1 & phenophase_status == 1, 1, 0))
  
  # Calculate accuracy rate, number of false positives/negatives and add
  # to phenoforecast dataframe
  pf$n_obs[i] <- nrow(obs)
  pf$n_positive[i] <- sum(obs$phenophase_status)
  pf$n_pred0_obs0[i] <- sum(obs$pred0_obs0)
  pf$n_pred1_obs1[i] <- sum(obs$pred1_obs1)
  pf$n_pred1_obs0[i] <- sum(obs$pred1_obs0)
  pf$n_pred0_obs1[i] <- sum(obs$pred0_obs1)
}

pf <- pf %>%
  mutate(n_correct = n_pred0_obs0 + n_pred1_obs1) %>%
  mutate(accuracy = round(n_correct/n_obs, 4))

# Write to file
write.csv(pf, 
          paste0("output/phenoforecast-validation-", yr, ".csv"),
          row.names = FALSE)

# Create bar graph
pf_fig <- pf %>%
  filter(n_obs > 0) %>%
  select(species, phenophase, n_obs, n_positive, n_correct, accuracy,
         n_pred1_obs0, n_pred0_obs1, n_pred0_obs0, n_pred1_obs1) %>%
  mutate(group = paste0(species, ": ", phenophase)) %>%
  arrange(desc(n_obs), desc(n_correct))
group_in_order <- pf_fig$group

pf_figl <- pf_fig %>%
  select(group, accuracy, contains("n_pred")) %>%
  pivot_longer(cols = contains("n_pred"),
               names_to = "type",
               values_to = "count") %>%
  mutate(group = factor(group, levels = group_in_order)) %>%
  mutate(type = factor(type,
                       levels = c("n_pred0_obs1",
                                  "n_pred1_obs0",
                                  "n_pred1_obs1",
                                  "n_pred0_obs0")))

pf_summaries <- pf_figl %>%
  group_by(group) %>%
  summarize(accuracy = round(accuracy[1] *100),
            total_count = sum(count)) %>%
  mutate(acc_label = paste0(accuracy, "%"),
         y = total_count + 10)

barchart <- ggplot(pf_figl, aes(x = group, y = count, fill = type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("gray60", "orange",
                               "forestgreen", "steelblue2"),
                    labels = c("Predicted no, reported yes",
                               "Predicted yes, reported no",
                               "Predicted and reported yes",
                               "Predicted and reported no")) +
  geom_text(data = pf_summaries, 
            aes(x = group, y = y, label = acc_label, fill = NULL), 
            vjust = 0, hjust = 0.50, size = 9/.pt) +
  labs(y = "Count") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.8),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

ggsave("output/phenoforecast-validation-2025-fig.png",
       barchart,
       dpi = 600,
       width = 10,
       height = 5,
       units = "in")
