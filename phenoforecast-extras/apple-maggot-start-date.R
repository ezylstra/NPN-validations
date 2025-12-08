
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

# Create table to hold results for apple maggot with Jan or Mar start date
ap <- pf %>%
  filter(species == "apple maggot")
ap <- rbind(ap, ap)
ap$gdd_start[2] <- "2025-01-01"

# Download observations 
obs_orig <- npn_download_status_data(
  request_source = "erinz", 
  years = yr, 
  species_ids = ap$species_id[1]
) %>% data.frame()

# What do NN data look like?
table(obs_orig$phenophase_description, obs_orig$phenophase_status)
# No observations of apple maggots in 2025 (all status == 0)
# Phenohases: adults, adults feeding, dead adults, dead larvae, egg laying,
# individuals in a trap, larvae, larvae feeding, mating

for (i in 1:nrow(ap)) {

  # Filter for phenophase of interest and start date for GDD accumulations and
  # remove any unknown status observations (status = -1)
  obs <- obs_orig %>%
    mutate(obsdate = ymd(observation_date)) %>%
    filter(phenophase_description == ap$npn_phenophase[i]) %>%
    filter(obsdate >= ap$gdd_start[i]) %>%
    filter(phenophase_status >= 0)

  # If there are no observations of live individuals this year, print message
  # but continue with calculations
  if (sum(obs$phenophase_status) == 0) {
    message("There are no observations of live ", ap$species[i], " ", 
            str_to_lower(ap$npn_phenophase[i]), " in ", yr)
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
      thresh_low = ap$gdd_base[i],
      thresh_up = ap$gdd_upper_threshold[i], # NA and NULL produce same results
      method = ap$dd_method[i],
      cumulative = FALSE
    )
    minmax$gdd <- gdd
    minmax <- minmax %>%
      mutate(agdd_working = ifelse(obsdate < ap$gdd_start[i], 0, gdd)) %>%
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
  last_temp_day <- ymd(str_remove(last(names(tmin)), "tmin_"))
  obs <- obs %>%
    filter(obsdate <= last_temp_day)
  obs <- obs %>%
    filter(!is.na(agdd))
  
  # Predict whether species-phenophase should be present based on GDD model
  obs <- obs %>%
    mutate(pred = ifelse(agdd >= ap$pred_lower_bound[i] &
                           agdd <= ap$pred_upper_bound[i],
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
  ap$n_obs[i] <- nrow(obs)
  ap$n_positive[i] <- sum(obs$phenophase_status)
  ap$n_pred0_obs0[i] <- sum(obs$pred0_obs0)
  ap$n_pred1_obs1[i] <- sum(obs$pred1_obs1)
  ap$n_pred1_obs0[i] <- sum(obs$pred1_obs0)
  ap$n_pred0_obs1[i] <- sum(obs$pred0_obs1)
  
  # Save obs/temps so we can look at it later:
  assign(paste0("obs", i), obs)
  assign(paste0("temps", i), minmax_all)
    
}
  
ap <- ap %>%
  mutate(n_correct = n_pred0_obs0 + n_pred1_obs1) %>%
  mutate(accuracy = round(n_correct/n_obs, 4))

# Look at predicted presence dates for different AGDD accumulation start dates
temps2 <- temps2 %>% rename(agdd_jan = agdd)
tempsj <- left_join(temps1, 
                    select(temps2, obsdate, agdd_jan), 
                    by = "obsdate") %>%
  mutate(pred = ifelse(agdd >= ap$pred_lower_bound[1] &
                         agdd <= ap$pred_upper_bound[1],
                       1, 0),
         pred_jan = ifelse(agdd_jan >= ap$pred_lower_bound[1] &
                             agdd_jan <= ap$pred_upper_bound[1],
                           1, 0))
tempsj[100:200,] %>% data.frame()
# If AGDD accumulation start date is Jan 1, predicted yes dates for site 47379
# start 8 days earlier (5/6 instead of 5/14) and ends 4 days earlier (6/21
# instead of 6/25).

obs1 %>%
  group_by(pred) %>%
  summarize(earliest = min(observation_date),
            latest = max(observation_date),
            ndates = n_distinct(observation_date)) %>%
  data.frame()
obs2 %>%
  group_by(pred) %>%
  summarize(earliest = min(observation_date),
            latest = max(observation_date),
            ndates = n_distinct(observation_date)) %>%
  data.frame()

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
