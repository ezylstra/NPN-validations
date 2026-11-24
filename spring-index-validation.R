# Validating Spring Indices: How well do the Spring Index predictions line up
# with Nature's Notebook observations?
# 24 Nov 2026

# This script is based in large part on previous work done by T. Crimmins

# Steps
# Download individual phenometrics data for cloned lilac, honeysuckle species
# Download predicted spring index values (based on PRISM)
# Compare predicted dates of first leaf/bloom with observed dates
# Merge these data with historical data

# Changes from previous workflows:
# Only using PRISM data for current year (no more NCEP)
# Remove observations after day 225 (~mid Aug). 
  # Latest predicted date in historical data is 206 (and 95% of them before day 150)

library(rnpn)
library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)

# Today
today <- Sys.Date()
# Year
yr <- year(today)

# Load csv with info on species, phenophase
sim <- read.csv("spring-index-matrix.csv")
sim <- sim %>%
  separate_wider_delim(npn_six_layer, 
                       delim = "_", 
                       names = c("spp", "php", "climate"),
                       cols_remove = FALSE) %>%
  select(-climate)

# For 2026 only, extract rasters we need from zipped folder
rast_names <- unique(paste0(sim$spp, "_", sim$php))
zip_names <- unzip("six-2026.zip", list = TRUE)$Name
unzip("six-2026.zip",
      files = str_subset(zip_names, str_flatten(rast_names, collapse = "|")),
      exdir = "six-2026",
      junkpaths = TRUE)

# Loop through each combination of species and phenophase 
# (Previously, we also looped through climate data source PRISM/NCEP, but 
# starting in 2026, we're only using PRISM)
for (i in 1:nrow(sim)) {

  # Get NN observations
  cat(paste0("Evaluating data for ", sim$species[i], ", ", sim$phenophase[i]))
  
  # Download observations 
  obs <- npn_download_individual_phenometrics(
    request_source = "erinz",
    years = yr,
    species_ids = sim$species_id[i],
    phenophase_ids = sim$npn_php_id[i]
  ) %>% data.frame()
  
  # Keep only the first observation of each plant that year
  obs <- obs %>%
    arrange(individual_id, first_yes_doy) %>%
    distinct(individual_id, .keep_all = TRUE) %>%
    data.frame()
  
  # Identify name of geospatial layer with spring index values
  lyr_name <- paste0("si-x:", sim$npn_six_layer[i])
  
  # Identify date for geospatial layer download
  geodate <- paste0(yr, "-09-08")
  
  # Load SpatRaster
  # six <- npn_download_geospatial(
  #   coverage_id = lyr_name, 
  #   date = geodate
  # )
  
  # For 2026, extract raster from those originally in zip folder
  six <- rast(paste0("six-2026/", sim$spp[i], "_", sim$php[i], "_",
                     str_remove_all(geodate, "-"), ".tif"))

  # Get dataframe with unique observation locations and convert to a SpatVector 
  # with the same projection as the geospatial layer
  sites <- distinct(obs, site_id, latitude, longitude) %>%
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  sites <- terra::project(sites, crs(six))
  
  # Extract spring index value from the raster cell each site is located in.
  # (Previously, we had used a buffer around each location, but I'm not sure
  # that was helpful given the size of PRISM cells and the fact that predictions
  # from adjacent cells can differ a lot due to elevation.)
  six_vals <- terra::extract(six, sites, bind = TRUE)
  sites <- data.frame(six_vals)
  colnames(sites)[2] <- "predicted"
  
  # Attach predicted values to observations
  obs <- obs %>%
    left_join(sites, by = "site_id") %>%
    mutate(preddiff = first_yes_doy - predicted)
  
  # Remove any observations without predicted dates (likely because they are
  # outside the US where we don't have climate data)
  obs <- obs %>%
    filter(!is.na(predicted))
  
  # Simplify and format dataframe
  obs <- obs %>%
    select(-contains("last_yes")) %>%
    select(-numdays_until_next_no)

  if (i == 1) {
    allobs <- obs
  } else {
    allobs <- rbind(allobs, obs)
  }
}

# Add "IndividualID-Phenophase-ID-Year" column
allobs <- allobs %>%
  mutate(indivID_ppID_year = paste(individual_id, 
                                   phenophase_id, 
                                   first_yes_year,
                                   sep = "-")) %>%
  mutate(Phenophase = ifelse(phenophase_id == 373, "Leaf", "Bloom"))

# Bring in table w/PRISM values for previous years (1981+)
output_files <- list.files("output", full.names = TRUE) 
historical_file <- str_subset(output_files, "PRISM-NCEP-obs_2025")
historical <- read.csv(historical_file)

# Remove NCEP predictions (and rows that have no PRISM pred/preddiff)
historical <- historical %>%
  select(-contains("NCEP")) %>%
  filter(!is.na(PRISMObsPredDiff))

# Quick check that PRISM diff values calculated correctly:
historical <- historical %>%
  mutate(prismdiff = first_yes_doy - PRISMPred)

table(historical$PRISMObsPredDiff - historical$prismdiff) # No!
table(abs(historical$PRISMObsPredDiff - historical$prismdiff) > 1) # 25
# There's a problem. Values within 1 could be the result of rounding error
# But there 25 values that differ by 5 or more (up to 287!)

# Did these problems occur in one year or more often?
historical  <- historical %>%
  mutate(prob = if_else(abs(PRISMObsPredDiff - prismdiff) > 1, 1, 0))
historical %>%
  group_by(first_yes_year) %>%
  summarize(nobs = n(),
            nprob = sum(prob)) %>%
  data.frame()
# All but one in 2019. One issue 2022.

# For now, will remove these rows from historical dataframe
historical <- historical %>%
  filter(prob == 0) %>%
  select(-c(prismdiff, prob))

# Rename columns with PRISM pred/prediff
historical <- historical %>%
  rename(predicted = PRISMPred,
         preddiff = PRISMObsPredDiff)

# Remove observations past day 225 in current year and historical data
historical <- historical %>%
  filter(first_yes_doy <= 225)
allobs <- allobs %>%
  filter(first_yes_doy <= 225)

# Arrange columns in new data to match historical data
allobs <- allobs %>%
  select(colnames(historical))

# Join historical with current data
merged <- rbind(historical, allobs)

# Few checks:
ggplot(merged) +
  geom_histogram(aes(x = preddiff), binwidth = 1)
quantile(merged$preddiff, 
         probs = c(0, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
# On average, 2 days earlier than predicted (preddiff = -2)
# 95% of observations were within 30 days
# 80% of observations within 16 days

# Are the ones that are really off because they didn't have a prior no?
summary(merged$preddiff)
merged %>%
  filter(!is.na(numdays_since_prior_no)) %>%
  pull(preddiff) %>%
  summary()
# Doesn't look like it. 

# Write to file
write.csv(merged,
          paste0("output/PRISM-Predictions-Allyears-", 
                 str_remove_all(Sys.Date(), "-"), ".csv"),
          row.names = FALSE)


# Some old stuff below....
# -------------------------------------------------------------------------#
# An aside: Explore differences between last year's data in two datasets
allobs24 <- filter(allobs, first_yes_year == 2024)
dim(allobs24) # 405 obs
summary(allobs24$PRISMPred) # 4 NAs 
filter(allobs24, is.na(PRISMPred)) # All NAs are Canada locations, so ok
filter(historical, first_yes_year == 2024 & state %in% c("BC", "ON")) 
# And those same individuals are in the historical dataset

# Identify plants-php from 2024 with predictions
prism24 <- allobs24 %>%
  filter(!is.na(PRISMPred)) %>%
  pull(indivID_ppID_year)
length(prism24) # 401

historical24 <- historical %>%
  filter(first_yes_year == 2024 & !state %in% c("BC", "ON"))
dim(historical24) # 395
setdiff(historical24$indivID_ppID_year, prism24) # all plants in historical are in new data

# There are 6 plant-php in new data that aren't in historical -- maybe data added later?
filter(allobs24, !indivID_ppID_year %in% historical$indivID_ppID_year)

# Now, see whether for all the rest, other column data match
allobs24_compare <- allobs24 %>%
  filter(!state %in% c("BC", "ON")) %>%
  filter(indivID_ppID_year %in% historical$indivID_ppID_year) %>%
  select(-contains("Pred")) %>%
  arrange(indivID_ppID_year, first_yes_doy)
historical24_compare <- historical24 %>%
  select(-contains("Pred")) %>%
  arrange(indivID_ppID_year, first_yes_doy)

all.equal(allobs24_compare, historical24_compare)
# Slight differences in site_id, lat, lon, elevation???

which(allobs24_compare$site_id != historical24_compare$site_id)
which(allobs24_compare$latitude != historical24_compare$latitude)
which(allobs24_compare$longitude != historical24_compare$longitude)
# Two rows: 193, 194

which(allobs24_compare$elevation_in_meters != historical24_compare$elevation_in_meters)
# Most elevations don't match, but that's likely because Jeff edited this
# field recently in the database....

allobs24_compare[193:194, ]
historical24_compare[193:194, ]
weird_sites <- npn_stations(state_code = "CT") %>% 
  filter(station_id %in% c(8210, 2540)) %>%
  data.frame()
lpls <- npn_groups()
lpls %>% filter(network_id == 958)

# These are two sites within the Weir Farm NHS, and somehow they 
# ended up with duplicate plant IDs.
# Won't worry about this for now, but it's weird...
# To merge these datasets, match everything except elevation. 
# -------------------------------------------------------------------------#