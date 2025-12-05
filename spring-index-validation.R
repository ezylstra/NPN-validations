# Validating Spring Indices: How well do the Spring Index predictions line up
# with Nature's Notebook observations?
# 5 Dec 2025

# This script is based in large part on previous work done by T. Crimmins

# Steps
# Download individual phenometrics data for cloned lilac, honeysuckle species
# Download predicted spring index values (for current year with NCEP; for
# previous year with PRISM)
# Compare predicted dates of first leaf/bloom with observed dates
# Merge these data with historical data

library(rnpn)
library(lubridate)
library(stringr)
library(dplyr)
library(terra)

# Load csv with info on species, phenophase, gdd parameters
sim <- read.csv("spring-index-matrix.csv")

# Loop through each combination of species, phenophase, climate data source
for (i in 1:nrow(sim)) {

  # Get NN observations
  cat(paste0("Evaluating data for ", sim$species[i], ", ", sim$phenophase[i], 
             ", ", sim$climate_source[i]))
  
  # Download observations 
  obs <- npn_download_individual_phenometrics(
    request_source = "erinz",
    years = sim$year[i],
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
  
  # Identify date for geospatial layer download (PRISM will be Jan 1, and 
  # NCEP need to be date later in the current year)
  geodate <- ifelse(sim$climate_source[i] == "PRISM",
                    paste0(sim$year[i], "-01-01"),
                    paste0(sim$year[i], "-07-01")) # Error if I try much later than this
  
  # Load SpatRaster
  six <- npn_download_geospatial(
    coverage_id = lyr_name, 
    date = geodate
  )
  
  # Get dataframe with unique observation locations and convert to a SpatVector 
  # with the same projection as the geospatial layer
  sites <- distinct(obs, site_id, latitude, longitude) %>%
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  sites <- terra::project(sites, crs(six))
  
  # Create polygons layer with a X-m circular buffer around each location
  radius <- ifelse(sim$climate_source[i] == "PRISM", 3000, 2000)
  buffers <- terra::buffer(sites, width = radius)
  
  # Average spring index values within radius of each plant observed
  six_mns <- terra::extract(six, buffers, "mean", na.rm = TRUE, bind = TRUE)
  sites <- data.frame(six_mns)
  colnames(sites)[2] <- "predicted"
  
  # Attach predicted values to observations
  obs <- obs %>%
    left_join(sites, by = "site_id") %>%
    mutate(preddiff = first_yes_doy - predicted)
  
  # Simplify and format dataframe
  obs <- obs %>%
    select(-contains("last_yes")) %>%
    select(-numdays_until_next_no)
  if (sim$climate_source[i] == "PRISM") {
    obs <- obs %>%
      mutate(NCEPPred = NA,
             NCEPObsPredDiff = NA) %>%
      rename(PRISMPred = predicted,
             PRISMObsPredDiff = preddiff)
  } else {
    obs <- obs %>%
      rename(NCEPPred = predicted,
             NCEPObsPredDiff = preddiff) %>%
      mutate(PRISMPred = NA,
             PRISMObsPredDiff = NA)
  }

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

# Bring in NCEP & PRISM table w/values for previous years (1981+)
historical <- read.csv("PRISM-NCEP-obs_12-18-24.csv")

# Remove first column with row numbers
historical <- historical[, -1]

# Arrange columns in new data to match historical data
allobs <- allobs %>%
  select(colnames(historical))

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
    
# Create function to get max value or return NA if all values in vector are NA
max_na <- function(x) {
  if(sum(is.na(x)) == length(x)) {
    y <- NA
  } else {
    y <- max(x, na.rm = TRUE)
  }
  return(y)
}

# Join historical with current data
merged <- rbind(historical, allobs) %>%
  group_by(across(-c(NCEPPred, NCEPObsPredDiff, 
                     PRISMPred, PRISMObsPredDiff,
                     elevation_in_meters))) %>%
  summarize(NCEPPred = max_na(NCEPPred),
            NCEPObsPredDiff = max_na(NCEPObsPredDiff),
            PRISMPred = max_na(PRISMPred),
            PRISMObsPredDiff = max(PRISMObsPredDiff),
            elevation_in_meters = elevation_in_meters[1],
            .groups = "keep") %>%
  data.frame()

# A few impossible values in the historical dataset, so removing them
merged <- merged %>%
  mutate(NCEPObsPredDiff = ifelse(NCEPPred < 0 | NCEPPred > 366, NA, NCEPObsPredDiff),
         NCEPPred = ifelse(NCEPPred < 0 | NCEPPred > 366, NA, NCEPPred),
         PRISMObsPredDiff = ifelse(PRISMPred < 0 | PRISMPred > 366, NA, PRISMObsPredDiff),
         PRISMPred = ifelse(PRISMPred < 0 | PRISMPred > 366, NA, PRISMPred))

# Checks -- 
# count(filter(merged, first_yes_year %in% 2016:2023),
#       !is.na(NCEPPred), !is.na(PRISMPred))
# count(filter(merged, first_yes_year == 2024),
#       !is.na(NCEPPred), !is.na(PRISMPred))
# count(filter(merged, first_yes_year == 2025),
#       !is.na(NCEPPred), !is.na(PRISMPred))

# Put columns back in order 
merged <- merged %>%
  select(colnames(historical))

# Write to file
write.csv(merged,
          paste0("output/PRISM-NCEP-obs_", Sys.Date(), ".csv"),
          row.names = FALSE)
