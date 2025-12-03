# Testing degday package
# 3 Dec 2025

library(degday)
library(rnpn)
library(terra)
library(dplyr)
library(stringr)

# This is not a great test of the package since I can't download an AGDD time
# series from rnpn based on PRISM data and the only temperature rasters I have
# available are from PRISM, but it should at least give me a general sense of 
# things...

# Download AGDD time series using rnpn (with end date before Jun 2025 to avoid
# errors)
ds_npn <- npn_get_custom_agdd_time_series(
  method = "double-sine",
  start_date = "2025-03-01",
  end_date = "2025-04-01",
  base_temp = 50,
  climate_data_source = "NCEP", 
  temp_unit = "Fahrenheit", 
  lat = 40.0,
  lon = -100.0, 
  upper_threshold = 90
)
sa_npn <- npn_get_custom_agdd_time_series(
  method = "simple",
  start_date = "2025-03-01",
  end_date = "2025-04-01",
  base_temp = 50,
  climate_data_source = "NCEP", 
  temp_unit = "Fahrenheit", 
  lat = 40.0,
  lon = -100.0, 
  upper_threshold = 90
)

# Read in PRISM data from rasters I downloaded previously
tmin <- readRDS("C:/Users/erin/Documents/PRISMData/tmin-2025.rds")
tmax <- readRDS("C:/Users/erin/Documents/PRISMData/tmax-2025.rds")

# Extract layers
dates <- seq.Date("2025-03-01", "2025-04-01")
tminsub <- tmin[[paste0("tmin_", str_remove_all(dates, "-"))]]
tmaxsub <- tmax[[paste0("tmax_", str_remove_all(dates, "-"))]]

# Extract temperature data for specific location and convert to deg F
loc <- vect(data.frame(lat = 40, lon = -100), crs = "epsg:4269")
tmins <- terra::extract(tminsub, loc, ID = FALSE)
tmaxs <- terra::extract(tmaxsub, loc, ID = FALSE)
temps <- data.frame(tmin = t(tmins), tmax = t(tmaxs)) %>%
  mutate(tmin = tmin * 9/5 + 32,
         tmax = tmax * 9/5 + 32)

# Use degday package to calculate GDD and AGDD values for double-sine method
ds_gdd <- dd_dbl_sine(
  daily_min = temps$tmin,
  daily_max = temps$tmax,
  thresh_low = 50,
  thresh_up = 90,
  cumulative = FALSE
)
ds_agdd <- dd_dbl_sine(
  daily_min = temps$tmin,
  daily_max = temps$tmax,
  thresh_low = 50,
  thresh_up = 90,
  cumulative = TRUE
)

# Compare output from double-sine methods
ds <- data.frame(date = as.character(dates), 
                 gdd_degday = ds_gdd,
                 agdd_degday = ds_agdd)
ds <- left_join(ds_npn, ds, by = "date")
cor(ds[, c("gdd", "gdd_degday")]) # 0.9978

# Use degday package to calculate GDD and AGDD values for simple method
sa_gdd <- dd_simp_avg(
  daily_min = temps$tmin,
  daily_max = temps$tmax,
  thresh_low = 50,
  thresh_up = 90,
  cumulative = FALSE
)
sa_agdd <- dd_simp_avg(
  daily_min = temps$tmin,
  daily_max = temps$tmax,
  thresh_low = 50,
  thresh_up = 90,
  cumulative = TRUE
)

# Compare output from simple averaging methods
sa <- data.frame(date = as.character(dates), 
                 gdd_degday = sa_gdd,
                 agdd_degday = sa_agdd)
sa <- left_join(sa_npn, sa, by = "date")
cor(sa[, c("gdd", "gdd_degday")]) # 0.9973

# Correlation between GDDs based on different methods?
cor(sa$gdd, ds$gdd) # 0.8871
cor(sa$gdd_degday, ds$gdd_degday) # 0.8819
plot(sa$gdd, ds$gdd) # Lots more 0 GDD values using simple averaging

# All of this suggests that the NPN and the degday package are calculating 
# GDD values in the same ways. Differences are probably all due to differences
# between min/max temperatures based on NCEP and PRISM datasets.
