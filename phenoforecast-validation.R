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
library(dplyr)
library(terra)
library(caret)

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
# if they haven't already been downloaded or they need to be overwritten. 
# We'll get data from PRISM at 4-km resolution.
if (!file.exists(tmin_rast) | replace_prism) {
  tmp_variable <- "tmin"
  source("prism-download.R")
}
if (!file.exists(tmax_rast) | replace_prism) {
  tmp_variable <- "tmax"
  source("prism-download.R")
}




