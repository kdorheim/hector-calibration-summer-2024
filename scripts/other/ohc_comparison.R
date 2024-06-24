# Script to compare ocean heat content data from Hector with observations
# Author: Peter Scully
# Date: 6/21/24


### Constants and Imports ###

# Importing libraries
library(hector)
library(ggplot2)

# Setting up file paths
COMP_DATA_DIR <- file.path(here::here(), "comparison_data")
SCRIPTS_DIR <- file.path(here::here(), "scripts")
RESULTS_DIR <- file.path(here::here(), "results")

OHC_FILE <- file.path(COMP_DATA_DIR, "OHC_ensemble_Kuhlbrodt_etal_2022.csv")
YEARS <- 1957:2014
BASELINE_START <- 2005
BASELINE_END <- 2014

INI_FILE <- system.file("input/hector_ssp245.ini", package = "hector")
PARAMS <- c(BETA(), Q10_RH(), DIFFUSIVITY(), ECS(), AERO_SCALE())

OUTPUT <- file.path(RESULTS_DIR, "ohc_plot.jpeg")


source(file.path(SCRIPTS_DIR, "major_functions.R"))


### Functions ###


# calc_ohc - function to calculate the ocean heat content for comparison with
#            data from Kuhlbrodt et al
#
# args:
#   params      - vector of Hector params to modify. if NULL, does default run
#   vals        - vector of values to use for those Hector parameters
#   include_unc - boolean indicating whether to include upper/lower bounds on
#                 values (default: FALSE)
#
# returns: Modified Hector data frame with year-by-year OHC values relative to 
#          2005-2014. Note that new data frame will only contain OHC values, not
#          heat flux values
calc_ohc <- function(params, vals, include_unc = F) {
  hect_data <- run_hector(INI_FILE, 
                          params = params, 
                          vals = vals, 
                          yrs = YEARS, 
                          vars = HEAT_FLUX(),
                          include_unc = include_unc)
  
  # Converting heat fluxes to OHC changes by year
  hect_data$value <- hect_data$value * OCEAN_AREA * W_TO_ZJ
  
  # Converting OHC changes by year to total OHC (relative to 1957)
  hect_data$value <- cumsum(hect_data$value)
  hect_data$variable <- "OHC"
  
  # Making OHC relative to 2005-2014 average
  hect_data <- rel_to_interval(data = hect_data,
                               var = "OHC",
                               start = BASELINE_START,
                               end = BASELINE_END)
  
  return(hect_data)
}


# get_ohc_data - function to get ocean heat content data
#
# args:
#   file - path to historical OHC data file
#   scenario - name of scenario being run (default: "historical")
#   include_unc - boolean indicating whether to include uncertainty data
#                 (default: F)
#
# returns: Hector-style data frame with OHC data
get_ohc_data <- function(file, scenario = "historical", include_unc = F) {
  
  # Reading in only OHC data
  ohc_data <- read.table(file, 
                         skip = 2, 
                         sep = ",",
                         colClasses = c("numeric", "NULL", "NULL", "NULL",
                                        "NULL", "NULL", "NULL", "numeric",
                                        "numeric"))
  
  # Fixing table formatting
  ohc_data <- na.omit(ohc_data)
  colnames(ohc_data) <- c("year", "value", "unc")
  
  # Getting rid of non-integer years
  ohc_data$year <- ohc_data$year - 0.5
  
  # Adding in new columns to match Hector data frames
  ohc_data$scenario <- scenario
  ohc_data$variable <- "OHC"
  ohc_data$units <- "ZJ"
  
  # Adding in confidence interval (if applicable)
  if (include_unc) {
    ohc_data$lower <- ohc_data$value - ohc_data$unc
    ohc_data$upper <- ohc_data$value + ohc_data$unc
  }
  
  # Getting rid of raw uncertainty column
  ohc_data$unc <- NULL
  
  return(ohc_data)
}

obs_data <- get_ohc_data(OHC_FILE, include_unc = T)
default_data <- calc_ohc(NULL, NULL, include_unc= T)
reg_box_data <- calc_ohc(PARAMS, c(0.57, 1.76, 2.38, 2.95, 0.49), include_unc= T)
low_diff_data <- calc_ohc(PARAMS, c(0.57, 1.76, 1.1, 2.95, 0.49), include_unc= T)

default_data$scenario <- "Hector - Default"
reg_box_data$scenario <- "Hector - All Params, Reg Box"
low_diff_data$scenario <- "Hector - All Params, Reg Box, Diff = 1.1"

comb_data <- rbind(obs_data, default_data, reg_box_data, low_diff_data)

ggplot(data = comb_data, aes(x = year, y = value, color = scenario)) + 
  geom_ribbon(data = 
                filter(comb_data, scenario == "historical" & variable == "OHC"),
              aes(ymin = lower, ymax = upper),
              fill = 'aquamarine1',
              color = NA) +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  ggtitle("Comparing Parameterizations")
ggsave(OUTPUT, width = 15)