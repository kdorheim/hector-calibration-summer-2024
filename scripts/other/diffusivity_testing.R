# Script to get diffusivity data from 20 runs with different diffusivity values
# Goal: See how different diffusivity values (using both the Brown et al 2024
# and Dorheim et al 2024 ranges) affect Hector's historical predictions
# Author: Peter Scully
# Date: 6/19/24

### Constants and Imports ###

# Importing libraries
library(hector)
library(ggplot2)

# Setting up file paths
COMP_DATA_DIR <- file.path(here::here(), "comparison_data")
SCRIPTS_DIR <- file.path(here::here(), "scripts")
RESULTS_DIR <- file.path(here::here(), "results")

CO2_PATH <- file.path(COMP_DATA_DIR,
                      "Supplementary_Table_UoM_GHGConcentrations-1-1-0_annualmeans_v23March2017.csv")
TEMP_PATH <-
  file.path(COMP_DATA_DIR,
            "HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")

INI_FILE <- system.file("input/hector_ssp245.ini", package = "hector")

OUTPUT <- file.path(RESULTS_DIR, "diff_plots.jpeg")


source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
co2_data <- get_co2_data(CO2_PATH)
temp_data <- get_temp_data(TEMP_PATH)
temp_data <- filter(temp_data, year <= 2014)

obs_data <- rbind(co2_data, temp_data)
obs_data$category <- "base"

### Running Hector ###
all_data <- run_hector(ini_file = INI_FILE, 
                           params = NULL,
                           vals = NULL,
                           yrs = 1750:2014, 
                           vars = c(GMST(), CONCENTRATIONS_CO2()))
all_data$scenario <- "default"
all_data$category <- "base"

# This is probably pretty inefficient
for (diff in c(seq(0.8, 1.5, by =0.1), seq(2, 2.6, by=0.1))) {
  curr_data <- run_hector(ini_file = INI_FILE,
                          params = DIFFUSIVITY(),
                          vals = diff,
                          yrs = 1750:2014,
                          vars = c(GMST(), CONCENTRATIONS_CO2()))
  curr_data$scenario <- paste("Diff =", diff)
  if (diff < 2) {
    curr_data$category <- "lower range"
  } else {
    curr_data$category <- "upper range"
  }
  all_data <- rbind(all_data, curr_data)
}

# Saving diff data
save(all_data, file = file.path(RESULTS_DIR, "diff_sens.rda"))

# Plotting
comb_data <- rbind(obs_data, all_data)

ggplot(data = comb_data, aes(x = year, y = value, color = scenario)) + 
  geom_line(aes(linetype = category)) +
  facet_wrap(~ variable, scales = "free") +
  ggtitle("Comparing Diffusivity Values")
ggsave(OUTPUT, width = 15)