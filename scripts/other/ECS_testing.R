# Script to get ECS data from 20 runs with different ECS values to see how
# changing it would affect Hector's historical predictions
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

OUTPUT <- file.path(RESULTS_DIR, "ecs_plots.jpeg")


source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
co2_data <- get_co2_data(CO2_PATH)
temp_data <- get_temp_data(TEMP_PATH)
temp_data <- filter(temp_data, year <= 2014)

obs_data <- rbind(co2_data, temp_data)

### Running Hector ###
all_data <- run_hector(ini_file = INI_FILE, 
                           params = NULL,
                           vals = NULL,
                           yrs = 1750:2014, 
                           vars = c(GMST(), CONCENTRATIONS_CO2()))
all_data$scenario <- "default"

# This is probably pretty inefficient
for (ecs in seq(2, 5, by=0.2)) {
  curr_data <- run_hector(ini_file = INI_FILE,
                          params = ECS(),
                          vals = ecs,
                          yrs = 1750:2014,
                          vars = c(GMST(), CONCENTRATIONS_CO2()))
  curr_data$scenario <- paste("ECS =", ecs)
  all_data <- rbind(all_data, curr_data)
}

# Saving ECS data
save(all_data, file = file.path(RESULTS_DIR, "ecs_sens.rda"))

# Plotting
comb_data <- rbind(obs_data, all_data)

ggplot(data = comb_data, aes(x = year, y = value, color = scenario)) + 
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  ggtitle("Comparing ECS Values")
ggsave(OUTPUT, width = 15)