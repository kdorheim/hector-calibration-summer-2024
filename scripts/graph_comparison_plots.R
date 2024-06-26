# Script to compare plots that use different hector variables
# Author: Peter Scully
# Date: 6/13/24

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
PARAMS <- c(BETA(), Q10_RH(), DIFFUSIVITY(), ECS(), AERO_SCALE())

OUTPUT <- file.path(RESULTS_DIR, "alpha_comparison_plots.jpeg")


source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
co2_data <- get_co2_data(CO2_PATH)
co2_data$lower <- co2_data$value
co2_data$upper <- co2_data$value

temp_data <- get_temp_data(TEMP_PATH, include_unc = T)
temp_data <- filter(temp_data, year <= 2014)

obs_data <- rbind(co2_data, temp_data)

### Running Hector ###
default_data <- run_hector(ini_file = INI_FILE, 
                           params = NULL,
                           vals = NULL,
                           yrs = 1750:2014, 
                           vars = c(GMST(), CONCENTRATIONS_CO2()))
default_data$scenario <- "Hector - Default Fit"


nmse_data <- run_hector(ini_file = INI_FILE,
                        params = PARAMS,
                        vals = c(0.268, 2.64, 2.2, 3, 1),
                        yrs = 1750:2014,
                        vars = c(GMST(), CONCENTRATIONS_CO2()))
nmse_data$scenario <- "Hector - NMSE w/ unc Fit"

nmse_bb_data <- run_hector(ini_file = INI_FILE,
                        params = PARAMS,
                        vals = c(0.028, 1.76, 2.6, 3, 1),
                        yrs = 1750:2014,
                        vars = c(GMST(), CONCENTRATIONS_CO2()))
nmse_bb_data$scenario <- "Hector - NMSE w/ unc, Big Box Fit"

ecs_data <- run_hector(ini_file = INI_FILE,
                                  params = PARAMS,
                                  vals = c(0.268, 1.95, 2.6, 3.97, 1),
                                  yrs = 1750:2014,
                                  vars = c(GMST(), CONCENTRATIONS_CO2()))
ecs_data$scenario <- "Hector - NMSE w/ unc, Fit w/ S"

ecs_bb_data <- run_hector(ini_file = INI_FILE,
                       params = PARAMS,
                       vals = c(0.006, 1, 2.6, 3.16, 1),
                       yrs = 1750:2014,
                       vars = c(GMST(), CONCENTRATIONS_CO2()))
ecs_bb_data$scenario <- "Hector - NMSE w/ unc, Big Box Fit w/ S"

alpha_data <- run_hector(ini_file = INI_FILE,
                       params = PARAMS,
                       vals = c(0.57, 1.76, 2.38, 2.96, 0.492),
                       yrs = 1750:2014,
                       vars = c(GMST(), CONCENTRATIONS_CO2()))
alpha_data$scenario <- "Hector - NMSE w/ unc, Fit w/ S, Alpha"

alpha_bb_data <- run_hector(ini_file = INI_FILE,
                            params = PARAMS,
                            vals = c(0.502, 0.99, 2, 2.88, 0.5),
                            yrs = 1750:2014,
                            vars = c(GMST(), CONCENTRATIONS_CO2()))
alpha_bb_data$scenario <- "Hector - NMSE w/ unc, Big Box Fit w/ S, Alpha"

alpha_ohc_data <- run_hector(ini_file = INI_FILE,
                         params = PARAMS,
                         vals = c(0.65, 1.76, 1.04, 2.33, 0.438),
                         yrs = 1750:2014,
                         vars = c(GMST(), CONCENTRATIONS_CO2()))
alpha_ohc_data$scenario <- "Hector - NMSE w/ unc (incl. OHC), Fit w/ S, Alpha"

#hector_data <- rbind(default_data, nmse_data, nmse_bb_data, ecs_data, ecs_bb_data, alpha_data, alpha_bb_data)
hector_data <- rbind(default_data, alpha_data, alpha_ohc_data)
hector_data$lower <- hector_data$value
hector_data$upper <- hector_data$value

comb_data <- rbind(obs_data, hector_data)

ggplot(data = comb_data, aes(x = year, y = value, color = scenario)) + 
  geom_ribbon(data = 
               filter(comb_data, scenario == "historical" & variable == GMST()),
              aes(ymin = lower, ymax = upper),
              fill = 'orchid1',
              color = NA) +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  ggtitle("Comparing Parameterizations")
ggsave(OUTPUT, width = 16)