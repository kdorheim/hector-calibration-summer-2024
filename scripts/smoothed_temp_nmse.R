# Script for experiments #6-#8
# Script to run smoothed temperature runs
# To keep running mean of a different number of observations, change ROLL_WIDTH
# Author: Peter Scully
# Date: 6/17/24

### Constants and Imports ###

# Importing libraries
library(hector)
library(zoo)

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

# Setting other constants
PARAMS <- c(BETA(), Q10_RH(), DIFFUSIVITY())
ROLL_WIDTH <- 10

OUTPUT <- file.path(RESULTS_DIR, 
                    paste("08_smooth_temps_nmse_", ROLL_WIDTH, "yr_big_box.txt", sep=""))



source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
co2_data <- get_co2_data(CO2_PATH)
temp_data <- get_temp_data(TEMP_PATH)

# Smoothing temperatures
smoothed_temps <- rollmean(temp_data$value, 
                           k = ROLL_WIDTH, 
                           align = "center",
                           fill = NA)
smoothed_temp_data <- temp_data
smoothed_temp_data$variable <- "Smooth T"
smoothed_temp_data$value <- smoothed_temps

obs_data <- rbind(co2_data, temp_data, smoothed_temp_data)

### Calling optim ###
best_pars <- run_optim(obs_data = obs_data,
                       ini_file = INI_FILE,
                       params = PARAMS,
                       lower = c(0, 2.2 - 0.44 * 3, 2.3 - 0.1 * 3),
                       upper = c(0.5 + 0.232 * 3, 2.2 + 0.44 * 3, 2.3 + 0.1 * 3),
                       yrs = 1750:2014,
                       vars = c(GMST(), CONCENTRATIONS_CO2()),
                       error_fn = smooth_T_CO2_nmse,
                       method = "L-BFGS-B",
                       output_file = OUTPUT)

### Outputting individual MSEs ###
hector_data <- run_hector(ini_file = INI_FILE, 
                          params = PARAMS, 
                          vals = best_pars, 
                          yrs = 1750:2014, 
                          vars = c(GMST(), CONCENTRATIONS_CO2()))

T_mse <- get_var_mse(obs_data = obs_data, 
                     hector_data = hector_data, 
                     var = GMST(), 
                     yrs = 1850:2014)
CO2_mse <- get_var_mse(obs_data = obs_data, 
                       hector_data = hector_data, 
                       var = CONCENTRATIONS_CO2(), 
                       yrs = c(1750, 1850:2014))

write_metric("CO2 MSE:", CO2_mse, OUTPUT)
write_metric("T MSE:  ", T_mse, OUTPUT)
write_metric("RMSE:   ", sqrt(mean(CO2_mse, T_mse)), OUTPUT) # not 100% sure this is how we want to calculate this
write("", OUTPUT, append = TRUE)

### Outputting table metrics ###
calc_table_metrics(PARAMS, best_pars, OUTPUT)