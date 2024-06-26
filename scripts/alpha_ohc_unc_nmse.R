# Script for experiment 12
# Script to use normalized T, OHC and CO2 MSEs while accounting for T, OHC unc
# Also includes ECS and alpha as params to optimize over
# Uses OHC range from Matilda table rather than manuscript script
# Author: Peter Scully
# Date: 6/25/24

### Constants and Imports ###

# Importing libraries
library(hector)

# Setting up file paths
COMP_DATA_DIR <- file.path(here::here(), "comparison_data")
SCRIPTS_DIR <- file.path(here::here(), "scripts")
RESULTS_DIR <- file.path(here::here(), "results")

CO2_PATH <- file.path(COMP_DATA_DIR,
                      "Supplementary_Table_UoM_GHGConcentrations-1-1-0_annualmeans_v23March2017.csv")
TEMP_PATH <-
  file.path(COMP_DATA_DIR,
            "HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")

OHC_PATH <- file.path(COMP_DATA_DIR, "OHC_ensemble_Kuhlbrodt_etal_2022.csv")


INI_FILE <- system.file("input/hector_ssp245.ini", package = "hector")
PARAMS <- c(BETA(), Q10_RH(), DIFFUSIVITY(), ECS(), AERO_SCALE())

OUTPUT <- file.path(RESULTS_DIR, "12_alpha_ecs_ohc_unc_nmse.txt")


source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
co2_data <- get_co2_data(CO2_PATH, include_unc = TRUE)
temp_data <- get_temp_data(TEMP_PATH, include_unc = TRUE)
ohc_data <- get_ohc_data(OHC_PATH, include_unc = T)
obs_data <- rbind(co2_data, temp_data, ohc_data)

### Calling optim ###
best_pars <- run_optim(obs_data = obs_data,
                       ini_file = INI_FILE,
                       params = PARAMS,
                       lower = c(0.5 - 0.232, 2.2 - 0.44, 1.16 - 0.118, 2, 0),
                       upper = c(0.5 + 0.232, 2.2 + 0.44, 1.16 + 0.118, 5, 3),
                       yrs = 1750:2014,
                       vars = c(GMST(), CONCENTRATIONS_CO2(), HEAT_FLUX()),
                       error_fn = mean_T_CO2_OHC_nmse_unc,
                       include_unc = T,
                       method = "L-BFGS-B",
                       output_file = OUTPUT)

### Outputting individual MSEs ###
hector_data <- run_hector(ini_file = INI_FILE, 
                          params = PARAMS, 
                          vals = best_pars, 
                          yrs = 1750:2014, 
                          vars = c(GMST(), CONCENTRATIONS_CO2(), HEAT_FLUX()))

T_mse <- get_var_mse(obs_data = obs_data, 
                     hector_data = hector_data, 
                     var = GMST(), 
                     yrs = 1850:2014)
CO2_mse <- get_var_mse(obs_data = obs_data, 
                       hector_data = hector_data, 
                       var = CONCENTRATIONS_CO2(), 
                       yrs = c(1750, 1850:2014))
T_mse_unc <- get_var_mse_unc(obs_data = obs_data,
                             hector_data = hector_data,
                             var = GMST(),
                             yrs = 1850:2014,
                             mse_fn = mse_unc)
OHC_mse_unc <- get_var_mse_unc(obs_data = obs_data, 
                                hector_data = hector_data, 
                                var = "OHC", 
                                yrs = 1957:2014,
                                mse_fn = mse_unc)

write_metric("CO2 MSE:", CO2_mse, OUTPUT)
write_metric("T MSE:  ", T_mse, OUTPUT)
write_metric("RMSE:   ", sqrt(mean(CO2_mse, T_mse)), OUTPUT) # not 100% sure this is how we want to calculate this
write_metric("T MSE accounting for unc:", T_mse_unc, OUTPUT)
write_metric("OHC MSE accounting for unc:", OHC_mse_unc, OUTPUT)
write("", OUTPUT, append = TRUE)

### Outputting table metrics ###
calc_table_metrics(PARAMS, best_pars, OUTPUT)