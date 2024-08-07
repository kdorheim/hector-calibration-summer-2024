# Script to calculate OHC RMSE for all runs (since we didn't calculate it
# automatically for most runs)
# Author: Peter Scully
# Date: 7/1/24

### Constants and Imports ###

# Importing libraries
library(hector)
library(ggplot2)
theme_set(theme_bw(base_size = 20))

# Setting up file paths
COMP_DATA_DIR <- file.path(here::here(), "comparison_data")
SCRIPTS_DIR <- file.path(here::here(), "scripts")
RESULTS_DIR <- file.path(here::here(), "results")

OHC_PATH <- file.path(COMP_DATA_DIR, "OHC_ensemble_Kuhlbrodt_etal_2022.csv")

INI_FILE <- system.file("input/hector_ssp245.ini", package = "hector")
PARAMS <- c(BETA(), Q10_RH(), DIFFUSIVITY(), ECS(), AERO_SCALE())

OUTPUT <- file.path(RESULTS_DIR, "OHC_RMSE.txt")


source(file.path(SCRIPTS_DIR, "major_functions.R"))

### Getting observational data ###
obs_data <- get_ohc_data(OHC_PATH, include_unc = T)

### Running Hector ###
# Default (and initial smoothing) Results [Exp. 1-4]
default_data <- run_hector(ini_file = INI_FILE, 
                           params = NULL,
                           vals = NULL,
                           yrs = 1750:2014, 
                           vars = HEAT_FLUX())
default_data$scenario <- "Hector - Default"

# NMSE 3-Parameter Results [Exp.5-9]
exp5_9A <- run_hector(ini_file = INI_FILE,
                      params = PARAMS,
                      vals = c(0.268, 2.64, 2.2, 3, 1),
                      yrs = 1750:2014,
                      vars = HEAT_FLUX())
exp5_9A$scenario <- "Hector - NMSE"

exp5B <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0, 1.5, 2.6, 3, 1),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp5B$scenario <- "Hector - NMSE \nBig Box"

exp6B <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0, 1.58, 2.6, 3, 1),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp6B$scenario <- "Hector - NMSE, Smoothing (k = 3) \nBig Box"

exp8B <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0, 1.95, 2.6, 3, 1),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp8B$scenario <- "Hector - NMSE, Smoothing (k = 10) \nBig Box"

exp9B <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0.028, 1.76, 2.6, 3, 1),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp9B$scenario <- "Hector - NMSE w/ unc \nBig Box"


# Optimizing S, Alpha [Exp. 10-11]
exp10A <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.268, 1.95, 2.6, 3.97, 1),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp10A$scenario <- "Hector - NMSE w/ unc \nTuning S"

exp10B <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.006, 1, 2.6, 3.16, 1),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp10B$scenario <- "Hector - NMSE w/ unc \nBig Box, Tuning S"

exp11A <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.57, 1.76, 2.38, 2.96, 0.492),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp11A$scenario <- "Hector - NMSE w/ unc \nTuning S, Alpha"

exp11B <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.502, 0.99, 2, 2.88, 0.5),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp11B$scenario <- "Hector - NMSE w/ unc \nBig Box, Tuning S, Alpha"

# Optimizing for OHC & Further Refinements [Exp. 12-16]
exp12 <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0.65, 1.76, 1.04, 2.33, 0.438),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp12$scenario <- "Hector - NMSE w/ unc, incl. OHC \nTuning S, Alpha"

exp13 <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0.53, 2.31, 1.04, 2.83, 1.405),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp13$scenario <- "Hector - MVSSE, incl. OHC \nTuning S, Alpha"

exp14A <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.732, 1.76, 1.04, 3, 0.613),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp14A$scenario <- "Hector - NMSE w/ unc, incl. OHC \nTuning Alpha"

exp14B <- run_hector(ini_file = INI_FILE,
                     params = PARAMS,
                     vals = c(0.904, 0.88, 0.806, 3, 0.46),
                     yrs = 1750:2014,
                     vars = HEAT_FLUX())
exp14B$scenario <- "Hector - NMSE w/ unc, incl. OHC \nBig Box, Tuning Alpha"

exp15 <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0.57, 2.49, 1.06, 3.14, 1.08),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp15$scenario <- "Hector - MAE w/ unc, incl. OHC \nTuning S, Alpha"

exp16 <- run_hector(ini_file = INI_FILE,
                    params = PARAMS,
                    vals = c(0.59, 1.76, 1.04, 2.17, 0.411),
                    yrs = 1750:2014,
                    vars = HEAT_FLUX())
exp16$scenario <- "Hector - NMAE w/ unc, incl. OHC \nTuning S, Alpha"

# Calculating OHC RMSE
all_exp <- list(default_data, 
             exp5_9A, exp5B, exp6B, exp8B, exp9B,  # NMSEs
             exp10A, exp10B,                       # Add S
             exp11A, exp11B,                       # Add alpha
             exp12,                                # Add OHC, Mat Diff
             exp13,                                # Try MVSSE
             exp14A, exp14B,                       # Try remove S
             exp15, exp16)                         # Try MAE/NMAE

# TODO: clean up this output sorry
sapply(all_exp, get_var_mse_unc, 
       obs_data = obs_data, var = "OHC", yrs = 1957:2014, mse_fn = mse_unc)

