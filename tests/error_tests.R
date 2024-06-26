# File to test the functions in scripts/error_functions.R

# loading functions
source(here::here("scripts", "error_functions.R"))

library(assertthat)
library(dplyr)

# Confirming that mse works as expected
mse_tests <- function() {
  assert_that(mse(0, 0) == 0)
  assert_that(mse(c(0, 1, 2), c(1, 2, 1)) == 1)
  assert_that(mse(c(0, 1, 2), c(0, 1, -1)) == 3)
  assert_that(mse(c(-2, 1, 2), c(0, 3, 0)) == 4)
}

# Trying to use get_var_mse
get_var_mse_tests <- function() {
  # First, we need to create data frames to pass to get_var_mse:
  obs_data <- data.frame(year = 2000:2015, 
                         variable = "GMST",
                         value = 1:16)
  hect_data <- data.frame(year = 2000:2015, 
                         variable = "GMST",
                         value = seq(0, 30, 2))
  
  # Adding in another variable to our data frames
  more_data <- data.frame(year = 2000:2015,
                          variable = "CO2_concentrations",
                          value = 400:415)
  
  obs_data <- rbind(obs_data, more_data)
  hect_data <- rbind(hect_data, more_data)
  
  # Confirming get_var_mse works as intended
  assert_that(get_var_mse(obs_data, 
                          hect_data, 
                          var = "CO2_concentrations", 
                          yrs = 2000:2015) == 0)
  
  assert_that(get_var_mse(obs_data, 
                          hect_data, 
                          var = "CO2_concentrations", 
                          yrs = 2000) == 0)
  
  assert_that(get_var_mse(obs_data, 
                          hect_data, 
                          var = "GMST", 
                          yrs = 2000) == 1)
  
  assert_that(get_var_mse(obs_data, 
                          hect_data, 
                          var = "GMST", 
                          yrs = 2000:2001) == 0.5)
  
  assert_that(get_var_mse(obs_data, 
                          hect_data, 
                          var = "GMST", 
                          yrs = 2000:2004) == 3)
}

# Trying to use mse_unc
mse_unc_tests <- function() {
  assert_that(mse_unc(15, 5, 10) == 0)
  assert_that(mse_unc(15, 5, 0)  == 25)
  assert_that(mse_unc(15, 5, 20) == 25)
  assert_that(mse_unc(c(15, 20, 25),
                      c(5,  10, 15),
                      c(10, 7,  28)) == 6)
  assert_that(mse_unc(c(15, 20, 25),
                      c(5,  10, 15),
                      c(NA, 7,  28)) == 9)
}

# Trying to use nmse_unc
nmse_unc_tests <- function() {
  assert_that(nmse_unc(10, 15, 5, 10) == 0)
  assert_that(nmse_unc(10, 15, 5, 0)  == .25)
  assert_that(nmse_unc(10, 15, 5, 20) == .25)
  assert_that(nmse_unc(c(10, 15, 20),
                       c(15, 20, 25),
                       c(5,  10, 15),
                       c(10, 7,  28)) == 6 / 725)
  assert_that(nmse_unc(c(NA, 15, 20),
                       c(NA, 20, 25),
                       c(NA, 10, 15),
                       c(10, 7,  28)) == 9 / 625)
}


# Calling all testing functions
mse_tests()
get_var_mse_tests()
mse_unc_tests()
nmse_unc_tests()