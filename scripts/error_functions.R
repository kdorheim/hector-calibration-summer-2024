# Title: Error functions
# Author: Peter Scully
# Date: 6/11/24


########################
### Helper Functions ###
########################


# mse - function to find mean squared error between two vectors
#
# args: 
#   x - observed values
#   y - predicted values
#
# returns: numeric vector length 1 containing the MSE between predicted and
#          observed values
mse <- function(x, y) {
  SE <- (x - y)^2
  return(mean(SE[!is.na(SE)]))
}

# mse_unc_helper - function to find SE between an individual prediction and the
#                  confidence interval for observations
#
# arg: input_vec - vector of upper bound, lower bound, and prediction
#   
#
# returns: (upper - prediction)^2 if prediction > upper
#          (lower - prediction)^2 if prediction < lower
#          0                      otherwise
mse_unc_helper <- function(input_vec) {
  upper <- input_vec[1]
  lower <- input_vec[2]
  y     <- input_vec[3]
  
  if (is.na(y) | is.na(upper) | is.na(lower)) {
    return(NA)
  }
  
  if (y > upper) {
    return((upper - y)^2)
  } else if (y < lower) {
    return((lower - y)^2)
  }
  return(0)
}

# mse_unc - function to find MSE between predictions and confidence interval for
#           observations (SE = 0 if predictions are within the CI)
#
# args: 
#   x_upper - upper bound for observed values
#   x_lower - lower bound for observed values
#   y       - predicted values
#   x       - dummy variable added to ensure compatibility with other functions
#
# returns: numeric vector length 1 containing the MSE between predicted and
#          observed values
mse_unc <- function(x_upper, x_lower, y, x = NULL) {
  
  # Combining data into a matrix to call apply on
  # Uses helper function defined above
  input_matrix <- matrix(c(x_upper, x_lower, y), ncol = 3)
  SE <- apply(input_matrix, 1, mse_unc_helper)
  
  return(mean(SE[!is.na(SE)]))
}

# rmse - function to find root mean squared error between two vectors
#
# args: 
#   x - observed values
#   y - predicted values
#
# returns: numeric vector length 1 containing the RMSE between predicted and
#          observed values
rmse <- function(x, y) {
  return(sqrt(mse(x, y)))
}

# nmse - function to find normalized mean squared error between two vectors
#        (NMSE = MSE / x^2)
#
# args: 
#   x - observed values
#   y - predicted values
#
# returns: numeric vector length 1 containing the NMSE between predicted and
#          observed values
nmse <- function(x, y) {
  return(mse(x, y) / sum((x[!is.na(x)])^2))
}

# nmse_unc - function to find normalized mean squared error between two vectors
#            while counting anything within the given confidence interval as
#            having an SE of 0
#            (NMSE = MSE / x^2)
#
# args: 
#   x       - observed values
#   x_upper - upper bound on observed values
#   x_lower - lower bound on observed values
#   y       - predicted values
#
# returns: numeric vector length 1 containing the NMSE between predicted and
#          observed values
nmse_unc <- function(x, x_upper, x_lower, y) {
  return(mse_unc(x_upper = x_upper, x_lower = x_lower, y = y) / 
           sum((x[!is.na(x)])^2))
}

# mvsse - function to find the mean variance-standardized squared error
#
# args: 
#   x  - observed values
#   sd - standard deviation of observed values
#   y  - predicted values
#
# returns: numeric vector length 1 containing the NMSE between predicted and
#          observed values
mvsse <- function(x, sd, y) {
  vsse <- (x - y)^2 / (sd)^2
  return(mean(vsse))
}


#######################################
### Single-Variable Error Functions ###
#######################################

# get_var_mse: function to find MSE between observed and predicted data for
#              a given variable
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#   var         - variable name
#   yrs         - vector of years for finding MSE
#   mse_fn      - function to calculate mse with (default = mse)
#
# Returns: MSE between predicted and observed data for var
get_var_mse <- function(obs_data, hector_data, var, yrs, mse_fn = mse) {
  obs_vec <- 
    filter(obs_data, year %in% yrs & variable == var)$value
  
  hector_vec <-
    filter(hector_data, year %in% yrs & variable == var)$value
  
  return(mse_fn(obs_vec, hector_vec))
}


# get_var_mse_unc: function to find NMSE between observed and predicted data
#                  for a given variable while accounting for uncertainty in
#                  observations
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#   var         - variable name
#   yrs         - vector of years for finding MSE
#   mse_fn      - function to calculate mse with 
#                 (function must take x, y, AND upper/lower bounds)
#
# Returns: MSE between predicted and observed data for var
get_var_mse_unc <- function(obs_data, hector_data, var, yrs, mse_fn) {
  x       <- filter(obs_data, year %in% yrs & variable == var)$value
  x_upper <- filter(obs_data, year %in% yrs & variable == var)$upper
  x_lower <- filter(obs_data, year %in% yrs & variable == var)$lower
  y       <- filter(hector_data, year %in% yrs & variable == var)$value
  
  return(mse_fn(x = x, x_upper = x_upper, x_lower = x_lower, y = y))
}


# get_var_mvsse: function to find MVSSE between observed and predicted data
#                  for a given variable
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#   var         - variable name
#   yrs         - vector of years for finding MSE
#
# Returns: MSE between predicted and observed data for var
#
# Note: Assumes observed data contains symmetric upper and lower bounds 1 SD 
#       away from actual value
get_var_mvsse <- function(obs_data, hector_data, var, yrs, mse_fn) {
  x       <- filter(obs_data, year %in% yrs & variable == var)$value
  x_upper <- filter(obs_data, year %in% yrs & variable == var)$upper
  sd      <- x_upper - x
  y       <- filter(hector_data, year %in% yrs & variable == var)$value
  
  return(mvsse(x = x, sd = sd, y = y))
}


#######################################
### Error Functions for Optim Usage ###
#######################################


# mean_T_CO2_mse: function to find the mean of the temperature and CO2 MSEs 
#                 between observed and predicted data for a given variable
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: MSE between predicted and observed data for var
mean_T_CO2_mse <- function(obs_data, hector_data) {
  T_mse <- get_var_mse(obs_data = obs_data, 
                       hector_data = hector_data, 
                       var = GMST(), 
                       yrs = 1850:2014)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014))
  
  return(mean(c(T_mse, CO2_mse)))
}

# mean_T_CO2_nmse: function to find the mean of the temperature and CO2 NMSEs 
#                  between observed and predicted data
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: Average NMSE between predicted and observed data
mean_T_CO2_nmse <- function(obs_data, hector_data) {
  T_mse <- get_var_mse(obs_data = obs_data, 
                       hector_data = hector_data, 
                       var = GMST(), 
                       yrs = 1850:2014,
                       mse_fn = nmse)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014),
                         mse_fn = nmse)  
  return(mean(c(T_mse, CO2_mse)))
}

# mean_T_CO2_nmse_unc: function to find the mean of the temperature and CO2 
#                      NMSEs between observed and predicted data while
#                      accounting for temperature uncertainty
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: Average NMSE between predicted and observed data
mean_T_CO2_nmse_unc <- function(obs_data, hector_data) {
  T_mse <- get_var_mse_unc(obs_data = obs_data, 
                           hector_data = hector_data, 
                           var = GMST(), 
                           yrs = 1850:2014,
                           mse_fn = nmse_unc)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014),
                         mse_fn = nmse)  
  return(mean(c(T_mse, CO2_mse)))
}

# mean_T_CO2_OHC_nmse_unc: function to find the mean of the temperature, CO2, 
#                          and OHC NMSEs between observed and predicted data
#                          while accounting for temperature and OHC uncertainty
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: Average NMSE between predicted and observed data
mean_T_CO2_OHC_nmse_unc <- function(obs_data, hector_data) {
  T_mse <- get_var_mse_unc(obs_data = obs_data, 
                           hector_data = hector_data, 
                           var = GMST(), 
                           yrs = 1850:2014,
                           mse_fn = nmse_unc)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014),
                         mse_fn = nmse)  
  OHC_mse <- get_var_mse_unc(obs_data = obs_data,
                             hector_data = hector_data,
                             var = "OHC",
                             yrs = 1957:2014,
                             mse_fn = nmse_unc)
  return(mean(c(T_mse, CO2_mse, OHC_mse)))
}

# smooth_T_CO2_mse: function to find the mean of smoothed temperature and CO2
#                   MSEs between observed & predicted data for a given variable
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: MSE between predicted and observed data for var
smooth_T_CO2_mse <- function(obs_data, hector_data) {
  
  # Getting data frame with just smoothed data that works with get_var_mse
  smooth_data <- filter(obs_data, variable == "Smooth T")
  smooth_data$variable <- GMST()
  
  # Getting  mses
  T_mse <- get_var_mse(obs_data = smooth_data, 
                       hector_data = hector_data, 
                       var = GMST(), 
                       yrs = 1850:2014)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014))  
  return(mean(c(T_mse, CO2_mse)))
}


# smooth_T_CO2_nmse: function to find mean of smoothed temperature and CO2 NMSEs
#                    between observed & predicted data
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: NMSE between predicted and observed data for var
smooth_T_CO2_nmse <- function(obs_data, hector_data) {
  
  # Getting data frame with just smoothed data that works with get_var_mse
  smooth_data <- filter(obs_data, variable == "Smooth T")
  smooth_data$variable <- GMST()
  
  # Getting  mses
  T_mse <- get_var_mse(obs_data = smooth_data, 
                       hector_data = hector_data, 
                       var = GMST(), 
                       yrs = 1850:2014,
                       mse_fn = nmse)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014),
                         mse_fn = nmse)  
  return(mean(c(T_mse, CO2_mse)))
}


# smooth_T_CO2_nmse_unc: function to find mean of smoothed temperature and CO2 
#                        NMSEs between observed & predicted data while
#                        accounting for uncertainty in temperature
#
# args: 
#   obs_data    - data frame of observed data formatted like Hector data frame
#   hector_data - data frame outputted by Hector
#
# Returns: NMSE between predicted and observed data
smooth_T_CO2_nmse_unc <- function(obs_data, hector_data) {
  
  # Getting data frame with just smoothed data that works with get_var_mse
  smooth_data <- filter(obs_data, variable == "Smooth T")
  smooth_data$variable <- GMST()
  
  # Getting  mses
  T_mse <- get_var_mse_unc(obs_data = smooth_data, 
                           hector_data = hector_data, 
                           var = GMST(), 
                           yrs = 1850:2014,
                           mse_fn = nmse_unc)
  CO2_mse <- get_var_mse(obs_data = obs_data, 
                         hector_data = hector_data, 
                         var = CONCENTRATIONS_CO2(), 
                         yrs = c(1750, 1850:2014),
                         mse_fn = nmse)  
  return(mean(c(T_mse, CO2_mse)))
}