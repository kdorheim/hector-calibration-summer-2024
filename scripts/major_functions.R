# Title: Major functions
# Author: Peter Scully
# Date: 6/11/24


### Imports and Constants ###
library(hector)
library(dplyr)
source(here::here("scripts", "error_functions.R"))


# File Paths
IDEAL_RUNS <- here::here("idealized_inputs")
INPUT_TCRE <- system.file("input/hector_ssp585.ini", package = "hector")
INPUT_TCR <- file.path(IDEAL_RUNS, "hector_1pctCO2.ini")

# Ocean heat content constants
OCEAN_AREA <- 5100656e8 * (1 - 0.29) # The total area of the ocean
W_TO_ZJ <- 3.155693e-14              # Watts to ZJ

# Other Numeric Constants
CO2_OBS_SD <- 0.12

# Constants for TCRE and ERF calculations
EMISSIONS_VARS <- c(FFI_EMISSIONS(), LUC_EMISSIONS())
UPTAKE_VARS <- c(LUC_UPTAKE(), DACCS_UPTAKE())
TCRE_VARS <- c(GLOBAL_TAS(),
               FFI_EMISSIONS(), LUC_EMISSIONS(),
               LUC_UPTAKE(), DACCS_UPTAKE())

AEROSOL_RF_VARS <- c(RF_BC(), RF_OC(), RF_NH3(), RF_SO2(), RF_ACI())
NONGHG_RF_VARS <- c(AEROSOL_RF_VARS, RF_VOL(), RF_ALBEDO(), RF_MISC())

################################
### Comparison Data Fetchers ###
################################

# get_co2_data - function to get the CO2 data from the appropriate file
#
# args:
#   file - path to historical CO2 data file
#   scenario - name of scenario being run (default: "historical")
#   include_unc - boolean indicating whether to include upper/lower bound cols
#                 (default: FALSE)
#
# returns: Hector-style data frame with CO2 data
#
# Note: The current data set does not support including uncertainty values, so
#       if uncertainty is included, upper/lower bound columns both equal the
#       observed value
get_co2_data <- function(file, scenario = "historical", include_unc = F) {
  
  # Reading in only CO2 data
  co2_data <- read.table(file, 
                         skip = 23, 
                         sep = ",",
                         colClasses = c("numeric", "numeric", "NULL", "NULL",
                                        "NULL", "NULL", "NULL", "NULL"))
  
  # Fixing table formatting
  co2_data <- na.omit(co2_data)
  colnames(co2_data) <- c("year", "value")
  
  # Adding in new columns to match Hector data frames
  co2_data$scenario <- scenario
  co2_data$variable <- CONCENTRATIONS_CO2()
  co2_data$units <- " ppmv CO2"
  
  if (include_unc) {
    co2_data$upper <- co2_data$value + CO2_OBS_SD
    co2_data$lower <- co2_data$value - CO2_OBS_SD
  }
  
  return(co2_data)
}

# get_temp_data - function to get the temperature data from the appropriate file
#
# args:
#   file - path to historical temperature data file
#   scenario - name of scenario being run (default: "historical")
#   include_unc - boolean indicating whether to include upper/lower bounds on
#                 values (default: FALSE)
#
# returns: Hector-style data frame with temperature data
get_temp_data <- function(file, scenario = "historical", include_unc = F) {
  temp_data <- read.csv(file)
  colnames(temp_data) <- c("year", "value", "lower", "upper")
  
  # Adding in new columns to match Hector data frames
  temp_data$scenario <- scenario
  temp_data$variable <- GMST()
  temp_data$units <- "degC"
  
  # Removing upper and lower bounds, if desired
  if (!include_unc) {
    temp_data$lower <- NULL
    temp_data$upper <- NULL
  }
  
  return(temp_data)
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


############################
### Data Frame Functions ###
############################


# get_var_change - function to find the change in a variable between a start
#                  and end date
#
# args:
#   data  - data frame outputted by fetchvars
#   var   - string containing the Hector variable to find the change in
#   start - start year
#   end   - end year
#
# returns: single value indicating the change in variable from start year to
#          end year
get_var_change <- function(data, var, start, end) {
  initial_val <- filter(data, variable == var & year == start)$value
  final_val   <- filter(data, variable == var & year == end)$value
  return(final_val - initial_val)
}


# get_interval_avg - function to find the average value of a Hector variable
#                    over a given time interval
#
# args:
#   data  - data frame outputted by run_hector (or fetchvars)
#   var   - string containing the Hector variable to average
#   start - start year for finding average
#   end   - end year for finding average
#
# returns: average value of var across the provided interval
get_interval_avg <- function(data, var, start, end) {
  data %>%
    subset(variable == var) %>%
    subset(start <= year & year <= end) %>%
    .$value -> interval_vals
  return(mean(interval_vals))
}


# rel_to_val - function to adjust values of a Hector variable to be relative
#              to a new value
#
# args:
#   data      - data frame outputted by fetchvars
#   var       - string containing the Hector variable to adjust
#   benchmark - value to make variable relative to
#
# returns: data frame containing the new values for var relative to benchmark
rel_to_val <- function(data, var, benchmark) {
  # Updating all of the values for var that are in data
  data %>%
    mutate(value = ifelse(variable == var, value - benchmark, value)) -> data
  return(data)
}


# rel_to_interval - function to normalize values of a Hector variable to a
#                   reference period
#
# args:
#   data  - data frame outputted by fetchvars
#   var   - string containing the Hector variable to adjust
#   start - start year for reference period
#   end   - end year for reference period
#
# returns: data frame containing the normalized values for that variable
rel_to_interval <- function(data, var, start, end) {
  benchmark <- get_interval_avg(data, var, start, end)
  return(rel_to_val(data, var, benchmark))
}


# sum_vars - function to sum the values of several Hector variables at each
#            timestep. adds new rows to the provided data frame containing the
#            summed values at each time step
#
# args:
#   data - data frame outputted by fetchvars
#   vars - vector of Hector variables to sum together
#   name - name to call the new variable containing the sums of the values in
#           vars
#   unit - units to assign to the new variable
#   yrs  - vector of years to sum these variables for
#
# returns: original data frame with additional rows containing the sums of the
#          values for the inputted variables
#
# Note: Currently, all columns other than year, variable, and value will be set
#       equal to the same values as the last row in the original data frame
# Note: This function is currently inefficient and should ideally be reworked to
#       take advantage of dplyr functions like group_by
sum_vars <- function(data, vars, name, unit, yrs) {
  
  # Setting up iterator to iterate through new rows of data frame
  curr_row <- nrow(data) + 1
  for (yr in yrs) {
    # Adding new row to data frame
    data[curr_row,] <- data[curr_row - 1,]
    data[curr_row,]$year <- yr
    data[curr_row,]$variable <- name
    data[curr_row,]$units <- unit
    
    # Summing the values across all variables
    data[curr_row,]$value <-
      sum(filter(data, year == yr & variable %in% vars)$value)
    curr_row <- curr_row + 1
  }
  return(data)
}


###############################
### Optim-Related Functions ###
###############################

# run_hector - function to run Hector using a given ini file and set of params
#
# args:
#   ini_file    - path to the ini file to use to run Hector
#   params      - vector of Hector params to modify. if NULL, does default run
#   vals        - vector of values to use for those Hector parameters
#   yrs         - year range to get Hector data from
#   vars        - Hector variables to get data on
#   include_unc - boolean indicating whether to include upper/lower bounds on
#                 values (default: FALSE)
#
# returns: data frame containing the variables' values for the given date range
#
# note: assumes params and vals are the same length and that each paraneter is
#       stored at the same index as its corresponding value
run_hector <- function(ini_file, params, vals, yrs, vars, include_unc = F) {
  core <- newcore(ini_file)
  
  # Setting parameter values
  if (!is.null(params)) {
    for (i in 1:length(params)) {
      setvar(core = core, 
             dates = NA, 
             var = params[i], 
             values = vals[i], 
             unit = getunits(params[i]))
    }
  }
  reset(core)
  
  # Running core and fetching data
  run(core)
  data <- fetchvars(core, yrs, vars = vars)
  shutdown(core)
  
  # Rescaling temperatures (if applicable)
  if (GMST() %in% vars) {
    data <- rel_to_interval(data = data, var = GMST(), start = 1961, end = 1990)
  }
  
  # Measuring ocean heat content (if applicable)
  if (HEAT_FLUX() %in% vars) {
    
    ### Converting heat fluxes to cumulative heat content
    
    # Getting heat flux data
    ohc_data <- filter(data, variable == HEAT_FLUX() & year %in% 1957:2014)
    
    # Converting heat flux to OHC change by year
    ohc_data$value <- ohc_data$value * OCEAN_AREA * W_TO_ZJ
    
    # Converting OHC changes by year to total OHC (relative to 1957)
    ohc_data$value <- cumsum(ohc_data$value)
    ohc_data$variable <- "OHC"
    
    # Making OHC relative to 2005-2014 average
    ohc_data <- rel_to_interval(data = ohc_data,
                                 var = "OHC",
                                 start = 2005,
                                 end = 2014)
    
    # Adding ohc data to Hector data frame
    data <- rbind(data, ohc_data)
  }
  
  # Adding in upper and lower bounds (if applicable)
  if (include_unc) {
    data$upper <- data$value
    data$lower <- data$value
  }
  
  return(data)
}


# write_metric - writes the value for a metric to the output file
#
# args:
#   name - name of metric
#   val  - value of metric
#
# returns: nothing
#
# Note: The value of the metric is rounded to 3 sig figs when outputted
write_metric <- function(name, val, file) {
  write(paste(name, signif(val, 3)), file = file, append = TRUE)
}


# objective_fn - flexible objective function that can be used to calculate the
#                error between observed and predicted data for a Hector run with
#                the given values for a set of parameters
#
# args:
#   obs_data    - Hector-style data frame with observed data
#   ini_file    - path to the ini file to use to run Hector
#   params      - vector of Hector parameters to modify
#   pars        - vector of values to use for those Hector parameters
#   yrs         - year range to get Hector data from
#   vars        - Hector variables to get data on
#   error_fn    - function to calculate error between observed/predicted vals
#   include_unc - boolean indicating whether upper/lower bounds needed for 
#                 error_fn (default: FALSE)
#
# returns: error between observed and predicted data calculated via error_fn
#
# note: uses an error function from error_functions.R
# note: assumes error_fn takes obs_data and hector_data as inputs and outputs a
#       real-valued error
objective_fn <- function(obs_data, ini_file, params, par, yrs, vars, error_fn, 
                         include_unc = F) {
  # Running Hector
  hector_data <- run_hector(ini_file = ini_file,
                            params = params,
                            vals = par,
                            yrs = yrs,
                            vars = vars,
                            include_unc = include_unc)
  
  # Returning error from model run
  return(error_fn(obs_data, hector_data))
}

# run_optim - function to run optim with an objective function based on a given
#             set of Hector parameters and an error function. Also outputs 
#             results of optim call
#
# args:
#   obj_fn      - objective function to use
#   obs_data    - Hector-style data frame with observed data
#   ini_file    - path to the ini file to use to run Hector
#   params      - vector of Hector parameters to modify
#   par         - vector of initial values for params. Default is NULL, meaning
#                 that Hector default parameters will be used as initial values
#   sd          - vector of standard derivs for par. Only used if using LBFGS-B.
#                 If NULL, use provided lower and upper bounds
#   lower       - lower bound for using LBFGS-B
#   upper       - upper bound for using LBFGS-B
#   yrs         - year range to get Hector data from
#   vars        - Hector variables to get data on
#   error_fn    - function to calculate error between observed/predicted vals
#   method      - optimization method to use. Default: Nelder-Mead
#   include_unc - boolean indicating whether upper/lower bounds needed for 
#                 error_fn (default: FALSE)
#   output_file - path to file to append table to
#
# returns: Optimal parameters. Also outputs optimal parameters and objective 
#          function value to given output file
#
# note: uses an error function from error_functions.R
run_optim <- function(obs_data, ini_file, params, par = NULL, sd = NULL, 
                      lower = NULL, upper = NULL, yrs, vars, error_fn, 
                      include_unc = F, method = "Nelder-Mead", output_file) {
  # Creating vector of default parameters
  if (is.null(par)) {
    default_core <- newcore(ini_file)
    par <- fetchvars(default_core, dates = NA, vars = params)$value
    shutdown(default_core)
  }
  
  # Applying optim
  if (method == "Nelder-Mead") {
    optim_output <- optim(par = par, 
                          fn = objective_fn, 
                          obs_data = obs_data, 
                          ini_file = ini_file, 
                          params = params,
                          yrs = yrs,
                          vars = vars, 
                          error_fn = error_fn, 
                          include_unc = include_unc)
  } else if (method == "L-BFGS-B") {
    if (!is.null(sd)) {
      lower <- par - sd
      upper <- par + sd
    }
    optim_output <- optim(par = par, 
                          fn = objective_fn, 
                          obs_data = obs_data, 
                          ini_file = ini_file, 
                          params = params,
                          yrs = yrs,
                          vars = vars, 
                          error_fn = error_fn, 
                          include_unc = include_unc,
                          method = method,
                          lower = lower,
                          upper = upper)
  }
  
  # Extract vals from optim output
  best_pars <- optim_output$par
  min_error <- optim_output$value
  
  # Output optimal parameterization
  write.table(data.frame(parameters = params, values = best_pars),
              file = output_file,
              append = TRUE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
  write("", output_file, append = TRUE)
  
  # Output value of objective function
  write_metric("Objective Function Value:", min_error, output_file)
  
  # Output counts
  write_metric("Counts:", optim_output$counts, output_file)
  
  # Output convergence
  write_metric("Convergence: ", optim_output$convergence, output_file)
  write(paste("Messages:", optim_output$message), output_file, append = TRUE)
  write("", output_file, append = TRUE)
  
  return(best_pars)
}

##############################
### OTHER OUTPUT FUNCTIONS ###
##############################

# calc_table_metrics - function to calculate the key values, historical metrics,
#                      and future predictions for a given set of Hector params
#
# args:
#   params      - vector of Hector params to modify. if NULL, does default run
#   vals        - vector of values to use for those Hector parameters
#   output_file - path to file to append table to
#
# returns: Nothing, but outputs table to given output file
calc_table_metrics <- function(params, vals, output_file) {
  
  ### KEY METRICS ###
  
  ## Finding TCRE ##
  
  tcre_data <- run_hector(ini_file = INPUT_TCRE,
                          params = params,
                          vals = vals,
                          yrs = 1750:2300,
                          vars = TCRE_VARS)
  
  # Combining emissions data
  tcre_data %>%
    subset(variable %in% EMISSIONS_VARS) %>%
    group_by(year) %>%
    summarize(all_emissions = sum(value)) %>%
    .$all_emissions -> tcre_emissions
  
  total_emissions <- cumsum(tcre_emissions)
  
  tcre_data %>%
    subset(variable %in% UPTAKE_VARS) %>%
    group_by(year) %>%
    summarize(all_uptake = sum(value)) %>%
    .$all_uptake -> tcre_uptake
  
  total_uptake <- cumsum(tcre_uptake)
  
  net_emissions <- total_emissions - total_uptake
  
  # Getting temperature data
  tcre_temps <- filter(tcre_data, variable == GLOBAL_TAS())$value
  
  # Calculating TCRE via linear regression
  tcre_reg <- lm(tcre_temps ~ net_emissions)
  tcre <- tcre_reg$coefficients[2] * 1000
  
  
  ## Finding TCR ##
  
  # Running Hector and extracting relevant data
  tcr_data <- run_hector(ini_file = INPUT_TCR,
                         params   = params,
                         vals     = vals,
                         yrs      = 1800:2000,
                         vars     = c(GLOBAL_TAS(), CONCENTRATIONS_CO2()))
  
  tcr_co2   <- filter(tcr_data, variable == CONCENTRATIONS_CO2())$value
  tcr_temps <- filter(tcr_data, variable == GLOBAL_TAS())$value
  
  # Finding TCR from linear fit
  tcr_reg <- lm(tcr_temps ~ tcr_co2)
  
  # Getting initial CO2
  core <- newcore(INPUT_TCR)
  initial_co2 <- fetchvars(core, dates = NA, vars = PREINDUSTRIAL_CO2())$value
  shutdown(core)
  
  tcr <- tcr_reg$coefficients[1] + tcr_reg$coefficients[2] * initial_co2 * 2
  
  
  ### HISTORICAL WARMING AND ERF ###
  
  ## Initial Hector Run ##
  
  # Setting up variables
  hist_yrs <- 1750:2019
  hist_vars <- c(GLOBAL_TAS(),               # For GSAT warming
                 HEAT_FLUX(),                # For ocean heat content change
                 RF_CH4(),                   # For Methane ERF
                 NONGHG_RF_VARS, RF_TOTAL()  # For other ERF calcs
  )
  hist_file <- system.file("input/hector_ssp245.ini", package = "hector")
  
  # Running Hector
  hist_data <- run_hector(ini_file = hist_file,
                          params = params,
                          vals = vals, 
                          yrs = hist_yrs, 
                          vars = hist_vars)
  
  
  ## Finding GSAT Warming ##
  
  # Normalizing temperatures
  hist_data <- rel_to_interval(data  = hist_data,
                               var   = GLOBAL_TAS(),
                               start = 1850,
                               end   = 1900)
  
  # Getting 1995-2014 avg
  GSAT_warming <- get_interval_avg(data  = hist_data,
                                   var   = GLOBAL_TAS(),
                                   start = 1995,
                                   end   = 2014)
  
  
  ## Finding ocean heat content change ##
  
  # Summing ocean heat (by taking average and multiplying by number of yrs)
  avg_flux <- get_interval_avg(data  = hist_data,
                               var   = HEAT_FLUX(),
                               start = 1971,
                               end   = 2018)
  total_flux <- avg_flux * length(1971:2018)
  
  # Converting flux to heat content change
  ocean_heat_content_change <- total_flux * OCEAN_AREA * W_TO_ZJ
  
  
  ###Finding total aerosol ERF ##
  
  # Getting total aerosol forcing for relevant years
  hist_data <- sum_vars(data = hist_data,
                        vars = AEROSOL_RF_VARS,
                        name = "tot_aer_ERF",
                        unit = "w/m^2",
                        yrs  = 2005:2015)
  
  # Getting average forcing from 2005-2015
  tot_aer_erf <- get_interval_avg(data  = hist_data,
                                  var   = "tot_aer_ERF",
                                  start = 2005,
                                  end   = 2015)
  
  ## Finding WMGHG ERF ##
  
  # Getting total non-GHG ERF
  hist_data <- sum_vars(data = hist_data,
                        vars = NONGHG_RF_VARS,
                        name = "non_ghg_ERF",
                        unit = "W/m^2",
                        yrs  = 2019)
  
  # Subtracting non-GHG ERF from total RF
  wmghg_erf <- filter(hist_data, variable == RF_TOTAL() & year == 2019)$value -
    filter(hist_data, variable == "non_ghg_ERF" & year == 2019)$value
  
  
  ## Finding methane ERF ##
  methane_erf <- filter(hist_data, variable == RF_CH4() & year == 2019)$value
  
  
  ### FUTURE WARMING ###
  
  # Setting up variables
  future_yrs <- 1995:2100
  future_vars <- GLOBAL_TAS()
  
  # Getting the names of each scenario file
  scenarios <- c("119", "126", "245", "370", "585")
  scenario_names <- paste("input/hector_ssp", scenarios, ".ini", sep = "")
  scenario_files <- system.file(scenario_names, package = "hector")
  
  # Setting up results data frame
  future_results <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(future_results) <- c("scenario", "start", "end", "GSAT")
  
  for (scen_counter in 1:length(scenario_files)) {
    # Getting data
    scen_data <- run_hector(ini_file = scenario_files[scen_counter],
                            params = params,
                            vals = vals,
                            yrs = future_yrs,
                            vars = future_vars)
    
    # Normalizing temperatures
    scen_data <- rel_to_interval(data = scen_data,
                                 var = GLOBAL_TAS(),
                                 start = 1995,
                                 end = 2014)
    
    # Finding the averages for each time interval
    start_yrs <- c(2021, 2041, 2081)
    end_yrs <- start_yrs + 19
    
    # Getting the necessary information for each time interval
    for (i in 1:length(start_yrs)) {
      
      # Finding average temperatures for the given time interval
      temps <- get_interval_avg(data = scen_data,
                                var = GLOBAL_TAS(),
                                start = start_yrs[i],
                                end = end_yrs[i])
      
      # Adding the scenario, start and end dates, and avg temperature to the
      # data frame
      future_results[nrow(future_results) + 1, ] <-
        list(paste("ssp", scenarios[scen_counter], sep = ""),
             start_yrs[i],
             end_yrs[i],
             signif(temps, 3))
    }
  }
  
  
  ### OUTPUTTING RESULTS ###
  write("***Key Metrics***", file = output_file, append = TRUE)
  write_metric("TCRE:", tcre, output_file)
  write_metric("TCR: ", tcr, output_file)
  write("", file = output_file, append = TRUE)
  
  write("***Historical Warming and ERF***", file = output_file, append = TRUE)
  write_metric("GSAT Warming:             ", GSAT_warming, output_file)
  write_metric("Ocean Heat Content Change:", ocean_heat_content_change, 
               output_file)
  write_metric("Total Aerosol ERF:        ", tot_aer_erf, output_file)
  write_metric("WMGHG ERF:                ", wmghg_erf, output_file)
  write_metric("Methane ERF:              ", methane_erf, output_file)
  write("", file = output_file, append = TRUE)
  
  write("***Future Warming***", file = output_file, append = TRUE)
  write.table(future_results,
              file = output_file,
              append = TRUE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
}