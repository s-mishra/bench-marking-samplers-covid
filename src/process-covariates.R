process_covariates <- function(regions,
                               df_data,
                               ifr_by_regions,
                               num_days_sim,
                               formula_full,
                               formula_macro,
                               formula_micro,
                               create_features = FALSE) {
  # Read in serial interval
  serial_interval <- read_csv(paste(here(),
                                    'data/serial-interval.csv',
                                    sep = "/"))
  serial_interval <-
    c(serial_interval$fit, rep(1e-17, num_days_sim - 100))
  # infection to onset
  mean1 <- 5.1
  cv1 <- 0.86
  # onset to death
  mean2 <- 17.8
  cv2 <- 0.45
  x1 <- rgammaAlt(1e6, mean1, cv1) # infection-to-onset distribution
  x2 <- rgammaAlt(1e6, mean2, cv2) # onset-to-death distribution
  f_cached <-
    ecdf(x1 + x2) # empirical cumulative distribtion function
  
  dates <- list()
  reported_cases <- list()
  reported_deaths <- list()
  stan_data <- list(
    M = length(regions), # Number of regions
    N0 = 6, # Number of days in seeding
    N = NULL,  # Number of time points with data
    N2 = NULL,  # Number of time points with data plus forecast
    cases = NULL,  # daily cases
    deaths =  NULL, # daily deaths
    f = NULL, # Hazard times survival
    X = NULL,  # Covariates for full pool
    X_partial_macro = NULL, # Covariates for partial pool macro
    X_partial_micro = NULL, # Covariates for partial pool micro
    P = NULL, # Number of covariates for full pool
    P_partial_macro = NULL, # Number of covariates for partial pool macro
    P_partial_micro = NULL, # Number of covariates for partial pool micro
    SI = serial_interval,  # Serial interval fit
    EpidemicStart = NULL,  # Date to start epidemic in each state
    pop = NULL,  #  region population
    Q = NULL,  # Number of macro regions
    Region = NULL,  # Macro region index for each state
    W = NULL, # Number of weeks for weekly effects
    week_index = NULL # Week index for each state
  )
  
  covariate_list <- list()
  covariate_list_macro <- list()
  covariate_list_micro <- list()
  k = 1
  for (region in regions) {
    print(region)
    # Selects correct IFR for each country
    ifr_region <- ifr_by_regions$ifr[ifr_by_regions$country == region]
    # Subsets data by state
    data_region <- filter(df_data, country_name == region)
    # Sorts data into time order
    data_region <- arrange(data_region, Date)
    # Selects index when first case
    index <- which(data_region$Cases > 0)[1]
    index1 <-  which(cumsum(data_region$Deaths) >= 10)[1]
    index2 <- index1 - 30 # was 30 days before 10 cumulative deaths
    print(
      sprintf(
        "First non-zero cases is on day %d, and 30 days before 10 deaths is day %d",
        index,
        index2
      )
    )
    
    # Works out state start point
    data_region <- data_region[index2:nrow(data_region), ]
    
    # Work out how much padding need at end
    N <- length(data_region$Cases)
    forecast_length <- num_days_sim - N[[1]]
    print(sprintf(
      "%s has %d days of data. Forcasting %d days",
      region,
      N,
      forecast_length
    ))
    # IFR is the overall probability of dying given infection
    convolution = function(u)
      (ifr_region * f_cached(u))
    
    f = rep(0, num_days_sim) # f is the probability of dying on day i given infection
    f[1] = (convolution(1.5) - convolution(0))
    for (i in 2:num_days_sim) {
      f[i] = (convolution(i + .5) - convolution(i - .5))
    }
    
    deaths <- c(data_region$Deaths, rep(-1, forecast_length))
    cases <- c(data_region$Cases, rep(-1, forecast_length))
    
    covariates <- select(data_region, -c(Date, country_name, Cases, Deaths))
    # padding last row for simulated days
    covariates <- bind_rows(covariates, covariates %>% slice(rep(n(), each = num_days_sim-N)))
    if(create_features){
      features <- model.matrix(formula_full, covariates)
      features_macro <- model.matrix(formula_macro, covariates)
      features_micro <- model.matrix(formula_micro, covariates)
    } else{
      features <- model.matrix(as.formula("~-1"), covariates)
      features_macro <- model.matrix(as.formula("~-1"), covariates)
      features_micro <- model.matrix(as.formula("~-1"), covariates)
    }
    # padding features
    covariate_list[[k]] <- features
    covariate_list_macro[[k]] <- features_macro
    covariate_list_micro[[k]] <- features_micro
    k <- k + 1
    
    ## Append data to stan data
    stan_data$EpidemicStart <-
      c(stan_data$EpidemicStart, 31) # 30 days before 10 cumulative deaths
    stan_data$pop <- c(stan_data$pop, ifr_by_regions$pop[ifr_by_regions$country == region])
    stan_data$f <- cbind(stan_data$f, f)
    stan_data$deaths <- cbind(stan_data$deaths, deaths)
    stan_data$cases <- cbind(stan_data$cases, cases)
    stan_data$N2 <- num_days_sim
    stan_data$N <- c(stan_data$N, N)
    stan_data$region <-
      c(stan_data$region, ifr_by_regions$macro.region[ifr_by_regions$country ==
                                                        region])
    
    # Saves other data for each state
    dates[[region]] <- data_region$Date
    reported_cases[[region]] <- data_region$Cases
    reported_deaths[[region]] <- data_region$Deaths
  }
  
  stan_data$P = dim(features)[2]
  stan_data$P_partial_macro = dim(features_macro)[2]
  stan_data$P_partial_micro = dim(features_micro)[2]
  if (stan_data$P == 0) {
    stan_data$X = array(0, dim = c(stan_data$M , stan_data$N2, 1))
  } else{
    stan_data$X = array(NA, dim = c(stan_data$M , stan_data$N2 , stan_data$P))
  }
  if (stan_data$P_partial_macro == 0) {
    stan_data$X_partial_macro = array(0, dim = c(stan_data$M , stan_data$N2, 1))
  } else{
    stan_data$X_partial_macro = array(NA,
                                      dim = c(stan_data$M , stan_data$N2 , stan_data$P_partial_macro))
  }
  if (stan_data$P_partial_micro == 0) {
    stan_data$X_partial_micro = array(0, dim = c(stan_data$M , stan_data$N2, 1))
  } else{
    stan_data$X_partial_micro = array(NA,
                                      dim = c(stan_data$M , stan_data$N2 , stan_data$P_partial_micro))
  }
  
  for (i in 1:stan_data$M) {
    if (stan_data$P != 0)
      stan_data$X[i, , ] = covariate_list[[i]]
    if (stan_data$P_partial_macro != 0)
      stan_data$X_partial_macro[i, , ] = covariate_list_macro[[i]]
    if (stan_data$P_partial_micro != 0)
      stan_data$X_partial_micro[i, , ] = covariate_list_micro[[i]]
  }
  if (stan_data$P == 0)
    stan_data$P = 1
  if (stan_data$P_partial_macro == 0)
    stan_data$P_partial_macro = 1
  if (stan_data$P_partial_micro == 0)
    stan_data$P_partial_micro = 1
  stan_data$Q <- max(stan_data$region)
  # Special case when we don't have macro regions
  if (is.na(stan_data$Q) || stan_data$Q == "") {
    stan_data$Q <- 1
    stan_data$P_partial_macro = 1
    stan_data$X_partial_macro = array(0, dim = c(stan_data$M , stan_data$N2, 1))
    stan_data$Region = rep(1, stan_data$M)
  }
  stan_data$W <- ceiling(stan_data$N2 / 7)
  stan_data$week_index <- matrix(1, stan_data$M, stan_data$N2)
  for (state.i in 1:stan_data$M) {
    stan_data$week_index[state.i, ] <-
      rep(2:(stan_data$W + 1), each = 7)[1:stan_data$N2]
    last_ar_week = which(dates[[state.i]] == max(df_data$Date) - 21)
    stan_data$week_index[state.i, last_ar_week:ncol(stan_data$week_index)] <-
      stan_data$week_index[state.i, last_ar_week]
  }
  if(length(stan_data$N) == 1) {
    stan_data$N = as.array(stan_data$N)
    stan_data$pop = as.array(stan_data$pop)
    stan_data$EpidemicStart = as.array(stan_data$EpidemicStart)
    stan_data$Region = as.array(stan_data$Region)
  }
  return(
    list(
      "stan_data" = stan_data,
      "dates" = dates,
      "reported_cases" = reported_cases,
      "reported_deaths" = reported_deaths
    )
  )
}
