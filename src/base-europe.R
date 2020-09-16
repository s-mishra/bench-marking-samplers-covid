library(rstan)
library(tidyverse)
library(EnvStats)
library(scales)
library(stringr)
library(abind)
library(optparse)
library(utils)
library(httr)
library(zoo)
library(here)
library(ggplot2)
library(bayesplot)
library(tidybayes)

source(paste(here(),'src/process-covariates.R', sep='/'))
source(paste(here(),'src/plot-3-pannel.R', sep='/'))
source(paste(here(),"src/covariates-effect-sizes.R", sep='/'))
source(paste(here(),"src/get-data.R", sep='/'))
# read data from delve, uncomment next line if you want new new data
# get_europe_data()
df_data <- readRDS(paste(here(), 'data/europe-data.rds', sep = '/'))

# making sure that we are reading only required data
df_data <- 
  df_data %>%
  select(DATE, country_name, starts_with('npi'), cases_new, deaths_new)
# Commandline options and parsing
option_list <- list(
  make_option(c("--debug", "-D"),action="store_true", default=FALSE,help="Whether running in DEBUG MODE [default \"%default\"]"),
  make_option(c("--full", "-F"),action="store_true", default=FALSE,help="Whether running in FULL MODE [default \"%default\"]. If both Debug and Full are passed code will run in FULL mode."),
  make_option(c("--stanFile"),action="store", default="base",help="Which Stan file to use [default \"%default\"]"),
  make_option(c("--nchains"),action="store", type="integer", default=4,help="Number of Chains [default \"%default\"]"),
  make_option(c("--iter"),action="store", type="integer", default=1500,help="Number of iterations [default \"%default\"]"),
  make_option(c("--thin"),action="store", type="integer", default=1,help="Amount of thinning of results [default \"%default\"]"),
  make_option(c("--formula_full"),action="store", default="~ -1 + Stringency",  help="Features to be used for full pooling [default \"%default\"]"),
  make_option(c("--formula_macro"),action="store", default="~ -1",  help="Features to be used for macro level pooling [default \"%default\"]"),
  make_option(c("--formula_micro"),action="store", default = "~ -1 + Stringency",  help="Features to be used for micro level pooling [default \"%default\"]")
)
opt <- parse_args(OptionParser(option_list=option_list))
StanModel <- opt$stanFile
formula_full <- opt$formula_full
formula_macro <- opt$formula_macro
formula_micro <- opt$formula_micro
DEBUG <- opt$debug
FULL <- opt$full
cat(sprintf("Running:\nStanModel = %s\nformula_full = %s\nformula_macro = %s\nformula_micro = %s\nDebug: %s\nFull: %s\n",
            StanModel, formula_full, formula_macro, formula_micro, DEBUG, FULL))
formula_full <- as.formula(formula_full)
formula_macro <- as.formula(formula_macro)
formula_micro <- as.formula(formula_micro)

# Read which countries to use
countries.df <- read_csv(paste(here(), 'data/regions-europe.csv', sep = '/'))
regions <- countries.df$Regions
# Read ifr 
ifr_by_regions <- read_csv(paste(here(),'data/ifr-pop-europe.csv',sep = '/'))
# selecting countries that we need
df_data <- 
  df_data %>%
  filter(country_name %in% regions) %>%
  rename(Date = DATE, Cases = cases_new, Deaths = deaths_new, Stringency = npi_stringency_index) %>%
  filter(Date <= as.Date('2020-05-05')) %>%
  mutate(Stringency = Stringency / 100) %>%
  na.locf(Stringency = na.locf(Stringency)) %>% # for days we don't have stringency just copy old
  mutate(Cases = case_when( Cases <= 0 ~ 0, TRUE ~ Cases), Deaths = case_when( Deaths <= 0 ~ 0, TRUE ~ Deaths)) # just make neagtive data as zero
# num of days to forecast
forecast <- 0
# Maximum number of days to simulate
num_days_sim <- (max(df_data$Date) - min(df_data$Date) + 1 + forecast)[[1]]


processed_data <- process_covariates(regions = regions,
                                     df_data = df_data,
                                     ifr_by_regions = ifr_by_regions, 
                                     num_days_sim = num_days_sim,
                                     formula_full = formula_full,
                                     formula_macro = formula_macro,
                                     formula_micro = formula_micro,
                                     create_features = TRUE)
stan_data = processed_data$stan_data
dates <- processed_data$dates
reported_deaths <- processed_data$reported_deaths
reported_cases <- processed_data$reported_cases

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stan_file_name = paste0('stan-models/',StanModel,'.stan')
m = stan_model(paste(here(),stan_file_name,sep='/'))

JOBID = Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID = as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s",JOBID))

if(FULL) {
  fit = sampling(m,data=stan_data,iter=opt$iter,warmup=500,chains=opt$nchains,thin=opt$thin,control = list(adapt_delta = 0.95, max_treedepth = 20))
} else if (DEBUG) {
  fit = sampling(m,data=stan_data,iter=100,warmup=50,chains=2)
} else { 
  fit = sampling(m,data=stan_data,iter=800,warmup=400,chains=2,thin=1,control = list(adapt_delta = 0.90, max_treedepth = 15))
}   
out = rstan::extract(fit)
estimated_cases_raw = out$prediction
estimated_deaths_raw = out$E_deaths

filename <- paste0(StanModel,'-',JOBID)

data_save <- list('fit'= fit, 'dates'= dates, 'reported_cases' = reported_cases, 'reported_deaths'= reported_deaths, 
                  'regions' = regions, 'estimated_cases_raw' = estimated_cases_raw, 'estimated_deaths_raw' = estimated_deaths_raw,
                  'formula_full'= formula_full, 'formula_macro'= formula_macro,'formula_micro'= formula_micro, 'stan_data'= stan_data, 
                  'observations'=df_data,'JOBID'= JOBID)
save_file = paste0('results/',StanModel,'-',JOBID,'-stanfit.rds') 
saveRDS(data_save, paste(here(),save_file,sep='/'))


make_three_pannel_plot(filename)
make_covariate_effect_size_plot(filename)
