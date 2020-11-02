library(here)
data = readRDS(here("src/epidemia_fixed_seed.rds"))
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(here("stan-models/epidemia_fixed_seed.stan"))

fit = sampling(m, data=data, iter=400, seed=7, control = list(max_treedepth=15, adapt_delta=0.95),thin=2)