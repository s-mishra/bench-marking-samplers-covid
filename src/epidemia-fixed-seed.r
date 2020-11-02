library(here)
data = readRDS(here("src/epidemia_fixed_seed.rds"))
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(here("stan-models/epidemia_fixed_seed.stan"))

fit = sampling(m, data=data, iter=4000, control = list(max_treedepth=20, adapt_delta=0.99),thin=4,seed=1)

