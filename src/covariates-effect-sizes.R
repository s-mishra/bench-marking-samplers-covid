library(ggplot2)
library(scales)
library(bayesplot)
library(matrixStats)
library(cowplot)
library(here)
library(tidybayes)

make_covariate_effect_size_plot <- function(filename) {
  print(sprintf("loading: %s/results/%s-stanfit.rds",here(),filename))
  fit_file = paste0('results/', filename,'-stanfit.rds') 
  data_save <- readRDS(paste(here(), fit_file, sep='/'))
  out <- rstan::extract(data_save$fit)
  alpha_index  = 1:length(dim(out$alpha)[2])
  g <- fit %>%
    spread_draws(alpha[alpha_index]) %>%
    median_qi(estimate = (1 - 2 *plogis(-alpha)), .width = c(.95, .66)) %>%
    ggplot(aes(y = alpha_index, x = estimate, xmin = .lower, xmax = .upper)) +
    geom_pointinterval()
  figure_file <- paste0("covars-effect-size-", filename, ".png")
  save_plot(filename = paste(here(),"figures",figure_file, sep='/' ), 
            g, base_width = 14)
}
