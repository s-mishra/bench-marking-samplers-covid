library(ggplot2)
library(tidybayes)
library(rstan)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(lubridate)
library(here)
#---------------------------------------------------------------------------
make_three_pannel_plot <- function(filename){
  
  print(sprintf("loading: %s/results/%s-stanfit.rds",here(),filename))
  fit_file = paste0('results/', filename,'-stanfit.rds') 
  data_save <- readRDS(paste(here(), fit_file, sep='/'))
  out <- rstan::extract(data_save$fit)
  dates <- data_save$dates
  
  for(i in 1:length(data_save$regions)){
    print(sprintf('Plotting %s', data_save$regions[[i]]))
    N <- length(data_save$dates[[i]])
    region <- data_save$regions[[i]]
    predicted_cases <- colMeans(data_save$estimated_cases_raw [,1:N,i])
    predicted_cases_li <- colQuantiles(data_save$estimated_cases_raw [,1:N,i], probs=.025)
    predicted_cases_ui <- colQuantiles(data_save$estimated_cases_raw [,1:N,i], probs=.975)
    predicted_cases_li2 <- colQuantiles(data_save$estimated_cases_raw [,1:N,i], probs=.25)
    predicted_cases_ui2 <- colQuantiles(data_save$estimated_cases_raw [,1:N,i], probs=.75)
    
    
    estimated_deaths <- colMeans(data_save$estimated_deaths_raw[,1:N,i])
    estimated_deaths_li <- colQuantiles(data_save$estimated_deaths_raw[,1:N,i], probs=.025)
    estimated_deaths_ui <- colQuantiles(data_save$estimated_deaths_raw[,1:N,i], probs=.975)
    estimated_deaths_li2 <- colQuantiles(data_save$estimated_deaths_raw[,1:N,i], probs=.25)
    estimated_deaths_ui2 <- colQuantiles(data_save$estimated_deaths_raw[,1:N,i], probs=.75)
    
    rt <- colMeans(out$Rt_adj[,1:N,i])
    rt_li <- colQuantiles(out$Rt_adj[,1:N,i],probs=.025)
    rt_ui <- colQuantiles(out$Rt_adj[,1:N,i],probs=.975)
    rt_li2 <- colQuantiles(out$Rt_adj[,1:N,i],probs=.25)
    rt_ui2 <- colQuantiles(out$Rt_adj[,1:N,i],probs=.75)
    
    data_region <- data.frame("time" = as_date(as.character(dates[[i]])),
                              "region" = rep(region, length(dates[[i]])),
                              "reported_cases" = data_save$reported_cases[[i]], 
                              "predicted_cases" = predicted_cases,
                              "predicted_min" = predicted_cases_li,
                              "predicted_max" = predicted_cases_ui,
                              "predicted_min2" = predicted_cases_li2,
                              "predicted_max2" = predicted_cases_ui2,
                              "deaths" = data_save$reported_deaths[[i]],
                              "estimated_deaths" = estimated_deaths,
                              "death_min" = estimated_deaths_li,
                              "death_max"= estimated_deaths_ui,
                              "death_min2" = estimated_deaths_li2,
                              "death_max2"= estimated_deaths_ui2,
                              "rt" = rt,
                              "rt_min" = rt_li,
                              "rt_max" = rt_ui,
                              "rt_min2" = rt_li2,
                              "rt_max2" = rt_ui2)
    
    make_plots(data_region = data_region, 
               filename2 = filename,
               region = region)
    
  }
}

#---------------------------------------------------------------------------
make_plots <- function(data_region, filename2, region){
  
  data_cases_95 <- data.frame(data_region$time, data_region$predicted_min, 
                              data_region$predicted_max)
  names(data_cases_95) <- c("time", "cases_min", "cases_max")
  data_cases_95$key <- rep("nintyfive", length(data_cases_95$time))
  data_cases_50 <- data.frame(data_region$time, data_region$predicted_min2, 
                              data_region$predicted_max2)
  names(data_cases_50) <- c("time", "cases_min", "cases_max")
  data_cases_50$key <- rep("fifty", length(data_cases_50$time))
  data_cases <- rbind(data_cases_95, data_cases_50)
  levels(data_cases$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data_region) +
    geom_bar(data = data_region, aes(x = time, y = reported_cases), 
             fill = "coral4", stat='identity', alpha=0.5) + 
    geom_ribbon(data = data_cases, 
                aes(x = time, ymin = cases_min, ymax = cases_max, fill = key)) +
    xlab("") +
    ylab("Daily number of infections\n") +
    scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b")) + 
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55), 
                                 alpha("deepskyblue4", 0.45))) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "None") + ggtitle(region) +
    guides(fill=guide_legend(ncol=1))
  
  data_deaths_95 <- data.frame(data_region$time, data_region$death_min, 
                               data_region$death_max)
  names(data_deaths_95) <- c("time", "death_min", "death_max")
  data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
  data_deaths_50 <- data.frame(data_region$time, data_region$death_min2, 
                               data_region$death_max2)
  names(data_deaths_50) <- c("time", "death_min", "death_max")
  data_deaths_50$key <- rep("fifty", length(data_deaths_50$time))
  data_deaths <- rbind(data_deaths_95, data_deaths_50)
  levels(data_deaths$key) <- c("ninetyfive", "fifty")
  
  
  p2 <-   ggplot(data_region, aes(x = time)) +
    geom_bar(data = data_region, aes(y = deaths, fill = "reported"),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_ribbon(
      data = data_deaths,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b")) +
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55), 
                                 alpha("deepskyblue4", 0.45))) + 
    ylab("Daily number of deaths\n") + 
    xlab("") +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "None") + 
    guides(fill=guide_legend(ncol=1))
  
  
  # Plotting interventions
  data_rt_95 <- data.frame(data_region$time, 
                           data_region$rt_min, data_region$rt_max)
  names(data_rt_95) <- c("time", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
  data_rt_50 <- data.frame(data_region$time, data_region$rt_min2, 
                           data_region$rt_max2)
  names(data_rt_50) <- c("time", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$time))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  p3 <- ggplot(data_region) +
    geom_ribbon(data = data_rt, aes(x = time, ymin = rt_min, ymax = rt_max, 
                                    group = key,
                                    fill = key)) +
    geom_hline(yintercept = 1, color = 'black', size = 0.1) + 
    xlab("") +
    ylab(expression(R[t])) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) + 
    scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b"), 
                 limits = c(data_region$time[1], 
                            data_region$time[length(data_region$time)])) + 
    scale_y_continuous(expand = expansion(mult=c(0,0.1))) + 
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="right")
  
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
  figure_file <- paste0(region, "_three_pannel_", filename2, ".png")
  save_plot(filename = paste(here(),"figures",figure_file, sep='/' ), 
            p, base_width = 14)
}
