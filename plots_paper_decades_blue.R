rm(list = ls())

library(loo)
library(dplyr)
library(reshape2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(grid)
library(rstan)
library(tidyverse)
library(reshape2)
library(bayesplot)
library(loo)
library(pracma)
library(cowplot)
library(grid)
library(gridExtra)

source('fun/additional_functions.R')


res <- readRDS('res/mod_decsVEEV 2012 Pi-Pi 010.RDS')

plot_my_results <- function(res, max_lambda=0.3, max_prev=1, horizontal = TRUE){
  
  
  
  dat    <- res$dat
  fit    <- res$res$fit
  virus     <- dat$virus[1]
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  tsur      <- dat$tsur[1]
  extract_foi_summary <- function (res)
  {
    dat  <- res$dat
    N       <- sum(dat$total)
    RealYexpo <- res$res$RealYexpo
    foi_chain <- data.frame(rstan::extract(fit, 'foi')[[1]])
    foi_summary <- data.frame(sapply(foi_chain, function(i) c(quantile(i, c(0.5, 0.025, 0.975)))))
    foi_summary_yexpo <- matrix(NA, ncol = 50, nrow = 3)
    
    foi_summary_yexpo[, 1:10]  <- foi_summary[,1]
    foi_summary_yexpo[, 11:20] <- foi_summary[,2]
    foi_summary_yexpo[, 21:30] <- foi_summary[,3]
    foi_summary_yexpo[, 31:40] <- foi_summary[,4]
    foi_summary_yexpo[, 41:50] <- foi_summary[,5]
    foi_summary_y <- as.data.frame(foi_summary_yexpo)
    names(foi_summary_y) <- (tsur-49):tsur
    foi_summary_y$metric <- c('median', 'lower', 'upper') 
    foi_summary_yo <- melt(foi_summary_y, id='metric', variable.name = 'year') %>% spread(metric, value) 
    foi_summary_yo$year <- (as.numeric(as.character(foi_summary_yo$year)))
    return(foi_summary_yo)
    
  }
  
  foi_overt <- extract_foi_summary(res)
  
  size_text <- 25
  space_years <- 5
  breaks_x <- round(seq(1960, max(foi_overt$year), by = space_years),0) 
  foi_overt$model <- 'decades'
  g1 <- 
    foi_overt %>% filter(year > 1960) %>%
    ggplot(aes(x=year, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = '#6baed6', alpha = .2) +
    geom_line(aes(y = median), color = 'blue') +
    scale_fill_brewer("FOI", palette = "Set1") +
    scale_color_brewer("FOI", palette = "Set1") +
    coord_cartesian(ylim = c(0, max_lambda)) + 
    theme_bw() + 
    theme(legend.position = 'none') +
    ylab ('FOI') +
    theme(
      legend.position="none", 
      text = element_text(size=size_text),
      plot.title = element_text(size = size_text*.8, face = "bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    scale_x_continuous(breaks = breaks_x)  +
    ggtitle (paste(PlaceName))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  foi_fit <- data.frame(rstan::extract(fit, 'foi')[[1]])
  iterf = NROW(foi_fit)
  lambdaYexpo <- as.matrix(foi_fit)
  # lambdaYexpoExpanded <- matrix(NA, ncol = 60, nrow = NROW(lambdaYexpo))
  # lambdaYexpoExpanded[, 1:10]  <- lambdaYexpo[,1]
  # lambdaYexpoExpanded[, 11:20] <- lambdaYexpo[,2]
  # lambdaYexpoExpanded[, 21:30] <- lambdaYexpo[,3]
  # lambdaYexpoExpanded[, 31:40] <- lambdaYexpo[,4]
  # lambdaYexpoExpanded[, 41:50] <- lambdaYexpo[,5]
  

  min_age         <- 1
  max_age         <- 60
  yexpo_expanded <- min_age:max_age
  
  lambdaYexpoExpanded <- matrix(NA, nrow = iterf, ncol= max_age)
  Ymax <- max_age
  lambdaYexpoExpanded[,(Ymax-9):Ymax] <- lambdaYexpo[,5]
  lambdaYexpoExpanded[,(Ymax-19):(Ymax-10)] <- lambdaYexpo[,4]
  lambdaYexpoExpanded[,(Ymax-29):(Ymax-20)] <- lambdaYexpo[,3]
  lambdaYexpoExpanded[,(Ymax-39):(Ymax-30)] <- lambdaYexpo[,2]
  lambdaYexpoExpanded[,1:(Ymax-40)] <- lambdaYexpo[,1]
  
  ExposureMatrixFit <- get_exposure_matrix_expanded(yexpo_expanded, min_age ,max_age)
  
  PrevP <- matrix(NA, nrow = iterf, ncol = max(yexpo_expanded))
  for (i in 1:iterf)
  {
    PrevP[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpoExpanded[i,])
    
  } 
  
  PPP <- matrix(NA, ncol = 3, nrow = max_age)
  for (j in seq_along(min_age:max_age)) 
  {
    PPP[j,] <- quantile(PrevP[,j], c(.025, .5, .975))
    
  }
  
  
  PPP        <- as.data.frame(PPP) %>% mutate(age = min_age:max_age,
                                              model = 'time_varying')
  
  colnames(PPP) <- c("L", "M", "U", "age", 'model')
  
  PPP_both <- PPP
  
  
  g2 <-
    ggplot() +
    geom_errorbar(data=dat, 
                  aes(x=age_mean_f, ymin=prev_obs_lower, ymax = prev_obs_upper),
                  color = 'grey', width = .1) +
    geom_ribbon(data=PPP_both, aes(x=age, ymin = L, ymax = U), fill = '#6baed6', alpha = .2) +
    geom_line  (data=PPP_both, aes(x=age, y= M), color ='blue') + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    geom_point (data=dat, aes(x=age_mean_f, y=prev_obs, size = total), 
                pch = 21,fill = 'white', color = 'black') +
    xlab("Age, years") +
    ylab("seropositivity") +
    theme_bw() + coord_cartesian(ylim = c(0,max_prev), xlim = c(3,50)) +
    theme(axis.text.x = element_text(angle= 0)) +
    theme(
      legend.position="none", 
      text = element_text(size=size_text),
      plot.title = element_text(size = size_text*.8, face = "bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  
  if(horizontal == TRUE){num_rows = 1} else {num_rows = 2}
  
  p_res <- plot_grid(g1, g2,
                     nrow = num_rows, 
                     align="hv",
                     rel_heights = c(0.9, 1))
  
  
  
  
  
  
  return(p_res)
  
  
}

VEEV2012_P1_05 <- readRDS('res/mod_decsVEEV 2012 Pi-Pi 005.RDS')
VEEV2012_P1_10 <- readRDS('res/mod_decsVEEV 2012 Pi-Pi 010.RDS')
VEEV2012_P2_05 <- readRDS('res/mod_decsVEEV 2012 Merca 005.RDS')
VEEV2012_P2_10 <- readRDS('res/mod_decsVEEV 2012 Merca 010.RDS')
VEEV2012_P3_05 <- readRDS('res/mod_decsVEEV 2012 Tamar 005.RDS')
VEEV2012_P3_10 <- readRDS('res/mod_decsVEEV 2012 Tamar 010.RDS')
VEEV2012_P4_05 <- readRDS('res/mod_decsVEEV 2012 Real 005.RDS')
VEEV2012_P4_10 <- readRDS('res/mod_decsVEEV 2012 Real 010.RDS')
VEEV2012_P5_05 <- readRDS('res/mod_decsVEEV 2012 Aruza 005.RDS')
VEEV2012_P5_10 <- readRDS('res/mod_decsVEEV 2012 Aruza 010.RDS')
VEEV2017_P6_05 <- readRDS('res/mod_decsVEEV 2017 Mogue 005.RDS')
VEEV2017_P6_10 <- readRDS('res/mod_decsVEEV 2017 Mogue 010.RDS')

MADV2012_P1_05 <- readRDS('res/mod_decsMADV 2012 Pi-Pi 005.RDS')
MADV2012_P1_10 <- readRDS('res/mod_decsMADV 2012 Pi-Pi 010.RDS')
MADV2012_P2_05 <- readRDS('res/mod_decsMADV 2012 Merca 005.RDS')
MADV2012_P2_10 <- readRDS('res/mod_decsMADV 2012 Merca 010.RDS')
MADV2012_P3_05 <- readRDS('res/mod_decsMADV 2012 Tamar 005.RDS')
MADV2012_P3_10 <- readRDS('res/mod_decsMADV 2012 Tamar 010.RDS')
MADV2012_P4_05 <- readRDS('res/mod_decsMADV 2012 Real 005.RDS')
MADV2012_P4_10 <- readRDS('res/mod_decsMADV 2012 Real 010.RDS')
MADV2012_P5_05 <- readRDS('res/mod_decsMADV 2012 Aruza 005.RDS')
MADV2012_P5_10 <- readRDS('res/mod_decsMADV 2012 Aruza 010.RDS')
MADV2017_P6_05 <- readRDS('res/mod_decsMADV 2017 Mogue 005.RDS')
MADV2017_P6_10 <- readRDS('res/mod_decsMADV 2017 Mogue 010.RDS')



UNAV2017_P6_05 <- readRDS('res/UNAV 2017 Mogue 005.RDS')
UNAV2017_P6_10 <- readRDS('res/UNAV 2017 Mogue 010.RDS')



png("res/foi_decs_VEEV.png",
    width = 350 *2, height = 350 *6)

title1=textGrob("VEEV", gp=gpar(fontface="bold"))
title2=textGrob("MADV", gp=gpar(fontface="bold"))

gridExtra::grid.arrange(
  plot_my_results(VEEV2012_P1_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P2_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P3_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P4_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P5_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2017_P6_10, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  nrow = 6,
  top= title1
)

dev.off()


png("res/foi_decs_VEEV_5y.png",
    width = 350 *2, height = 350 *6)

title1=textGrob("VEEV", gp=gpar(fontface="bold"))
title2=textGrob("MADV", gp=gpar(fontface="bold"))

gridExtra::grid.arrange(
  plot_my_results(VEEV2012_P1_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P2_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P3_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P4_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2012_P5_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  plot_my_results(VEEV2017_P6_05, max_lambda = 0.3, max_prev = 1, horizontal = TRUE),
  nrow = 6,
  top= title1
)

dev.off()


png("res/foi_decs_MADV_5y.png",
    width = 350 *2, height = 350 *6)

gridExtra::grid.arrange(
  plot_my_results(MADV2012_P1_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P2_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P3_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P4_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P5_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2017_P6_05, max_lambda = 0.1, max_prev = 0.6, horizontal = TRUE),
  nrow = 6,
  top= textGrob("MADV 5y", gp=gpar(fontface="bold"))
)

dev.off()








png("res/foi_decs_MADV.png",
    width = 350 *2, height = 350 *6)

gridExtra::grid.arrange(
  plot_my_results(MADV2012_P1_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P2_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P3_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P4_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2012_P5_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  plot_my_results(MADV2017_P6_10, max_lambda = 0.07, max_prev = 0.6, horizontal = TRUE),
  nrow = 6,
  top= title2
)

dev.off()




extract_ll_DIC <- function(res,model_name)
{
  dat <- res$dat
  fit <- res$res$fit
  loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
  yexpo   <- make_yexpo(dat)
  yexpo <- yexpo[-length(yexpo)]
  ExposureMatrixFit <- get_exposure_matrix(dat, yexpo)
  foi <- data.frame(rstan::extract(fit, 'foi')[[1]])
  lambdaTM <- (as.numeric(sapply(foi, function(i) c(quantile(i, 0.5)))))
  
  if(model_name == 'decades'){
    lambdaTM_expanded <- matrix(NA, nrow = 1, ncol=max(yexpo))
    Ymax <- max(yexpo)
    lambdaTM_expanded[(Ymax-9):Ymax] <- lambdaTM[5]
    lambdaTM_expanded[(Ymax-19):(Ymax-10)] <- lambdaTM[4]
    lambdaTM_expanded[(Ymax-29):(Ymax-20)] <- lambdaTM[3]
    lambdaTM_expanded[(Ymax-39):(Ymax-30)] <- lambdaTM[2]
    lambdaTM_expanded[1:(Ymax-40)] <- lambdaTM[1]
    
    lambdaTM <- lambdaTM_expanded
    
  }
  
  
  PrevM    <-  1-(exp(- ExposureMatrixFit %*% as.numeric(lambdaTM)))
  llnew    <-  sum(dbinom(as.matrix(dat$pos), as.matrix(dat$total), PrevM, log = TRUE))
  
  # then, it calculates DIC (using the meadian values) and makes the DIC table
  ll <- loo_fit$estimates[1,1]
  dev_a    <- -2 * ll
  dev_b    <- -2 * llnew
  Pd       <- dev_a - dev_b
  DIC      <- Pd + dev_a
  
  
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  virus    <- dat$virus[1]
  res_DIC <- data.frame(Place = PlaceName, 
                        virus = dat$virus[1], 
                        sample_size = sum(dat$total), 
                        log_lik = ll, DIC = DIC,
                        model   = model_name)
  
  return(res_DIC)
  
}

extract_ll_DIC(res = VEEV2012_P1_05, model_name=
                 'decades')
