rm(list = ls())

library(loo)
library(dplyr)
library(reshape2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(grid)

source('fun/additional_functions.R')


plot_my_results <- function (res, max_lambda, max_prev, horizontal = TRUE) {
  
  model1 = res$fit1 
  model2 = res$fit3
  fit1    <- model1$fit
  fit2    <- model2$fit
  dat    <- res$dat
  # if(dat$tsur[1] == 2012) {virus  <- paste0(dat$virus[1],', 2012')}
  # if(dat$tsur[1] == 2017) {virus  <- paste0(dat$virus[1],', 2017')}
  virus     <- dat$virus[1]
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  RealYexpo <- model1$RealYexpo
  
  
  
  
  extract_foi_summary <- function (res, model, fit, dat, RealYexpo)
  {
    N       <- sum(dat$total)
    foi_chain <- data.frame(rstan::extract(fit, 'foi')[[1]])
    foi_summary <- data.frame(sapply(foi_chain, function(i) c(quantile(i, c(0.5, 0.025, 0.975)))))
    colnames(foi_summary) <- RealYexpo
    foi_summary$metric <- c('median', 'lower', 'upper') 
    foi_summary <- melt(foi_summary, id='metric', variable.name = 'year') %>% spread(metric, value) 
    foi_summary$year <- (as.numeric(as.character(foi_summary$year)))
    return(foi_summary)
    
  }
  
  foi1 <- extract_foi_summary(res, model1, fit1, dat, RealYexpo) %>%
    mutate (model = 'constant')
  foi2 <- extract_foi_summary(res, model2, fit2, dat, RealYexpo) %>%
    mutate (model = 'time_varying')
  foi <- rbind(foi1, foi2)
  size_text <- 25
  space_years <- 5
  breaks_x <- round(seq(1960, max(foi_overt$year), by = space_years),0) 
  foi_overt$model <- 'decades'
  g1 <- 
    foi_overt %>% filter(year > 1960) %>%
    ggplot(aes(x=year, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = .2) +
    geom_line(aes(y = median,  color = model)) +
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
  
  
  foi_fit1 <- data.frame(rstan::extract(fit, 'foi')[[1]])
  iterf = NROW(foi_fit1)
  lambdaYexpo1 <- as.matrix(foi_fit1)
  ExposureMatrixFit <- get_exposure_matrix(dat, yexpo)[, -1]
  
  PrevP1 <- matrix(NA, nrow = iterf, ncol = length(dat$age_mean_f))
  PrevP2 <- matrix(NA, nrow = iterf, ncol = length(dat$age_mean_f))
  
  for (i in 1:iterf)
  {
    PrevP1[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpo1[i,])
    PrevP2[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpo2[i,])
  } 
  
  PPP1 <- matrix(NA, ncol = 3, nrow = length(dat$age_mean_f))
  PPP2 <- matrix(NA, ncol = 3, nrow = length(dat$age_mean_f))
  for (j in seq_along(dat$age_mean_f)) 
  {
    PPP1[j,] <- quantile(PrevP1[,j], c(.025, .5, .975))
    PPP2[j,] <- quantile(PrevP2[,j], c(.025, .5, .975))
  }
  
  
  PPP1        <- as.data.frame(PPP1) %>% mutate(age = dat$age_mean_f,
                                                model = 'constant')
  PPP2        <- as.data.frame(PPP2) %>% mutate(age = dat$age_mean_f,
                                                model = 'time_varying')
  
  PPP_both <- rbind(PPP1, PPP2)
  colnames(PPP_both) <- c("L", "M", "U", "age", 'model')
  
  
  age_mean_fix <- (dat$age_min + dat$age_max)/2
  PPP_both$age_mean_fix <- age_mean_fix 
  dat$age_mean_fix      <- age_mean_fix
  
  g2 <-
    ggplot() +
    geom_errorbar(data=dat, 
                  aes(x=age_mean_fix, ymin=prev_obs_lower, ymax = prev_obs_upper),
                  color = 'grey', width = .1) +
    geom_ribbon(data=PPP_both, aes(x=age_mean_fix, ymin = L, ymax = U, fill = model), alpha = .2) +
    geom_line  (data=PPP_both, aes(x=age_mean_fix, y= M, color = model)) + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
      geom_point (data=dat, aes(x=age_mean_fix, y=prev_obs, size = total), 
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
                     rel_heights = c(0.9, 1)
  )
  
  
  return(p_res)
}

VEEV2012_P1_05 <- readRDS('res/VEEV 2012 Pi-Pi 005.RDS')
VEEV2012_P1_10 <- readRDS('res/VEEV 2012 Pi-Pi 010.RDS')
VEEV2012_P2_05 <- readRDS('res/VEEV 2012 Merca 005.RDS')
VEEV2012_P2_10 <- readRDS('res/VEEV 2012 Merca 010.RDS')
VEEV2012_P3_05 <- readRDS('res/VEEV 2012 Tamar 005.RDS')
VEEV2012_P3_10 <- readRDS('res/VEEV 2012 Tamar 010.RDS')
VEEV2012_P4_05 <- readRDS('res/VEEV 2012 Real 005.RDS')
VEEV2012_P4_10 <- readRDS('res/VEEV 2012 Real 010.RDS')
VEEV2012_P5_05 <- readRDS('res/VEEV 2012 Aruza 005.RDS')
VEEV2012_P5_10 <- readRDS('res/VEEV 2012 Aruza 010.RDS')
VEEV2017_P6_05 <- readRDS('res/VEEV 2017 Mogue 005.RDS')
VEEV2017_P6_10 <- readRDS('res/VEEV 2017 Mogue 010.RDS')

MADV2012_P1_05 <- readRDS('res/MADV 2012 Pi-Pi 005.RDS')
MADV2012_P1_10 <- readRDS('res/MADV 2012 Pi-Pi 010.RDS')
MADV2012_P2_05 <- readRDS('res/MADV 2012 Merca 005.RDS')
MADV2012_P2_10 <- readRDS('res/MADV 2012 Merca 010.RDS')
MADV2012_P3_05 <- readRDS('res/MADV 2012 Tamar 005.RDS')
MADV2012_P3_10 <- readRDS('res/MADV 2012 Tamar 010.RDS')
MADV2012_P4_05 <- readRDS('res/MADV 2012 Real 005.RDS')
MADV2012_P4_10 <- readRDS('res/MADV 2012 Real 010.RDS')
MADV2012_P5_05 <- readRDS('res/MADV 2012 Aruza 005.RDS')
MADV2012_P5_10 <- readRDS('res/MADV 2012 Aruza 010.RDS')
MADV2017_P6_05 <- readRDS('res/MADV 2017 Mogue 005.RDS')
MADV2017_P6_10 <- readRDS('res/MADV 2017 Mogue 010.RDS')



UNAV2017_P6_05 <- readRDS('res/UNAV 2017 Mogue 005.RDS')
UNAV2017_P6_10 <- readRDS('res/UNAV 2017 Mogue 010.RDS')


results_model <- MADV2017_P6_10
name_file <- paste0('res/', results_model$dat$virus[1], ' ', results_model$dat$PlaceName[1], ".png")
png(name_file, width = 350 *2, height = 350 *1)
plot_my_results(results_model, max_lambda = 0.3, max_prev = 1, horizontal = TRUE)
dev.off()



png("res/best_foi.png",
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


gridExtra::grid.arrange(
  plot_my_results(MADV2012_P1_10, max_lambda = 0.07, max_prev = 0.6),
  plot_my_results(MADV2012_P2_10, max_lambda = 0.07, max_prev = 0.6),
  plot_my_results(MADV2012_P3_10, max_lambda = 0.07, max_prev = 0.6),
  plot_my_results(MADV2012_P4_10, max_lambda = 0.07, max_prev = 0.6),
  plot_my_results(MADV2012_P5_10, max_lambda = 0.07, max_prev = 0.6),
  plot_my_results(MADV2017_P6_10, max_lambda = 0.07, max_prev = 0.6),
  nrow = 1,
  top= title2
)

dev.off()




extract_metrics <- function(model, res, fit)
{
  DIC <- (loo(fit, save_psis = TRUE, 'logLikelihood')[1])$estimates[,1]
  virus <- res$dat$virus[1]
  PlaceName  <- res$dat$PlaceName[1]
  YearSurvey <- res$dat$tsur[1]
  SampleSize <- sum(res$dat$total)
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  res_DIC <- data.frame(virus = virus, Place = PlaceName, Year_survey = YearSurvey,
                        sample_size =SampleSize,
                        ll = DIC[1], DIC = DIC[3], model = model)
  row.names(res_DIC) <- NULL
  return(res_DIC)

}


res_DIC <-rbind(
        extract_metrics('constant',     VEEV2012_P1_10, VEEV2012_P1_10$fit1$fit),
        extract_metrics('time-varying', VEEV2012_P1_10, VEEV2012_P1_10$fit3$fit),
        extract_metrics('constant',     VEEV2012_P2_10, VEEV2012_P2_10$fit1$fit),
        extract_metrics('time-varying', VEEV2012_P2_10, VEEV2012_P2_10$fit3$fit),
        extract_metrics('constant',     VEEV2012_P3_10, VEEV2012_P3_10$fit1$fit),
        extract_metrics('time-varying', VEEV2012_P3_10, VEEV2012_P3_10$fit3$fit),
        extract_metrics('constant',     VEEV2012_P4_10, VEEV2012_P4_10$fit1$fit),
        extract_metrics('time-varying', VEEV2012_P4_10, VEEV2012_P4_10$fit3$fit),
        extract_metrics('constant',     VEEV2012_P5_10, VEEV2012_P5_10$fit1$fit),
        extract_metrics('time-varying', VEEV2012_P5_10, VEEV2012_P5_10$fit3$fit),
        extract_metrics('constant',     VEEV2017_P6_10, VEEV2017_P6_10$fit1$fit),
        extract_metrics('time-varying', VEEV2017_P6_10, VEEV2017_P6_10$fit3$fit),
        
        extract_metrics('constant',     MADV2012_P1_10, MADV2012_P1_10$fit1$fit),
        extract_metrics('time-varying', MADV2012_P1_10, MADV2012_P1_10$fit3$fit),
        extract_metrics('constant',     MADV2012_P2_10, MADV2012_P2_10$fit1$fit),
        extract_metrics('time-varying', MADV2012_P2_10, MADV2012_P2_10$fit3$fit),
        extract_metrics('constant',     MADV2012_P3_10, MADV2012_P3_10$fit1$fit),
        extract_metrics('time-varying', MADV2012_P3_10, MADV2012_P3_10$fit3$fit),
        extract_metrics('constant',     MADV2012_P4_10, MADV2012_P4_10$fit1$fit),
        extract_metrics('time-varying', MADV2012_P4_10, MADV2012_P4_10$fit3$fit),
        extract_metrics('constant',     MADV2012_P5_10, MADV2012_P5_10$fit1$fit),
        extract_metrics('time-varying', MADV2012_P5_10, MADV2012_P5_10$fit3$fit),
        extract_metrics('constant',     MADV2017_P6_10, MADV2017_P6_10$fit1$fit),
        extract_metrics('time-varying', MADV2017_P6_10, MADV2017_P6_10$fit3$fit)
        
)

write_csv(res_DIC, 'res/res_DIC.csv')

# 
# library(ggpubr)
# empty_plot <- ggplot(dat = data.frame(x= 1:10, y = 1:10, 
#                                       model = c('constant', 'time varying'))) +
# geom_line(aes(x, y, color = model)) +
# geom_ribbon(aes(x, ymin = y, ymax=y, fill = model), alpha = 0.2) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_color_brewer(palette = "Set1") +
#   theme(legend.position="right",
#         legend.text=element_text(size=25))
# legend_plots <- get_legend(empty_plot) 
# 
# 
# png('res/plots.png', width = 480*3.5, height = 480 * 1.5)
# gridExtra::grid.arrange(plot_my_results(MADV2012, max_lambda = 0.05, max_prev = 0.6), 
#                         plot_my_results(VEEV2012, max_lambda = 0.05, max_prev = 0.6),
#                         plot_my_results(MADV2017, max_lambda = 0.05, max_prev = 0.6),
#                         plot_my_results(VEEV2017, max_lambda = 0.05, max_prev = 0.6),
#                         plot_my_results(UNAV2017, max_lambda = 0.05, max_prev = 0.6),
#                         as_ggplot(legend_plots),
#                         nrow = 1)
# dev.off()
# 
# 
# extract_metrics <- function(model, res, fit)
# {
#   DIC <- (loo(fit, save_psis = TRUE, 'logLikelihood')[1])$estimates[,1]
#   virus <- paste(res$dat$virus[1], res$dat$tsur[1])
#   res_DIC <- data.frame(virus = virus, ll = DIC[1], DIC = DIC[3], model = model)
#   row.names(res_DIC) <- NULL
#   return(res_DIC)
#   
# }
# 
# 
# res_DIC <-rbind(
#       extract_metrics('constant',     MADV2012, MADV2012$fit1$fit),
#       extract_metrics('time-varying', MADV2012, MADV2012$fit2$fit),
#       extract_metrics('constant',     VEEV2012, VEEV2012$fit1$fit),
#       extract_metrics('time-varying', VEEV2012, VEEV2012$fit2$fit),
#       
#       extract_metrics('constant',     MADV2017, MADV2017$fit1$fit),
#       extract_metrics('time-varying', MADV2017, MADV2017$fit2$fit),
#       extract_metrics('constant',     VEEV2017, VEEV2017$fit1$fit),
#       extract_metrics('time-varying', VEEV2017, VEEV2017$fit2$fit),
#       extract_metrics('constant',     UNAV2017, UNAV2017$fit1$fit),
#       extract_metrics('time-varying', UNAV2017, UNAV2017$fit2$fit)
#       )
#       
# 
# write_csv(res_DIC, 'res/res_DIC.csv')
# 
# # extract_ll_DIC <- function(res, fit, model) 
# # {
# #   dat <- res$dat
# #   loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
# #   ExposureMatrixFit <- get_exposure_matrix(dat)[, -1]
# #   foi <- data.frame(rstan::extract(fit, 'foi')[[1]])
# #   lambdaTM <- (as.numeric(sapply(foi, function(i) c(quantile(i, 0.5)))))
# #   
# #   PrevM    <-  1-(exp(- ExposureMatrixFit %*% lambdaTM))
# #   llnew    <-  sum(dbinom(as.matrix(dat$pos), as.matrix(dat$total), PrevM, log = TRUE))
# #   
# #   # then, it calculates DIC (using the meadian values) and makes the DIC table
# #   ll <- loo_fit$estimates[1,1]
# #   dev_a    <- -2 * ll
# #   dev_b    <- -2 * llnew
# #   Pd       <- dev_a - dev_b
# #   DIC      <- Pd + dev_a
# #   
# #   virus    <- dat$virus[1]
# #   model    <- 
# #     
# #     res_DIC <- data.frame(log_lik = ll, DIC = DIC, virus, model)
# #   
# #   return(res_DIC)
# #   
# # }
# # 
# # extract_ll_DIC(MADV2012, MADV2012$fit1$fit, model = 'constant')
# # extract_ll_DIC(MADV2012, MADV2012$fit2$fit, model = 'time-varying')
# # extract_ll_DIC(VEEV2012, VEEV2012$fit1$fit, model = 'constant')
# # extract_ll_DIC(MADV2012, MADV2012$fit2$fit, model = 'time-varying')
# 
# 
# 
# ip <- as.data.frame(installed.packages()[,c(1,3:4)])
# rownames(ip) <- NULL
# ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
# print(ip, row.names=FALSE)
# 
# stan_version()
