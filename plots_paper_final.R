rm(list = ls())

library(loo)
library(dplyr)
library(reshape2)
library(tidyverse)
library(cowplot)
source('fun/additional_functions.R')


plot_my_results <- function (res, max_lambda, max_prev) {
  
  model1 = res$fit1 
  model2 = res$fit2
  fit1    <- model1$fit
  fit2    <- model2$fit
  dat    <- res$dat
  if(dat$tsur[1] == 2012) {virus  <- paste0(dat$virus[1],', 2012')}
  if(dat$tsur[1] == 2017) {virus  <- paste0(dat$virus[1],', 2017')}
  
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
  space_years <- 1 
  if(length(model1$yexpo)>15 ) {space_years <- 2}
  if(length(model1$yexpo)>30 ) {space_years <- 5}
  if(length(model1$yexpo)>50 ) {space_years <- 8}
  
  
  breaks_x <- round(seq(1960, max(foi$year), by = space_years),0) 
  breaks_x <- round(seq(min(breaks_x) + (max(foi$year) - max(breaks_x)) , max(foi$year), by = space_years),0)
  g1 <- 
    foi %>% filter(year > 1960) %>%
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
    ggtitle (virus)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  foi_fit1 <- data.frame(rstan::extract(fit1, 'foi')[[1]])
  foi_fit2 <- data.frame(rstan::extract(fit2, 'foi')[[1]])
  iterf = NROW(foi_fit1)
  lambdaYexpo1 <- as.matrix(foi_fit1)
  lambdaYexpo2 <- as.matrix(foi_fit2)
  ExposureMatrixFit <- get_exposure_matrix(dat)[, -1]
  
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
  
  

  
  g2 <-
    ggplot() +
    geom_errorbar(data=dat, 
                  aes(x=age_mean_f, ymin=prev_obs_lower, ymax = prev_obs_upper),
                  color = 'grey', width = .1) +
    geom_ribbon(data=PPP_both, aes(x=age, ymin = L, ymax = U, fill = model), alpha = .2) +
    geom_line  (data=PPP_both, aes(x=age, y= M, color = model)) + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
      geom_point (data=dat, aes(x=age_mean_f, y=prev_obs, size = total), 
                  pch = 21,fill = 'white', color = 'black') +
    xlab("Age, years") +
    ylab("seropositivity") +
    theme_bw() + coord_cartesian(ylim = c(0,max_prev)) +
    theme(axis.text.x = element_text(angle= 0)) +
    theme(
      legend.position="none", 
      text = element_text(size=size_text),
      plot.title = element_text(size = size_text*.8, face = "bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p_res <- plot_grid(g1, g2,
                     nrow = 2, 
                     align="hv",
                     rel_heights = c(0.9, 1)
  )
  
  
  return(p_res)
}


# MADV2012 <- readRDS('res/MADV 2012 005.RDS')
# VEEV2012 <- readRDS('res/VEEV 2012 005.RDS')
# MADV2017 <- readRDS('res/MADV 2017 005.RDS')
# UNAV2017 <- readRDS('res/UNA 2017 005.RDS')
# VEEV2017 <- readRDS('res/VEEV 2017 005.RDS')

# # MADV2012 <- readRDS('res/MADV 2012 010.RDS')
# # VEEV2012 <- readRDS('res/VEEV 2012 010.RDS')
# MADV2017 <- readRDS('res/MADV 2017 010.RDS')
# UNAV2017 <- readRDS('res/UNA 2017 010.RDS')
# VEEV2017 <- readRDS('res/VEEV 2017 010.RDS')


MADV2012_P5_05 <- readRDS('res/MADV 2012 P5 005.RDS')
MADV2012_P5_10 <- readRDS('res/MADV 2012 P5 010.RDS')
MADV2012_P6_05 <- readRDS('res/MADV 2012 P6 005.RDS')
MADV2012_P6_10 <- readRDS('res/MADV 2012 P6 010.RDS')
MADV2012_P7_05 <- readRDS('res/MADV 2012 P7 005.RDS')
MADV2012_P7_10 <- readRDS('res/MADV 2012 P7 010.RDS')
MADV2012_P8_05 <- readRDS('res/MADV 2012 P8 005.RDS')
MADV2012_P8_10 <- readRDS('res/MADV 2012 P8 010.RDS')
MADV2012_P9_05 <- readRDS('res/MADV 2012 P9 005.RDS')
MADV2012_P9_10 <- readRDS('res/MADV 2012 P9 010.RDS')

plot_my_results(MADV2012_P9_10, max_lambda = 0.05, max_prev = 0.6)

gridExtra::grid.arrange(
  plot_my_results(MADV2012_P5_10, max_lambda = 0.05, max_prev = 0.6),
  plot_my_results(MADV2012_P6_10, max_lambda = 0.05, max_prev = 0.6),
  plot_my_results(MADV2012_P7_10, max_lambda = 0.05, max_prev = 0.6),
  plot_my_results(MADV2012_P8_10, max_lambda = 0.05, max_prev = 0.6),
  plot_my_results(MADV2012_P9_10, max_lambda = 0.05, max_prev = 0.6),
  nrow = 1
)



VEEV2012_P5_05 <- readRDS('res/VEEV 2012 P5 005.RDS')
VEEV2012_P5_10 <- readRDS('res/VEEV 2012 P5 010.RDS')
VEEV2012_P6_05 <- readRDS('res/VEEV 2012 P6 005.RDS')
VEEV2012_P6_10 <- readRDS('res/VEEV 2012 P6 010.RDS')
VEEV2012_P7_05 <- readRDS('res/VEEV 2012 P7 005.RDS')
VEEV2012_P7_10 <- readRDS('res/VEEV 2012 P7 010.RDS')
VEEV2012_P8_05 <- readRDS('res/VEEV 2012 P8 005.RDS')
VEEV2012_P8_10 <- readRDS('res/VEEV 2012 P8 010.RDS')
VEEV2012_P9_05 <- readRDS('res/VEEV 2012 P9 005.RDS')
VEEV2012_P9_10 <- readRDS('res/VEEV 2012 P9 010.RDS')

plot_my_results(VEEV2012_P6_10, max_lambda = 0.05, max_prev = 0.6)

gridExtra::grid.arrange(
  plot_my_results(VEEV2012_P5_10, max_lambda = 0.1, max_prev = 1),
  plot_my_results(VEEV2012_P6_10, max_lambda = 0.1, max_prev = 1),
  plot_my_results(VEEV2012_P7_05, max_lambda = 0.3, max_prev = 1),
  plot_my_results(VEEV2012_P8_10, max_lambda = 0.1, max_prev = 1),
  plot_my_results(VEEV2012_P9_10, max_lambda = 0.1, max_prev = 1),
  nrow = 1
)









library(ggpubr)
empty_plot <- ggplot(dat = data.frame(x= 1:10, y = 1:10, 
                                      model = c('constant', 'time varying'))) +
geom_line(aes(x, y, color = model)) +
geom_ribbon(aes(x, ymin = y, ymax=y, fill = model), alpha = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position="right",
        legend.text=element_text(size=25))
legend_plots <- get_legend(empty_plot) 


png('res/plots.png', width = 480*3.5, height = 480 * 1.5)
gridExtra::grid.arrange(plot_my_results(MADV2012, max_lambda = 0.05, max_prev = 0.6), 
                        plot_my_results(VEEV2012, max_lambda = 0.05, max_prev = 0.6),
                        plot_my_results(MADV2017, max_lambda = 0.05, max_prev = 0.6),
                        plot_my_results(VEEV2017, max_lambda = 0.05, max_prev = 0.6),
                        plot_my_results(UNAV2017, max_lambda = 0.05, max_prev = 0.6),
                        as_ggplot(legend_plots),
                        nrow = 1)
dev.off()


extract_metrics <- function(model, res, fit)
{
  DIC <- (loo(fit, save_psis = TRUE, 'logLikelihood')[1])$estimates[,1]
  virus <- paste(res$dat$virus[1], res$dat$tsur[1])
  res_DIC <- data.frame(virus = virus, ll = DIC[1], DIC = DIC[3], model = model)
  row.names(res_DIC) <- NULL
  return(res_DIC)
  
}


res_DIC <-rbind(
      extract_metrics('constant',     MADV2012, MADV2012$fit1$fit),
      extract_metrics('time-varying', MADV2012, MADV2012$fit2$fit),
      extract_metrics('constant',     VEEV2012, VEEV2012$fit1$fit),
      extract_metrics('time-varying', VEEV2012, VEEV2012$fit2$fit),
      
      extract_metrics('constant',     MADV2017, MADV2017$fit1$fit),
      extract_metrics('time-varying', MADV2017, MADV2017$fit2$fit),
      extract_metrics('constant',     VEEV2017, VEEV2017$fit1$fit),
      extract_metrics('time-varying', VEEV2017, VEEV2017$fit2$fit),
      extract_metrics('constant',     UNAV2017, UNAV2017$fit1$fit),
      extract_metrics('time-varying', UNAV2017, UNAV2017$fit2$fit)
      )
      

write_csv(res_DIC, 'res/res_DIC.csv')

# extract_ll_DIC <- function(res, fit, model) 
# {
#   dat <- res$dat
#   loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
#   ExposureMatrixFit <- get_exposure_matrix(dat)[, -1]
#   foi <- data.frame(rstan::extract(fit, 'foi')[[1]])
#   lambdaTM <- (as.numeric(sapply(foi, function(i) c(quantile(i, 0.5)))))
#   
#   PrevM    <-  1-(exp(- ExposureMatrixFit %*% lambdaTM))
#   llnew    <-  sum(dbinom(as.matrix(dat$pos), as.matrix(dat$total), PrevM, log = TRUE))
#   
#   # then, it calculates DIC (using the meadian values) and makes the DIC table
#   ll <- loo_fit$estimates[1,1]
#   dev_a    <- -2 * ll
#   dev_b    <- -2 * llnew
#   Pd       <- dev_a - dev_b
#   DIC      <- Pd + dev_a
#   
#   virus    <- dat$virus[1]
#   model    <- 
#     
#     res_DIC <- data.frame(log_lik = ll, DIC = DIC, virus, model)
#   
#   return(res_DIC)
#   
# }
# 
# extract_ll_DIC(MADV2012, MADV2012$fit1$fit, model = 'constant')
# extract_ll_DIC(MADV2012, MADV2012$fit2$fit, model = 'time-varying')
# extract_ll_DIC(VEEV2012, VEEV2012$fit1$fit, model = 'constant')
# extract_ll_DIC(MADV2012, MADV2012$fit2$fit, model = 'time-varying')



ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
print(ip, row.names=FALSE)

stan_version()
