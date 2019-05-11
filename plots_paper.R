rm(list = ls())

library(loo)
library(dplyr)
library(reshape2)
library(tidyverse)
library(cowplot)
source('fun/additional_functions.R')


plot_my_results <- function (res, model, max_lambda) {
  
  fit    <- model$fit
  dat    <- res$dat
  RealYexpo <- model$RealYexpo
  
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
  
  foi <- extract_foi_summary(res, model, fit, dat, RealYexpo) 
  size_text <- 15
  space_years <- 1 
  if(length(model$yexpo)>15 ) {space_years <- 2}
  if(length(model$yexpo)>30 ) {space_years <- 3}
  if(length(model$yexpo)>50 ) {space_years <- 4}
  
  
  breaks_x <- round(seq(1960, max(foi$year), by = space_years),0) 
  breaks_x <- round(seq(min(breaks_x) + (max(foi$year) - max(breaks_x)) , max(foi$year), by = space_years),0)
  g1 <- 
    foi %>% filter(year > 1960) %>%
    ggplot(aes(x=year, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
    geom_line() +
    scale_fill_brewer("FOI", palette = "Dark2") +
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
    scale_x_continuous(breaks = breaks_x)  
  
  
  foi_chain <- data.frame(rstan::extract(fit, 'foi')[[1]])
  iterf = NROW(foi_chain)
  lambdaYexpo <- as.matrix(foi_chain)
  ExposureMatrixFit <- get_exposure_matrix(dat)[, -1]
  
  PrevP <- matrix(NA, nrow = iterf, ncol = length(dat$age_mean_f))
  
  for (i in 1:iterf)
  {
    PrevP[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpo[i,])
  } 
  
  PPP <- matrix(NA, ncol = 3, nrow = length(dat$age_mean_f))
  for (j in seq_along(dat$age_mean_f)) 
  {
    PPP[j,] <- quantile(PrevP[,j], c(.025, .5, .975))
  }
  
  PPP        <- as.data.frame(PPP)
  PPP$age    <- dat$age_mean_f
  names(PPP) <- c("L", "M", "U", "age")
  
  
  
  g2 <-
    ggplot() +
    geom_errorbar(data=dat, 
                  aes(x=age_mean_f, ymin=prev_obs_lower, ymax = prev_obs_upper),
                  color = 'grey', width = .1) +
    geom_point (data=dat, aes(x=age_mean_f, y=prev_obs)) +
    geom_ribbon(data=PPP, aes(x=age, ymin = L, ymax = U), alpha = .2) +
    geom_line  (data=PPP, aes(x=age, y= M)) + 
    xlab("Age, years") +
    ylab("Seropositivity") +
    theme_bw() + coord_cartesian(ylim = c(0,1)) +
    theme(axis.text.x = element_text(angle= 0)) +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(size=10)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p_res <- plot_grid(g1, g2,
                     nrow = 1, 
                     rel_heights = c(1, 1)
  )
  
  
  return(p_res)
}


MADV2012 <- readRDS('res/MADV 2012 005.RDS')
VEEV2012 <- readRDS('res/VEEV 2012 005.RDS')


MADV2017 <- readRDS('res/MADV 2012 005.RDS')
UNAV2017 <- readRDS('res/UNA 2017 005.RDS')
VEEV2017 <- readRDS('res/VEEV 2017 010.RDS')


res <- VEEV2012
VEEV2012_cons <- plot_my_results(res, model = res$fit1, max_lambda = 0.1)
VEEV2012_vary <- plot_my_results(res, model = res$fit2, max_lambda = 0.1)

res <- VEEV2017
VEEV2017_cons <- plot_my_results(res, model = res$fit1, max_lambda = 0.1)
VEEV2017_vary <- plot_my_results(res, model = res$fit2, max_lambda = 0.1)

pdf('res/plots3.pdf')
plot_grid(VEEV2012_cons, VEEV2012_vary, nrow = 2)
plot_grid(VEEV2017_cons, VEEV2017_vary, nrow = 2)
dev.off()


