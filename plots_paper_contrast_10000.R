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



plot_my_results <- function(res, 
                            max_lambda=0.3  , 
                            max_prev=1 , 
                            horizontal  = TRUE
                            ){
  
  
  
  dat    <- res$dat
  mod1    <- res$mod1
  mod2    <- res$mod2
  
  virus     <- dat$virus[1]
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  tsur      <- dat$tsur[1]
  extract_foi_summary <- function (mod, dat)
  {

    fit <- mod$fit
    N       <- sum(dat$total)
    RealYexpo <- mod$RealYexpo
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
  
  foi_years1 <- extract_foi_summary(mod1, dat) %>% mutate(model = 'constant')
  foi_years2 <- extract_foi_summary(mod2, dat) %>% mutate(model = 'time-varying')
  
  foi_overt <- rbind(foi_years1, foi_years2)
  size_text <- 25
  space_years <- 10
  breaks_x <- round(seq(1960, max(foi_overt$year), by = space_years),0) 

  g1 <- 
    foi_overt %>% filter(year > 1960) %>%
    ggplot(aes(x=year, y=median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = .4) +
    geom_line(aes(y = median, colour = model)) +
    scale_fill_brewer("FOI", palette = "Set1") +
    scale_color_brewer("FOI", palette = "Set1") +
    coord_cartesian(ylim = c(0, max_lambda), xlim = c(1960, 2017)) + 
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
  
  

  
  
  get_lambda_yexpo_expanded <- function (fit, 
                                         model,
                                         min_age = 1, 
                                         max_age = 60)
  {
    
    foi_fit <- data.frame(rstan::extract(fit, 'foi')[[1]])
    lambdaYexpo <- as.matrix(foi_fit)
    iterf = NROW(foi_fit)
    lambdaYexpoExpanded <- matrix(NA, nrow = iterf, ncol= max_age)
    Ymax <- max_age
    
    if(model == 'time_varying') {


    lambdaYexpoExpanded[,(Ymax-9):Ymax] <- lambdaYexpo[,5]
    lambdaYexpoExpanded[,(Ymax-19):(Ymax-10)] <- lambdaYexpo[,4]
    lambdaYexpoExpanded[,(Ymax-29):(Ymax-20)] <- lambdaYexpo[,3]
    lambdaYexpoExpanded[,(Ymax-39):(Ymax-30)] <- lambdaYexpo[,2]
    lambdaYexpoExpanded[,1:(Ymax-40)] <- lambdaYexpo[,1]
    } else {
      lambdaYexpoExpanded[,1:Ymax] <- lambdaYexpo[,1]
    }
    
    
    return(lambdaYexpoExpanded)
    
  }
  
  

  
  min_age <-1
  max_age <-60
  yexpo_expanded <- min_age:max_age
  lambdaYexpoExpanded1 <- get_lambda_yexpo_expanded(mod1$fit, 'constant')
  lambdaYexpoExpanded2 <- get_lambda_yexpo_expanded(mod2$fit, 'time_varying')

  ExposureMatrixFit <- get_exposure_matrix_expanded(yexpo_expanded, min_age ,max_age)
  

  iterf <- NROW(lambdaYexpoExpanded2)
  PrevP1 <- matrix(NA, nrow = iterf, ncol = max(yexpo_expanded))
  PrevP2 <- matrix(NA, nrow = iterf, ncol = max(yexpo_expanded))
  for (i in 1:iterf)
  {
    PrevP1[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpoExpanded1[i,])
    PrevP2[i,] <- 1 - exp(- ExposureMatrixFit %*% lambdaYexpoExpanded2[i,])
    
  } 
  
  PPP1 <- matrix(NA, ncol = 3, nrow = max_age)
  PPP2 <- PPP1
  for (j in seq_along(min_age:max_age)) 
  {
    PPP1[j,] <- quantile(PrevP1[,j], c(.025, .5, .975))
    PPP2[j,] <- quantile(PrevP2[,j], c(.025, .5, .975))
  }
  
  
  PPP1        <- as.data.frame(PPP1) %>% mutate(age = min_age:max_age,
                                              model = 'constant')
  PPP2        <- as.data.frame(PPP2) %>% mutate(age = min_age:max_age,
                                               model = 'time_varying')
  
  PPP_both <- rbind(PPP1, PPP2)
  colnames(PPP_both) <- c("L", "M", "U", "age", 'model')
  
  
  g2 <-
    ggplot() +
    geom_errorbar(data=dat, 
                  aes(x=age_mean_f, ymin=prev_obs_lower, ymax = prev_obs_upper),
                  color = 'grey', width = .1) +
    geom_ribbon(data=PPP_both, aes(x=age, ymin = L, ymax = U, fill = model), alpha = 0.3) +
    geom_line  (data=PPP_both, aes(x=age, y= M, color = model)) + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    geom_point (data=dat, aes(x=age_mean_f, y=prev_obs, size = total), 
                pch = 21,fill = 'white', color = 'black') +
    xlab("Age, years") +
    ylab("seropositivity") +
    theme_bw() + coord_cartesian(ylim = c(0,max_prev), xlim = c(1,60)) +
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



VP1_05 <- readRDS('res_final_10000/VEEV 2012 Pi-Pi 005.RDS')
VP2_05 <- readRDS('res_final_10000/VEEV 2012 Merca 005.RDS')
VP3_05 <- readRDS('res_final_10000/VEEV 2012 Tamar 005.RDS')
VP4_05 <- readRDS('res_final_10000/VEEV 2012 Real 005.RDS')
VP5_05 <- readRDS('res_final_10000/VEEV 2012 Aruza 005.RDS')
VP6_05 <- readRDS('res_final_10000/VEEV 2017 Mogue 005.RDS')

MP1_05 <- readRDS('res_final_10000/MADV 2012 Pi-Pi 005.RDS')
MP2_05 <- readRDS('res_final_10000/MADV 2012 Merca 005.RDS')
MP3_05 <- readRDS('res_final_10000/MADV 2012 Tamar 005.RDS')
MP4_05 <- readRDS('res_final_10000/MADV 2012 Real 005.RDS')
MP5_05 <- readRDS('res_final_10000/MADV 2012 Aruza 005.RDS')
MP6_05 <- readRDS('res_final_10000/MADV 2017 Mogue 005.RDS')


plot_my_results(VP1_05, max_lambda = 0.3, max_prev = 1, horizontal = FALSE)


title1=textGrob("VEEV", gp=gpar(fontface="bold"))
title2=textGrob("MADV", gp=gpar(fontface="bold"))



png("res_final_10000/foi_VEEV_5y.png",
    width = 350 *2, height = 350 *6)

title1=textGrob("VEEV", gp=gpar(fontface="bold"))
title2=textGrob("MADV", gp=gpar(fontface="bold"))

gridExtra::grid.arrange(
  plot_my_results(VP1_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(VP2_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(VP3_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(VP4_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(VP5_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(VP6_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  nrow = 6,
  top= title1
)

dev.off()


png("res_final_10000/foi_MADV_5y.png",
    width = 350 *2, height = 350 *6)

gridExtra::grid.arrange(
  plot_my_results(MP1_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(MP2_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(MP3_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(MP4_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(MP5_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  plot_my_results(MP6_05, max_lambda = 0.2, max_prev = 1, horizontal = TRUE),
  nrow = 6,
  top= textGrob("MADV 5y", gp=gpar(fontface="bold"))
)

dev.off()










get_DIC <- function(mod,model_name, dat)
{
  fit <- mod$fit
  loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
  yexpo   <- make_yexpo(dat)
  yexpo <- yexpo[-length(yexpo)]
  ExposureMatrixFit <- get_exposure_matrix(dat, yexpo)
  foi <- data.frame(rstan::extract(fit, 'foi')[[1]])
  # lambdaTM <- (as.numeric(sapply(foi, function(i) c(quantile(i, 0.5)))))
  lambdaTM <- (as.numeric(sapply(foi, function(i) c(mean(i)))))
  
  if(model_name == 'time_varying'){
    lambdaTM_expanded <- matrix(NA, nrow = 1, ncol=max(yexpo))
    Ymax <- max(yexpo)
    lambdaTM_expanded[(Ymax-9):Ymax] <- lambdaTM[5]
    lambdaTM_expanded[(Ymax-19):(Ymax-10)] <- lambdaTM[4]
    lambdaTM_expanded[(Ymax-29):(Ymax-20)] <- lambdaTM[3]
    lambdaTM_expanded[(Ymax-39):(Ymax-30)] <- lambdaTM[2]
    lambdaTM_expanded[1:(Ymax-40)] <- lambdaTM[1]
    
    lambdaTM <- lambdaTM_expanded
    
  }
  
  
  PrevM          <-  1-(exp(- ExposureMatrixFit %*% as.numeric(lambdaTM)))
  ll_best_params <-  sum(dbinom(as.matrix(dat$pos), as.matrix(dat$total), PrevM, log = TRUE))
  
  ll_mean <- loo_fit$estimates[1,1]
  Dbar =  -2 * ll_mean          #------> Dbar is the posterior mean of the deviance 
  Dhat =  -2 * ll_best_params   #------> Dhat point estimate of the deviance obtained by substituting in the posterior means 
  pD    =  Dbar - Dhat          #------> Effective number of parameters
  DIC   = Dbar + pD  
  
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  virus    <- dat$virus[1]
  res_DIC <- data.frame(Place = PlaceName, 
                        virus = dat$virus[1], 
                        sample_size = sum(dat$total), 
                        log_lik = ll_mean, DIC = DIC,
                        n_eff   = pD,
                        model   = model_name)
  
  return(res_DIC)
  
}




results <- 
rbind(
get_DIC(VP1_05$mod1, model_name ='constant',     VP1_05$dat),
get_DIC(VP1_05$mod2, model_name ='time_varying', VP1_05$dat),
get_DIC(VP2_05$mod1, model_name ='constant',     VP2_05$dat),
get_DIC(VP2_05$mod2, model_name ='time_varying', VP2_05$dat),
get_DIC(VP3_05$mod1, model_name ='constant',     VP3_05$dat),
get_DIC(VP3_05$mod2, model_name ='time_varying', VP3_05$dat),
get_DIC(VP4_05$mod1, model_name ='constant',     VP4_05$dat),
get_DIC(VP4_05$mod2, model_name ='time_varying', VP4_05$dat),
get_DIC(VP5_05$mod1, model_name ='constant',     VP5_05$dat),
get_DIC(VP5_05$mod2, model_name ='time_varying', VP5_05$dat),
get_DIC(VP6_05$mod1, model_name ='constant',     VP6_05$dat),
get_DIC(VP6_05$mod2, model_name ='time_varying', VP6_05$dat),



get_DIC(MP1_05$mod1, model_name ='constant',     MP1_05$dat),
get_DIC(MP1_05$mod2, model_name ='time_varying', MP1_05$dat),
get_DIC(MP2_05$mod1, model_name ='constant',     MP2_05$dat),
get_DIC(MP2_05$mod2, model_name ='time_varying', MP2_05$dat),
get_DIC(MP3_05$mod1, model_name ='constant',     MP3_05$dat),
get_DIC(MP3_05$mod2, model_name ='time_varying', MP3_05$dat),
get_DIC(MP4_05$mod1, model_name ='constant',     MP4_05$dat),
get_DIC(MP4_05$mod2, model_name ='time_varying', MP4_05$dat),
get_DIC(MP5_05$mod1, model_name ='constant',     MP5_05$dat),
get_DIC(MP5_05$mod2, model_name ='time_varying', MP5_05$dat),
get_DIC(MP6_05$mod1, model_name ='constant',     MP6_05$dat),
get_DIC(MP6_05$mod2, model_name ='time_varying', MP6_05$dat)

)


write_csv(results, 'res_final_10000/DIC.csv')



fit <- MP5_05$mod2$fit
traceplot(fit, 'lp__')


diagnostic_plots <- function (mod, virus)
{
  
 dat <- mod$dat
  
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  fit1 <- mod$mod1$fit
  fit2 <- mod$mod2$fit
  

  title_plot1 <- paste(virus, PlaceName, '(constant model)')
  title_plot2 <- paste(virus, PlaceName, '(time_varying model)')
  
  p1 <- traceplot(fit1, 'lp__', inc_warmup = TRUE, size = 0.3) + ggtitle(title_plot1) 
  p2 <- traceplot(fit2, 'lp__', inc_warmup = TRUE, size = 0.3) + ggtitle(title_plot2)
  
  pars1 <- traceplot(fit1, 'foi[1]', inc_warmup = TRUE, size = 0.3, nrow = 5) +
    ggtitle('params') +  ylab ('foi')
  pars2 <- traceplot(fit2, 'foi'   , inc_warmup = TRUE, size = 0.3, nrow = 2)  +
    ggtitle('params')
  
  rhats1  <- rhat(fit1, 'foi[1]')
  rhats2  <- rhat(fit2, 'foi')
 
  color_scheme_set("viridisC")
  rh1 <- mcmc_rhat(rhats1) + ggtitle('R^')
  rh2 <- mcmc_rhat(rhats2) + ggtitle('R^')
  
  grid.arrange(p1, pars1, rh1, 
               p2, pars2, rh2, 
               nrow = 2, widths = c(1, 1.5, 0.5))
 
  

}

pdf('res_final_10000/diagnostics_VEEV.pdf',  height = 8, width = 20)
diagnostic_plots(VP1_05, 'VEEV')
diagnostic_plots(VP2_05, 'VEEV')
diagnostic_plots(VP3_05, 'VEEV')
diagnostic_plots(VP4_05, 'VEEV')
diagnostic_plots(VP5_05, 'VEEV')
diagnostic_plots(VP6_05, 'VEEV')

dev.off()


pdf('res_final_10000/diagnostics_MADV.pdf',  height = 8, width = 20)
diagnostic_plots(MP1_05, 'MADV')
diagnostic_plots(MP2_05, 'MADV')
diagnostic_plots(MP3_05, 'MADV')
diagnostic_plots(MP4_05, 'MADV')
diagnostic_plots(MP5_05, 'MADV')
diagnostic_plots(MP6_05, 'MADV')

dev.off()
