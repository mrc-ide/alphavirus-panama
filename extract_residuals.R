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



plot_residuals <- function (model, res) 
{
  
  if (model == 'constant') {
    fit <- res$mod1$fit 
  }
  
  if (model == 'time-varying') {
    fit <- res$mod2$fit 
  }
  
  dat <- res$dat
  PlaceName <- dat$PlaceName[1]
  if (PlaceName== 'Pi-Pi')  {PlaceName = 'Pirries & Pijivasal'}
  if (PlaceName== 'Tamar')  {PlaceName = 'Tamarindo'}
  if (PlaceName== 'Merca')  {PlaceName = 'Mercadeo'}
  if (PlaceName== 'Real')   {PlaceName = 'El Real'}
  
  
  
  
  Npos_sim <- rstan::extract(fit, 'Npos_sim')[[1]]
  P_sim <- rstan::extract(fit, 'P_sim')[[1]]
  colnames(P_sim) <- dat$age_mean_f
  P_sim <- as.data.frame(P_sim)
  P_sim$iteration <- seq_along(P_sim[,1])
  P_sim <- P_sim %>% 
    melt(id.vars='iteration') %>%
    mutate(source='simulated') %>% 
    mutate(variable=as.character(variable))
  dat$source <- 'actual'
  dat$iteration <- NA
  dat$prev_obs <- dat$counts / dat$total
  
  
  dat1 <- dat %>% 
    mutate(variable=as.character(age_mean_f)) %>%
    select(variable, source, iteration, prev_obs, age_min, age_max) %>%
    rename(value=prev_obs) %>% 
    bind_rows(P_sim) %>% 
    mutate(variable=as.numeric(variable))
  
  dat2 <- dat1 %>% 
    filter(source=="actual") %>% 
    group_by(variable) %>% 
    summarise(value=mean(value))
  
  dat3 <- dat1 %>% 
    filter(source=="simulated") %>% 
    left_join(dat2, by='variable') %>% 
    mutate(resid=value.y-value.x) 
  
  dat$age_minc <- as.character(dat$age_min)
  dat$age_maxc <- as.character(dat$age_max)
  dat$age_minc[dat$age_min <10] <- paste0(0, dat$age_min[dat$age_min <10])
  dat$age_maxc[dat$age_max <10] <- paste0(0, dat$age_max[dat$age_max <10])
  dat$age_class <- paste0(dat$age_minc, '-', dat$age_maxc)
  dat1$pobs_lower <- NA
  dat1$pobs_upper <- NA
  
  
  dat1$age_class <-NA
  dat1$pobs_lower <-NA
  dat1$pobs_upper <-NA
  dat1$total <-NA
  
  for (i in unique(dat1$variable)) {
    dat1$age_class[dat1$variable==i] <- dat$age_class[dat$age_mean_f==i]
    dat1$pobs_lower[dat1$variable==i] <- dat$prev_obs_lower[dat$age_mean_f==i]
    dat1$pobs_upper[dat1$variable==i] <- dat$prev_obs_upper[dat$age_mean_f==i]
    dat1$total[dat1$variable==i]      <- dat$total[dat$age_mean_f==i]
  }
  
  
  
  g2 <-
    ggplot(data=filter(dat1, source=="simulated"), aes(x=age_class, y=value)) +
    geom_boxplot(aes(group=age_class), fill = 'gray', alpha = .4, outlier.shape = NA) +
    geom_errorbar(data=filter(dat1, source=="actual"), 
                  aes(x=age_class, ymin = pobs_lower, ymax = pobs_upper), 
                  colour ='red', width = .02) +
    geom_point(data=filter(dat1, source=="actual"), 
               aes(x=age_class, y=value, size = total), pch = 21, fill="red") +
    xlab("Age, years") +
    ylab("Seropositivity") +
    theme_bw() + coord_cartesian(ylim = c(0,1)) +
    theme(axis.text.x = element_text(angle=45)) +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(size=10)) +
    ggtitle (paste(PlaceName, model))
  
  
  g3 <- 
    dat3 %>% 
    select(variable, iteration, resid) %>%
    ggplot(aes(x=as.factor(variable), y=resid, group=as.factor(iteration))) +
    geom_line(alpha=0.4) +
    geom_hline(yintercept = 0, linetype=2, colour="orange") + xlab('') +
    theme_bw () + 
    theme(axis.text.x = element_blank()) +
    theme(plot.title = element_text(size=10)) +
    ggtitle (paste(PlaceName, model))
  
  
  pt <- plot_grid (g2, g3, nrow = 1)
  
  return (pt)
  
}

#================== Residuals VEEV


VP1_05 <- readRDS('res_final_10000/VEEV 2012 Pi-Pi 005.RDS')
VP2_05 <- readRDS('res_final_10000/VEEV 2012 Merca 005.RDS')
VP3_05 <- readRDS('res_final_10000/VEEV 2012 Tamar 005.RDS')
VP4_05 <- readRDS('res_final_10000/VEEV 2012 Real 005.RDS')
VP5_05 <- readRDS('res_final_10000/VEEV 2012 Aruza 005.RDS')
VP6_05 <- readRDS('res_final_10000/VEEV 2017 Mogue 005.RDS')



VP1_A <- plot_residuals(model = 'constant', VP1_05)
VP1_B <- plot_residuals(model = 'time-varying', VP1_05)
VP2_A <- plot_residuals(model = 'constant', VP2_05)
VP2_B <- plot_residuals(model = 'time-varying', VP2_05)
VP3_A <- plot_residuals(model = 'constant', VP3_05)
VP3_B <- plot_residuals(model = 'time-varying', VP3_05)
VP4_A <- plot_residuals(model = 'constant', VP4_05)
VP4_B <- plot_residuals(model = 'time-varying', VP4_05)
VP5_A <- plot_residuals(model = 'constant', VP5_05)
VP5_B <- plot_residuals(model = 'time-varying', VP5_05)
VP6_A <- plot_residuals(model = 'constant', VP6_05)
VP6_B <- plot_residuals(model = 'time-varying', VP6_05)


png('VEEP_residuals.png', width = 300 *3, height = 480 * 3)
plot_grid (VP1_A, VP1_B,
           VP2_A, VP2_B,
           VP3_A, VP3_B,
           VP4_A, VP4_B,
           VP5_A, VP5_B,
           VP6_A, VP6_B,
           nrow = 6
           )
dev.off()
  

#================== Residuals MADV 

MP1_05 <- readRDS('res_final_10000/MADV 2012 Pi-Pi 005.RDS')
MP2_05 <- readRDS('res_final_10000/MADV 2012 Merca 005.RDS')
MP3_05 <- readRDS('res_final_10000/MADV 2012 Tamar 005.RDS')
MP4_05 <- readRDS('res_final_10000/MADV 2012 Real 005.RDS')
MP5_05 <- readRDS('res_final_10000/MADV 2012 Aruza 005.RDS')
MP6_05 <- readRDS('res_final_10000/MADV 2017 Mogue 005.RDS')




MP1_A <- plot_residuals(model = 'constant', MP1_05)
MP1_B <- plot_residuals(model = 'time-varying', MP1_05)
MP2_A <- plot_residuals(model = 'constant', MP2_05)
MP2_B <- plot_residuals(model = 'time-varying', MP2_05)
MP3_A <- plot_residuals(model = 'constant', MP3_05)
MP3_B <- plot_residuals(model = 'time-varying', MP3_05)
MP4_A <- plot_residuals(model = 'constant', MP4_05)
MP4_B <- plot_residuals(model = 'time-varying', MP4_05)
MP5_A <- plot_residuals(model = 'constant', MP5_05)
MP5_B <- plot_residuals(model = 'time-varying', MP5_05)
MP6_A <- plot_residuals(model = 'constant', MP6_05)
MP6_B <- plot_residuals(model = 'time-varying', MP6_05)


png('MADV_residuals.png', width = 300 *3, height = 480 * 3)
plot_grid (MP1_A, MP1_B,
           MP2_A, MP2_B,
           MP3_A, MP3_B,
           MP4_A, MP4_B,
           MP5_A, MP5_B,
           MP6_A, MP6_B,
           nrow = 6
)
dev.off()



