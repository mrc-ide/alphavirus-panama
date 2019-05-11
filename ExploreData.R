
# Variabkles
rm(list= ls())

library(haven)
library(dplyr)
library(reshape2)
library(Hmisc)
library(ggplot2)



dat <- read_dta("data/Mogue_2017_vf_243n_copia.dta") 
dat$veevprnt <- as.numeric(dat$veevprnt)


una <- data.frame(id = dat$idmuestra, 
                  age = dat$edad, 
                  results = as.character(dat$una), 
                  screen = dat$unascreen, 
                  prnt = dat$unaprnt)


una <- data.frame(id = dat$idmuestra, 
                  age = dat$edad, 
                  results = as.character(dat$veevuna), 
                  screen = dat$veevprnt_2, 
                  prnt = dat$unaprnt)

virus <- una
ggplot(data = virus, aes(x = (prnt))) +
  # geom_histogram(binwidth = 20) +
  geom_density(aes(colour = results, y=..density..)) +
  # facet_wrap(~una) +
  coord_cartesian(xlim = c(0,20))
  


rm(list=ls())
dat2  <- read_excel("data/HH simple 103014_MADV_VEEV_Zulma.xlsx")
igm <- read_excel("data/ELISA_IgM_Mogue_2017.xlsx")
igm <- read_excel("data/ELISA_IgM_Mogue_2017.xlsx")
mogue17 <- read_dta("data/Mogue_2017_vf_243n_copia.dta")
prntmogue <- read_dta( 'data/prnt_titers_edad_Mogue_cross_sectional .dta')


plot(prntmogue$edad, prntmogue$mayprnt)
plot(prntmogue$edad, prntmogue$unaprnt)
