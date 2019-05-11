
# Variabkles

library(haven)
library(dplyr)
library(reshape2)

rm(list= ls())

dat <- read_dta("data/Mogue_2017_vf_243n_copia.dta")

dat$edad

dat$age <- as.numeric(dat$edad)
dat$result <- dat$una
dat$age_class <- dat$edadcat
dat$sample <- 1
hist(dat$age)

dat$age_class_new <- NA
dat$age_class_new[dat$age  < 10] <- 10/2
dat$age_class_new[dat$age >= 10 & dat$age < 20] <- 15
dat$age_class_new[dat$age >= 20 & dat$age < 30] <- 25
dat$age_class_new[dat$age >= 30 & dat$age < 40] <- 35
dat$age_class_new[dat$age >= 40 & dat$age < 50] <- 45
dat$age_class_new[dat$age >= 50] <- 55


dd1 <- aggregate(result ~ age_class_new, data = dat, FUN = sum)
dd2 <- aggregate(sample ~ age_class_new, data = dat, FUN = sum)
dd <- merge(dd1, dd2, by ='age_class_new')

plot(dd$age_class_new, dd$result, type = 'o')

dd$age_min <- NA
dd$age_min[dd$age_class_new == 5] <- 1
dd$age_min[dd$age_class_new == 15] <- 10
dd$age_min[dd$age_class_new == 25] <- 20
dd$age_min[dd$age_class_new == 35] <- 30
dd$age_min[dd$age_class_new == 45] <- 40
dd$age_min[dd$age_class_new == 55] <- 50



dd$age_max <- NA
dd$age_max[dd$age_class_new == 5] <- 10
dd$age_max[dd$age_class_new == 15] <- 20
dd$age_max[dd$age_class_new == 25] <- 30
dd$age_max[dd$age_class_new == 35] <- 40
dd$age_max[dd$age_class_new == 45] <- 50
dd$age_max[dd$age_class_new == 55] <- 60 # BE aware this is NOT the max age!!


dd$survey    <- 'XX1'
dd$country   <- 'Panama'
dd$Code      <- '001'
dd$Ref       <- 'Carrasco et al. 2018'
dd$Admin0    <- 'PAN'
dd$PlaceName <- 'XXX'
dd$zone      <- 'rural'
dd$Total     <- dd$sample
dd$total     <- dd$sample
dd$n_pos     <- dd$result
dd$Neg       <- dd$total - dd$n_pos
dd$neg       <- dd$total - dd$n_pos
dd$dataset_id<- '0011'
dd$age_mean  <- dd$age_class_new
dd$netspec   <- 1
dd$netspec   <- 1 
dd$age_mean_f <-  dd$age_mean
dd$tSur   <- 2017
dd$tSur    <- 2017


saveRDS(dd, 'data/Panama.RDS')
