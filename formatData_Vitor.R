
# Variabkles

library(haven)
library(dplyr)
library(reshape2)
library(Hmisc)

rm(list= ls())


formatingData <- function(dat, virus){


  if (virus == 'MADV') {dat$result <- dat$eeepos}
  if (virus == 'VEEV') {dat$result <- dat$veepos}
  
  dat <- filter(dat, !is.na(result))
  dat$sample <- 1
  
  
  # 10 years age class
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
  # plot(dd$age_class_new, dd$result/dd$sample, type = 'o')
  d_10y <- dd ;rm(dd1, dd2, dd)
  
  
  # Single age class
  dd1 <- aggregate(result ~ age, data = dat, FUN = sum)
  dd2 <- aggregate(sample ~ age, data = dat, FUN = sum)
  dd <- merge(dd1, dd2, by ='age')
  # plot(dd$age, dd$result/dd$sample, type = 'p')
  
  d_01y <- dd ;rm(dd1, dd2, dd)
  
  #Two Years
  dat$age_class_new <- NA
  dat$age_class_new[dat$age  < 3] <- 2
  dat$age_class_new[dat$age >= 3 & dat$age < 5] <- 4
  dat$age_class_new[dat$age >= 5 & dat$age < 7] <- 6
  dat$age_class_new[dat$age >= 7 & dat$age < 9] <- 8
  dat$age_class_new[dat$age >= 9 & dat$age < 11] <-10
  dat$age_class_new[dat$age >= 11 & dat$age < 13] <- 12
  dat$age_class_new[dat$age >= 13 & dat$age < 15] <- 14
  dat$age_class_new[dat$age >= 15 & dat$age < 17] <- 16
  dat$age_class_new[dat$age >= 17 & dat$age < 19] <- 18
  dat$age_class_new[dat$age >= 19 & dat$age < 21] <- 20
  dat$age_class_new[dat$age >= 21 & dat$age < 23] <- 22
  dat$age_class_new[dat$age >= 23 & dat$age < 25] <- 24
  dat$age_class_new[dat$age >= 25 & dat$age < 27] <- 26
  dat$age_class_new[dat$age >= 27 & dat$age < 29] <- 28
  dat$age_class_new[dat$age >= 29 & dat$age < 31] <- 30
  dat$age_class_new[dat$age >= 31 & dat$age < 33] <- 32
  dat$age_class_new[dat$age >= 33 & dat$age < 35] <- 34
  dat$age_class_new[dat$age >= 35 & dat$age < 37] <- 36
  dat$age_class_new[dat$age >= 37 & dat$age < 39] <- 38
  dat$age_class_new[dat$age >= 39 & dat$age < 41] <- 40
  dat$age_class_new[dat$age >= 41 & dat$age < 43] <- 42
  dat$age_class_new[dat$age >= 43 & dat$age < 45] <- 44
  dat$age_class_new[dat$age >= 45 & dat$age < 47] <- 46
  dat$age_class_new[dat$age >= 47 & dat$age < 49] <- 48
  dat$age_class_new[dat$age >= 50] <- 50
  
  dd1 <- aggregate(result ~ age_class_new, data = dat, FUN = sum)
  dd2 <- aggregate(sample ~ age_class_new, data = dat, FUN = sum)
  dd <- merge(dd1, dd2, by ='age_class_new')
  d_02y <- dd ;rm(dd1, dd2, dd)
  
  
  
  #Five Years
  dat$age_class_new <- NA
  dat$age_class_new[dat$age  < 5] <- 3
  dat$age_class_new[dat$age >= 5  & dat$age < 10] <- 7
  dat$age_class_new[dat$age >= 10 & dat$age < 15] <- 12
  dat$age_class_new[dat$age >= 15 & dat$age < 20] <- 17
  dat$age_class_new[dat$age >= 20 & dat$age < 25] <- 22
  dat$age_class_new[dat$age >= 25 & dat$age < 30] <- 27
  dat$age_class_new[dat$age >= 30 & dat$age < 35] <- 32
  dat$age_class_new[dat$age >= 35 & dat$age < 40] <- 37
  dat$age_class_new[dat$age >= 40 & dat$age < 45] <- 42
  dat$age_class_new[dat$age >= 45 & dat$age < 50] <- 47
    dat$age_class_new[dat$age >= 50] <- 52
  
  dd1 <- aggregate(result ~ age_class_new, data = dat, FUN = sum)
  dd2 <- aggregate(sample ~ age_class_new, data = dat, FUN = sum)
  dd <- merge(dd1, dd2, by ='age_class_new')
  d_05y <- dd ;rm(dd1, dd2, dd)
  
  
  variables <- c("dataset_id", "age_min","age_max", "result", "sample" )
  d_01y$age_min <- d_01y$age
  d_01y$age_max <- d_01y$age
  d_01y$dataset_id    <- '001'
  d_01y <- d_01y[,variables]
  
  d_02y$age_min <- d_02y$age_class_new-1
  d_02y$age_max <- d_02y$age_class_new
  d_02y$dataset_id    <- '002'
  d_02y <- d_02y[,variables]
  
  
  
  d_05y$age_min <- d_05y$age_class_new-2
  d_05y$age_max[1]<- d_05y$age_class_new[1] + 2
  d_05y$age_max[2:length(d_05y$age_min)]<- d_05y$age_class_new[2:length(d_05y$age_min)] + 3
  d_05y$dataset_id    <- '005'
  d_05y <- d_05y[,variables]
  
  d_10y$age_min <- d_10y$age_class_new-4
  d_10y$age_max <- d_10y$age_class_new+5
  d_10y$dataset_id    <- '010'
  d_10y <- d_10y[,variables]
  
  dd <- rbind(d_01y, d_02y, d_05y, d_10y)
  

  dd$country   <- 'Panama'
  dd$Ref       <- 'Vitor et al. 2016'
  dd$Admin0    <- 'PAN'
  dd$zone      <- 'rural'
  dd$ISO       <- 'PAN'
  dd$Total     <- dd$sample
  dd$total     <- dd$sample
  dd$n_pos     <- dd$result
  dd$Neg       <- dd$total - dd$n_pos
  dd$neg       <- dd$total - dd$n_pos
  dd$netspec   <- 1
  dd$netsens   <- 1 
  dd$age_mean_f <-  floor((dd$age_min + dd$age_max)/2)
  dd$tSur   <- 2012
  dd$tSur    <- 2012
  dd$pos <- dd$n_pos
  dd$survey <- dd$dataset_id
  dd$PlaceName <- virus
  
  dd$birth_year <- dd$tSur - dd$age_min
  dd$age_mean_f <- dd$age_min
  dd$tsur <- dd$tSur
  
  conf <- data.frame(Hmisc::binconf(dd$n_pos, dd$total,method="exact"))
  prev <- cbind(dd, conf)
  prev$country <- prev$Admin0
  dd$prev_obs    <- dd$pos/dd$total
  dd$birth_year  <- dd$tSur - dd$age_mean_f
  dd$netsens <- 1
  dd$netspec <- 1
  dd$prev_obs_lower <- prev$Lower
  dd$prev_obs_upper <- prev$Upper
  dd$virus <- virus
  
  dd$dataset_id <- paste(dd$dataset_id, virus)
  dd$survey <- dd$dataset_id
  
  return(dd)
  
}

  

dat  <- read_excel("data/HH simple 103014_MADV_VEEV_Zulma.xlsx")

MADV2012 <- formatingData(dat, 'MADV')
VEEV2012 <- formatingData(dat, 'VEEV')

datf <- rbind(MADV2012, VEEV2012)
saveRDS(datf, 'data/Panama2012.RDS')
