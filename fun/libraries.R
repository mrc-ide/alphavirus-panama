
library(rstan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(bayesplot)
library(ggthemr)
library(loo)
library(gridExtra)
library(grid)
library(gtable)
library(knitr)
library(kableExtra)
library(rstanarm)



rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
