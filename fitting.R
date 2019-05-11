
rm(list=ls())
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
source('fun/FitAndPlot_function.R')

# model0 <- stan_model("models/constant_foi.stan")
# model1 <- stan_model("models/continuous_foi.stan")
# model2 <- stan_model("./models/continuous_foi_unsmoothed.stan")
# model3 <- stan_model("./models/laplace_foi.stan")
# model4 <- stan_model("./models/continuous_foi_student_t.stan")
# model5 <- stan_model("./models/continuous_foi_student_t_age_dependent_faster.stan")
model6 <- stan_model("./models/continuous_foi_salje2016.stan")


# saveRDS(model0, 'models/model0.RDS')
# saveRDS(model1, 'models/model1.RDS')
# saveRDS(model2, 'models/model2.RDS')
# saveRDS(model3, 'models/model3.RDS')
# saveRDS(model4, 'models/model4.RDS')
# saveRDS(model5, 'models/model5.RDS')
# saveRDS(model6, 'models/model6.RDS')

# ---- Models 
model0 <-readRDS('models/model0.RDS')
model1 <-readRDS('models/model1.RDS')
model2 <-readRDS('models/model2.RDS')
model3 <-readRDS('models/model3.RDS')
model4 <-readRDS('models/model4.RDS')
# model5 <-readRDS('models/model5.RDS')
model6 <-readRDS('models/model6.RDS')

# ---- Datasets
dat0 <- readRDS('data/Panama.RDS') %>% mutate(counts = pos)
datasets <- unique(dat0$survey)

s = 1
dat <- filter(dat0, survey == datasets[s])

fit0  <- fFitModel(model0, dat)
fit1  <- fFitModel(model1, dat)
fit2  <- fFitModel(model2, dat)
fit3  <- fFitModel(model3, dat)
fit4  <- fFitModel(model4, dat)
# fit5  <- fFitModel(model5, dat)
# fit6  <- fFitModel(model6, dat)

pp0 <- fPlotModel(fit0, dat, 'model 0', 'constant foi')
pp1 <- fPlotModel(fit1, dat, 'model 1', 'continuous')
pp2 <- fPlotModel(fit2, dat, 'model 2', 'contin. unsmoothed')
pp3 <- fPlotModel(fit3, dat, 'model 3', 'laplace')
pp4 <- fPlotModel(fit4, dat, 'model 4', 'Student t ')
# pp5 <- fPlotModel(fit5, dat, 'model 5', 'Student t Age Dep')
# pp6 <- fPlotModel(fit6, dat, 'model 6', 'Salje')

res_survey <- list(dat  = dat,
                   fit0 = fit0,
                   fit1 = fit1,
                   fit2 = fit2,
                   fit3 = fit3,
                   fit4 = fit4,
                   # fit5 = fit5,
                   # fit6 = fit6,
                   pp0  = pp0,
                   pp1  = pp1,
                   pp2  = pp2,
                   pp3  = pp3,
                   pp4  = pp4
                   # pp5  = pp5,
                   # pp6  = pp6 
)

pdf(paste('survey', s, '.pdf'))
grid.arrange(pp0, pp1, pp2, pp3, pp4, nrow = 5)
dev.off()
