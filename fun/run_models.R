rm(list=ls())

source('fun/libraries.R')
source('fun/set_up_model_function.R')
set.seed(123)


dat <- readRDS('data/chik_serop_IND 1965(S095).RDS')
sm0 <- set_up_function(model      = 'm0',
                      duration   = NA,
                      n_iter     = 500,
                      burnin     = 100,
                      n_thin     = 1,
                      dat        = dat)

fit_m0 <- sampling(sm0$stan_file, data=sm0$stan_data,iter = sm0$n_iter)

source('fun/make_report_function.R')
source('fun/utility.R')
make_report_model(fit_m0, sm0)




