
# --------- Get exposure matrix
get_exposure_matrix <- function(dat)
{
  age_class <- dat$age_mean_f
  yexpo     <- seq_along(min(dat$birth_year):dat$tsur[1]) 
  ly  <- length(yexpo)
  exposure       <- matrix(0,nrow=length(age_class), ncol=ly)   
  for (k in 1:length(age_class)) exposure[k,(ly-age_class[k]+1):ly] <-  1 
  exposure_output <- exposure
  return(exposure_output)
  
}

make_yexpo <- function(dat) {
  yexpo <- (seq_along(min(dat$birth_year):dat$tsur[1]))
  }


get_posterior_summary <- function(results_chain)
{
  res <- sapply(results_chain, 
                function(i) c(quantile(i, c(0.5, 0.025, 0.975))))
  row.names(res) <- c('Median', 'Lower', 'Upper')
  return(res)
}


make_lambda_yexpo <- function(model,
                              nbreaks, 
                              params_to_fit,
                              yexpo,
                              duration_outbreak)
{
  
  
  ly <- length(yexpo)
  
  if (model == 'constant')
  {
    lambda_yexpo    <- rep(params_to_fit$lambdac, length(yexpo))
  }
  
  
  if (model == 'epidemic') {
    lambda_yexpo    <- rep(0.000000001, ly)
    if (nbreaks >= 1)
    {
      
      lambda1_x_duration <- rep(params_to_fit$lambda1, duration_outbreak)
      period_outbreak1   <- params_to_fit$year1:(params_to_fit$year1 + duration_outbreak - 1)
      lambda_yexpo [period_outbreak1] <- lambda1_x_duration
    }
    if (nbreaks >= 2)
    {
      lambda2_x_duration <- rep(params_to_fit$lambda2, duration_outbreak)
      period_outbreak2   <- params_to_fit$year2:(params_to_fit$year2 + duration_outbreak - 1)
      lambda_yexpo [period_outbreak2] <- lambda2_x_duration
    }
    
    if (nbreaks >= 3)
    {
      
      
      lambda3_x_duration <- rep(params_to_fit$lambda3, duration_outbreak)
      period_outbreak3   <- params_to_fit$year3:(params_to_fit$year3 + duration_outbreak - 1)
      lambda_yexpo [period_outbreak3] <- lambda3_x_duration
    }
  }
  
  if (model == 'inter-endemic') {
    if (nbreaks == 1) {
      period1 <- length(1: params_to_fit$year1)
      period2 <- length((params_to_fit$year1 + 1): ly)
      lambda1_x_period1 <- rep(params_to_fit$lambda1, period1)
      lambda2_x_period2 <- rep(params_to_fit$lambda2, period2)
      lambda_yexpo <- c(lambda1_x_period1, lambda2_x_period2)
    } 
    
    if (nbreaks == 2) {
      period1 <- length(1:params_to_fit$year1)
      period2 <- length((params_to_fit$year1 + 1): params_to_fit$year2)
      period3 <- length((params_to_fit$year2 + 1): ly)
      lambda1_x_period1 <- rep(params_to_fit$lambda1, period1)
      lambda2_x_period2 <- rep(params_to_fit$lambda2, period2)
      lambda3_x_period3 <- rep(params_to_fit$lambda3, period3)
      lambda_yexpo <- c(lambda1_x_period1, lambda2_x_period2, lambda3_x_period3)
    }
    
  }
  
  if(length(lambda_yexpo) > ly){
    print('lambda yexpo longer than yexpo')
    lambda_yexpo <- lambda_yexpo[1:ly] # CAREFULL!!! A really bad solution!!
    
  }
  
  
  return(lambda_yexpo)
  
}



obtain_prevalence_extended <- function(dato, exposure, ly, nbreaks, lambdaYexpo) {
  
  ly            <- length(yexpo)
  
  if(nbreaks == 0 & ly < 100) 
  {
    ly = 100
    lambdaYexpo <- matrix(lambdaYexpo, nrow = length(lambdaYexpo), ncol = ly )
  }
  
  new_dat       <- data.frame(age = 1:ly) 
  exposure_new  <- matrix(0,nrow=length(new_dat$age), ncol=ly)   
  for (k in 1:length(new_dat$age)) exposure_new[k,(ly - new_dat$age[k]+1):ly] <-  1 
  
  olders <- data.frame(age = (ly+1):99) 
  exposure_olders <- matrix(1,nrow=length(olders$age), ncol=ly)   
  
  exposure_total <- rbind(exposure_new, exposure_olders)
  new_age_clases <- c(new_dat$age,   olders$age)
  
  iterf <- nrow(lambdaYexpo)
  PrevPn <- matrix(NA, nrow = iterf, ncol = length(new_age_clases))
  for (i in 1:iterf){
    PrevPn[i,] <- 1 - exp( - exposure_total %*% lambdaYexpo[i,])}
  new_PPP <- matrix(NA, ncol = 3, nrow = length(new_age_clases))
  for (j in seq_along(new_age_clases)){
    new_PPP[j,] <- quantile(PrevPn[,j], c(.025, .5, .975))}
  new_PPP        <- as.data.frame(new_PPP)
  new_PPP$age    <- new_age_clases
  names(new_PPP) <- c("L", "M", "U", "age")
  
  return(new_PPP)
  
}

thin_chain <- function(results_chain, thin = 10)
{
  results_chain <- results_chain[seq(1, nrow(results_chain), thin),]
  return (results_chain)
}



#  Streaming Stan Fit file
cleanObject <- function(object, pars){
  pars <- c(pars,'lp__')
  nn <- paste0('^',pars,'(\\[|$)',collapse="|")
  ids <-  grep(nn,  object@sim$fnames_oi)
  ids.2 <- which(names(object@par_dims) %in% pars)
  for(i in 1:4){
    a <- attributes(object@sim$samples[[i]])
    x <- object@sim$samples[[i]][ids]
    for(j in c('names','inits','mean_pars'))
      a[[j]] <- a[[j]][ids]
    attributes(x) <- a
    object@sim$samples[[i]] <- x
  }
  object@par_dims <- object@par_dims[ids.2]
  object@sim$dims_oi <-   object@sim$dims_oi[ids.2]  
  object@sim$pars_oi<- object@sim$pars_oi[ids.2]
  object@sim$fnames_oi <-  object@sim$fnames_oi[ids]
  object@sim$n_flatnames <- length(object@sim$fnames_oi)
  
  return(object)
  
}


sub_sample <- function(stanfit, n, keep_warmup = TRUE) {
  sim <- stanfit@sim
  samp <- sim$samples
  W <- sim$warmup
  I <- sim$iter
  sel <- c(if (keep_warmup) 1:W, sample((W + 1):I, size = n))
  subsamp <- lapply(samp, function(chain_samp) {
    lapply(chain_samp, function(x) x[sel])
  })
  stanfit@sim$samples <- subsamp
  stanfit
} 

                             