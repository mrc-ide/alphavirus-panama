
#------------------------- Fit Function

fFitModel <- function(model, dat){
  
  yexpo <- make_yexpo(dat)
  yexpo <- yexpo[-length(yexpo)]
  RealYexpo <- (min(dat$birth_year):dat$tsur[1])[-1]
  ExposureMatrix <- get_exposure_matrix(dat)[, -1]
  Nobs <- nrow(dat)
  
  stan_data <- list(
    Nobs   = nrow(dat),
    Npos   = dat$counts,
    Ntotal = dat$total,
    Age    = dat$age_mean_f,
    Ymax   = max(yexpo),
    Dur    = 1,
    AgeExpoMatrix = ExposureMatrix)
  
  fInit <- function(){
    list(foi=rep(0.01, max(yexpo)))
  } 
  
  # suppressMessages(library(rstan, warn.conflicts = FALSE, quietly=TRUE))
  
  fit <- sampling(model, data=stan_data, iter=200, chains=4, init=fInit, 
                  verbose=FALSE,refresh=0)
  
  res <- list(fit=fit,
                  stan_data = stan_data,
                  RealYexpo = RealYexpo,
                  yexpo     = yexpo)
  
  return(res)
  
}



#------------------------- Plot Function

fPlotModel <- function(res, dat, model, model_details) {
  survey  <- dat$survey[1]
  fit <- res$fit
  
  if(class(fit@sim$samples)  != "NULL") 
    
  {
    country <- paste(dat$country[1])
    Ref     <- dat$Ref[1]
    N       <- sum(dat$total)
    loo_fit <- loo(fit, save_psis = TRUE, 'logLikelihood')
    info <- data.frame(c(model, model_details, round(loo_fit$estimates[,1],3), survey, country, Ref, N))
    colnames(info) <- NULL
    rownames(info) <- c('model', 'details', rownames(loo_fit$estimates), 'survey', 'country', 'Ref', 'N')

    blank <- data.frame(x= 0:1, y = 0:1)
    
    g0 <- 
      ggplot(blank, aes(x, y)) +
      geom_blank() + ylab('') + xlab ('') + theme_void() +
      annotation_custom(tableGrob(info)) 
    
    
    RealYexpo <- res$RealYexpo
    foi_est <- colMeans(rstan::extract(fit, 'foi')[[1]])
    foi_chain <- data.frame(rstan::extract(fit, 'foi')[[1]])
    foi_summary <- data.frame(sapply(foi_chain, function(i) c(quantile(i, c(0.5, 0.025, 0.975)))))
    colnames(foi_summary) <- RealYexpo
    foi_summary$metric <- c('median', 'lower', 'upper') 
    foi_summary <- melt(foi_summary, id='metric', variable.name = 'year') %>% spread(metric, value) 
    foi_summary$year <- (as.numeric(as.character(foi_summary$year)))
    
    
    
    
    
    max_lambda <- 0.2
    if(any(foi_summary$median > 0.2)) {max_lambda <- 0.5}
    if(any(foi_summary$median > 0.5)) {max_lambda <- 1.0}
    if(any(foi_summary$median > 1.0)) {max_lambda <- 1.5}
    
    space_years <- 1 
    if(length(res$yexpo)>15 ) {space_years <- 2}
    if(length(res$yexpo)>30 ) {space_years <- 3}
    if(length(res$yexpo)>50 ) {space_years <- 4}
    
    g1 <- 
      foi_summary %>%
      ggplot(aes(x=year, y=median)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#7fcdbb") +
      geom_line() +
      scale_color_brewer("FOI", palette = "Dark2") +
      coord_cartesian(ylim = c(0, max_lambda)) + theme_bw() + theme(legend.position = 'none') +
      ylab ('FOI')+
      theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
      scale_x_continuous(breaks = round(seq(min(foi_summary$year), max(foi_summary$year), by = space_years),0)) 
    
    
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
      geom_boxplot(aes(group=age_class), outlier.shape = NA) +
      geom_errorbar(data=filter(dat1, source=="actual"), aes(x=age_class, ymin = pobs_lower, ymax = pobs_upper), 
                    colour ='red', width = .05) +
      geom_point(data=filter(dat1, source=="actual"), aes(x=age_class, y=value, size = total), pch = 21, colour="red", fill = 'white') +
      xlab("Age, years") +
      ylab("Seropositivity") +
      theme_bw() + coord_cartesian(ylim = c(0,1)) +
      theme(axis.text.x = element_text(angle=45)) +
      theme(legend.position = 'none') +
      theme(plot.title = element_text(size=10))
    
    
    g3 <- 
      dat3 %>% 
      select(variable, iteration, resid) %>%
      ggplot(aes(x=as.factor(variable), y=resid, group=as.factor(iteration))) +
      geom_line(alpha=0.4) +
      geom_hline(yintercept = 0, linetype=2, colour="orange") + xlab('') +
      theme_bw () + theme(axis.text.x = element_blank()) +
      theme(plot.title = element_text(size=10))
    
    

    
  
    } else 
   
    {
      print_warning <- 'errors'
      df <- data.frame() 
      g0 <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x = 4, y = 5, label = print_warning) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
        ylab('') + xlab('') 
      g1 <- g0
      g2 <- g0
      g3 <- g0
      g0 <- g0 +   labs(subtitle = model)+
        theme(plot.title = element_text(size=10))
      
      
    }
  
  
  p_all <- plot_grid(g0, g1, g2, g3,
                     nrow = 4, 
                     rel_heights = c(1, 1.3, 1, 1)
                     ) #, align = "hv", 
  
  return(p_all)
}
