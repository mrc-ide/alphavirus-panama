
check_plots <- function(res_total, rres, variable, yexpo, dat) {
  

  varr <- res_total[,grep(variable, names(res_total))]
  Dat <- data.frame(varr) 
  colnames(Dat) <- 1:4
  Dat$iter <- 1:nrow(Dat)
  Dat <- melt (Dat, id = 'iter', variable.name = 'chain')
  c1 <- ggplot(Dat) + 
    geom_line(aes(y=value, x = iter, colour=chain), size = .1) +
    theme(legend.position = 'none') + xlab('') + ylab (variable)
  
  
  if (variable == 'lambdac') {bin_hist = 0.001} else {
    bin_hist = 0.03
  }
  
  rr <- data.frame(var = rres[, variable])
  c2 <-  ggplot(rr, aes(var)) +
    geom_histogram(binwidth = bin_hist, 
                   colour = 'black', size = .05) + 
    xlab(variable)
  
  var_years <- c('year1', 'year2', 'year3')
  var_lambdas <- c('lambda1', 'lambda2', 'lambda3')
  
  if (variable %in%  var_years) {
    c2 <- c2 + coord_cartesian(xlim = c(min(yexpo), max(yexpo)))
  } 
  
  
  if (variable %in%  var_lambdas) {
    c2 <- c2 + coord_cartesian(xlim = c(0, 1.5))
  } 
  
  
  if (variable == 'timelag') {
    c2 <- c2 + coord_cartesian(xlim = c(0, 50))
  } 
  
  
  pfinal <- gridExtra::grid.arrange(c1, c2, nrow = 1)
  
  return(pfinal) 
}
