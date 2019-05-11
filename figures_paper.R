#=============================== REPORTING
#         Force-of-infection Una virus
#          (By Individual sero-surveys)
#
# ============================================ 


rm(list=ls())

source('fun/plot_ress_paper.R')
source('fun/libraries.R')
library(cowplot)




plot_fitted <- function(i){
  
  # Data
  
  datt      <- readRDS('data/Panama.RDS')
  datasets  <- sort(unique(datt$dataset_id))
  dataset    <- datasets[i]
  result_date = 'res_17.12.2018'
  report_date = 'report_17.12.2018'
  
  
  nbreaks <- 0
  dat1 <- readRDS(paste0('results/', result_date,'/',dataset,'_model_',nbreaks, '.RDS'))
  
  nbreaks <- 1
  dat2 <- readRDS(paste0('results/', result_date,'/',dataset,'_model_',nbreaks, '.RDS'))
  
  virus= dat1$dat_imputed$virus[1]
  res1 <- processing_res (dat1)
  res2 <- processing_res (dat2)
  
  p1 <- plot_lambda (res1, res2) 
  p2 <- plot_prevalence (res1, res2) 
  
  pf <- plot_grid(p1,
                  p2,
                  ncol=1, align="v", rel_heights=c(0.8, 1))
  
  res_virus <- list(virus= virus,
                    plot = pf,
                    DIC_model1  = res1$DIC,
                    DIC_model2  = res2$DIC)
  
  
  return(res_virus)
  
}




#### =================
MADV <- plot_fitted(i = 10)
UNA  <- plot_fitted(i = 11)
VEEV <- plot_fitted(i = 12)

grid.arrange(
  arrangeGrob(UNA$plot, MADV$plot, VEEV$plot,

              nrow = 1,
              ncol = 3
  )
)


library(ggpubr)
ggarrange(UNA$plot, 
          MADV$plot, 
          VEEV$plot, 
          # labels = c("A", "B", "C"),
          heights = c(1, 1, 1),
          ncol = 3, nrow = 1)


# 
# #### =================
# MADV <- plot_fitted(i = 1)
# UNA  <- plot_fitted(i = 2)
# VEEV <- plot_fitted(i = 3)
# 
# grid.arrange(
#   arrangeGrob(UNA$plot, MADV$plot, VEEV$plot, 
#               
#               nrow = 1,
#               ncol = 3
#   )
# )
# 
# 
# #### =================
# MADV <- plot_fitted(i = 4)
# UNA  <- plot_fitted(i = 5)
# VEEV <- plot_fitted(i = 6)
# 
# grid.arrange(
#   arrangeGrob(UNA$plot, MADV$plot, VEEV$plot, 
#               
#               nrow = 1,
#               ncol = 3
#   )
# )
# 
# #### =================
# MADV <- plot_fitted(i = 7)
# UNA  <- plot_fitted(i = 8)
# VEEV <- plot_fitted(i = 9)
# 
# grid.arrange(
#   arrangeGrob(UNA$plot, MADV$plot, VEEV$plot, 
#               
#               nrow = 1,
#               ncol = 3
#   )
# )



