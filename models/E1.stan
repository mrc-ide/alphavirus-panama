
functions {
  
    int bin_search (real x, int min_val, int max_val){

    int range = (max_val - min_val+1)/2; 
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        range =  (range+1)/2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }
  


row_vector lambda_yexpo_fun (real lambda1,
                             real year1,
                             int Ymax,
                             int Dur){

 row_vector[Ymax] lambda_yexpo;
 real lambda0 = 1e-10;
 real year_rounded = round(year1);
 int YearOutbreak;
 int Ymin = 1;
 int Duration = Dur;
 
 YearOutbreak = bin_search(year_rounded, Ymin, Ymax);
  
 for (i in 1:Ymax) {lambda_yexpo[i] = lambda0;}
 
 for (d in 1:Duration) {
  lambda_yexpo [YearOutbreak + Duration -1] = lambda1;
  }

  return lambda_yexpo;
  }
}

data {
     int<lower=0> Nobs;
     int Npos [Nobs];
     int Ntotal [Nobs];
     int Age [Nobs];
     int <lower=1>Ymin;
     int <lower=1>Ymax;
     int <lower=1>Dur;
     matrix[Nobs,Ymax] AgeExpoMatrix;
}

parameters {
           real<lower=0,upper=2> lambda1;
           real<lower=1,upper=Ymax-Dur> year1;
}

transformed parameters {
  real P [Nobs];
  row_vector [Ymax] lambda_yexpo;
  real ScalerDotProduct [Nobs];

  lambda_yexpo = lambda_yexpo_fun (lambda1, year1,Ymax, Dur);
  

 for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product (AgeExpoMatrix[i,], lambda_yexpo) ;
   P[i] = 1 - exp(-ScalerDotProduct[i]);
 }


}

model {
  lambda1 ~ uniform (0, 2);
  year1   ~ uniform (1, Ymax-Dur);
  
  for (i in 1:Nobs){
  Npos[i] ~ binomial(Ntotal[i], P[i]) ;
  }
  
}


generated quantities{
  vector[Nobs] Npos_sim;
  vector[Nobs] P_sim;
  vector[Nobs] log_lik;
  real sum_log_lik;
  
  for(i in 1:Nobs){
    Npos_sim[i] = binomial_rng(Ntotal[i], P[i]);
    P_sim[i] = Npos_sim[i] / Ntotal[i];
  }
  
  for(i in 1:Nobs){
    log_lik[i] = binomial_logit_lpmf(Npos[i] | Ntotal[i], P[i]);
  }
  
  sum_log_lik = sum (log_lik);
  
}








