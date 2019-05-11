functions {

row_vector lambda_yexpo_fun (real lambda0, int Ymax) {

row_vector[Ymax] lambda_yexpo;

for (i in 1:Ymax) {lambda_yexpo[i] = lambda0;}

return lambda_yexpo;

 }
 
}
  

data {
     int<lower=0> Nobs;
     int Npos [Nobs];
     int Ntotal [Nobs];
     int Age [Nobs];
     int <lower=1>Ymax;
     int <lower=1>Dur;
     matrix[Nobs,Ymax] AgeExpoMatrix;
     }

parameters {
           real<lower=0,upper=5> lambda0;
           }

transformed parameters {
  real P [Nobs];
  row_vector [Ymax] lambda_yexpo;
  real ScalerDotProduct [Nobs];

  lambda_yexpo = lambda_yexpo_fun (lambda0, Ymax);

 for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product (AgeExpoMatrix[i,], lambda_yexpo) ;
   P[i] = 1 - exp(-ScalerDotProduct[i]);
   }

}

model {
  lambda0 ~ uniform (0, 2);
  
  for (i in 1:Nobs){
  Npos[i] ~ binomial(Ntotal[i], P[i]);
  }
  
}



generated quantities{
  vector[Nobs] Npos_sim;
  vector[Nobs] P_sim;
  vector[Nobs] logLikelihood;
  for(i in 1:Nobs){
    Npos_sim[i] = binomial_rng(Ntotal[i], P[i]);
    P_sim[i] = Npos_sim[i] / Ntotal[i];
    logLikelihood[i] = binomial_lpmf(Npos[i] | Ntotal[i], P[i]);
  }
}



