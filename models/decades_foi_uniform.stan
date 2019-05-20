functions {
row_vector lambda_yexpo_fun (row_vector foi,
                             int Ymax){

 row_vector[Ymax] lambda_yexpo;
  int StartDec5 = Ymax - 10;
  int StartDec4 = Ymax - 20;
  int StartDec3 = Ymax - 30;
  int StartDec2 = Ymax - 40;
  int StartDec1 = Ymax - 50;
  int EndDec4   = StartDec5-1;
  int EndDec3   = StartDec4-1;
  int EndDec2   = StartDec3-1;
  int EndDec1   = StartDec2-1;
  
   for (j in 1: EndDec1) lambda_yexpo[j] = foi[1];
   for (k in StartDec2: EndDec2) lambda_yexpo[k] = foi[2];
   for (l in StartDec3: EndDec3) lambda_yexpo[l] = foi[3];
   for (m in StartDec4: EndDec4) lambda_yexpo[m] = foi[4];
   for (n in StartDec5: Ymax) lambda_yexpo[n] = foi[5];


  return lambda_yexpo;
  }
}



data {
     int<lower=0> Nobs;
     int Npos[Nobs];
     int Ntotal[Nobs];
     int Age[Nobs];
     int <lower=1>Ymax;
     matrix[Nobs, Ymax] AgeExpoMatrix;
}

parameters {
   
   row_vector<lower=0>[5]foi;
   // real<lower=0> sigma;
   // real<lower=0> nu;
}

transformed parameters {
  real P[Nobs];
  real ScalerDotProduct[Nobs];
  row_vector<lower=0>[Ymax] foi_yexpo;


foi_yexpo = lambda_yexpo_fun (foi,Ymax);

 for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product(AgeExpoMatrix[i,], foi_yexpo);
   P[i] = 1 - exp(-ScalerDotProduct[i]);
 }
}

model 
{
  for (i in 1:Nobs)
    Npos[i] ~ binomial(Ntotal[i], P[i]) ;
      // sigma ~ cauchy(0, 1);
      //    nu ~ cauchy(0, 1);
  
  for(i in 2:5)
    // foi[i] ~ student_t(nu, foi[i - 1], sigma);
    // foi[1] ~ normal(0, 1);
    foi[i] ~ uniform (0, 2);
 

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

