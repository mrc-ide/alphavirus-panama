data {
     int<lower=0> Nobs;
     int Npos[Nobs];
     int Ntotal[Nobs];
     int Age[Nobs];
     int <lower=1>Ymax;
     matrix[Nobs, Ymax] AgeExpoMatrix;
}

parameters {
   row_vector<lower=0, upper=2>[Ymax] alpha;
   simplex[2] theta[2];
}

transformed parameters {
  real P[Nobs];
  row_vector<lower=0, upper=2>[Ymax] foi;
  real ScalerDotProduct[Nobs];
  matrix[Ymax, 2] lp;
  real acc[2];
  
  // forward algorithm to determine log_p | t=1,2,..,i
  lp[1] = rep_row_vector(-log(2), 2);
  for(i in 2:Ymax){
    for(k in 1:2){
      for(j in 1:2)
        acc[j] = lp[i - 1, j] + log(theta[j, k]);
      lp[i, k] = log_sum_exp(acc);
    }
  }
  
  // Viterbi to determine state sequence
  {
    real log_p_y_star;
    int y_star[Ymax];
    int back_ptr[Ymax, 2];
    real best_logp[Ymax, 2];
    real best_total_logp;
    for (k in 1:2)
      best_logp[1, k] = 0.0;
    for (t in 2:Ymax) {
      for (k in 1:2) {
        best_logp[t, k] = negative_infinity();
        for (j in 1:2) {
          real logp;
          logp = best_logp[t-1, j] + log(theta[j, k]);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_y_star = max(best_logp[Ymax]);
    for (k in 1:2)
      if (best_logp[Ymax, k] == log_p_y_star)
        y_star[Ymax] = k;
    for (t in 1:(Ymax - 1))
      y_star[Ymax - t] = back_ptr[Ymax - t + 1, y_star[Ymax - t + 1]];
    
    // print(y_star);
    for(i in 1:Ymax)
      foi[i] = alpha[i] * (y_star[i] - 1);
  }
    
  
 for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product(AgeExpoMatrix[i,], foi);
   P[i] = 1 - exp(-ScalerDotProduct[i]);
 }
}

model {
  for (i in 1:Nobs)
    Npos[i] ~ binomial(Ntotal[i], P[i]) ;
  target += log_sum_exp(lp[Ymax]);
}

generated quantities{
  int y_star[Ymax];
  {
    real log_p_y_star;
    int back_ptr[Ymax, 2];
    real best_logp[Ymax, 2];
    real best_total_logp;
    for (k in 1:2)
      best_logp[1, k] = 0.0;
    for (t in 2:Ymax) {
      for (k in 1:2) {
        best_logp[t, k] = negative_infinity();
        for (j in 1:2) {
          real logp;
          logp = best_logp[t-1, j] + log(theta[j, k]);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_y_star = max(best_logp[Ymax]);
    for (k in 1:2)
      if (best_logp[Ymax, k] == log_p_y_star)
        y_star[Ymax] = k;
    for (t in 1:(Ymax - 1))
      y_star[Ymax - t] = back_ptr[Ymax - t + 1, y_star[Ymax - t + 1]];
  }
}
