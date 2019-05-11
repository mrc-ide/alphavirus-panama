functions{
  real[] short_solve(real s0, real i0, real lambda, real mu){
    real s = exp(-lambda) * s0;
    real i = exp(-(lambda + mu)) * (-i0 * lambda * exp(lambda) - exp(lambda) * s0 * lambda + exp(mu) * s0 * lambda + exp(lambda) * i0 * mu) / (mu - lambda);
    real y[2];
    y[1] = s;
    y[2] = i;
    return(y);
  }
}


data {
   int<lower=0> Nobs;
   int Npos[Nobs];
   int Ntotal[Nobs];
   int Age[Nobs];
   int <lower=1>Ymax;
   matrix[Nobs, Ymax] AgeExpoMatrix;
   int birthdays[Nobs];
}


parameters {
   real<lower=0> foi[Ymax];
   real<lower=0> sigma;
   real<lower=0> nu;
   real<lower=0> mu;
}

transformed parameters {
  real P[Nobs];
  real ytemp[2];

  for(i in 1:Nobs){
    real s0 = 1;
    real i0 = 0;
    for(t in 1:Ymax){
      if(t >= birthdays[i]){
        ytemp = short_solve(s0, i0, foi[t], mu);
        s0 = ytemp[1];
        i0 = ytemp[2];
      }
    }
    P[i] = ytemp[2] / sum(ytemp);
  }
}

model {
  for (i in 1:Nobs)
    Npos[i] ~ binomial(Ntotal[i], P[i]) ;
  sigma ~ cauchy(0, 1);
  nu ~ cauchy(0, 1);
  mu ~ cauchy(0, 1);
  
  for(i in 2:Ymax)
    foi[i] ~ student_t(nu, foi[i - 1], sigma);
  foi[1] ~ normal(0, 1);
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