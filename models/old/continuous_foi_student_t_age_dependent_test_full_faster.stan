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
   real<lower=0> foi[Ymax];
   real<lower=0> mu;
}

model {
}

generated quantities{
  real P[Nobs];
  int Npos_sim[Nobs];
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
    Npos_sim[i] = binomial_rng(Ntotal[i], P[i]);
  }
}