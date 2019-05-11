functions{
  int indexer(real t, int tmax){
    int j;
    for(i in 0:tmax){
      if((t >= i) && (t < (i + 1))){
        j = i + 1;
        break;
      }
    }
    return(j);
  } 
  
  real[] si_derivative(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    real dydt[2];
    real s = y[1];
    real i = y[2];
    int tmax = x_i[1];
    int index_temp = indexer(t, tmax);
    int birth = x_i[2];
    real mu = theta[1];
    if(t < birth){
      dydt[1] = 0.0;
      dydt[2] = 0.0;
    }else{
      real lambda_temp = theta[index_temp + 1];
      dydt[1] = -lambda_temp * s;
      dydt[2] = lambda_temp * s - mu * i;
    }
    return dydt;
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
   int nt;
   real ts[nt];
   real<lower=0> foi[Ymax];
   real<lower=0> sigma;
   real<lower=0> nu;
   real<lower=0> mu;
}

transformed data{
  real y0[2];
  real t0 = 0.0;
  int x_i[Nobs, 2];
  real x_r[0];
  for(i in 1:Nobs){
    x_i[i, 1] = Ymax;
    x_i[i, 2] = birthdays[i];
  }
  y0[1] = 1.0;
  y0[2] = 0.0;
}

model {
}

generated quantities{
  real P[Nobs];
  real y_hat[Nobs, nt, 2];
  real theta[Ymax + 1];
  int Npos_sim[Nobs];
  theta[1] = mu;
  for(i in 1:Ymax)
    theta[i + 1] = foi[i];

  for(i in 1:Nobs){
    y_hat[i] = integrate_ode_rk45(si_derivative, y0, t0, ts, theta, x_r, x_i[i]);
    P[i] = y_hat[i, nt, 2] / sum(y_hat[i, nt]);
    Npos_sim[i] = binomial_rng(Ntotal[i], P[i]) ;
  }
}