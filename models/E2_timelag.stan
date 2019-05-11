
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
                             real lambda2,
                             real year1,
                             real timelag,
                             int Ymax,
                             int Dur){

 row_vector[Ymax] lambda_yexpo;
 real base_foi = 0.000000001;
 real year1_rounded = round(year1);
 real timelag_rounded = round(timelag);
 int YearOutbreak1;
 int YearOutbreak2;
 int TimeLag;
 int Ymin = 1;
  
 for (i in 1:Ymax) {lambda_yexpo[i] = base_foi;}
 YearOutbreak1 = bin_search(year1_rounded, Ymin, Ymax);
 TimeLag = bin_search(timelag_rounded, Ymin, Ymax);
 YearOutbreak2 = YearOutbreak1 + TimeLag;
lambda_yexpo [YearOutbreak1] = lambda1; 
lambda_yexpo [YearOutbreak2] = lambda2; 
if(Dur > 1) {
lambda_yexpo [YearOutbreak1+1] = lambda1; 
lambda_yexpo [YearOutbreak2+1] = lambda2; 
}

if(Dur > 2) {
lambda_yexpo [YearOutbreak1+2] = lambda1; 
lambda_yexpo [YearOutbreak2+2] = lambda2; }

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
           real<lower=0,upper=2> lambda1;
           real<lower=0,upper=2> lambda2;
           real<lower=1,upper=Ymax> year1;
           real<lower=1 + Dur,upper=Ymax-year1-Dur> timelag;
}

transformed parameters {
  real P [Nobs];
  row_vector [Ymax] lambda_yexpo;
  real ScalerDotProduct [Nobs];

  lambda_yexpo = lambda_yexpo_fun (lambda1, lambda2,
                                    year1, timelag,Ymax,Dur);
  

 for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product (AgeExpoMatrix[i,], lambda_yexpo) ;
   P[i] = 1 - exp(-ScalerDotProduct[i]);
 }


}

model {
  lambda1 ~ uniform (0, 2);
  lambda2 ~ uniform (0, 2);
  year1   ~ uniform (1, Ymax - Dur);
  timelag ~ uniform (1 + Dur, Ymax-year1-Dur);
  
  for (i in 1:Nobs){
  Npos[i] ~ binomial(Ntotal[i], P[i]) ;
  }
  
}


