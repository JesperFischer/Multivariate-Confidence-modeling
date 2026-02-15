
data {
  int<lower=0> N;
  vector[N] RT;

  vector[N] X;
  real minRT;

}

transformed data{
  int P = 4;
}

parameters {
  vector[P] gm;
  real<lower=0, upper = minRT> rt_ndt;

}

transformed parameters{


  real a = (gm[1]);
  real b = (gm[2]);
  real c = gm[3];
  real rt_prec = exp(gm[4]);


}
model {

  gm[1] ~ normal(0,3); //global mean of beta
  gm[2] ~ normal(0,3); //global mean of beta
  gm[3] ~ normal(0,3); //global mean of beta
  gm[4] ~ normal(0,3); //global mean of beta


  rt_ndt ~ normal(0.3,0.05);

  target += lognormal_lpdf(RT - rt_ndt | a*X^2+b*X+c, rt_prec);


}


