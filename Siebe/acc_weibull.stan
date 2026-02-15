functions {


  real psycho_ACC(real x, real alpha, real beta, real lapse, real gamma){
    return (gamma + (1-gamma-lapse) * (1 - exp(-(x/alpha)^beta)));
   }
}


data {
  int<lower=0> N;

  array[N] int ACC;

  vector[N] X;

}

transformed data{
  int P = 3;

}

parameters {
  vector[P] gm;


}

transformed parameters{


  real alpha = (gm[1]);
  real beta = (gm[2]);
  real lapse = inv_logit(gm[3]) / 2;
  real gamma = 0.5;



  vector[N] theta;

  for (n in 1:N) {
  theta[n] = psycho_ACC(X[n], exp(alpha), exp(beta), lapse,gamma);
}


}
model {
  gm[1] ~ normal(-2,1); //global mean of beta
  gm[2] ~ normal(1,1); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta


  for(n in 1:N){
    target += bernoulli_lpmf(ACC[n] | theta[n]);
  }

}

generated quantities {


  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


  for(n in 1:N){
    log_lik_bin[n] = binomial_lpmf(ACC[n] | 1, theta[n]);
    log_lik[n] = log_lik_bin[n];
  }



}
