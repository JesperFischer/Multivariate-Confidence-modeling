functions {


  real psycho_ACC(real x, real alpha, real beta, real lapse){
    return (lapse + (1-2*lapse) * inv_logit(beta * (x - alpha)));
   }


 real get_prob_cor(real theta, real x){
  if(x > 0){
    return(theta);
  }else if(x < 0){
    return(1-theta);
  }else{
    return(0);
  }

}
}


data {
  int<lower=0> N;

  array[N] int binom_y;

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
  real lapse = gm[3];



  vector[N] theta;

  for (n in 1:N) {
  theta[n] = get_prob_cor(psycho_ACC(X[n], (alpha), exp(beta), inv_logit(lapse)/ 2), X[n]);

}


}
model {
  gm[1] ~ normal(0,0.5); //global mean of beta
  gm[2] ~ normal(1,2); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta

  for(n in 1:N){
    target += bernoulli_lpmf(binom_y[n] | theta[n]);
  }

}

generated quantities {


  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


  for(n in 1:N){
    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, theta[n]);
    log_lik[n] = log_lik_bin[n];
  }



}
