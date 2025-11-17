functions {


  real psycho(real x, real alpha, real beta, real lapse){
    return (lapse + (1-2*lapse) * inv_logit(beta * (x - alpha)));
   }

}


data {
  int<lower=0> N;
  int<lower=0> S;
  array[N] int S_id;

  array[N] int binom_y;

  vector[N] X;

}

transformed data{
  int P = 3;
}

parameters {
  vector[P] gm;
  vector<lower=0>[P] tau_u;
  cholesky_factor_corr[P] L_u;    // Between participant cholesky decomposition
  matrix[P, S] z_expo;    // Participant deviation from the group means


}

transformed parameters{

  // Extracting individual deviations for each subject for each parameter
  matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

  matrix[S, P] param;

  for(p in 1:P){
    param[,p]= gm[p] + indi_dif[,p];
  }

  vector[S] alpha = (param[,1]);
  vector[S] beta = (param[,2]);
  vector[S] lapse = inv_logit(param[,3]) / 2;

  vector[N] theta;

  profile("likelihood") {
  for (n in 1:N) {
  theta[n] = psycho(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]]), lapse[S_id[n]]);
  }
  }

}
model {
  gm[1] ~ normal(0,10); //global mean of beta
  gm[2] ~ normal(-2,2); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta



  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(10 , 3);
  tau_u[2] ~ normal(0 , 3);
  tau_u[3] ~ normal(0 , 3);

  L_u ~ lkj_corr_cholesky(2);

  for(n in 1:N){
    binom_y[n] ~ bernoulli(theta[n]);
  }




}

generated quantities {

  matrix[P,P] correlation_matrix = L_u * L_u';

  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


  for(n in 1:N){
    log_lik_bin[n] = bernoulli_lpmf(binom_y[n] | theta[n]);
    log_lik[n] = log_lik_bin[n];
  }



}
