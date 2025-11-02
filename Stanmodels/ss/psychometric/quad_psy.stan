functions {


  real psycho_ACC(real x, real alpha, real beta, real lapse){
    return (lapse + (1-2*lapse) * inv_logit(beta * (x - alpha)));

   }



  real entropy(real p){
    return(-p * log(p) - (1-p) * log(1-p));
  }

  real gauss_copula_cholesky_lpdf(matrix u, matrix L) {
    array[rows(u)] row_vector[cols(u)] q;
    for (n in 1:rows(u)) {
      q[n] = inv_Phi(u[n]);
    }

    return multi_normal_cholesky_lpdf(q | rep_row_vector(0, cols(L)), L)
            - std_normal_lpdf(to_vector(to_matrix(q)));
  }

   vector gauss_copula_cholesky_per_row(matrix u, matrix L) {
    int N = rows(u);
    int D = cols(u);
    array[N] row_vector[D] q;
    vector[N] loglik;

    for (n in 1:N) {
        q[n,] = inv_Phi(u[n,]);
        loglik[n] = multi_normal_cholesky_lpdf(to_row_vector(q[n,]) |
                                                 rep_row_vector(0, D), L) - std_normal_lpdf(to_vector(to_matrix(q[n,])));
    }

    return loglik;
  }





  matrix uvar_bounds(array[] int binom_y, vector gm, vector tau_u,matrix L_u, matrix z_expo, array[] int S_id, vector X,
                     int is_upper) {

    int N = size(binom_y);
    matrix[N, 1] u_bounds;

    int S = cols(z_expo);
    int P = rows(z_expo);

    matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

    matrix[S, P] param;

    for(p in 1:P){
      param[,p]= gm[p] + indi_dif[,p];
    }


    vector[S] alpha = (param[,1]);
    vector[S] beta = (param[,2]);
    vector[S] lapse = inv_logit(param[,3]) / 2;


    for (n in 1:N) {
      real theta = get_prob_cor(psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]]), lapse[S_id[n]]), X[n]);
      if (is_upper == 0) {
        u_bounds[n, 1] = binom_y[n] == 0.0
                          ? 0.0 : binomial_cdf(binom_y[n] - 1 | 1, theta);
      } else {
        u_bounds[n, 1] = binomial_cdf(binom_y[n] | 1, theta);
      }
    }

    return u_bounds;
  }


  matrix uvar_bounds_conf(array[] int conf, vector gm, vector tau_u,matrix L_u, matrix z_expo, array[] int S_id, vector X, vector ACC,
                     int is_upper) {


    int S = cols(z_expo);
    int P = rows(z_expo);


    int N = size(conf);

    matrix[N, 1] u_bounds;

    matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

    matrix[S, P] param;

    for(p in 1:P){
      param[,p]= gm[p] + indi_dif[,p];
    }


    vector[S] alpha = (param[,1]);
    vector[S] beta = (param[,2]);
    vector[S] lapse = inv_logit(param[,3]) / 2;

    vector[S] meta_un = param[,7];
    vector[S] meta_bias = param[,8];

    for (n in 1:N) {
      real theta_conf = psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]] + meta_un[S_id[n]]), lapse[S_id[n]]);
      real theta = get_conf(ACC[n],theta_conf, X[n], alpha[S_id[n]]);

      if (is_upper == 0) {
        u_bounds[n, 1] = conf[n] == 0.0
                          ? 0.0 : binomial_cdf(conf[n] - 1 | 1, theta);
      } else {
        u_bounds[n, 1] = binomial_cdf(conf[n] | 1, theta);
      }
    }

    return u_bounds;
  }



  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    if(cutnum==1) {

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    // divide in half for the two cutpoints

    // don't forget the ordered transformation

      return   dirichlet_lpdf(p | alpha)
           + log_determinant(J) + cut2;
    } else {
      return(0);
    }
  }

    real get_conf(real ACC, real theta, real x, real alpha){
  if(ACC == 1 && x > alpha){
    return(theta);
  }else if(ACC == 1 && x < alpha){
    return(1-theta);
  }else if(ACC == 0 && x > alpha){
    return(1-theta);
  }else if(ACC == 0 && x < alpha){
    return(theta);
  }else{
    return(0);
  }
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
  int<lower=0> S;
  array[N] int S_id;

  array[S] int starts;
  array[S] int ends;

  array[N] int binom_y;
  vector[N] RT;
  array[N] int Conf;

  vector[N] X;

  vector[S] minRT;

  vector[N] ACC;

  array[S] int t_p_s;



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
  theta[n] = get_prob_cor(psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]]), lapse[S_id[n]]), X[n]);

  }
  }


}

model {

  gm[1] ~ normal(0,20); //global mean of beta
  gm[2] ~ normal(-2,3); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta

  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(5 , 10);
  tau_u[2:3] ~ normal(0 , 3);

  L_u ~ lkj_corr_cholesky(2);




  for (n in 1:N) {
    target += binomial_lpmf(binom_y[n] | 1, theta[n]);
  }


}

generated quantities {

  matrix[P,P] correlation_matrix = L_u * L_u';

  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


  for (n in 1:N) {
    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, theta[n]);

    log_lik[n] = log_lik_bin[n];
  }


}
