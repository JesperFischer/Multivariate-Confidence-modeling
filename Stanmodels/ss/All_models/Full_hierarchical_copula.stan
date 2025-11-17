functions {


  // psychometric function to get the probability of responding (1)
  real psycho_ACC(real x, real alpha, real beta, real lapse){
    return (lapse + (1-2*lapse) * inv_logit(beta * (x - alpha)));
   }

  // entropy of a binary decision variable
  real entropy(real p){
    return(-p * log(p) - (1-p) * log(1-p));
  }

  // the log probability density function for the copula
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


  // the bounds for binary variables in copulas based on data augmentation
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


 // cummulative ordered beta distribution density
  real ord_beta_reg_cdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

    real p0 = 1-inv_logit(mu - thresh[1]);

    real p_m = (inv_logit(mu - thresh[1])-inv_logit(mu - thresh[2]))  * beta_cdf(y | exp(log_inv_logit(mu) + log(phi)), exp(log1m_inv_logit(mu) + log(phi)));



    if (y < 0) {
      return 0;
    } else if (y == 0) {
      return p0;
    } else if (y == 1) {
      return 1-(1e-12);
    } else {
      return (p0 + p_m);
    }
  }

 // ordered beta distribution density
  real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu  - thresh[2]);
    } else {
      return log_diff_exp(log_inv_logit(mu - thresh[1]), log_inv_logit(mu - thresh[2])) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }


  // priors for the cutpoints so that they are ordered
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


 // translate the probability of being correct into the two curves of being correct or being wrong
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

  // get the probability of being correct from probability of responding "(1)"
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
  // number of data points
  int<lower=0> N;
  // number of subjects
  int<lower=0> S;
  // subject identifier i.e. (1,1,1,1,1,1,..,2,2,2,2,2)
  array[N] int S_id;

  // which row does a new subject start
  array[S] int starts;
  // which row does a new subject end
  array[S] int ends;

  // binary responses (here correct or incorrect)
  array[N] int binom_y;
  //Response times
  vector[N] RT;
  // confidence ratings
  vector[N] Conf;

  // stimulus value
  vector[N] X;

  // the minimum response time
  vector[S] minRT;

  // the accuracy of the subject at trial.
  vector[N] ACC;

}

transformed data{
  // number of parameters for each subject:
  int P = 9;
}

parameters {
  // number of group mean parameters
  vector[P] gm;
  // between subject variances
  vector<lower=0>[P] tau_u;
  // correlation matrix (correlation between subject level parameters)
  cholesky_factor_corr[P] L_u;    // Between participant cholesky decomposition
  //The random effect i.e. subject deviations from the group mean.
  matrix[P, S] z_expo;    // Participant deviation from the group means

  // the lower and upper bound on the binary variables (cummulative density of the bernoulli) ((used for the copula))
  matrix<
    lower=uvar_bounds(binom_y, gm, tau_u,L_u,z_expo, S_id,X, 0),
    upper=uvar_bounds(binom_y, gm, tau_u,L_u,z_expo, S_id,X, 1)
  >[N, 1] u;


  // each subject get a residual copula correlation matrix (there are 3 correlations, between B,RT and C)
  array[S] cholesky_factor_corr[3] rho_chol;

  // each subject gets a cutpoint for responding 0
  vector[S] c0;
  // each subject gets a cutpoint for responding 1
  vector[S] c11;
    // each subject gets a non-decision time parameter that is between 0 and their lowest response time. As the non-decision time can't be less than the lowest response time.
  vector<lower=0, upper = minRT>[S] rt_ndt;

}



transformed parameters{

  // Extracting individual deviations for each subject for each parameter
  matrix[S, P] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';

  // defining matrix for all subject level parameters
  matrix[S, P] param;

  // filling this matrix with group level and indicidual difference:
  for(p in 1:P){
    param[,p]= gm[p] + indi_dif[,p];
  }

  // extracting the subject level paramters
  //psychometric function parametesr
  vector[S] alpha = (param[,1]);
  vector[S] beta = (param[,2]);
  vector[S] lapse = inv_logit(param[,3]) / 2;

  // response time distribution
  vector[S] rt_int = param[,4];
  vector[S] rt_slope = param[,5];
  vector[S] rt_prec = exp(param[,6]);

  //confidence parameters.
  vector[S] conf_prec = exp(param[,7]);
  vector[S] meta_un = param[,8];
  vector[S] meta_bias = param[,9];
  vector[S] rt_stim = param[,10];



 //trial by trial predictions of the model for:
 //entropy
  vector[N] entropy_t;
  //mean confidence
  vector[N] conf_mu;
  //probability of responding 1
  vector[N] theta;
  // probability of responding 1 (from the perspective of confidence), i.e. with meta cognitive noise:
  vector[N] theta_conf;


  // this is filling up these above trial by trial predictions with the model.
  // it loops through all trials and puts in the model prediction at that trial for that subject notice the S_id[n] for all subject level parameters.
  for (n in 1:N) {
  theta[n] = psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]]), lapse[S_id[n]]);

  entropy_t[n] = entropy(psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]]), lapse[S_id[n]]));

  theta_conf[n] = psycho_ACC(X[n], (alpha[S_id[n]]), exp(beta[S_id[n]] + meta_un[S_id[n]]), lapse[S_id[n]]);

  conf_mu[n] = get_conf(ACC[n],theta_conf[n], X[n], alpha[S_id[n]]);
  }

}
model {

  // Priors

  //psychometric functiuon parameters
  gm[1] ~ normal(0,5);
  gm[2] ~ normal(-2,2);
  gm[3] ~ normal(-4,2);

  // RT and confidence parameters
  gm[4:9] ~ normal(0,2);

  //subject level devbiations from the group mean are standard normally distributioned as prior
  to_vector(z_expo) ~ std_normal();

  // between subject varability parameters for each of the same parametesr from above:
  tau_u[1] ~ normal(3 , 3);
  tau_u[2] ~ normal(0 , 3);
  tau_u[3] ~ normal(0 , 3);
  tau_u[4:9] ~ normal(0 , 2);

  // prior on the correlation matrix between subject level parameters:
  L_u ~ lkj_corr_cholesky(2);

  // prior on the non-decision time.
  rt_ndt ~ normal(0.3,0.05);


  // calculating the uniform variables for the copula bit, by using the cummulative marginal distribution to get a uniform variable (based on the model), for each of the three response variables (B,RT,C)
  matrix[N, 3] u_mix;
  for (n in 1:N) {
    u_mix[n, 1] = u[n,1];

    u_mix[n, 2] = lognormal_cdf(RT[n] - rt_ndt[S_id[n]] | rt_int[S_id[n]] + rt_slope[S_id[n]] * entropy_t[n] , rt_prec[S_id[n]]);

    u_mix[n, 3] = ord_beta_reg_cdf(Conf[n] | logit(conf_mu[n]) + meta_bias[S_id[n]], conf_prec[S_id[n]], c0[S_id[n]], c11[S_id[n]]);

    // here we then add the data (RT and Conf[n]) to the model that we have (i.e. telling stan that this is what need to be optimized)
    // first the likeliood for the response times
    target += lognormal_lpdf(RT[n] - rt_ndt[S_id[n]] | rt_int[S_id[n]] + rt_slope[S_id[n]] * entropy_t[n], rt_prec[S_id[n]]);

    // then the likelihood for the confidnnce ratings (notice the meta noise is adding in logit space)
    target += ord_beta_reg_lpdf(Conf[n] | logit(conf_mu[n])+ meta_bias[S_id[n]], conf_prec[S_id[n]], c0[S_id[n]], c11[S_id[n]]);

  }

  // next we need to also account for the copula (and in this also the bianry responses)
  for(s in 1:S){
    // this is the prior for the cutpoints (they are a bit special as they need to be ordered)
    c0[s] ~ induced_dirichlet([1,10,1]', 0, 1, c0[s], c11[s]);
    c11[s] ~ induced_dirichlet([1,10,1]', 0, 2, c0[s], c11[s]);

    //prior for the correlation between measures (i.e. copula correlation)
    rho_chol[s] ~ lkj_corr_cholesky(2);

    //the likelihood for the copula part (and herein the binary responses.)
    u_mix[starts[s]:ends[s],] ~ gauss_copula_cholesky(rho_chol[s]);
  }

}


// all of this are just extract computations for convinence afters (all this can be calculated from the above in R afterwards.)
generated quantities {

 vector[S] c1 = c0 + exp(c11);
  vector[S] rho_p_rt;
  vector[S] rho_p_conf;
  vector[S] rho_rt_conf;

  matrix[P,P] correlation_matrix = L_u * L_u';

  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik_rt = rep_vector(0,N);
  vector[N] log_lik_conf = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);



 matrix[N, 3] u_mixx;
  for (n in 1:N) {
    u_mixx[n, 1] = u[n,1];

    u_mixx[n, 2] = lognormal_cdf(RT[n] - rt_ndt[S_id[n]] | rt_int[S_id[n]] + rt_slope[S_id[n]] * entropy_t[n], rt_prec[S_id[n]]);

    u_mixx[n, 3] = ord_beta_reg_cdf(Conf[n] | logit(conf_mu[n])+ meta_bias[S_id[n]], conf_prec[S_id[n]], c0[S_id[n]], c11[S_id[n]]);
  }

  vector[N] log_lik_cop;
  int pos;
  pos = 1;
  for (s in 1:S) {
    int n_s = ends[s] - starts[s] + 1;
    vector[n_s] log_lik_s;

    log_lik_s = gauss_copula_cholesky_per_row(u_mixx[starts[s]:ends[s], ], rho_chol[s]);

    // store results in the big vector
    log_lik_cop[pos:(pos + n_s - 1)] = log_lik_s;

    pos += n_s;
  }

  for(s in 1:S){

    rho_p_rt[s] = multiply_lower_tri_self_transpose(rho_chol[s])[1, 2];
    rho_p_conf[s] = multiply_lower_tri_self_transpose(rho_chol[s])[1, 3];
    rho_rt_conf[s] = multiply_lower_tri_self_transpose(rho_chol[s])[2, 3];

  }

  for(n in 1:N){
    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, get_prob_cor(theta[n], X[n]));
    log_lik_rt[n] = lognormal_lpdf(RT[n] - rt_ndt[S_id[n]] | rt_int[S_id[n]] + rt_slope[S_id[n]] * entropy_t[n], rt_prec[S_id[n]]);
    log_lik_conf[n] = ord_beta_reg_lpdf(Conf[n] | logit(conf_mu[n]) + meta_bias[S_id[n]], conf_prec[S_id[n]], c0[S_id[n]], c11[S_id[n]]);
    log_lik[n] = log_lik_bin[n] + log_lik_rt[n] + log_lik_conf[n] + log_lik_cop[n];
  }



}
