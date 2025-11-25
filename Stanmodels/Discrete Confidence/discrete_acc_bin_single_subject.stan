functions {
  real psycho_ACC(real x, real alpha, real beta, real lapse){
    return (0.5+0.5*((1-2*lapse) * inv_logit(beta * (x - alpha))));
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
  // Accuracy bounds
  matrix uvar_bounds(array[] int binom_y, real alpha, real beta, real lapse, vector X, int is_upper) {
    int N = size(binom_y);
    matrix[N, 1] u_bounds;
    for (n in 1:N) {
      real theta = psycho_ACC(X[n], alpha, exp(beta), inv_logit(lapse) / 2);
      if (is_upper == 0) {
        u_bounds[n, 1] = binom_y[n] == 0 ? 0.0 : binomial_cdf(binom_y[n] - 1 | 1, theta);
      } else {
        u_bounds[n, 1] = binomial_cdf(binom_y[n] | 1, theta);
      }
    }
    return u_bounds;
  }
  // Confidence bounds
  matrix uvar_bounds_conf(array[] int conf, real alpha, real beta, real lapse, real meta_un, real meta_bias, vector cutpoints, vector X, vector ACC, int K, int is_upper) {
    int N = size(conf);
    matrix[N, 1] u_bounds;
    for (n in 1:N) {
      real theta_conf = psycho_ACC(X[n], alpha, exp(beta + meta_un), inv_logit(lapse) / 2);
      real conf_mu = theta_conf; // Or use get_conf if you want
      real eta = logit(conf_mu) + meta_bias;
      if (is_upper == 0) {
        if (conf[n] == 1) {
          u_bounds[n, 1] = 0.0;
        } else {
          u_bounds[n, 1] = inv_logit(cutpoints[conf[n] - 1] - eta);
        }
      } else {
        if (conf[n] == K) {
          u_bounds[n, 1] = 1.0;
        } else {
          u_bounds[n, 1] = inv_logit(cutpoints[conf[n]] - eta);
        }
      }
    }
    return u_bounds;
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  array[N] int binom_y;
  vector[N] RT;
  array[N] int Conf;
  vector[N] X;
  vector[N] ACC;
}

transformed data{
  vector[K-1] cutmeans = linspaced_vector(K-1,-3,3);
  real cutsd = 4.0 / (K-1);
}

parameters {
  real alpha;
  real beta;
  real lapse;
  real rt_int;
  real rt_slope;
  real rt_prec;
  real meta_un;
  real meta_bias;
  ordered[K-1] cutpoints;
  real<lower=0, upper=min(RT)> rt_ndt;
  cholesky_factor_corr[3] rho_chol;

  matrix<
    lower=uvar_bounds(binom_y, alpha, beta, lapse, X, 0),
    upper=uvar_bounds(binom_y, alpha, beta, lapse, X, 1)
  >[N, 1] u;

  matrix<
    lower=uvar_bounds_conf(Conf, alpha, beta, lapse, meta_un, meta_bias, cutpoints, X, ACC, K, 0),
    upper=uvar_bounds_conf(Conf, alpha, beta, lapse, meta_un, meta_bias, cutpoints, X, ACC, K, 1)
  >[N, 1] u_conf;
}

transformed parameters {
  vector[N] theta;
  vector[N] entropy_t;
  vector[N] theta_conf;
  vector[N] conf_mu;
  for (n in 1:N) {
    theta[n] = psycho_ACC(X[n], alpha, exp(beta), inv_logit(lapse) / 2);
    entropy_t[n] = entropy(theta[n]);
    theta_conf[n] = psycho_ACC(X[n], alpha, exp(beta + meta_un), inv_logit(lapse) / 2);
    conf_mu[n] = theta_conf[n];
  }
}

model {
  alpha ~ normal(0, 1);
  beta ~ normal(2, 2);
  lapse ~ normal(-4, 2);
  rt_int ~ normal(0, 1);
  rt_slope ~ normal(0, 2);
  rt_prec ~ normal(0, 1);
  meta_un ~ normal(0, 1);
  meta_bias ~ normal(0, 1);
  cutpoints ~ normal(cutmeans,cutsd);
  rt_ndt ~ normal(0.3, 0.05);
  rho_chol ~ lkj_corr_cholesky(2);

  matrix[N, 3] u_mix;
  for (n in 1:N) {
    u_mix[n, 1] = u[n, 1];
    u_mix[n, 2] = lognormal_cdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n], exp(rt_prec));
    u_mix[n, 3] = u_conf[n, 1];
    target += lognormal_lpdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n], exp(rt_prec));
  }
  u_mix ~ gauss_copula_cholesky(rho_chol);
}

generated quantities {
  matrix[3, 3] rho = multiply_lower_tri_self_transpose(rho_chol);
}
