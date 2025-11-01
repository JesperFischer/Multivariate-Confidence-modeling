# Helper functions
psycho = function(x, alpha, beta, lapse) {
  lapse + (1 - 2 * lapse) * (brms::inv_logit_scaled(beta * (x - alpha)))
}

entropy = function(p) {
  -p * log(p) - (1-p) * log(1-p)
}

get_conf = function(ACC, theta, x, alpha) {
  # ACC=1 means correct response
  # If x > alpha and ACC=1 (correct), return theta
  # If x > alpha and ACC=0 (incorrect), return 1-theta
  # If x < alpha and ACC=1 (correct), return 1-theta
  # If x < alpha and ACC=0 (incorrect), return theta
  ifelse(ACC == 1 & x > alpha, theta,
         ifelse(ACC == 1 & x < alpha, 1 - theta,
                ifelse(ACC == 0 & x > alpha, 1 - theta,
                       ifelse(ACC == 0 & x < alpha, theta, 0.5))))
}

get_prob_cor = function(theta, x) {
  # If x > 0, P(correct) = theta
  # If x < 0, P(correct) = 1-theta
  ifelse(x > 0, theta,
         ifelse(x < 0, 1 - theta, 0.5))
}


# Helper function for ordered beta inverse CDF
qordbeta = function(p, mu, phi, c0, c11) {
  cutpoints = c(c0, c0 + exp(c11))

  # Probability of 0
  p0 = 1 - brms::inv_logit_scaled(mu - cutpoints[1])
  # Probability of 1
  p1 = brms::inv_logit_scaled(mu - cutpoints[2])
  # Probability of middle region
  pm = brms::inv_logit_scaled(mu - cutpoints[1]) - brms::inv_logit_scaled(mu - cutpoints[2])

  # Beta parameters for middle region
  mu_beta = brms::inv_logit_scaled(mu)
  alpha_beta = mu_beta * phi
  beta_beta = (1 - mu_beta) * phi

  # Inverse CDF
  ifelse(p < p0, 0,
         ifelse(p > (1 - p1), 1,
                qbeta((p - p0) / pm, alpha_beta, beta_beta)))
}
