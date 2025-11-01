# Plotting functions for quadmodel.stan and discrete confidence models

library(tidyverse)
library(posterior)
library(brms)
library(copula)
library(future.apply)

# Source utility functions
source("Analysis/Functions/Correct_models/utility.R")


# ============================================================================
# QUAD MODEL FUNCTIONS (Binary Confidence)
# ============================================================================

# Function to generate posterior predictive samples for quad model
Get_predictive_quad = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 10000 * 1024^2

  # Parameters in the quad model
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec", "rt_ndt",
                 "meta_un", "meta_bias",
                 "rho_p_rt", "rho_p_conf", "rho_rt_conf")

  # Extract draws
  df_param = as_draws_df(fit$draws(parameters)) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable", "subject"),
            regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE)

  # Set up parallel processing
  plan(multisession, workers = workers)
  options(future.globals.maxSize = memory)

  subjects <- unique(df$subject)
  draws <- 1:n_draws

  # Only use the number of draws that the user wants:
  dfq = df_param %>% filter(draw %in% draws)

  # Function to generate predictions per subject and draw
  pred_list <- future_lapply(subjects, function(s) {
    lapply(draws, function(d) {

      # Extract parameters for this subject and draw
      params <- dfq %>%
        filter(subject == s, draw == d) %>%
        select(variable, value) %>%
        pivot_wider(names_from = "variable", values_from = "value")

      data = df %>% filter(subject == s)
      x = data$X
      n_trials = length(x)

      # Get probability correct for each trial
      prob_faster = psycho(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse) / 2)
      prob_cor = get_prob_cor(prob_faster, x)

      # Entropy for RT model
      entropy_t = entropy(prob_faster)

      # Theta for confidence (with meta_un)
      prob_faster_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un),
                                brms::inv_logit_scaled(params$lapse) / 2)

      # Create correlation matrix from rho parameters
      copula_obj = normalCopula(param = c(params$rho_p_rt,
                                         params$rho_p_conf,
                                         params$rho_rt_conf),
                               dim = 3, dispstr = "un")

      u_samples = rCopula(n_trials, copula_obj)

      # Transform uniforms to marginal distributions
      # 1. Accuracy
      ACC_pred = rbinom(n_trials, 1, prob_cor)

      # 2. RT (lognormal)
      rt_mu = params$rt_int + params$rt_slope * entropy_t
      RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = params$rt_prec) + params$rt_ndt

      # Calculate expected RT mean
      rt_mu_expected = exp(rt_mu + params$rt_prec^2 / 2) + params$rt_ndt

      # 3. Binary Confidence (probability of high confidence given accuracy)
      prob_cor_conf = get_conf(ACC_pred, prob_faster_conf, x, params$alpha)

      # Apply meta_bias in logit space
      conf_mu_correct = brms::inv_logit_scaled(qlogis(prob_cor_conf) + params$meta_bias)

      # Sample binary confidence
      conf_pred_binary = rbinom(n_trials, 1, conf_mu_correct)

      predictions = data.frame(
        X = x,
        prob_cor = prob_cor,
        ACC_pred = ACC_pred,
        prob_faster = prob_faster,
        RT_pred = RT_pred,
        rt_mu = rt_mu_expected,
        conf_pred_binary = conf_pred_binary,
        conf_mu_correct = conf_mu_correct,
        draw = d,
        subject = s
      )

      return(predictions)
    })
  }, future.seed = TRUE)

  # Flatten nested list and create a tidy long dataframe
  return(map_dfr(pred_list, bind_rows))
}


# Function to plot psychometric curves for quad model
Plot_psychometric_quad = function(predictions, df, bin = 7, draws = F) {

  if(!draws){
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(p_correct = mean(Correct),
                se = sqrt(p_correct * (1 - p_correct) / n()),
                .groups = "drop")

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, draw) %>%
      summarize(p_correct = mean(prob_cor), .groups = "drop")

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_correct,
                          ymin = p_correct - 2*se, ymax = p_correct + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_correct, group = draw),
                alpha = 0.2, color = "red") +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(psychometric_plot)

  } else {
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(p_correct = mean(Correct),
                se = sqrt(p_correct * (1 - p_correct) / n()),
                .groups = "drop")

    # Predicted data summary with quantiles
    pred_summary <- predictions %>%
      group_by(X, subject) %>%
      summarize(p_correct = mean(prob_cor),
                q5 = quantile(prob_cor, 0.05),
                q10 = quantile(prob_cor, 0.1),
                q20 = quantile(prob_cor, 0.2),
                q95 = quantile(prob_cor, 0.95),
                q90 = quantile(prob_cor, 0.90),
                q80 = quantile(prob_cor, 0.80),
                .groups = "drop")

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_correct,
                          ymin = p_correct - 2*se, ymax = p_correct + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_correct),
                alpha = 0.2, color = "red") +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q5, ymax = q95), alpha = 0.1) +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q10, ymax = q90), alpha = 0.3) +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q20, ymax = q80), alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(psychometric_plot)
  }
}


# Function to plot RT for quad model
Plot_RT_quad = function(predictions, df, bin = 7, draws = F) {

  if(!draws){
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(RT_mean = mean(RT),
                se = sd(RT) / sqrt(n()),
                .groups = "drop")

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, draw) %>%
      summarize(RT_mean = mean(rt_mu), .groups = "drop")

    # Plot
    rt_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                     aes(x = X_mid, y = RT_mean,
                         ymin = RT_mean - 2*se, ymax = RT_mean + 2*se),
                     size = 0.3) +
      geom_line(data = pred_summary,
               aes(x = X, y = RT_mean, group = draw),
               alpha = 0.2, color = "blue") +
      facet_wrap(~subject, scales = "free") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Response Time",
           title = "RT by stimulus strength") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(rt_plot)
  } else {
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(RT_mean = mean(RT),
                se = sd(RT) / sqrt(n()),
                .groups = "drop")

    # Predicted data summary with quantiles
    pred_summary <- predictions %>%
      group_by(X, subject) %>%
      summarize(RT_mean = mean(rt_mu),
                q5 = quantile(rt_mu, 0.05),
                q10 = quantile(rt_mu, 0.1),
                q20 = quantile(rt_mu, 0.2),
                q95 = quantile(rt_mu, 0.95),
                q90 = quantile(rt_mu, 0.90),
                q80 = quantile(rt_mu, 0.80),
                .groups = "drop")

    # Plot
    rt_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = RT_mean,
                          ymin = RT_mean - 2*se, ymax = RT_mean + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = RT_mean),
                alpha = 0.2, color = "blue") +
      geom_ribbon(data = pred_summary, aes(x = X, y = RT_mean, ymin = q5, ymax = q95),
                  alpha = 0.1, fill = "blue") +
      geom_ribbon(data = pred_summary, aes(x = X, y = RT_mean, ymin = q10, ymax = q90),
                  alpha = 0.3, fill = "blue") +
      geom_ribbon(data = pred_summary, aes(x = X, y = RT_mean, ymin = q20, ymax = q80),
                  alpha = 0.5, fill = "blue") +
      facet_wrap(~subject, scales = "free") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Response Time",
           title = "RT by stimulus strength") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(rt_plot)
  }
}


# Function to plot binary confidence by accuracy for quad model
Plot_Conf_quad = function(predictions, df, bin = 7, draws = F) {

  if(!draws){
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject, Correct) %>%
      summarize(Conf_mean = mean(Confidence),  # Binary confidence
                se = sqrt(Conf_mean * (1 - Conf_mean) / n()),
                .groups = "drop") %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred, draw) %>%
      summarize(Conf_mean = mean(conf_mu_correct), .groups = "drop") %>%
      rename(Correct = ACC_pred) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                     aes(x = X_mid, y = Conf_mean,
                         ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                         color = Correct),
                     size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
               aes(x = X, y = Conf_mean, group = interaction(draw, Correct),
                   color = Correct),
               alpha = 0.2) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(High Confidence)",
           title = "Binary confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top") +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(conf_plot)
  } else {
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject, Correct) %>%
      summarize(Conf_mean = mean(Confidence),
                se = sqrt(Conf_mean * (1 - Conf_mean) / n()),
                .groups = "drop") %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary with quantiles
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred) %>%
      summarize(Conf_mean = mean(conf_mu_correct),
                q5 = quantile(conf_mu_correct, 0.05),
                q10 = quantile(conf_mu_correct, 0.1),
                q20 = quantile(conf_mu_correct, 0.2),
                q95 = quantile(conf_mu_correct, 0.95),
                q90 = quantile(conf_mu_correct, 0.90),
                q80 = quantile(conf_mu_correct, 0.80),
                .groups = "drop") %>%
      rename(Correct = ACC_pred) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = Conf_mean,
                          ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                          color = Correct),
                      size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
                aes(x = X, y = Conf_mean, group = Correct,
                    color = Correct),
                alpha = 0.2) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q5, ymax = q95, fill = Correct),
                  alpha = 0.1) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q10, ymax = q90, fill = Correct),
                  alpha = 0.3) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q20, ymax = q80, fill = Correct),
                  alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(High Confidence)",
           title = "Binary confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top") +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(conf_plot)
  }
}


# Combined plotting function for quad model
Plot_all_quad = function(fit, df, n_draws = 50, bin = 7, draws = F) {

  # Generate predictions
  cat("Generating predictions for quad model...\n")
  predictions = Get_predictive_quad(fit, df, n_draws)

  # Create plots
  cat("Creating psychometric plot...\n")
  p1 = Plot_psychometric_quad(predictions, df, bin, draws)

  cat("Creating RT plot...\n")
  p2 = Plot_RT_quad(predictions, df, bin, draws)

  cat("Creating binary confidence plot...\n")
  p3 = Plot_Conf_quad(predictions, df, bin, draws)

  cat("Done!\n")

  return(list(
    predictions = predictions,
    psychometric = p1,
    RT = p2,
    confidence = p3
  ))
}


# ============================================================================
# DISCRETE CONFIDENCE MODEL FUNCTIONS
# ============================================================================

# Function to generate posterior predictive samples for discrete confidence model
Get_predictive_discrete = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 10000 * 1024^2

  # Parameters in the discrete confidence model
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec", "rt_ndt",
                 "conf_ACC", "conf_entropy", "conf_entropy_ACC",
                 "meta_un",
                 "rho_p_rt", "rho_p_conf", "rho_rt_conf")

  # Extract draws
  df_param = as_draws_df(fit$draws(parameters)) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable", "subject"),
            regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE)

  # Extract cutpoints (these may have different dimensions per subject)
  cutpoints_raw = as_draws_df(fit$draws("cutpoints")) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable", "subject", "cutpoint"),
            regex = "([a-zA-Z0-9_]+)\\[(\\d+),(\\d+)\\]", convert = TRUE)

  # Set up parallel processing
  plan(multisession, workers = workers)
  options(future.globals.maxSize = memory)

  subjects <- unique(df$subject)
  draws <- 1:n_draws

  # Only use the number of draws that the user wants:
  dfq = df_param %>% filter(draw %in% draws)
  cutq = cutpoints_raw %>% filter(draw %in% draws)

  # Function to generate predictions per subject and draw
  pred_list <- future_lapply(subjects, function(s) {
    lapply(draws, function(d) {

      # Extract parameters for this subject and draw
      params <- dfq %>%
        filter(subject == s, draw == d) %>%
        select(variable, value) %>%
        pivot_wider(names_from = "variable", values_from = "value")

      # Extract cutpoints for this subject and draw
      cuts <- cutq %>%
        filter(subject == s, draw == d) %>%
        arrange(cutpoint) %>%
        pull(value)

      data = df %>% filter(subject == s)
      x = data$X
      n_trials = length(x)

      # Get probability correct for each trial
      prob_faster = psycho_ACC(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse) / 2)
      prob_cor = prob_faster  # For ACC models, this is already P(correct)

      # Entropy for RT and confidence models
      entropy_t = entropy(prob_faster)

      # Theta for confidence (with meta_un)
      prob_faster_conf = psycho_ACC(x, params$alpha, exp(params$beta + params$meta_un),
                                    brms::inv_logit_scaled(params$lapse) / 2)

      # Create correlation matrix from rho parameters
      copula_obj = normalCopula(param = c(params$rho_p_rt,
                                         params$rho_p_conf,
                                         params$rho_rt_conf),
                               dim = 3, dispstr = "un")

      u_samples = rCopula(n_trials, copula_obj)

      # Transform uniforms to marginal distributions
      # 1. Accuracy
      ACC_pred = rbinom(n_trials, 1, prob_cor)

      # 2. RT (lognormal)
      rt_mu = params$rt_int + params$rt_slope * entropy_t
      RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = params$rt_prec) + params$rt_ndt

      # Calculate expected RT mean
      rt_mu_expected = exp(rt_mu + params$rt_prec^2 / 2) + params$rt_ndt

      # 3. Discrete Confidence (ordered categorical)
      # Compute latent confidence variable
      conf_latent = params$conf_ACC * ACC_pred +
                    params$conf_entropy * entropy(prob_faster_conf) +
                    params$conf_entropy_ACC * ACC_pred * entropy(prob_faster_conf)

      # Convert to probabilities using cutpoints
      K = length(cuts) + 1  # Number of categories

      # For each trial, determine confidence category
      conf_pred_discrete = numeric(n_trials)
      conf_probs = matrix(0, nrow = n_trials, ncol = K)

      for (i in 1:n_trials) {
        # Compute cumulative probabilities
        cum_probs = c(0, brms::inv_logit_scaled(cuts - conf_latent[i]), 1)

        # Category probabilities
        for (k in 1:K) {
          conf_probs[i, k] = cum_probs[k+1] - cum_probs[k]
        }

        # Sample discrete confidence
        conf_pred_discrete[i] = sample(1:K, size = 1, prob = conf_probs[i, ])
      }

      # Mean confidence (expected value)
      conf_mean = rowSums(conf_probs * matrix(1:K, nrow = n_trials, ncol = K, byrow = TRUE))

      predictions = data.frame(
        X = x,
        prob_cor = prob_cor,
        ACC_pred = ACC_pred,
        RT_pred = RT_pred,
        rt_mu = rt_mu_expected,
        conf_pred_discrete = conf_pred_discrete,
        conf_mean = conf_mean,
        conf_latent = conf_latent,
        draw = d,
        subject = s
      )

      return(predictions)
    })
  }, future.seed = TRUE)

  # Flatten nested list and create a tidy long dataframe
  return(map_dfr(pred_list, bind_rows))
}


# Helper function for psycho_ACC used in discrete models
psycho_ACC = function(x, alpha, beta, lapse) {
  0.5 + 0.5 * ((1 - 2*lapse) * brms::inv_logit_scaled_scaled(beta * (x - alpha)))
}


# Function to plot psychometric curves for discrete model
Plot_psychometric_discrete = function(predictions, df, bin = 7, draws = F) {

  if(!draws){
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(p_correct = mean(Correct),
                se = sqrt(p_correct * (1 - p_correct) / n()),
                .groups = "drop")

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, draw) %>%
      summarize(p_correct = mean(prob_cor), .groups = "drop")

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_correct,
                          ymin = p_correct - 2*se, ymax = p_correct + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_correct, group = draw),
                alpha = 0.2, color = "red") +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(psychometric_plot)

  } else {
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject) %>%
      summarize(p_correct = mean(Correct),
                se = sqrt(p_correct * (1 - p_correct) / n()),
                .groups = "drop")

    # Predicted data summary with quantiles
    pred_summary <- predictions %>%
      group_by(X, subject) %>%
      summarize(p_correct = mean(prob_cor),
                q5 = quantile(prob_cor, 0.05),
                q10 = quantile(prob_cor, 0.1),
                q20 = quantile(prob_cor, 0.2),
                q95 = quantile(prob_cor, 0.95),
                q90 = quantile(prob_cor, 0.90),
                q80 = quantile(prob_cor, 0.80),
                .groups = "drop")

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_correct,
                          ymin = p_correct - 2*se, ymax = p_correct + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_correct),
                alpha = 0.2, color = "red") +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q5, ymax = q95), alpha = 0.1) +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q10, ymax = q90), alpha = 0.3) +
      geom_ribbon(data = pred_summary, aes(x = X, y = p_correct, ymin = q20, ymax = q80), alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(psychometric_plot)
  }
}


# Function to plot RT for discrete model (same as quad)
Plot_RT_discrete = Plot_RT_quad


# Function to plot discrete confidence by accuracy
Plot_Conf_discrete = function(predictions, df, bin = 7, draws = F) {

  if(!draws){
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject, Correct) %>%
      summarize(Conf_mean = mean(Conf),  # Discrete confidence levels
                se = sd(Conf) / sqrt(n()),
                .groups = "drop") %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred, draw) %>%
      summarize(Conf_mean = mean(conf_mean), .groups = "drop") %>%
      rename(Correct = ACC_pred) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                     aes(x = X_mid, y = Conf_mean,
                         ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                         color = Correct),
                     size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
               aes(x = X, y = Conf_mean, group = interaction(draw, Correct),
                   color = Correct),
               alpha = 0.2) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Mean Confidence",
           title = "Discrete confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top") +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(conf_plot)
  } else {
    df$subject = as.numeric(as.factor(df$subject))

    # Bin the X values
    bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
    bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Add binned X to data
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
             X_mid = bin_mids[as.integer(X_bin)])

    # Observed data summary
    obs_summary <- df %>%
      group_by(X_mid, subject, Correct) %>%
      summarize(Conf_mean = mean(Conf),
                se = sd(Conf) / sqrt(n()),
                .groups = "drop") %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary with quantiles
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred) %>%
      summarize(Conf_mean = mean(conf_mean),
                q5 = quantile(conf_mean, 0.05),
                q10 = quantile(conf_mean, 0.1),
                q20 = quantile(conf_mean, 0.2),
                q95 = quantile(conf_mean, 0.95),
                q90 = quantile(conf_mean, 0.90),
                q80 = quantile(conf_mean, 0.80),
                .groups = "drop") %>%
      rename(Correct = ACC_pred) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = Conf_mean,
                          ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                          color = Correct),
                      size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
                aes(x = X, y = Conf_mean, group = Correct,
                    color = Correct),
                alpha = 0.2) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q5, ymax = q95, fill = Correct),
                  alpha = 0.1) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q10, ymax = q90, fill = Correct),
                  alpha = 0.3) +
      geom_ribbon(data = pred_summary, aes(x = X, y = Conf_mean, ymin = q20, ymax = q80, fill = Correct),
                  alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Mean Confidence",
           title = "Discrete confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top") +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )

    return(conf_plot)
  }
}


# Combined plotting function for discrete confidence model
Plot_all_discrete = function(fit, df, n_draws = 50, bin = 7, draws = F) {

  # Generate predictions
  cat("Generating predictions for discrete confidence model...\n")
  predictions = Get_predictive_discrete(fit, df, n_draws)

  # Create plots
  cat("Creating psychometric plot...\n")
  p1 = Plot_psychometric_discrete(predictions, df, bin, draws)

  cat("Creating RT plot...\n")
  p2 = Plot_RT_discrete(predictions, df, bin, draws)

  cat("Creating discrete confidence plot...\n")
  p3 = Plot_Conf_discrete(predictions, df, bin, draws)

  cat("Done!\n")

  return(list(
    predictions = predictions,
    psychometric = p1,
    RT = p2,
    confidence = p3
  ))
}
