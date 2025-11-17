# Plotting functions for ACC_Bin_Hier_Correct.stan model

library(tidyverse)
library(posterior)
library(brms)
library(future.apply)
library(patchwork)

# Source utility functions
source(here::here("Analysis/Functions/Correct_models/utility.R"))

# Function to generate posterior predictive samples
Get_predictive_ACC_hier = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 25000 * 1024^2

  # Parameters in the model
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec", "rt_ndt",
                 "conf_prec", "meta_un", "meta_bias","rt_stim",
                 "c0", "c11",
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

      X_scaled = data$X_scaled

      # Get probability correct for each trial
      prob_faster = psycho(x, params$alpha, exp(params$beta), params$lapse)
      prob_cor = get_prob_cor(prob_faster, x)

      # Entropy for RT model
      entropy_t = entropy(prob_faster)

      # Theta for confidence (with meta_un)
      prob_faster_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un), params$lapse)




      # Build correlation matrix from subject parameters
      R = matrix(c(1, params$rho_p_rt, params$rho_p_conf,
                   params$rho_p_rt, 1, params$rho_rt_conf,
                   params$rho_p_conf, params$rho_rt_conf, 1), nrow = 3, byrow = TRUE)

      # Sample multivariate normals
      z_samples = MASS::mvrnorm(n = n_trials, mu = rep(0, 3), Sigma = R)

      # Transform to uniform [0,1] via standard normal CDF
      u_samples = pnorm(z_samples)  # N x 3
      # Transform uniforms to each marginal distribution


      # Transform uniforms to marginal distributions
      ACC_pred = rbinom(length(prob_cor),1,prob_cor)
      bin_pred = rbinom(length(prob_faster),1,prob_faster)


      # 2. RT (lognormal)
      rt_mu = params$rt_int + params$rt_slope * entropy_t + params$rt_stim * X_scaled
      RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = params$rt_prec) + params$rt_ndt

      # Calculate expected RT mean
      rt_mu_expected = exp(rt_mu + params$rt_prec^2 / 2) + params$rt_ndt


      # 3 Confidence mean (probability of getting it correct from confidence)
      prob_cor_conf = get_conf(ACC_pred, prob_faster_conf, x, params$alpha)


      # Apply meta_bias in logit space
      conf_mu_correct = brms::inv_logit_scaled(qlogis(prob_cor_conf) + params$meta_bias)


      conf_pred_correct = qordbeta(u_samples[, 3],
                                   mu = conf_mu_correct,
                                   phi = params$conf_prec,
                                   cutzero = params$c0,
                                   cutone = exp(params$c0) + params$c11)


      # Add observed data for residual calculation
      predictions = data.frame(
        X = x,
        prob_cor = prob_cor,
        ACC_pred = ACC_pred,
        prob_faster = prob_faster,
        bin_pred = bin_pred,
        RT_pred = RT_pred,
        rt_mu = rt_mu_expected,
        conf_pred_correct = conf_pred_correct,
        conf_mu_correct = conf_mu_correct,

        # Observed values
        Y_obs = data$Y,
        RT_obs = data$RT,
        Correct_obs = data$Correct,
        Confidence_obs = data$Confidence,

        # Residuals (observed - predicted)
        residual_Y = data$Correct - prob_cor,
        residual_RT = data$RT - rt_mu_expected,
        residual_Confidence = data$Confidence - conf_mu_correct,

        draw = d,
        subject = s
      )

      return(predictions)
    })
  }, future.seed = TRUE)

  # Flatten nested list and create a tidy long dataframe
  return(map_dfr(pred_list, bind_rows))
}


# Function to plot psychometric curves with predictions
Plot_psychometric_ACC_hier = function(predictions, df, bin = 7, draws = F) {

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
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)   +
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(psychometric_plot)

  }else{
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
      group_by(X, subject) %>%
      summarize(p_correct = mean(prob_cor), .groups = "drop",
                q5 = quantile(prob_cor, 0.05),
                q10 = quantile(prob_cor, 0.1),
                q20 = quantile(prob_cor, 0.2),
                q95 = quantile(prob_cor, 0.95),
                q90 = quantile(prob_cor, 0.90),
                q80 = quantile(prob_cor, 0.80))

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_correct,
                          ymin = p_correct - 2*se, ymax = p_correct + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_correct),
                alpha = 0.2, color = "red") +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_correct, ymin = q5, ymax = q95), alpha = 0.1) +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_correct, ymin = q10, ymax = q90), alpha = 0.3) +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_correct, ymin = q20, ymax = q80), alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)   +
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(psychometric_plot)


  }
}

# Function to plot psychometric curves with predictions
Plot_psychometric_resp_hier = function(predictions, df, bin = 7, draws = F) {

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
      summarize(p_faster = mean(resp ),
                se = sqrt(p_faster * (1 - p_faster) / n()),
                .groups = "drop")

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, draw) %>%
      summarize(p_faster = mean(prob_faster), .groups = "drop")

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_faster,
                          ymin = p_faster - 2*se, ymax = p_faster + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_faster, group = draw),
                alpha = 0.2, color = "red") +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)   +
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(psychometric_plot)

  }else{
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
      summarize(p_faster = mean(resp ),
                se = sqrt(p_faster * (1 - p_faster) / n()),
                .groups = "drop")


    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject) %>%
      summarize(p_faster = mean(prob_faster), .groups = "drop",
                q5 = quantile(prob_faster, 0.05),
                q10 = quantile(prob_faster, 0.1),
                q20 = quantile(prob_faster, 0.2),
                q95 = quantile(prob_faster, 0.95),
                q90 = quantile(prob_faster, 0.90),
                q80 = quantile(prob_faster, 0.80))

    # Plot
    psychometric_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = p_faster,
                          ymin = p_faster - 2*se, ymax = p_faster + 2*se),
                      size = 0.3) +
      geom_line(data = pred_summary,
                aes(x = X, y = p_faster),
                alpha = 0.2, color = "red") +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_faster, ymin = q5, ymax = q95), alpha = 0.1) +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_faster, ymin = q10, ymax = q90), alpha = 0.3) +
      geom_ribbon(data = pred_summary,aes(x = X, y = p_faster, ymin = q20, ymax = q80), alpha = 0.5) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "P(Correct)",
           title = "Psychometric curves") +
      geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)   +
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(psychometric_plot)


  }
}

# Function to plot RT patterns with predictions
Plot_RT_hier = function(predictions, df, bin = 7, draws = F) {

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
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)+
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(rt_plot)
  }else{
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
      geom_ribbon(data = pred_summary,aes(x = X, y = RT_mean, ymin = q5, ymax = q95), alpha = 0.1, fill = "blue") +
      geom_ribbon(data = pred_summary,aes(x = X, y = RT_mean, ymin = q10, ymax = q90), alpha = 0.3, fill = "blue") +
      geom_ribbon(data = pred_summary,aes(x = X, y = RT_mean, ymin = q20, ymax = q80), alpha = 0.5, fill = "blue") +
      facet_wrap(~subject, scales = "free") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Response Time",
           title = "RT by stimulus strength") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5)+
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

    return(rt_plot)



  }
}


# Function to plot confidence by accuracy with predictions
Plot_Conf_ACC_hier = function(predictions, df, bin = 7, draws = F) {

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
      summarize(Conf_mean = mean(Confidence),
                se = sd(Confidence) / sqrt(n()),
                .groups = "drop") %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred , draw) %>%
      summarize(Conf_mean = mean(conf_mu_correct ), .groups = "drop") %>%
      rename(Correct = ACC_pred) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))


    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                     aes(x = X_mid, y = Conf_mean,
                         ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                         color = as.factor(Correct)),
                     size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
               aes(x = X, y = Conf_mean, group = interaction(draw, Correct),
                   color = as.factor(Correct)),
               alpha = 0.2) +
      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Confidence",
           title = "Confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top")+
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )

  return(conf_plot)
  }else{

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
                se = sd(Confidence) / sqrt(n()),
                .groups = "drop") %>%       mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Predicted data summary
    pred_summary <- predictions %>%
      group_by(X, subject, ACC_pred ) %>%
      summarize(Conf_mean = mean(conf_mu_correct ),
                q5 = quantile(conf_mu_correct , 0.05),
                q10 = quantile(conf_mu_correct , 0.1),
                q20 = quantile(conf_mu_correct , 0.2),
                q95 = quantile(conf_mu_correct , 0.95),
                q90 = quantile(conf_mu_correct , 0.90),
                q80 = quantile(conf_mu_correct , 0.80),
                .groups = "drop") %>%
      rename(Correct = ACC_pred ) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Plot
    conf_plot <- ggplot() +
      geom_pointrange(data = obs_summary,
                      aes(x = X_mid, y = Conf_mean,
                          ymin = Conf_mean - 2*se, ymax = Conf_mean + 2*se,
                          color = (Correct)),
                      size = 0.3, position = position_dodge(width = 0.3)) +
      geom_line(data = pred_summary,
                aes(x = X, y = Conf_mean, group = Correct,
                    color = (Correct)),
                alpha = 0.2) +
      geom_ribbon(data = pred_summary,aes(x = X, y = Conf_mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
      geom_ribbon(data = pred_summary,aes(x = X, y = Conf_mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
      geom_ribbon(data = pred_summary,aes(x = X, y = Conf_mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +

      facet_wrap(~subject, scales = "free_x") +
      theme_classic(base_size = 16) +
      labs(x = "Stimulus strength (X)", y = "Confidence",
           title = "Confidence by accuracy and stimulus strength",
           color = "Correct") +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
      theme(legend.position = "top")+
      theme_classic(base_size = 16) +
      theme(
        strip.background = element_blank(),  # remove facet boxes
        strip.text = element_blank(),        # remove facet labels
        axis.text = element_blank(),         # remove x and y axis numbers
        axis.ticks = element_blank(),        # remove tick marks
        axis.title = element_blank()         # remove axis titles
      )
    return(conf_plot)


  }
}

Plot_residuals = function(predictions){

  # Aggregate residuals by subject and X for plotting
  residual_summary <- predictions %>%
    group_by(subject, X) %>%
    summarize(
      mean_resid_Y = mean(residual_Y),
      mean_resid_RT = mean(residual_RT),
      q5_resid_Y = quantile(residual_Y,0.05),
      q95_resid_Y = quantile(residual_Y,0.95),
      q5_resid_RT = quantile(residual_RT,0.05),
      q95_resid_RT = quantile(residual_RT,0.95),
      .groups = "drop"
    )

  # Plot 1: Accuracy/Type-1 residuals
  plot_resid_Y <- residual_summary %>%
    ggplot(aes(x = X, y = mean_resid_Y)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_Y,
                        ymax = q95_resid_RT),
                    size = 0.3, alpha = 0.7) +
    facet_wrap(~subject, scales = "free") +
    labs(y = "Accuracy Residual", x = "Stimulus strength (X)",
         title = "Subject-level Accuracy residuals") +
    scale_y_continuous("Accuracy Residual", breaks = c(-0.5,0,0.5),labels = c(-0.5,0,0.5))+
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "top"
    )
  # Plot 2: RT residuals
  plot_resid_RT <- residual_summary %>%
    ggplot(aes(x = X, y = mean_resid_RT)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_RT,
                        ymax = q95_resid_RT),
                    size = 0.3, alpha = 0.7) +
    facet_wrap(~subject, scales = "free") +
    labs(y = "RT Residual", x = "Stimulus strength (X)",
         title = "Subject-level RT residuals") +
    scale_y_continuous("RT Residual", breaks = c(-0.5,0,0.5),labels = c(-0.5,0,0.5))+
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "top"
    )


  # Plot 3: Confidence residuals (split by Correct)
  residual_conf_summary <- predictions %>%
    group_by(subject, X, Correct_obs) %>%
    summarize(
      mean_resid_Conf = mean(residual_Confidence),
      q5_resid_Conf = quantile(residual_Confidence,0.05),
      q95_resid_Conf = quantile(residual_Confidence,0.95),
      .groups = "drop"
    ) %>%
    mutate(Correct = ifelse(Correct_obs == 1, "Correct", "Incorrect"))

  plot_resid_Conf <- residual_conf_summary %>%
    ggplot(aes(x = X, y = mean_resid_Conf, color = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_Conf,
                        ymax = q95_resid_Conf),
                    size = 0.3, alpha = 0.7,
                    position = position_dodge(width = 0.5)) +
    facet_wrap(~subject, scales = "free") +
    labs(x = "Stimulus strength (X)",
         title = "Subject-level Confidence residuals",
         color = "Accuracy") +
    scale_y_continuous("Confidence Residual", breaks = c(-0.5,0,0.5),labels = c(-0.5,0,0.5))+
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "top"
    )


  # group
  residual_summary_group <- predictions %>%
    group_by(X) %>%
    summarize(
      mean_resid_Y = mean(residual_Y),
      mean_resid_RT = mean(residual_RT),
      q5_resid_Y = quantile(residual_Y,0.05),
      q95_resid_Y = quantile(residual_Y,0.95),
      q5_resid_RT = quantile(residual_RT,0.05),
      q95_resid_RT = quantile(residual_RT,0.95),
      .groups = "drop"
    )

  plot_resid_Y_group <- residual_summary_group %>%
    ggplot(aes(x = X, y = mean_resid_Y)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_Y,
                        ymax = q95_resid_RT),
                    size = 0.3, alpha = 0.7) +
    labs(y = "Accuracy Residual", x = "Stimulus strength (X)",
         title = "Subject-level Accuracy residuals") +
    scale_y_continuous("Accuracy Residual", breaks = c(-1,0,1),labels = c(-1,0,1))+
    theme_classic(base_size = 16) +
    theme(
      legend.position = "top"
    )

  # Plot 2: RT residuals
  plot_resid_RT_group <- residual_summary_group %>%
    ggplot(aes(x = X, y = mean_resid_RT)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_RT,
                        ymax = q95_resid_RT),
                    size = 0.3, alpha = 0.7) +
    labs(y = "RT Residual", x = "Stimulus strength (X)",
         title = "Subject-level RT residuals") +
    scale_y_continuous("RT Residual", breaks = c(-1,0,1),labels = c(-1,0,1))+
    theme_classic(base_size = 16) +
    theme(
      legend.position = "top"
    )



  # Plot 3: Confidence residuals (split by Correct)
  residual_conf_summary <- predictions %>%
    group_by(X, Correct_obs) %>%
    summarize(
      mean_resid_Conf = mean(residual_Confidence),
      q5_resid_Conf = quantile(residual_Confidence,0.05),
      q95_resid_Conf = quantile(residual_Confidence,0.95),
      .groups = "drop"
    ) %>%
    mutate(Correct = ifelse(Correct_obs == 1, "Correct", "Incorrect"))

  plot_resid_Conf_group <- residual_conf_summary %>%
    ggplot(aes(x = X, y = mean_resid_Conf, color = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = q5_resid_Conf,
                        ymax = q95_resid_Conf),
                    size = 0.3, alpha = 0.7,
                    position = position_dodge(width = 0.5)) +
    labs(x = "Stimulus strength (X)",
         title = "Subject-level Confidence residuals",
         color = "Accuracy") +
    scale_y_continuous("Confidence Residual", breaks = c(-0.5,0,0.5),labels = c(-0.5,0,0.5))+
    theme_classic(base_size = 16) +
    theme(
      legend.position = "top"
    )

  return(list(
    accuracy = plot_resid_Y,
    RT = plot_resid_RT,
    confidence = plot_resid_Conf,
    residual_summary = residual_summary,
    residual_conf_summary = residual_conf_summary
  ))
}


# Function to generate group-level predictions (average across subjects)
Get_predictive_group = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 25000 * 1024^2

  # Group-level parameters (from gm)
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec",
                 "conf_prec", "meta_un", "meta_bias","rt_stim")

  # Extract group means (gm)
  df_param = as_draws_df(fit$draws("gm")) %>%
    select(-contains(".")) %>%
    rename_with(~parameters) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw, names_to = "variable")

  # Get average values for subject-specific parameters
  constants = as_draws_df(fit$draws(c("rt_ndt", "c0", "c11",
                                      "rho_p_rt", "rho_p_conf", "rho_rt_conf"))) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable", "subject"),
            regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE) %>%
    group_by(variable) %>%
    summarize(mean = mean(value)) %>%
    pivot_wider(names_from = variable, values_from = mean)

  # Set up parallel processing
  plan(multisession, workers = workers)
  options(future.globals.maxSize = memory)

  draws <- 1:n_draws

  # Only use the number of draws that the user wants
  dfq = df_param %>% filter(draw %in% draws)

  # Function to get the draws
  pred_list <- future_lapply(draws, function(d) {

    # Extract parameter vectors for this draw
    params <- dfq %>%
      filter(draw == d) %>%
      select(variable, value) %>%
      pivot_wider(names_from = "variable", values_from = "value")

    # Generate X values
    x = seq(-40, 40, by = 0.5)
    n_trials = length(x)
    X_scaled = (x - mean(df$X)) / sd(df$X)

    # Get probability correct for each trial (using same approach as subject-level)
    prob_faster = psycho(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse) / 2)
    prob_cor = get_prob_cor(prob_faster, x)

    # Entropy for RT model
    entropy_t = entropy(prob_faster)

    # Theta for confidence (with meta_un)
    prob_faster_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un), brms::inv_logit_scaled(params$lapse) / 2)

    # Build correlation matrix from averaged copula parameters
    R = matrix(c(1, constants$rho_p_rt, constants$rho_p_conf,
                 constants$rho_p_rt, 1, constants$rho_rt_conf,
                 constants$rho_p_conf, constants$rho_rt_conf, 1),
               nrow = 3, byrow = TRUE)

    # Sample multivariate normals (same as subject-level code)
    z_samples = MASS::mvrnorm(n = n_trials, mu = rep(0, 3), Sigma = R)

    # Transform to uniform [0,1] via standard normal CDF
    u_samples = pnorm(z_samples)

    # Transform uniforms to marginal distributions
    ACC_pred = rbinom(length(prob_cor), 1, prob_cor)

    # 2. RT (lognormal)
    rt_mu = params$rt_int + params$rt_slope * entropy_t + params$rt_stim * X_scaled
    RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = exp(params$rt_prec)) + constants$rt_ndt

    # Calculate expected RT mean
    rt_mu_expected = exp(rt_mu + exp(params$rt_prec)^2 / 2) + constants$rt_ndt

    # 3 Confidence mean (probability of getting it correct from confidence)
    prob_cor_conf = get_conf(ACC_pred, prob_faster_conf, x, params$alpha)

    # Apply meta_bias in logit space
    conf_mu_correct = brms::inv_logit_scaled(qlogis(prob_cor_conf) + params$meta_bias)

    # Sample confidence values
    conf_pred_correct = qordbeta(u_samples[, 3],
                                mu = conf_mu_correct,
                                phi = exp(params$conf_prec),
                                cutzero = constants$c0,
                                cutone = exp(constants$c0) + constants$c11)

    predictions = data.frame(
      X = x,
      Correct = ACC_pred,  # 1 if correct, 0 if incorrect
      prob = prob_cor,
      RT_pred = RT_pred,
      conf_mu_actual = conf_mu_correct,
      rt_mu = rt_mu_expected,
      Confidence = conf_pred_correct,
      prob_faster = prob_faster,
      draw = d
    )

    return(predictions)
  }, future.seed = TRUE)

  # Flatten nested list and create a tidy long dataframe
  return(predictions = map_dfr(pred_list, bind_rows))
}


# Function to plot group-level predictions with ribbons
Plot_group_predictive_psycho = function(predictions, df) {

  # Prepare observed data
  dataq = bind_rows(
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        n = n(),
        k = sum(Y),
        mean = (1 + k) / (2 + n),
        q5  = qbeta(0.05, 1 + k, 1 + n - k),
        q10 = qbeta(0.10, 1 + k, 1 + n - k),
        q20 = qbeta(0.20, 1 + k, 1 + n - k),
        q80 = qbeta(0.80, 1 + k, 1 + n - k),
        q90 = qbeta(0.90, 1 + k, 1 + n - k),
        q95 = qbeta(0.95, 1 + k, 1 + n - k),
        .groups = "drop"
      ),

    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT),
                q5 = mean(RT) - 2 * (sd(RT) / sqrt(n())),
                q95 = mean(RT) + 2 * (sd(RT) / sqrt(n())),
                .groups = "drop"),

    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence),
                q5 = mean(Confidence) - 2 * (sd(Confidence) / sqrt(n())),
                q95 = mean(Confidence) + 2 * (sd(Confidence) / sqrt(n())),
                .groups = "drop")
  ) %>%
    filter(abs(X) < 25)

  # Prepare predicted data (using expected means)
  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob_faster ),
                q5 = quantile(prob_faster , 0.05),
                q10 = quantile(prob_faster , 0.1),
                q20 = quantile(prob_faster , 0.2),
                q95 = quantile(prob_faster , 0.95),
                q90 = quantile(prob_faster , 0.90),
                q80 = quantile(prob_faster , 0.80),
                .groups = "drop"),

    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(rt_mu),
                q5 = quantile(rt_mu, 0.05),
                q10 = quantile(rt_mu, 0.1),
                q20 = quantile(rt_mu, 0.2),
                q95 = quantile(rt_mu, 0.95),
                q90 = quantile(rt_mu, 0.90),
                q80 = quantile(rt_mu, 0.80),
                .groups = "drop"),

    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(conf_mu_actual),
                q5 = quantile(conf_mu_actual, 0.05),
                q10 = quantile(conf_mu_actual, 0.1),
                q20 = quantile(conf_mu_actual, 0.2),
                q95 = quantile(conf_mu_actual, 0.95),
                q90 = quantile(conf_mu_actual, 0.90),
                q80 = quantile(conf_mu_actual, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < 25) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))

  # Calculate trial-level residuals properly
  # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
  pred_mean_per_trial = predictions %>%
    group_by(X) %>%
    summarize(
      pred_prob_faster = mean(prob_faster),
      pred_rt = mean(rt_mu),
      .groups = "drop"
    )

  # Join predictions to actual trial-level data
  df_with_pred = df %>%
    filter(abs(X) < 25) %>%
    left_join(pred_mean_per_trial, by = "X") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

  # Confidence predictions need to be split by Correct
  pred_conf_per_trial = predictions %>%
    group_by(X, Correct) %>%
    summarize(pred_conf = mean(conf_mu_actual), .groups = "drop") %>%
    mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

  df_with_pred = df_with_pred %>%
    left_join(pred_conf_per_trial %>% select(X, Correct_label, pred_conf),
              by = c("X", "Correct_label"))

  # Calculate residuals at trial level, then aggregate
  residuals_data = bind_rows(
    df_with_pred %>%
      mutate(residual = Y - pred_prob_faster,
             name = "Type-1") %>%
      group_by(X) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop"),

    df_with_pred %>%
      mutate(residual = RT - pred_rt,
             name = "RT") %>%
      group_by(X) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop"),

    df_with_pred %>%
      mutate(residual = Confidence - pred_conf,
             name = "Confidence") %>%
      group_by(X, Correct_label) %>%
      summarize(residual_mean = mean(residual, na.rm = TRUE),
                residual_se = sd(residual, na.rm = TRUE) / sqrt(n()),
                name = first(name),
                .groups = "drop") %>%
      rename(Correct = Correct_label)
  )

  # Plot 1: Expected means (main plot)
  plot_mean = predictionsq_mean %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_pointrange(data = dataq, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 16) +
    labs(color = "Correct", fill = "Correct",
         title = "Group predictions (expected means)",
         y = "Value") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  # Plot 2: Residuals
  plot_residuals = residuals_data %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    geom_pointrange(aes(ymin = residual_mean - 2*residual_se,
                        ymax = residual_mean + 2*residual_se),
                    alpha = 0.5, size = 0.3) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 16) +
    labs(x = "Stimulus strength (X)", y = "Residual (Obs - Pred)") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "none")

  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))

  combined_plot

  # Prepare predicted data (using actual samples)
  predictionsq_preds = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob_faster ),
                q5 = quantile(prob_faster , 0.05),
                q10 = quantile(prob_faster , 0.1),
                q20 = quantile(prob_faster , 0.2),
                q95 = quantile(prob_faster , 0.95),
                q90 = quantile(prob_faster , 0.90),
                q80 = quantile(prob_faster , 0.80),
                .groups = "drop"),

    predictions %>%
      group_by(X) %>%
      summarize(name = "RT",
                mean = mean(RT_pred),
                q5 = quantile(RT_pred, 0.05),
                q10 = quantile(RT_pred, 0.1),
                q20 = quantile(RT_pred, 0.2),
                q95 = quantile(RT_pred, 0.95),
                q90 = quantile(RT_pred, 0.90),
                q80 = quantile(RT_pred, 0.80),
                .groups = "drop"),

    predictions %>%
      group_by(X, Correct) %>%
      summarize(name = "Confidence",
                mean = mean(Confidence),
                q5 = quantile(Confidence, 0.05),
                q10 = quantile(Confidence, 0.1),
                q20 = quantile(Confidence, 0.2),
                q95 = quantile(Confidence, 0.95),
                q90 = quantile(Confidence, 0.90),
                q80 = quantile(Confidence, 0.80),
                .groups = "drop")
  ) %>%
    filter(abs(X) < 25) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))


  # Plot 2: Actual predictions
  plot_preds = predictionsq_preds %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_pointrange(data = dataq, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = as.factor(Correct)),
                    shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free", ncol = 3) +
    theme_classic(base_size = 16) +
    labs(color = "Correct", fill = "Correct",
         title = "Group predictions (posterior predictive samples)") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")

  return(list(
    plot_combined = combined_plot,
    plot_mean = plot_mean,
    plot_residuals = plot_residuals,
    plot_preds = plot_preds
  ))
}


# Combined plotting function that returns all plots
Plot_all_ACC_hier = function(fit, df, n_draws = 50, bin = 7) {

  # Generate predictions
  cat("Generating predictions...\n")
  predictions = Get_predictive_ACC_hier(fit, df, n_draws)

  # Create plots
  cat("Creating psychometric plot...\n")
  p1 = Plot_psychometric_ACC_hier(predictions, df, bin)

  cat("Creating RT plot...\n")
  p2 = Plot_RT_ACC_hier(predictions, df, bin)

  cat("Creating confidence plot...\n")
  p3 = Plot_Conf_ACC_hier(predictions, df, bin)

  cat("Done!\n")

  return(list(
    predictions = predictions,
    psychometric = p1,
    RT = p2,
    confidence = p3
  ))
}
