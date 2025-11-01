# Plotting functions for ACC_Bin_Hier_Correct.stan model

library(tidyverse)
library(posterior)
library(brms)
library(copula)
library(future.apply)

# Function to generate posterior predictive samples
Get_predictive_ACC_hier = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 10000 * 1024^2

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

      # Get probability correct for each trial
      prob_faster = psycho(x, params$alpha, exp(params$beta), params$lapse)
      prob_cor = get_prob_cor(prob_faster, x)

      # Entropy for RT model
      entropy_t = entropy(prob_faster)

      # Theta for confidence (with meta_un)
      prob_faster_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un), params$lapse)




      # Create correlation matrix from rho parameters
      rho_matrix = matrix(c(1, params$rho_p_rt, params$rho_p_conf,
                           params$rho_p_rt, 1, params$rho_rt_conf,
                           params$rho_p_conf, params$rho_rt_conf, 1),
                         nrow = 3, ncol = 3)

      # Sample from Gaussian copula
      copula_obj = normalCopula(param = c(params$rho_p_rt,
                                         params$rho_p_conf,
                                         params$rho_rt_conf),
                               dim = 3, dispstr = "un")

      u_samples = rCopula(n_trials, copula_obj)

      # Transform uniforms to marginal distributions
      ACC_pred = rbinom(length(prob_cor),1,prob_cor)
      bin_pred = rbinom(length(prob_faster),1,prob_faster)


      # 2. RT (lognormal)
      rt_mu = params$rt_int + params$rt_slope * entropy_t + params$rt_stim * x
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
                                   c0 = params$c0,
                                   c11 = params$c11)

      predictions = data.frame(
        X = x,
        prob_cor = prob_cor,
        ACC_pred = ACC_pred,
        prob_faster =prob_faster,
        bin_pred = bin_pred,
        RT_pred = RT_pred,
        rt_mu = rt_mu_expected,

        conf_pred_correct = conf_pred_correct,
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


# Function to generate group-level predictions (average across subjects)
Get_predictive_group = function(fit, df, n_draws = 50) {

  df$subject = as.numeric(as.factor(df$subject))

  workers = 7
  memory = 10000 * 1024^2

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

    # Probability correct
    theta = psycho_ACC(x, params$alpha, exp(params$beta),
                      brms::inv_logit_scaled(params$lapse) / 2)
    prob_cor = get_prob_cor(theta, x)

    # Entropy for RT model
    entropy_t = entropy(theta)

    # Confidence theta (with meta_un)
    theta_conf = psycho(x, params$alpha, exp(params$beta + params$meta_un),
                           brms::inv_logit_scaled(params$lapse) / 2)

    # Build correlation matrix from averaged copula parameters
    cop = normalCopula(param = c(constants$rho_p_rt,
                                 constants$rho_p_conf,
                                 constants$rho_rt_conf),
                      dim = 3, dispstr = "un")

    # Generate correlated uniform samples using Gaussian copula
    u_samples = rCopula(length(x), cop)

    # Transform uniforms to each marginal distribution

    # 1. Decision (accuracy) - whether response was correct
    ACC_pred = as.numeric(u_samples[, 1] < prob_cor)

    # 2. RT (lognormal)
    rt_mu = params$rt_int + params$rt_slope * entropy_t + params$rt_stim * x
    RT_pred = qlnorm(u_samples[, 2], meanlog = rt_mu, sdlog = exp(params$rt_prec)) + constants$rt_ndt

    # Expected RT mean
    rt_mu_expected = exp(rt_mu + exp(params$rt_prec)^2 / 2) + constants$rt_ndt

    # 3. Confidence (ordered beta)
    # Get confidence mean based on whether response was correct (ACC_pred)
    conf_mu_raw = get_conf(ACC_pred, theta_conf, x, params$alpha)


    # Apply meta_bias in logit space
    conf_mu_biased = brms::inv_logit_scaled(qlogis(conf_mu_raw) + params$meta_bias)

    Conf_pred = qordbeta(u_samples[, 3],
                        mu = conf_mu_logit,
                        phi = params$conf_prec,
                        c0 = constants$c0,
                        c11 = constants$c11)

    predictions = data.frame(
      X = x,
      Correct = ACC_pred,  # 1 if correct, 0 if incorrect
      prob = prob_cor,
      RT_pred = RT_pred,
      conf_mu_actual = brms::inv_logit_scaled(conf_mu_logit),
      conf_mu_biased = conf_mu_biased,
      rt_mu = rt_mu_expected,
      Confidence = Conf_pred,
      theta = theta,
      draw = d
    )

    return(predictions)
  }, future.seed = TRUE)

  # Flatten nested list and create a tidy long dataframe
  return(predictions = map_dfr(pred_list, bind_rows))
}


# Function to plot group-level predictions with ribbons
Plot_group_predictive = function(predictions, df) {

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
                mean = mean(prob),
                q5 = quantile(prob, 0.05),
                q10 = quantile(prob, 0.1),
                q20 = quantile(prob, 0.2),
                q95 = quantile(prob, 0.95),
                q90 = quantile(prob, 0.90),
                q80 = quantile(prob, 0.80),
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

  # Plot 1: Expected means
  plot_mean = predictionsq_mean %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_pointrange(data = dataq, aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct),
                    shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free", ncol = 3) +
    theme_classic(base_size = 16) +
    labs(color = "Correct", fill = "Correct",
         title = "Group predictions (expected means)") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")

  plot_mean

  # Prepare predicted data (using actual samples)
  predictionsq_preds = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob),
                q5 = quantile(prob, 0.05),
                q10 = quantile(prob, 0.1),
                q20 = quantile(prob, 0.2),
                q95 = quantile(prob, 0.95),
                q90 = quantile(prob, 0.90),
                q80 = quantile(prob, 0.80),
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
  plot_preds
  return(list(plot_mean = plot_mean, plot_preds = plot_preds))
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
