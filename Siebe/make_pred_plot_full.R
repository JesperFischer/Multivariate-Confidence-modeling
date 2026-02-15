make_pred_plot_full = function(fit,data1,number){

  n_draws = 2000
  workers = 7
  memory = 15000 * 1024^2

  # Group-level parameters (from gm)
  parameters = c("alpha", "beta", "lapse",
                 "rt_int", "rt_slope", "rt_prec",
                 "conf_prec", "meta_un", "meta_bias")

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
    extract(name, into = c("variable"),
            regex = "([a-zA-Z0-9_]+)", convert = TRUE) %>%
    group_by(variable) %>%
    summarize(mean = mean(value)) %>%
    pivot_wider(names_from = variable, values_from = mean)


  library(future.apply)
  # Set up parallel processing
  plan(multisession, workers = 2)
  options(future.globals.maxSize = 5000 * 1024^2)

  draws <- 1:n_draws

  # Only use the number of draws that the user wants
  dfq = df_param %>% filter(draw %in% draws)

  source("~/Multivariate-Confidence-modeling/Analysis/Functions/Correct_models/utility.R")

  # Function to get the draws
  pred_list <- future_lapply(draws, function(d) {

    # Extract parameter vectors for this draw
    params <- dfq %>%
      filter(draw == d) %>%
      select(variable, value) %>%
      pivot_wider(names_from = "variable", values_from = "value")

    # Generate X values
    x = seq(-1, 1, by = 0.01)
    n_trials = length(x)

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
    rt_mu = params$rt_int + params$rt_slope * entropy_t
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

  predictions = map_dfr(pred_list, bind_rows)


  Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {

    cutoff = 2

    dataq = bind_rows(

      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X) %>%
        summarize(
          name = "Type-1",
          mean = (Y),
          .groups = "drop"
        ),
      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X) %>%
        summarize(name = "RT",
                  mean = (RT),
                  .groups = "drop"),

      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X, Correct) %>%
        summarize(name = "Confidence",
                  mean = (Confidence),
                  .groups = "drop")
    ) %>%
      filter(abs(X) < cutoff)

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
      filter(abs(X) < cutoff) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))

    # Calculate trial-level residuals properly
    # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
    pred_mean_per_trial = predictions %>%
      mutate(X = round(X, 2)) %>%
      group_by(X) %>%
      summarize(
        pred_prob_faster = mean(prob_faster),
        pred_rt = mean(rt_mu),
        .groups = "drop"
      )

    # Join predictions to actual trial-level data
    df_with_pred = df %>% select(Y,X,Confidence,RT,Correct) %>%
      mutate(X = round(X, 2)) %>%
      filter(abs(X) < cutoff) %>%
      left_join(pred_mean_per_trial, by = "X") %>%
      mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Confidence predictions need to be split by Correct
    pred_conf_per_trial = predictions %>%
      mutate(X = round(X, 2)) %>%
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
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop"),

      df_with_pred %>%
        mutate(residual = RT - pred_rt,
               name = "RT") %>%
        group_by(X) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop"),

      df_with_pred %>%
        mutate(residual = Confidence - pred_conf,
               name = "Confidence") %>%
        group_by(X, Correct_label) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop") %>%
        rename(Correct = Correct_label)
    )




    # Plot 1: Expected means (main plot)
    plot_mean = predictionsq_mean %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      # mutate(name = ifelse(name == "RT","Response time",ifelse(name == "Type-1","Binary choice","Confidence"))) %>%
      ggplot() +
      geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
      geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
      geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
      geom_point(data = dataq %>%
                   mutate(name = ifelse(name == "RT","Response time",
                                        ifelse(name == "Type-1","Binary choice","Confidence")),
                          name = factor(name, levels = c("Binary choice",
                                                         "Response time",
                                                         "Confidence")))
                 , aes(x = X, y = mean, fill = Correct),
                 shape = 21, color = "black", alpha = 0.5, size = 3) +
      geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      theme_classic(base_size = 20) +
      labs(color = "Correct", fill = "Correct",
           # title = "Group predictions (expected means)",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(
        legend.position = "none",
        legend.text = element_text(size = 20),      # text of legend items
            legend.title = element_text(size = 20),      # title of legend
            axis.title.x = element_blank(),
            axis.text.x = element_blank())


    # Plot 2: Residuals
    plot_residuals = residuals_data  %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
      geom_point(alpha = 0.5, size = 3) +
      # geom_smooth(method = "loess", se = F, alpha = 0.2) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      theme_classic(base_size = 20) +
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      labs(x = "Coherence", y = "(Obs - Pred)") +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
      theme(legend.position = "none")

    library(patchwork)
    # Combine plots using patchwork
    combined_plot = plot_mean / plot_residuals +
      plot_layout(heights = c(2, 1))

    combined_plot


    # ggsave(here::here("Siebe","Siebe_results_bin.tiff"),combined_plot ,dpi = 300,width = 30, height = 22, units = "cm")



    #################### Probability of responding crorectly:

    cutoff = 2

    dataq = bind_rows(

      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X) %>%
        summarize(
          name = "Type-1",
          mean = (ACC),
          .groups = "drop"
        ),
      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X) %>%
        summarize(name = "RT",
                  mean = (RT),
                  .groups = "drop"),

      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X, Correct) %>%
        summarize(name = "Confidence",
                  mean = (Confidence),
                  .groups = "drop")
    ) %>%
      filter(abs(X) < cutoff)

    # Prepare predicted data (using expected means)
    predictionsq_mean = bind_rows(
      predictions %>%
        group_by(X) %>%
        summarize(name = "Type-1",
                  mean = mean(prob    ),
                  q5 = quantile(prob    , 0.05),
                  q10 = quantile(prob , 0.1),
                  q20 = quantile(prob , 0.2),
                  q95 = quantile(prob , 0.95),
                  q90 = quantile(prob , 0.90),
                  q80 = quantile(prob , 0.80),
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
      filter(abs(X) < cutoff) %>%
      mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))

    # Calculate trial-level residuals properly
    # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
    pred_mean_per_trial = predictions %>%
      mutate(X = round(X, 2)) %>%
      group_by(X) %>%
      summarize(
        pred_prob_faster = mean(prob),
        pred_rt = mean(rt_mu),
        .groups = "drop"
      )

    # Join predictions to actual trial-level data
    df_with_pred = df %>% select(Y,X,Confidence,RT,Correct) %>%
      mutate(X = round(X, 2)) %>%
      filter(abs(X) < cutoff) %>%
      left_join(pred_mean_per_trial, by = "X") %>%
      mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Confidence predictions need to be split by Correct
    pred_conf_per_trial = predictions %>%
      mutate(X = round(X, 2)) %>%
      group_by(X, Correct) %>%
      summarize(pred_conf = mean(conf_mu_actual), .groups = "drop") %>%
      mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

    df_with_pred = df_with_pred %>%
      left_join(pred_conf_per_trial %>% select(X, Correct_label, pred_conf),
                by = c("X", "Correct_label"))

    # Calculate residuals at trial level, then aggregate
    residuals_data = bind_rows(
      df_with_pred %>%
        mutate(residual = Correct - pred_prob_faster,
               name = "Type-1") %>%
        group_by(X) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop"),

      df_with_pred %>%
        mutate(residual = RT - pred_rt,
               name = "RT") %>%
        group_by(X) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop"),

      df_with_pred %>%
        mutate(residual = Confidence - pred_conf,
               name = "Confidence") %>%
        group_by(X, Correct_label) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop") %>%
        rename(Correct = Correct_label)
    )




    # Plot 1: Expected means (main plot)
    plot_mean = predictionsq_mean %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      # mutate(name = ifelse(name == "RT","Response time",ifelse(name == "Type-1","Binary choice","Confidence"))) %>%
      ggplot() +
      geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
      geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
      geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
      geom_point(data = dataq %>%
                   mutate(name = ifelse(name == "RT","Response time",
                                        ifelse(name == "Type-1","Binary choice","Confidence")),
                          name = factor(name, levels = c("Binary choice",
                                                         "Response time",
                                                         "Confidence"))) %>%
                   mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct)))
                 , aes(x = X, y = mean, fill = Correct),
                 shape = 21, color = "black", alpha = 0.5, size = 3) +
      geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      theme_classic(base_size = 20) +
      labs(color = "Correct", fill = "Correct",
           # title = "Group predictions (expected means)",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(
        legend.position = "none",
        legend.text = element_text(size = 20),      # text of legend items
            legend.title = element_text(size = 20),      # title of legend
            axis.title.x = element_blank(),
            axis.text.x = element_blank())


    # Plot 2: Residuals
    plot_residuals = residuals_data  %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      mutate(Correct = ifelse(name == "Binary choice" & residual_mean  < 0, "Incorrect",ifelse(name == "Binary choice" & residual_mean  > 0, "Correct",Correct))) %>%
      ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
      geom_point(alpha = 0.5, size = 3) +
      # geom_smooth(method = "loess", se = F, alpha = 0.2) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      theme_classic(base_size = 20) +
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      labs(x = "Coherence", y = "(Obs - Pred)") +
      geom_vline(xintercept = 0, linetype = 2) +

      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
      theme(legend.position = "none")

    # Combine plots using patchwork
    combined_plot_acc = plot_mean / plot_residuals +
      plot_layout(heights = c(2, 1))

    combined_plot_acc
    # ggsave(here::here("Siebe","Siebe_results_bin.tiff"),combined_plot ,dpi = 300,width = 30, height = 22, units = "cm")


    #######################################################



    return(list(
      combined_plot_acc = combined_plot_acc,
      combined_plot = combined_plot
    ))
  }


  pps = Plot_group_predictive_psycho(predictions,data1,n_bins = NULL)

  return(pps)

  # ggsave(here::here("Siebe",paste0("Siebe_results_",number,"_.tiff")),pps$plot_combined ,dpi = 300,width = 18, height = 14, units = "cm")





}











make_pred_plot_bin = function(fit,data1,number){

  n_draws = 2000
  workers = 7
  memory = 15000 * 1024^2

  # Group-level parameters (from gm)
  parameters = c("alpha", "beta", "lapse")

  # Extract group means (gm)
  df_param = as_draws_df(fit$draws("gm")) %>%
    select(-contains(".")) %>%
    rename_with(~parameters) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw, names_to = "variable")



  library(future.apply)
  # Set up parallel processing
  plan(multisession, workers = 2)
  options(future.globals.maxSize = 5000 * 1024^2)

  draws <- 1:n_draws

  # Only use the number of draws that the user wants
  dfq = df_param %>% filter(draw %in% draws)

  source("~/Multivariate-Confidence-modeling/Analysis/Functions/Correct_models/utility.R")

  # Function to get the draws
  pred_list <- future_lapply(draws, function(d) {

    # Extract parameter vectors for this draw
    params <- dfq %>%
      filter(draw == d) %>%
      select(variable, value) %>%
      pivot_wider(names_from = "variable", values_from = "value")

    # Generate X values
    x = seq(-1, 1, by = 0.01)
    n_trials = length(x)

    # Get probability correct for each trial (using same approach as subject-level)
    prob_faster = psycho(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse) / 2)
    prob_cor = get_prob_cor(prob_faster, x)



    # Transform uniforms to marginal distributions
    ACC_pred = rbinom(length(prob_cor), 1, prob_cor)


    predictions = data.frame(
      X = x,
      Correct = ACC_pred,  # 1 if correct, 0 if incorrect
      prob = prob_cor,
      prob_faster = prob_faster,
      draw = d
    )

    return(predictions)
  }, future.seed = TRUE)

  predictions = map_dfr(pred_list, bind_rows)


  Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {

    cutoff = 2

    dataq = bind_rows(

      df %>%
        mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
        group_by(X) %>%
        summarize(
          name = "Type-1",
          mean = (Y),
          .groups = "drop"
        ) %>%
      filter(abs(X) < cutoff)
      )

    # Prepare predicted data (using expected means)
    predictionsq_mean =
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
                  .groups = "drop") %>% filter(abs(X) < cutoff)

    # Calculate trial-level residuals properly
    # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
    pred_mean_per_trial = predictions %>%
      mutate(X = round(X, 2)) %>%
      group_by(X) %>%
      summarize(
        pred_prob_faster = mean(prob_faster),
        .groups = "drop"
      )

    # Join predictions to actual trial-level data
    df_with_pred = df %>% select(Y,X,Confidence,RT,Correct) %>%
      mutate(X = round(X, 2)) %>%
      filter(abs(X) < cutoff) %>%
      left_join(pred_mean_per_trial, by = "X") %>%
      mutate(Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"))

    # Calculate residuals at trial level, then aggregate
    residuals_data = bind_rows(
      df_with_pred %>%
        mutate(residual = Y - pred_prob_faster,
               name = "Type-1") %>%
        group_by(X) %>%
        summarize(residual_mean = (residual),
                  name = first(name),
                  .groups = "drop")
    )




    # Plot 1: Expected means (main plot)
    plot_mean = predictionsq_mean %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      # mutate(name = ifelse(name == "RT","Response time",ifelse(name == "Type-1","Binary choice","Confidence"))) %>%
      ggplot() +
      geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95), alpha = 0.1) +
      geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90), alpha = 0.3) +
      geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80), alpha = 0.5) +
      geom_point(data = dataq %>%
                   mutate(name = ifelse(name == "RT","Response time",
                                        ifelse(name == "Type-1","Binary choice","Confidence")),
                          name = factor(name, levels = c("Binary choice",
                                                         "Response time",
                                                         "Confidence")))
                 , aes(x = X, y = mean),
                 shape = 21, color = "black", alpha = 0.5, size = 3) +
      geom_line(aes(x = X, y = mean), linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      theme_classic(base_size = 20) +
      labs(color = "Correct", fill = "Correct",
           # title = "Group predictions (expected means)",
           y = "Value") +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(legend.position = "top",
            legend.text = element_text(size = 20),      # text of legend items
            legend.title = element_text(size = 20),      # title of legend
            axis.title.x = element_blank(),
            axis.text.x = element_blank())


    # Plot 2: Residuals
    plot_residuals = residuals_data  %>%
      mutate(name = ifelse(name == "RT","Response time",
                           ifelse(name == "Type-1","Binary choice","Confidence")),
             name = factor(name, levels = c("Binary choice",
                                            "Response time",
                                            "Confidence"))) %>%
      ggplot(aes(x = X, y = residual_mean)) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
      geom_point(alpha = 0.5, size = 3) +
      # geom_smooth(method = "loess", se = F, alpha = 0.2) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      theme_classic(base_size = 20) +
      scale_color_manual(values = c("darkgreen","darkred","grey"))+
      scale_fill_manual(values = c("darkgreen","darkred","grey"))+
      labs(x = "Coherence", y = "(Obs - Pred)") +
      geom_vline(xintercept = 0, linetype = 2) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
      theme(legend.position = "none")

    library(patchwork)
    # Combine plots using patchwork
    combined_plot = plot_mean / plot_residuals +
      plot_layout(heights = c(2, 1))

    combined_plot

    return(list(
      plot_combined = combined_plot
    ))
  }


  pps = Plot_group_predictive_psycho(predictions,data1,n_bins = NULL)

  return(pps)

}
