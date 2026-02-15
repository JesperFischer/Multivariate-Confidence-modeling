library(tidyverse)

# Weibull psychometric function for accuracy
psycho_weibull <- function(x, alpha, beta, lapse, gamma = 0.5) {
  gamma + (1 - gamma - lapse) * (1 - exp(-(x/alpha)^beta))
}

# Logistic psychometric function for response (full psychometric)
psycho_logistic <- function(x, alpha, beta, lapse) {
  # alpha is the threshold (where p = 0.5), set to 0 for no vertical shift
  # beta is the slope
  # lapse is the lapse rate
  lapse + (1 - 2*lapse) / (1 + exp(-beta * (x - alpha)))
}

# Grid-based parameter variation for Weibull (accuracy)
simulate_parameter_grid <- function(alpha_values = c(0.05, 0.1, 0.2),
                                   beta_values = c(1, 2, 3),
                                   lapse = 0.05,
                                   gamma = 0.5,
                                   coherence_values = seq(0, 1, length.out = 50),
                                   trials_per_coherence = 50) {

  # Create parameter grid
  param_grid <- expand.grid(
    alpha = alpha_values,
    beta = beta_values
  )

  all_data <- list()

  for(i in 1:nrow(param_grid)) {
    alpha_i <- param_grid$alpha[i]
    beta_i <- param_grid$beta[i]

    curve_data <- data.frame()

    for(coh in coherence_values) {
      p_correct <- psycho_weibull(coh, alpha_i, beta_i, lapse, gamma)

      responses <- rbinom(trials_per_coherence, 1, p_correct)

      temp_df <- data.frame(
        curve_id = i,
        coherence = coh,
        trial = 1:trials_per_coherence,
        ACC = responses,
        p_correct = p_correct,
        alpha_true = alpha_i,
        beta_true = beta_i,
        lapse_true = lapse
      )

      curve_data <- rbind(curve_data, temp_df)
    }

    all_data[[i]] <- curve_data
  }

  simulated_data <- bind_rows(all_data)

  return(simulated_data)
}

# Grid-based parameter variation for logistic (response psychometric)
simulate_logistic_grid <- function(alpha_values = c(0),  # threshold (0 = no shift)
                                   beta_values = c(5, 10, 20),  # slope
                                   lapse = 0.05,
                                   stim_values = seq(-0.5, 0.5, length.out = 20),
                                   trials_per_stim = 50) {

  # Create parameter grid
  param_grid <- expand.grid(
    alpha = alpha_values,
    beta = beta_values
  )

  all_data <- list()

  for(i in 1:nrow(param_grid)) {
    alpha_i <- param_grid$alpha[i]
    beta_i <- param_grid$beta[i]

    curve_data <- data.frame()

    for(stim in stim_values) {
      p_resp1 <- psycho_logistic(stim, alpha_i, beta_i, lapse)

      responses <- rbinom(trials_per_stim, 1, p_resp1)

      temp_df <- data.frame(
        curve_id = i,
        stim = stim,
        trial = 1:trials_per_stim,
        resp = responses,
        p_resp1 = p_resp1,
        alpha_true = alpha_i,
        beta_true = beta_i,
        lapse_true = lapse
      )

      curve_data <- rbind(curve_data, temp_df)
    }

    all_data[[i]] <- curve_data
  }

  simulated_data <- bind_rows(all_data)

  return(simulated_data)
}

# Example with parameter grid
sim_data_grid <- simulate_parameter_grid(
  alpha_values = c(0.05, 0.1, 0.2,0.3,0.4,0.5),
  beta_values = c(5,6,7,8,9,10),
  lapse = 0.05
)

# Plot curves in a grid (alpha in rows, beta in columns)
sim_data_grid %>%
  group_by(curve_id, coherence, alpha_true, beta_true) %>%
  summarise(
    mean_acc = mean(ACC),
    true_p = first(p_correct),
    lapse_true = first(lapse_true),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = coherence, y = mean_acc)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(aes(y = true_p), color = "red", size = 1) +
  labs(
    title = "Simulated Weibull Curves - Parameter Grid",
    x = "Coherence",
    y = "P(Correct)",
  ) +
  facet_grid(alpha_true ~ beta_true,
             labeller = label_both) +
  theme_classic(base_size = 16)
