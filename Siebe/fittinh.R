library(tidyverse)
library(posterior)
library(bayesplot)

font = "sans"
font_size = 16
font_size_small = 12
axis_width = 1.5
tick_width = 1.5

text = ggplot2::theme(text = ggplot2::element_text(family = font, size = font_size))
theme = theme_classic()

#patch theme
patchtheme = ggplot2::theme(plot.tag = ggplot2::element_text(size = (font_size+4),
                                                             family = "sans",
                                                             face = "bold",
                                                             hjust = 0.5,
                                                             vjust = 0.5),
                            text=ggplot2::element_text(family=font))






data = read_csv("Siebe/Data_siebe.csv") %>%
  mutate(resp = ifelse(resp == "['left']",1,0),
                stim = ifelse(`dots direction` == 180, -coherence,coherence)) %>%
  mutate(ACC = ifelse(resp == 1 & stim > 0, 1,
               ifelse(resp == 0 & stim < 0, 1,0))) %>%
  filter(Trialtype == "Main") %>% filter(SR_conf != "None") %>% mutate(SR_conf = as.numeric(SR_conf))

data %>% ggplot(aes(x = stim, y = resp))+
  geom_point()


data %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = rts), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = stim, y = resp), shape = 21, fill = "darkgrey", col = "black", size = 2, alpha = 0.5)+
  ylab("Response = `Right`")+
  xlab("Coherence")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("binary_dots.tiff"),dpi = 300,width = 24, height = 10, units = "cm")



mod = cmdstanr::cmdstan_model(here::here("Siebe","psyscho.stan"))


datastan = list(N = nrow(data),
                binom_y = data$resp,
                X = data$stim)


fit_psycho = mod$sample(data = datastan,
                        parallel_chains = 4,
                        iter_warmup = 3000,
                        adapt_delta = 0.99)

psyco = function(x,alpha,beta,lapse){return(lapse + (1-2*lapse) * brms::inv_logit_scaled(beta * (x-alpha)))}

expected = as_draws_df(fit_psycho$draws(c("alpha","beta","lapse"))) %>% select(-contains(".")) %>%
  mutate(beta = exp(beta)) %>%
  mutate(x = list(seq(min(data$stim),max(data$stim),by = 0.01))) %>%
  unnest(x) %>%
  mutate(p = psyco(x,alpha,beta,lapse)) %>%
  group_by(x) %>% summarize(E_mean = mean(p),
                            E_q5 = quantile(p, 0.10),
                            E_q95 = quantile(p, 0.90)) %>%
  rename(stim = x)


data %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = rts), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = stim, y = resp), shape = 21, fill = "darkgrey", col = "black", size = 2, alpha = 0.5)+
  geom_ribbon(data = expected, aes(x = stim, y = E_mean, ymin = E_q5, ymax = E_q95), fill = "grey", col = "darkgrey", alpha = 0.5)+
  geom_line(data = expected, aes(x = stim, y = E_mean), col = "black")+
  ylab("Response = `Right`")+
  xlab("Coherence")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))


ggsave(here::here("binary_dots_psychometric.tiff"),dpi = 300,width = 24, height = 10, units = "cm")


inner_join(data %>% mutate(stim = round(stim,2)),expected%>% mutate(stim = round(stim,2))) %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = rts), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = E_mean, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 2)+
  ylab("Response time (s)")+
  xlab("P(Response = `Right`)")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  scale_x_reverse() +  # <- this flips the axis
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("RT_dots_expect.tiff"),dpi = 300,width = 10, height = 10, units = "cm")




modrt = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_binRT.stan"))


datastan = list(N = nrow(data),
                binom_y = data$cor,
                RT = data$RTdec,
                Conf = (data$SR_conf+1)/2,
                X = data$stim,
                minRT = min(data$RTdec),
                ACC = data$cor)


fit_rt = modrt$sample(data = datastan,
                      parallel_chains = 4,
                      adapt_delta = 0.90)



H = function(E){
  return(-E * log(E) - (1-E) * log(1-E));
}

murts = as_draws_df(fit_rt$draws(c("rt_int","rt_slope","rt_prec","rt_ndt"))) %>%
         mutate(E = list(seq(0.01,0.99,by = 0.01))) %>% unnest() %>%
  mutate(mu_pred_RT = exp(rt_int + rt_slope * H(E) + (rt_prec^2 / 2)) + rt_ndt) %>%
  group_by(E) %>%
  summarize(RT_mean = mean(mu_pred_RT),
            RT_q5 = quantile(mu_pred_RT,0.05),
            RT_q95 = quantile(mu_pred_RT,0.95))



expected = data.frame(fit_rt$summary("theta")) %>% mutate(trials = 1:n()) %>% filter(trials %in% 1:nrow(data))

data %>% mutate(E_mean = expected$mean,
              E_q5 = expected$q5,
              E_q95 = expected$q95) %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = E_mean, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 2)+
  geom_ribbon(data = murts,aes(x = E, y = RT_mean, ymin = RT_q5, ymax = RT_q95), fill = "darkgrey", col = "black", alpha = 0.5)+
  geom_line(data = murts,aes(x = E, y = RT_mean), col = "black")+
  ylab("Response time (s)")+
  xlab("P(Response = `Right`)")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  scale_x_reverse() +  # <- this flips the axis
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("dots_with_fit_RT.tiff"),dpi = 300,width = 10, height = 10, units = "cm")




data %>% mutate(E_mean = expected$mean,
                E_q5 = expected$q5,
                E_q95 = expected$q95) %>%
  ggplot()+
  geom_pointrange(aes(x = E_mean, xmin = E_q5, xmax = E_q95, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  # geom_point(aes(x = E_mean, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 2)+
  geom_ribbon(data = murts,aes(x = E, y = RT_mean, ymin = RT_q5, ymax = RT_q95), fill = "darkgrey", col = "black", alpha = 0.5)+
  geom_line(data = murts,aes(x = E, y = RT_mean), col = "black")+
  ylab("Response time (s)")+
  xlab("P(Response = `Right`)")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  scale_x_reverse() +  # <- this flips the axis
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("dots_with_fit_RT_w_errorbars.tiff"),dpi = 300,width = 10, height = 10, units = "cm")



data %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = rts), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = stim, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 2)+
  ylab("Response time (s)")+
  xlab("Coherence")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("RT_dots_coherence.tiff"),dpi = 300,width = 24, height = 10, units = "cm")




modrt = cmdstanr::cmdstan_model(here::here("Siebe","quad.stan"))


datastan = list(N = nrow(data),
                RT = data$RTdec,
                X = data$stim,
                minRT = min(data$RTdec))


fit_rt = modrt$sample(data = datastan,
                      parallel_chains = 4,
                      adapt_delta = 0.90)



murts = as_draws_df(fit_rt$draws(c("a","b","c","rt_prec","rt_ndt"))) %>%
  mutate(C = list(seq(-0.97,0.97,by = 0.01))) %>% unnest() %>%
  mutate(mu_pred_RT = exp(a*C^2 + b*C + c + (rt_prec^2 / 2)) + rt_ndt) %>%
  group_by(C) %>%
  summarize(RT_mean = mean(mu_pred_RT),
            RT_q5 = quantile(mu_pred_RT,0.05),
            RT_q95 = quantile(mu_pred_RT,0.95))


data %>%
  ggplot()+
  # geom_pointrange(aes(x = E_mean*100, xmin = E_q5*100, xmax = E_q95*100, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 0.35)+
  geom_point(aes(x = stim, y = RTdec), shape = 21, fill = "darkgrey", col = "black", size = 2)+
  geom_ribbon(data = murts,aes(x = C, y = RT_mean, ymin = RT_q5, ymax = RT_q95), fill = "darkgrey", col = "black", alpha = 0.5)+
  geom_line(data = murts,aes(x = C, y = RT_mean), col = "black")+
  ylab("Response time (s)")+
  xlab("Coherence")+
  # geom_smooth(aes(x = E_mean, y = predRT))+
  theme+text+
  # scale_x_reverse() +  # <- this flips the axis
  theme(legend.position = "none",
        legend.key.width = unit(0.85, 'cm'), #size of bars
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width))

ggsave(here::here("dots_with_fit_RT.tiff"),dpi = 300,width = 24, height = 10, units = "cm")







data1 = data %>% rename(X = stim, RT = RTdec, Correct = cor, Y = resp) %>%
  mutate(Confidence = (data$SR_conf+1)/2)





mod2 = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_bin_ss_correct.stan"))


datastan = list(N = nrow(data1),
                binom_y = data1$Correct,
                RT = data1$RT,
                Conf = data1$Confidence,
                X = data1$X,
                minRT = min(data1$RT),
                ACC = data1$Correct)


fit = mod2$sample(data = datastan,
                 parallel_chains = 4,
                 adapt_delta = 0.90)




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


# Flatten nested list and create a tidy long dataframe
predictions = map_dfr(pred_list, bind_rows)
library(patchwork)
Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {

  cutoff = 2

  # Prepare observed data
  if (!is.null(n_bins)) {
    # Create common bin boundaries based on the range of both datasets
    all_X <- c(df$X, predictions$X)
    X_range <- range(all_X, na.rm = TRUE)
    bin_breaks <- seq(X_range[1], X_range[2], length.out = n_bins + 1)

    # Calculate bin centers (midpoints)
    bin_centers <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

    # Bin the data using the common breaks and assign bin centers
    df <- df %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)

    # Bin the predictions using the same common breaks and assign same bin centers
    predictions <- predictions %>%
      mutate(X_bin = cut(X, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE),
             X = bin_centers[X_bin]) %>%
      select(-X_bin)
  }

  dataq = bind_rows(
    # df %>%
    #   mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
    #   group_by(X) %>%
    #   summarize(
    #     name = "Type-1",
    #     n = n(),
    #     k = sum(Y),
    #     mean = (1 + k) / (2 + n),
    #     q5  = qbeta(0.05, 1 + k, 1 + n - k),
    #     q10 = qbeta(0.10, 1 + k, 1 + n - k),
    #     q20 = qbeta(0.20, 1 + k, 1 + n - k),
    #     q80 = qbeta(0.80, 1 + k, 1 + n - k),
    #     q90 = qbeta(0.90, 1 + k, 1 + n - k),
    #     q95 = qbeta(0.95, 1 + k, 1 + n - k),
    #     .groups = "drop"
    #   ),

    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X) %>%
      summarize(
        name = "Type-1",
        mean = mean(ACC),
        q5 = mean(ACC) - 2* (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        q95 = mean(ACC) + 2 * (mean(ACC) * (1-mean(ACC)) / sqrt(n())),
        # n = n(),
        # k = sum(Y),
        # mean = (1 + k) / (2 + n),
        # q5  = qbeta(0.05, 1 + k, 1 + n - k),
        # q10 = qbeta(0.10, 1 + k, 1 + n - k),
        # q20 = qbeta(0.20, 1 + k, 1 + n - k),
        # q80 = qbeta(0.80, 1 + k, 1 + n - k),
        # q90 = qbeta(0.90, 1 + k, 1 + n - k),
        # q95 = qbeta(0.95, 1 + k, 1 + n - k),
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
    filter(abs(X) < cutoff)

  # Prepare predicted data (using expected means)
  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob ),
                q5 = quantile(prob , 0.05),
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
  df_with_pred = df %>%
    mutate(X = round(X, 2)) %>%
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
      mutate(residual = ACC - pred_prob_faster,
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
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct",
         # title = "Group predictions (expected means)",
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
    theme_classic(base_size = 14) +
    labs(x = "Stimulus strength (X)", y = "(Obs - Pred)") +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")

  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))

  combined_plot

  ggsave(here::here("Siebe","Siebe_results_ACC.tiff"),combined_plot ,dpi = 300,width = 18, height = 14, units = "cm")

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
    filter(abs(X) < cutoff) %>%
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
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct") +
    # ggtitle("Group predictions (posterior predictive samples)")+
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")

  return(list(
    plot_combined = combined_plot,
    plot_mean = plot_mean,
    plot_residuals = plot_residuals,
    plot_preds = plot_preds
  ))
}


pps = Plot_group_predictive_psycho(predictions,data1,n_bins = 10)

pps$plot_combined

ggsave(here::here("Siebe","Siebe_results_bin.tiff"),pps$plot_combined ,dpi = 300,width = 18, height = 14, units = "cm")


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
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")

  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))

  combined_plot


  ggsave(here::here("Siebe","Siebe_results_bin.tiff"),combined_plot ,dpi = 300,width = 30, height = 22, units = "cm")



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
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none")

  # Combine plots using patchwork
  combined_plot = plot_mean / plot_residuals +
    plot_layout(heights = c(2, 1))

  combined_plot
  ggsave(here::here("Siebe","Siebe_results_bin.tiff"),combined_plot ,dpi = 300,width = 30, height = 22, units = "cm")


  #######################################################



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
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))


  # Plot 2: Actual predictions
  plot_preds = predictionsq_preds %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_point(data = dataq, aes(x = X, y = mean, fill = as.factor(Correct)),
                    shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free", ncol = 3) +
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct") +
    # ggtitle("Group predictions (posterior predictive samples)")+
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")

  return(list(
    plot_combined = combined_plot,
    plot_mean = plot_mean,
    plot_residuals = plot_residuals,
    plot_preds = plot_preds
  ))
}


pps = Plot_group_predictive_psycho(predictions,data1,n_bins = NULL)

pps$plot_combined

ggsave(here::here("Siebe","Siebe_results_real.tiff"),pps$plot_combined ,dpi = 300,width = 18, height = 14, units = "cm")





# mod2 = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","crazy.stan"))
#
#
# datastan = list(N = nrow(data1),
#                 binom_y = data1$Correct,
#                 RT = data1$RT,
#                 Conf = data1$Confidence,
#                 X = data1$X,
#                 minRT = min(data1$RT),
#                 ACC = data1$Correct)
#
#
# fit_crazy = mod2$sample(data = datastan,
#                   parallel_chains = 4,
#                   adapt_delta = 0.90)




modrt = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_binRT.stan"))


datastan = list(N = nrow(data),
                binom_y = data$cor,
                RT = data$RTdec,
                Conf = (data$SR_conf+1)/2,
                X = data$stim,
                minRT = min(data$RTdec),
                ACC = data$cor)


fit_rt = modrt$sample(data = datastan,
                  parallel_chains = 4,
                  adapt_delta = 0.90)



loos_bin = loo::loo_compare(list(psychometric = fit_psycho$loo(),
                                 RT = fit_rt$loo("log_lik_bin"),
                                 Full = fit$loo("log_lik_bin")))

loos_bin



#font size and style for all plots for everything
font = "sans"
font_size = 16
font_size_small = 12
axis_width = 1.5

text = ggplot2::theme(text = ggplot2::element_text(family = font, size = font_size))
theme = cowplot::theme_cowplot()



rbind(
  as_draws_df(fit_psycho$draws("gm")) %>% select(-contains(".")) %>% mutate(model = "Bin"),
  as_draws_df(fit_rt$draws(paste0("gm[",1:3,"]")))%>% select(-contains(".")) %>% mutate(model = "RT"),
  as_draws_df(fit$draws(paste0("gm[",1:3,"]")))%>% select(-contains(".")) %>% mutate(model = "Conf")
) %>% pivot_longer(-model) %>%
  mutate(name = ifelse(name == "gm[1]", "Threshold",ifelse(name == "gm[2]", "Slope","Lapse"))) %>%
  ggplot(aes(x = value, fill = model))+geom_histogram(col = "black", position = "identity", alpha = 0.5)+
  scale_fill_manual(values = c("red","black","orange"))+
  xlab("Parameter value")+
  ylab("Count")+
  facet_wrap(~name,scales = "free")+
  text+theme

# Example: extract your draws and reshape
df <- as_draws_df(fit$draws(c("rho_p_rt","rho_rt_conf","rho_p_conf"))) %>%
  select(-contains(".")) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value")

# Plot
ggplot(df, aes(x = value, fill = parameter)) +
  geom_histogram(alpha = 0.7, color = "black", position = "identity", bins = 30) +
  scale_fill_manual(
    values = c("rho_p_rt" = "black",
               "rho_rt_conf" = "darksalmon",
               "rho_p_conf" = "darkblue")
  ) +
  theme_classic(base_size = 18) +      # makes all text bigger for slideshow
  theme(
    # legend.position = "none",
    legend.title = element_blank(),    # remove legend title
    legend.text = element_text(size = 18),  # bigger legend text
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  labs(x = "Correlation", y = "Count")+
  scale_x_continuous(limits = c(-1,0.55),breaks = seq(-1,0.6,by = 0.4),labels = seq(-1,0.6,by = 0.4))

ggsave(here::here("copula_correlations.tiff"),dpi = 300,width = 18, height = 10, units = "cm")


ddf = data.frame()

for(n in c(seq(10,nrow(data),by = 10), nrow(data))){

  print(n)

  mod = cmdstanr::cmdstan_model(here::here("Siebe","psyscho.stan"))


  datastan = list(N = nrow(data1[1:n,]),
                  binom_y = data1$Correct[1:n],
                  X = data1$X[1:n])


  fit_psycho_x = mod$sample(data = datastan,
                          parallel_chains = 4,
                          init = 0,
                          seed = 1997,
                          adapt_delta = 0.999)

 divs = fit_psycho_x$diagnostic_summary()
  dq = fit_psycho_x$summary(c("alpha","beta","lapse")) %>% mutate(trial = n, div = max(divs$num_divergent))

  ddf = rbind(ddf,dq)
}


ddf %>%
  # filter(div == 0) %>%
  select(variable, mean, q5, q95, trial) %>%
  ggplot(aes(x = trial, y = mean, ymin = q5, ymax = q95, color = variable)) +
  geom_pointrange() +
  facet_wrap(~ variable, scales = "free_y")



ddfrt = data.frame()
#
for(n in c(seq(10,nrow(data1),by = 10), nrow(data1))){

  print(n)

  modrt = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_binRT.stan"))


  datastan = list(N = nrow(data1[1:n,]),
                  binom_y = data1$Correct[1:n],
                  RT = data1$RT[1:n],
                  Conf = data1$Confidence[1:n],
                  X = data1$X[1:n],
                  minRT = min(data1$RT[1:n]),
                  ACC = data1$Correct[1:n])


  fit = modrt$sample(data = datastan,
                    parallel_chains = 4,
                    init = 0,
                    adapt_delta = 0.99)

  divs = fit$diagnostic_summary()

  dq = fit$summary(c("alpha","beta","lapse",
                     "rt_int","rt_slope","rt_prec","rt_ndt",
                     "rho_p_rt")) %>% mutate(trial = n, div = max(divs$num_divergent))

  ddfrt = rbind(ddfrt,dq)
}


ddfrt %>%
  filter(div == 0) %>%
  select(variable, mean, q5, q95, trial) %>%
  ggplot(aes(x = trial, y = mean, ymin = q5, ymax = q95, color = variable)) +
  geom_pointrange() +
  facet_wrap(~ variable, scales = "free")+
  ylab("Parameter value [90% HDI]") +
  xlab("N trials") +
  theme + text+
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank())


ggsave(here::here("Siebe","rt_all_params.tiff") ,dpi = 300,width = 24, height = 16, units = "cm")




ddfq = data.frame()

for(n in c(seq(10,nrow(data1),by = 10), nrow(data1))){

  print(n)

  mod2 = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_bin_ss_correct.stan"))

  datastan = list(N = nrow(data1[1:n,]),
                  binom_y = data1$Correct[1:n],
                  RT = data1$RT[1:n],
                  Conf = data1$Confidence[1:n],
                  X = data1$X[1:n],
                  minRT = min(data1$RT[1:n]),
                  ACC = data1$Correct[1:n])


  fit = mod2$sample(data = datastan,
                    parallel_chains = 4,
                    init = 0,
                    seed = 1997,
                    max_treedepth = 12,
                    adapt_delta = 0.99)


  divs = fit$diagnostic_summary()
  dq = fit$summary(c("alpha","beta","lapse",
                     "rt_int","rt_slope","rt_prec","rt_ndt",
                     "conf_prec","meta_un","meta_bias",
                     "rho_p_rt","rho_p_conf","rho_rt_conf")) %>% mutate(trial = n, div = max(divs$num_divergent))

  ddfq = rbind(ddfq,dq)
}


ddfq %>%
  # filter(trial > 10) %>%
  select(variable, mean, q5, q95, trial) %>%
  ggplot(aes(x = trial, y = mean, ymin = q5, ymax = q95, color = variable)) +
  geom_pointrange() +
  facet_wrap(~ variable, scales = "free")+
  ylab("Parameter value [90% HDI]") +
  xlab("N trials") +
  theme + text+
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

ggsave(here::here("Siebe","all_params.tiff") ,dpi = 300,width = 24, height = 16, units = "cm")



variable = c("Threshold","Slope","Lapse")
mean = c(mean(rnorm(10000,0,0.5)),mean(rnorm(10000,1,2)), mean(rnorm(10000,-4,2)))
q5 = c(quantile(rnorm(10000,0,0.5), 0.05),quantile(rnorm(10000,1,2), 0.05), quantile(rnorm(10000,-4,2), 0.05))
q95 = c(quantile(rnorm(10000,0,0.5), 0.95),quantile(rnorm(10000,1,2), 0.95), quantile(rnorm(10000,-4,2), 0.95))

prior_df = data.frame(variable) %>% mutate(mean = mean,q5 = q5, q95 = q95) %>% mutate(model = "Prior", trial = 0)

bind_rows(
  prior_df,
  rbind(ddf %>% mutate(model = "Binary Choice"),
        ddfq %>% mutate(model = "Copula Based model")) %>%
  filter(variable %in% c("Threshold", "Slope", "Lapse") |
           variable %in% c("alpha","beta","lapse")) %>%
  mutate(variable = ifelse(variable == "alpha","Threshold",ifelse(variable == "beta","Slope","Lapse"))) %>%
  mutate(
    variable = recode(variable,
                      alpha = "Threshold",
                      beta  = "Slope",
                      lapse = "Lapse"),
    variable = factor(variable, levels = c("Threshold", "Slope", "Lapse"))
  )) %>%
  mutate(
    variable = recode(variable,
                      alpha = "Threshold",
                      beta  = "Slope",
                      lapse = "Lapse"),
    variable = factor(variable, levels = c("Threshold", "Slope", "Lapse"))
  ) %>%
  select(variable, mean, q5, q95, trial,model) %>%
  ggplot(aes(x = trial, y = mean, ymin = q5, ymax = q95, color = model)) +
  geom_pointrange(position = position_dodge(width = 3), size = 0.5, alpha = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_color_manual(values = c("red","black", "blue"))+
  facet_wrap(~ variable, scales = "free", ncol = 3)+
  ylab("Parameter value [90% HDI]") +
  xlab("N trials") +
  theme + text+
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

ggsave(here::here("Siebe","Siebe_results_params.tiff") ,dpi = 300,width = 24, height = 16, units = "cm")

rbind(ddf %>% mutate(model = "Binary Choice"),
      ddfq %>% mutate(model = "Copula Based model")) %>%
  filter(variable %in% c("Threshold", "Slope", "Lapse") |
           variable %in% c("alpha","beta","lapse")) %>%
  mutate(variable = ifelse(variable == "alpha","Threshold",ifelse(variable == "beta","Slope","Lapse"))) %>%
  mutate(
    variable = recode(variable,
                      alpha = "Threshold",
                      beta  = "Slope",
                      lapse = "Lapse"),
    variable = factor(variable, levels = c("Threshold", "Slope", "Lapse"))
  ) %>%
  mutate(
    variable = recode(variable,
                      alpha = "Threshold",
                      beta  = "Slope",
                      lapse = "Lapse"),
    variable = factor(variable, levels = c("Threshold", "Slope", "Lapse"))) %>%
  # filter(trial > 10) %>%
  select(variable, mean, q5, q95, trial,model) %>% mutate(interval = q95-q5) %>%
  ggplot(aes(x = trial, y = interval, color = model)) +
  geom_point(position = position_dodge(width = 3), size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_color_manual(values = c("red","black","orange"))+
  facet_wrap(~ variable, scales = "free", ncol = 3)+
  ylab("90% HDI") +
  xlab("N trials") +
  theme + text+
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank())


ggsave(here::here("Siebe","Siebe_results_params2_35.tiff") ,dpi = 300,width = 24, height = 16, units = "cm")









## ff


ddfq2 = data.frame()

for(n in c(seq(10,nrow(data1),by = 10), nrow(data1))){

  print(n)

  mod2 = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","real","ACC_bin_ss_correct.stan"))

  datastan = list(N = nrow(data1[1:n,]),
                  binom_y = data1$Correct[1:n],
                  RT = data1$RT[1:n],
                  Conf = data1$Confidence[1:n],
                  X = data1$X[1:n],
                  minRT = min(data1$RT[1:n]),
                  ACC = data1$Correct[1:n])


  fit = mod2$sample(data = datastan,
                    parallel_chains = 4,
                    init = 0,
                    seed = 1997,
                    max_treedepth = 12,
                    adapt_delta = 0.99)

  qq = make_pred_plot_full(fit,data1[1:n,],1)

  ggsave(here::here("Siebe","comb_plot_99.tiff"),qq$combined_plot ,dpi = 300,width = 24, height = 16, units = "cm")

  divs = fit$diagnostic_summary()
  dq = fit$summary(c("alpha","beta","lapse",
                     "rt_int","rt_slope","rt_prec","rt_ndt",
                     "conf_prec","meta_un","meta_bias",
                     "rho_p_rt","rho_p_conf","rho_rt_conf")) %>% mutate(trial = n, div = max(divs$num_divergent))

  ddfq2 = rbind(ddfq2,dq)
}






ddf = data.frame()

for(n in c(seq(10,nrow(data),by = 10), nrow(data))){

  print(n)

  mod = cmdstanr::cmdstan_model(here::here("Siebe","psyscho.stan"))


  datastan = list(N = nrow(data1[1:n,]),
                  binom_y = data1$Correct[1:n],
                  X = data1$X[1:n])


  fit_psycho_x = mod$sample(data = datastan,
                            parallel_chains = 4,
                            init = 0,
                            seed = 1997,
                            adapt_delta = 0.999)

  qq = make_pred_plot_bin(fit_psycho_x,data1[1:n,],1)

  ggsave(here::here("Siebe","pure_plot_99.tiff"),qq$plot_combined ,dpi = 300,width = 10, height = 16, units = "cm")



  divs = fit_psycho_x$diagnostic_summary()
  dq = fit_psycho_x$summary(c("alpha","beta","lapse")) %>% mutate(trial = n, div = max(divs$num_divergent))

  ddf = rbind(ddf,dq)
}








## animation


Plot_group_predictive_psycho = function(predictions, df, n_bins = NULL) {

  df = df %>% mutate(trial = 1:n())

  dataq = bind_rows(
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,trial) %>%
      summarize(
        name = "Type-1",
        mean = (ACC),
        .groups = "drop"
      ),
    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X,trial) %>%
      summarize(name = "RT",
                mean = (RT),
                .groups = "drop"),

    df %>%
      mutate(Correct = ifelse(Correct == 1, "Correct", "Incorrect")) %>%
      group_by(X, Correct,trial) %>%
      summarize(name = "Confidence",
                mean = (Confidence),
                .groups = "drop")
  )

  predictionsq_mean = bind_rows(
    predictions %>%
      group_by(X) %>%
      summarize(name = "Type-1",
                mean = mean(prob    ),
                q5 = quantile(prob    , 0.05),
                q10 = quantile(prob    , 0.1),
                q20 = quantile(prob    , 0.2),
                q95 = quantile(prob    , 0.95),
                q90 = quantile(prob    , 0.90),
                q80 = quantile(prob    , 0.80),
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
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))

  # Calculate trial-level residuals properly
  # For Type-1: observed Y vs predicted prob_faster (across draws, use mean prediction)
  pred_mean_per_trial = predictions %>%
    mutate(X = round(X, 2)) %>%
    group_by(X) %>%
    summarize(
      pred_prob_ACC = mean(prob),
      pred_rt = mean(rt_mu),
      .groups = "drop"
    )

  # Join predictions to actual trial-level data
  df_with_pred = df %>% select(ACC,X,Confidence,RT,Correct) %>%
    mutate(X = round(X, 2)) %>%
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
      mutate(trial = 1:n(),
             residual = ACC - pred_prob_ACC,
             name = "Type-1") %>%
      group_by(X, trial) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),

    df_with_pred %>%
      mutate(trial = 1:n(),
             residual = RT - pred_rt,
             name = "RT") %>%
      group_by(X,trial) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop"),

    df_with_pred %>%
      mutate(trial = 1:n(),
             residual = Confidence - pred_conf,
             name = "Confidence") %>%
      group_by(trial,X, Correct_label) %>%
      summarize(residual_mean = (residual),
                name = first(name),
                .groups = "drop") %>%
      rename(Correct = Correct_label)
  )



  predictionsq_mean = predictionsq_mean%>%
    mutate(name = factor(name, levels = c("Type-1", "RT", "Confidence")))


  library(
    gganimate)
  # Plot 1: Expected means (main plot)
  plot_mean = dataq %>%
    filter(name != "RT") %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    # filter(trial < 10) %>%
    mutate(
      name = factor(name, levels = c("Binary choice", "RT", "Confidence")),
      delay = case_when(
        name == "Binary choice"    ~ 0,
        # name == "RT"        ~ 0.1,
        name == "Confidence" ~ 0.6
      ),
      anim_time = trial + delay
    ) %>%
    mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct))) %>%
    ggplot() +
    geom_ribbon(data  = predictionsq_mean %>% filter(name != "RT") %>%
                  mutate(name = ifelse(name == "RT","Response time",
                                       ifelse(name == "Type-1","Binary choice","Confidence")),
                         name = factor(name, levels = c("Binary choice",
                                                        "Response time",
                                                        "Confidence"))) %>%
                  mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct))),
                aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +

    geom_ribbon(data  = predictionsq_mean%>% filter(name != "RT") %>%
                  mutate(name = ifelse(name == "RT","Response time",
                                       ifelse(name == "Type-1","Binary choice","Confidence")),
                         name = factor(name, levels = c("Binary choice",
                                                        "Response time",
                                                        "Confidence"))) %>%
                  mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct))),
                aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +

    geom_ribbon(data  = predictionsq_mean%>% filter(name != "RT") %>%
                  mutate(name = ifelse(name == "RT","Response time",
                                       ifelse(name == "Type-1","Binary choice","Confidence")),
                         name = factor(name, levels = c("Binary choice",
                                                        "Response time",
                                                        "Confidence"))) %>%
                  mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct))),
                aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +

    geom_line(data  = predictionsq_mean%>% filter(name != "RT") %>%
                mutate(name = ifelse(name == "RT","Response time",
                                     ifelse(name == "Type-1","Binary choice","Confidence")),
                       name = factor(name, levels = c("Binary choice",
                                                      "Response time",
                                                      "Confidence"))) %>%
                mutate(Correct = ifelse(name == "Binary choice" & mean < 0.5, "Incorrect",ifelse(name == "Binary choice" & mean > 0.5, "Correct",Correct))),
              aes(x = X, y = mean, color = Correct), linewidth = 1) +

    geom_point(aes(x = X, y = mean, fill = Correct),
               shape = 21, color = "black", alpha = 0.4, size = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2))+
    scale_color_manual(values = c("darkgreen","darkred"))+
    scale_fill_manual(values = c("darkgreen","darkred"))+
    theme_classic(base_size = 6) +
    labs(color = "Correct", fill = "Correct", x = "Coherence") +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.line  = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text  = element_text(size = 6),   # or larger if you want them visually thicker
          axis.title.x = element_blank())+
    transition_states(anim_time)+
    shadow_mark()+
    ease_aes("quintic-in")

  # plot_mean
  # --- Animate ---
  anim <- animate(plot_mean,
                  nframes = 1000,
                  # nframes = 10,
                  fps = 20,
                  height = 21*300/4 / 2.54,        # px
                  width = 14*300/4 / 2.54,        # px
                  res = 300,           # keep consistent with PowerPoint
                  renderer = gifski_renderer(file = here::here("ANI_data.gif"), loop = FALSE))


  # Plot 2: Residuals
  plot_residuals = residuals_data %>%
    filter(name == "Type-1") %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    # filter(trial < 10) %>%
    mutate(
      name = factor(name, levels = c("Binary choice", "RT", "Confidence")),
      delay = case_when(
        name == "Binary choice"    ~ 0,
        # name == "RT"        ~ 0.1,
        name == "Confidence" ~ 0.6
      ),
      anim_time = trial + delay
    ) %>%
    mutate(Correct = ifelse(name == "Binary choice" & residual_mean < 0, "Incorrect",ifelse(name == "Binary choice" & residual_mean > 0, "Correct",Correct))) %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_point(alpha = 0.5, size = 0.6) +
    # geom_smooth(method = "loess", se = F, alpha = 0.2) +
    # facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 6) +
    labs(x = "Coherence", y = "(Obs - Pred)") +

    scale_color_manual(values = c("darkgreen","darkred"))+
    scale_fill_manual(values = c("darkgreen","darkred"))+
    scale_y_continuous(breaks = c(-1.0,-0.5,0,0.5), labels = c("-1.0","-0.5","0.0","0.5"))+
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line  = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text  = element_text(size = 6),   # or larger if you want them visually thicker
          axis.title.x = element_blank())+
    transition_states(anim_time)+
    shadow_mark()+
    ease_aes("quintic-in")

  # plot_mean
  # --- Animate ---
  anim <- animate(plot_residuals,
                  nframes = 1000,
                  # nframes = 10,
                  fps = 20,
                  height = 11*300/4 / 2.54,        # px
                  width = 11*300/4 / 2.54,       # px
                  res = 300,           # keep consistent with PowerPoint
                  renderer = gifski_renderer(file = here::here("ANI_resid_resp.gif"), loop = FALSE))



  # Plot 2: Residuals
  plot_residuals = residuals_data %>%
    filter(name == "Confidence") %>%
    mutate(name = ifelse(name == "RT","Response time",
                         ifelse(name == "Type-1","Binary choice","Confidence")),
           name = factor(name, levels = c("Binary choice",
                                          "Response time",
                                          "Confidence"))) %>%
    # filter(trial < 10) %>%
    mutate(
      name = factor(name, levels = c("Binary choice", "RT", "Confidence")),
      delay = case_when(
        name == "Binary choice"    ~ 0,
        # name == "RT"        ~ 0.1,
        name == "Confidence" ~ 0.6
      ),
      anim_time = trial + delay
    ) %>%
    mutate(Correct = ifelse(name == "Binary choice" & residual_mean < 0, "Incorrect",ifelse(name == "Binary choice" & residual_mean > 0, "Correct",Correct))) %>%
    # mutate(Correct = ifelse(name == "Type-1" & residual_mean > 0, "Correct",ifelse(name == "Type-1" & residual_mean < 0, "Incorrect",NA))) %>%
    ggplot(aes(x = X, y = residual_mean, color = Correct, fill = Correct)) +
    geom_point(alpha = 0.5, size = 0.6) +
    # geom_smooth(method = "loess", se = F, alpha = 0.2) +
    # facet_wrap(~name, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 6) +

    scale_color_manual(values = c("darkgreen","darkred"))+
    scale_fill_manual(values = c("darkgreen","darkred"))+
    labs(x = "Confidence", y = "(Obs - Pred)") +
    scale_y_continuous(breaks = c(-0.2,-0.1,0.1,0.2), labels = c("-0.2","-0.1","0.1","0.2"))+
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line  = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text  = element_text(size = 6),   # or larger if you want them visually thicker
          axis.title.x = element_blank())+
    transition_states(anim_time)+
    shadow_mark()+
    ease_aes("quintic-in")

  # plot_mean
  # --- Animate ---
  anim <- animate(plot_residuals,
                  nframes = 1000,
                  # nframes = 10,
                  fps = 20,
                  height = 11*300/4 / 2.54,        # px
                  width = 11*300/4 / 2.54,       # px
                  res = 300,           # keep consistent with PowerPoint
                  renderer = gifski_renderer(file = here::here("ANI_resid_conf.gif"), loop = FALSE))





  residuals_wide <- residuals_data %>%
    select(trial, name, residual_mean) %>%
    pivot_wider(names_from = name, values_from = residual_mean) %>%
    rename(Type1 = `Type-1`)

  facet_cols <- c("Type1", "RT", "Confidence")

  pairs_data <- lapply(combn(facet_cols, 2, simplify = FALSE), function(pair) {
    residuals_wide %>%
      select(trial, all_of(pair)) %>%        # keep trial
      rename(val1 = !!pair[1], val2 = !!pair[2]) %>%
      mutate(facet1 = pair[1], facet2 = pair[2])
  }) %>%
    bind_rows()

  facet_order <- c("Confidence", "Type1", "RT")



  # 3. Plot with ggplot
  cop_cor = pairs_data %>%
    # filter(trial < 10) %>%
    mutate(
      facet1 = factor(facet1, levels = facet_order),
      facet2 = factor(facet2, levels = facet_order)
    ) %>% filter(facet1 == "Type1", facet2 == "Confidence") %>%
    mutate(Correct = ifelse(facet1 == "Type1" & val1 > 0, "Correct", "Incorrect")) %>%
    ggplot(aes(x = val1, y = val2, color = Correct, fill = Correct)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5, linewidth = 0.6) +
    # geom_smooth(method = "lm", se = FALSE, color = "blue") +
    # stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    facet_wrap(facet1 ~ facet2, scales = "free", ncol = 1) +
    scale_x_continuous(breaks = c(-1.0,-0.5,0,0.5), labels = c("-1.0","-0.5","0.0","0.5"))+
    scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1,0.2), labels = c("-0.2","-0.1","0.0","0.1","0.2"))+
    theme_classic(base_size = 6) +
    scale_color_manual(values = c("darkgreen","darkred"))+
    scale_fill_manual(values = c("darkgreen","darkred"))+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.line  = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          axis.text  = element_text(size = 6),   # or larger if you want them visually thicker
          axis.title.x = element_blank())+
    transition_states(trial)+
    shadow_mark()+
    ease_aes("quintic-in")

  anim <- animate(cop_cor,
                  nframes = 1000,
                  # nframes = 10,
                  fps = 20,
                  height = 11*300/4 / 2.54,        # px
                  width = 11*300/4 / 2.54,       # px
                  res = 300,           # keep consistent with PowerPoint
                  renderer = gifski_renderer(file = here::here("ANI_copcor.gif"), loop = FALSE))



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
    filter(abs(X) < cutoff) %>%
    mutate(Correct = ifelse(Correct == 1, "Correct",ifelse(Correct == 0, "Incorrect",NA)))


  # Plot 2: Actual predictions
  plot_preds = predictionsq_preds %>%
    ggplot() +
    geom_ribbon(aes(x = X, y = mean, ymin = q5, ymax = q95, fill = Correct), alpha = 0.1) +
    geom_ribbon(aes(x = X, y = mean, ymin = q10, ymax = q90, fill = Correct), alpha = 0.3) +
    geom_ribbon(aes(x = X, y = mean, ymin = q20, ymax = q80, fill = Correct), alpha = 0.5) +
    geom_point(data = dataq, aes(x = X, y = mean, fill = as.factor(Correct)),
               shape = 21, color = "black", alpha = 0.5) +
    geom_line(aes(x = X, y = mean, color = Correct), linewidth = 1) +
    facet_wrap(~name, scales = "free", ncol = 3) +
    theme_classic(base_size = 14) +
    labs(color = "Correct", fill = "Correct") +
    # ggtitle("Group predictions (posterior predictive samples)")+
    geom_vline(xintercept = 0, linetype = 2) +
    theme(legend.position = "top")

  return(list(
    plot_combined = combined_plot,
    plot_mean = plot_mean,
    plot_residuals = plot_residuals,
    plot_preds = plot_preds
  ))
}

