

Get_predictive = function(fit,df,n_draws){

  df$subject = as.numeric(as.factor(df$subject))

  workers = 15
  memory = 5000 * 1024^2

  parameters = c("alpha","beta","lapse","rt_int","rt_slope","rt_stim","rt_ndt","rt_prec",
                 "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC","c0","c11","conf_prec",
                 "rho_p_rt","rho_p_conf","rho_rt_conf")

  df_param = as_draws_df(pure$draws(parameters)) %>%
    select(-contains(".")) %>%
    mutate(draw = 1:n()) %>%
    pivot_longer(-draw) %>%
    extract(name, into = c("variable", "subject"),
            regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE)



  # First we select the number of workers and then the memory:
  plan(multisession, workers = workers)
  options(future.globals.maxSize = memory)

  #then the number of subjects
  subjects <- unique(df$subject)
  # and draws per subject
  draws <- 1:n_draws

  # Only use the number of draws that the user wants:
  dfq = df_param %>% filter(draw %in% draws)

  # function to get the draws (goes through subjects and then the draws)
  pred_list <- future_lapply(subjects, function(s) {
    lapply(draws, function(d) {

      # extract parameter vectors (mu, phi, c0 and c1) for that draw and that subject
      params <- dfq %>%
        filter(subject == s, draw == d) %>%
        select(variable, value) %>%
        pivot_wider(names_from = "variable", values_from = "value")


      data = df %>% filter(subject == s)

      psycho_ACC = function(x,alpha,beta,lapse){
        lapse + (1 - 2 * lapse)*((tanh(beta*(x-alpha))/2) + 0.5)
      }
      entropy = function(p){
        -p * log(p) - (1-p) * log(1-p)
      }

      x = data$X
      prob = psycho_ACC(x, params$alpha, params$beta, params$lapse)

      bin_pred = rbinom(length(prob),1,prob)

      rt_mu = params$rt_int + params$rt_slope * entropy(prob) + params$rt_stim * x
      rt_ndt = params$rt_ndt

      # apply inverse logit and scale
      RT_pred = rlnorm(length(rt_mu),rt_mu,params$rt_prec)+rt_ndt



      conf_dat = data.frame(prob = prob) %>%
        rowwise() %>%
        mutate(ACC = list(c(0,1))) %>% unnest(ACC)

      conf_mu = params$conf_int +
        params$conf_ACC * conf_dat$ACC +
        params$conf_entropy * entropy(conf_dat$prob) +
        params$conf_entropy_ACC * conf_dat$ACC * entropy(conf_dat$prob)

      # Collect the responses from the ordered beta
      confidence_pred = data.frame()
      for(i in 1:length(conf_mu)){
        q = rordbeta(1, mu = brms::inv_logit_scaled(conf_mu[i]),
                     phi = params$conf_prec,
                     cutpoints = c(params$c0, params$c0 + exp(params$c11)))

        dat = data.frame(pred = q, ACC = conf_dat$ACC[i])

        confidence_pred = rbind(confidence_pred, dat)
      }

      confidence_pred =  confidence_pred %>% mutate(ACC = ifelse(ACC == 1, "Correct","Incorrect")) %>%
        pivot_wider(names_from = "ACC", values_from = pred, values_fn = mean) %>% unnest(cols = c(Correct, Incorrect))

      predictions = confidence_pred %>% mutate(bin_pred = bin_pred, RT_pred = RT_pred,
                                               x = x, draw = d, subject = s, prob = prob)

      return(predictions)
    })
  },future.seed=TRUE)

  # flatten nested list and create a tidy long dataframe
  return(predictions = map_dfr(pred_list, bind_rows))

}


Plot_bin = function(predictions,df){

  df$subject = as.numeric(as.factor(df$subject))

  sum = df %>% group_by(subject,X) %>% summarize(mean = mean(Y),
                                                      q5 = mean *(Y*(1-Y)) / sqrt(n()),
                                                      q95 = mean *(Y*(1-Y)) / sqrt(n())
  ) %>% mutate(draw = NA)

  psychometrics = predictions %>%
    ggplot() +
    geom_line(aes(x = X, y = prob, group = draw), col = "red") +
    geom_pointrange(data = sum, aes(x = X, y = mean, ymin = q5,ymax = q95, group = draw)) +
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )


  accuarcy = predictions %>% group_by(draw,subject) %>% summarize(p_correct = sum(bin_pred) / n()) %>%
    ggplot() +
    geom_col(data = df %>% mutate(draw = NA) %>% group_by(subject) %>% summarize(p_correct = sum(Correct) / n()),
             aes(x = subject, y = p_correct))+
    geom_point(aes(x = subject,y = p_correct, group = draw), position = position_dodge(width = 0.1), alpha = 0.5)+
    theme_classic(base_size = 16)

  return(list(psychometrics,accuarcy))

}


Plot_RT = function(predictions,df){

  df$subject = as.numeric(as.factor(df$subject))


  # predictions %>%
  #   ggplot() +
  #   geom_line(aes(x = X, y = RT_pred, group = draw), col = "red")+
  #   geom_point(data = df %>% mutate(draw = NA), aes(x = X, y = RT, group = draw))+
  #   facet_wrap(~subject, scales = "free") +
  #   theme_classic(base_size = 16) +
  #   theme(
  #     strip.background = element_blank(),  # remove facet boxes
  #     strip.text = element_blank(),        # remove facet labels
  #     axis.text = element_blank(),         # remove x and y axis numbers
  #     axis.ticks = element_blank(),        # remove tick marks
  #     axis.title = element_blank()         # remove axis titles
  #   )


  conditional_X = predictions %>% group_by(subject,X) %>%
    summarize(mean = mean(RT_pred),
              q5 = quantile(RT_pred,0.05),
              q10 = quantile(RT_pred,0.1),
              q20 = quantile(RT_pred,0.2),
              q95 = quantile(RT_pred,0.95),
              q90 = quantile(RT_pred,0.90),
              q80 = quantile(RT_pred,0.80),
    ) %>%
    ggplot() +
    geom_line(aes(x = X, y = mean), col = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
    geom_point(data = df, aes(x = X, y = RT))+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )


  aggregate_conditional_x = predictions %>% group_by(subject,X) %>%
    summarize(mean = mean(RT_pred),
              q5 = quantile(RT_pred,0.05),
              q10 = quantile(RT_pred,0.1),
              q20 = quantile(RT_pred,0.2),
              q95 = quantile(RT_pred,0.95),
              q90 = quantile(RT_pred,0.90),
              q80 = quantile(RT_pred,0.80),
    ) %>%
    ggplot() +
    geom_line(aes(x = X, y = mean), col = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
    geom_pointrange(data = df %>% group_by(X,subject) %>% summarize(mean = mean(RT), q5 = quantile(RT,0.05),q95 = quantile(RT,0.95)),
                    aes(x = X, y = mean,ymin = q5,ymax = q95))+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )


  Marginal = predictions %>%
    ggplot() +
    geom_histogram(data = df %>% mutate(draw = NA),
                   aes(x = RT, y = after_stat(density)),
                   position = "identity",col = "black", alpha = 0.5)+
    geom_density(aes(x = RT_pred, group = draw), alpha = 0.5)+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )


    return(list(conditional_X,aggregate_conditional_x, Marginal))

}



