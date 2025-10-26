

Get_predictive = function(fit,df,n_draws, model){

  if(model == "pure"){
    df$subject = as.numeric(as.factor(df$subject))

    workers = 15
    memory = 8000 * 1024^2

    parameters = c("alpha","beta","lapse","rt_int","rt_slope","rt_stim","rt_ndt","rt_prec",
                   "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC","c0","c11","conf_prec",
                   "rho_p_rt","rho_p_conf","rho_rt_conf")

    df_param = as_draws_df(fit$draws(parameters)) %>%
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
          lapse + (1 - 2 * lapse)*(brms::inv_logit_scaled(beta*(x-alpha)))
        }
        entropy = function(p){
          -p * log(p) - (1-p) * log(1-p)
        }

        x = data$X
        prob = psycho_ACC(x, params$alpha, exp(params$beta), params$lapse)

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
          pivot_wider(names_from = "ACC", values_from = pred,values_fn = list) %>% unnest(cols = c(Correct, Incorrect))

        predictions = confidence_pred %>% mutate(bin_pred = bin_pred, RT_pred = RT_pred,
                                                 X = x, draw = d, subject = s, prob = prob)

        return(predictions)
      })
    },future.seed=TRUE)

    # flatten nested list and create a tidy long dataframe
    return(predictions = map_dfr(pred_list, bind_rows))
  }
  if(model == "full"){
    df$subject = as.numeric(as.factor(df$subject))

    workers = 15
    memory = 10000 * 1024^2

    parameters = c("alpha","beta","lapse","rt_int","rt_slope","rt_stim","rt_ndt","rt_prec",
                   "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC","c0","c11","conf_prec",
                   "meta_un","rt_un",
                   "rho_p_rt","rho_p_conf","rho_rt_conf")

    df_param = as_draws_df(fit$draws(parameters)) %>%
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
          lapse + (1 - 2 * lapse)*(brms::inv_logit_scaled(beta*(x-alpha)))
        }
        entropy = function(p){
          -p * log(p) - (1-p) * log(1-p)
        }

        x = data$X
        prob = psycho_ACC(x, params$alpha, exp(params$beta), params$lapse)

        bin_pred = rbinom(length(prob),1,prob)

        rt_mu = params$rt_int + params$rt_slope * entropy(psycho_ACC(x, params$alpha, exp(params$beta + params$rt_un), params$lapse)) + params$rt_stim * x
        rt_ndt = params$rt_ndt

        # apply inverse logit and scale
        RT_pred = rlnorm(length(rt_mu),rt_mu,params$rt_prec)+rt_ndt


        conf_prob = entropy(psycho_ACC(x, params$alpha, exp(params$beta + params$meta_un), params$lapse))

        conf_dat = data.frame(prob = conf_prob) %>%
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
          pivot_wider(names_from = "ACC", values_from = pred,values_fn = list) %>% unnest(cols = c(Correct, Incorrect))

        predictions = confidence_pred %>% mutate(bin_pred = bin_pred, RT_pred = RT_pred,
                                                 X = x, draw = d, subject = s, prob = prob)

        return(predictions)
      })
    },future.seed=TRUE)

    # flatten nested list and create a tidy long dataframe
    return(predictions = map_dfr(pred_list, bind_rows))

  }

}


Plot_bin = function(predictions,df,ACC){

  df$subject = as.numeric(as.factor(df$subject))

  sum = df %>% group_by(subject,X) %>% summarize(mean = mean(Y),
                                                      q5 = mean *(Y*(1-Y)) / sqrt(n()),
                                                      q95 = mean *(Y*(1-Y)) / sqrt(n())
  ) %>% mutate(draw = NA)

  psychometrics = predictions %>%
    ggplot() +
    geom_pointrange(data = sum, aes(x = X, y = mean, ymin = q5,ymax = q95, group = draw), alpha = 0.5, shape = 21) +
    geom_line(aes(x = X, y = prob, group = draw), col = "red") +
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )

  if(ACC == T){
    accuarcy = predictions %>% group_by(draw,subject) %>% summarize(p_correct = sum(bin_pred) / n()) %>%
      ggplot() +
      geom_col(data = df %>% mutate(draw = NA) %>% group_by(subject) %>% summarize(p_correct = sum(Correct) / n()),
               aes(x = subject, y = p_correct))+
      geom_point(aes(x = subject,y = p_correct, group = draw), position = position_dodge(width = 0.1), alpha = 0.5)+
      theme_classic(base_size = 16)

    return(list(psychometrics,accuarcy))
  }

  if(ACC == F){
    resp = predictions %>% group_by(draw,subject) %>% summarize(p_resp = sum(bin_pred) / n()) %>%
      ggplot() +
      geom_col(data = df %>% mutate(draw = NA) %>% group_by(subject) %>% summarize(p_resp = sum(Y) / n()),
               aes(x = subject, y = p_resp))+
      geom_point(aes(x = subject,y = p_resp, group = draw), position = position_dodge(width = 0.1), alpha = 0.5)+
      theme_classic(base_size = 16)

    return(list(psychometrics,resp))
  }

}


Plot_RT = function(predictions,df, bin = 7){

  df$subject = as.numeric(as.factor(df$subject))

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
    geom_ribbon(aes(x = X, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = X, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
    geom_point(data = df, aes(x = X, y = RT), shape = 21, alpha = 0.5)+
    geom_line(aes(x = X, y = mean), col = "red")+
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
    geom_density(aes(x = RT_pred, group = draw), alpha = 0.5, col = "red")+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )



  # compute bin breaks and midpoints first
  bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
  bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

  # helper to add numeric midpoint
  df <- df %>%
    mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
           X_mid = bin_mids[as.integer(X_bin)])

  predictions <- predictions %>%
    mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
           X_mid = bin_mids[as.integer(X_bin)])



  conditional_X_binned = predictions %>% group_by(subject,X_mid) %>%
    summarize(mean = mean(RT_pred),
              q5 = quantile(RT_pred,0.05),
              q10 = quantile(RT_pred,0.1),
              q20 = quantile(RT_pred,0.2),
              q95 = quantile(RT_pred,0.95),
              q90 = quantile(RT_pred,0.90),
              q80 = quantile(RT_pred,0.80),
    ) %>%
    ggplot() +
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
    geom_pointrange(data = df %>% group_by(subject,X_mid) %>% summarize(median = median(RT), q5 = quantile(RT,0.05), q95 = quantile(RT,0.95)),
                    aes(x = X_mid, y = median, ymin = q5, ymax = q95),
                    shape = 21, alpha = 0.5)+
    geom_line(aes(x = X_mid, y = mean), col = "red")+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )





    return(list(conditional_X,conditional_X_binned, Marginal))

}


Plot_conf = function(predictions,df, bin = 7){

  df$subject = as.numeric(as.factor(df$subject))

  Marginal = predictions %>%
    pivot_longer(cols = c("Incorrect","Correct"), values_to = "Confidence",names_to = "Correct") %>%
    ggplot() +
    geom_histogram(data = df %>% mutate(draw = NA)%>% mutate(Correct = ifelse(Correct == 1, "Correct","Incorrect")),
                   aes(x = Confidence, y = after_stat(density), fill = Correct),
                   position = "identity", alpha = 0.5, col = "black")+
    geom_density(aes(x = Confidence, group = interaction(Correct,draw), col = Correct), alpha = 0.5)+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )



  # compute bin breaks and midpoints first
  bin_breaks <- seq(min(df$X, na.rm = TRUE), max(df$X, na.rm = TRUE), length.out = bin + 1)
  bin_mids <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

  # helper to add numeric midpoint
  df <- df %>%
    mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
           X_mid = bin_mids[as.integer(X_bin)])

  predictions <- predictions %>%
    mutate(X_bin = cut(X, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
           X_mid = bin_mids[as.integer(X_bin)])



  conditional_X_binned = predictions %>%
    pivot_longer(cols = c("Incorrect","Correct"), values_to = "Confidence",names_to = "Correct") %>%
    group_by(subject,X_mid,Correct) %>%
    summarize(mean = mean(Confidence),
              q5 = quantile(Confidence,0.05),
              q10 = quantile(Confidence,0.1),
              q20 = quantile(Confidence,0.2),
              q95 = quantile(Confidence,0.95),
              q90 = quantile(Confidence,0.90),
              q80 = quantile(Confidence,0.80),
    ) %>%
    ggplot() +
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q5,ymax = q95, fill = Correct),alpha = 0.1)+
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q10,ymax = q90, fill = Correct),alpha = 0.3)+
    geom_ribbon(aes(x = X_mid, y = mean,ymin = q20,ymax = q80, fill = Correct),alpha = 0.5)+
    geom_pointrange(data = df %>%
                      mutate(Correct = ifelse(Correct == 1, "Correct","Incorrect")) %>%
                      group_by(subject,X_mid,Correct) %>%
                      summarize(median = median(Confidence),
                                q5 = quantile(Confidence,0.05),
                                q95 = quantile(Confidence,0.95)),
                    aes(x = X_mid, y = median, ymin = q5, ymax = q95, fill = Correct), col = "black",shape = 21)+
    geom_line(aes(x = X_mid, y = mean, color = Correct))+
    facet_wrap(~subject, scales = "free") +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_blank(),  # remove facet boxes
      strip.text = element_blank(),        # remove facet labels
      axis.text = element_blank(),         # remove x and y axis numbers
      axis.ticks = element_blank(),        # remove tick marks
      axis.title = element_blank()         # remove axis titles
    )





  return(list(conditional_X_binned, Marginal))

}



Get_predictive_group = function(fit,df,n_draws, model){

  if(model == "pure"){
    df$subject = as.numeric(as.factor(df$subject))

    workers = 15
    memory = 8000 * 1024^2

    parameters = c("alpha","beta","lapse","rt_int","rt_slope","rt_prec","rt_stim",
                   "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC","conf_prec")

    df_param = as_draws_df(fit$draws("gm")) %>%
      select(-contains(".")) %>%
      rename_with(~parameters) %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(-draw, names_to = "variable")


    constants = as_draws_df(fit$draws(c("rt_ndt","c0","c11"))) %>%
      select(-contains(".")) %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(-draw) %>%
      extract(name, into = c("variable", "subject"),
              regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE) %>% group_by(variable) %>%
      summarize(mean = mean(value)) %>% pivot_wider(names_from = variable,values_from = mean)


    # First we select the number of workers and then the memory:
    plan(multisession, workers = workers)
    options(future.globals.maxSize = memory)


    # and draws per subject
    draws <- 1:n_draws

    # Only use the number of draws that the user wants:
    dfq = df_param %>% filter(draw %in% draws)

    # function to get the draws (goes through subjects and then the draws)
    pred_list <- future_lapply(draws, function(d) {

        # extract parameter vectors (mu, phi, c0 and c1) for that draw and that subject
        params <- dfq %>%
          filter(draw == d) %>%
          select(variable, value) %>%
          pivot_wider(names_from = "variable", values_from = "value")


        psycho_ACC = function(x,alpha,beta,lapse){
          lapse + (1 - 2 * lapse)*(brms::inv_logit_scaled(beta*(x-alpha)))
        }
        entropy = function(p){
          -p * log(p) - (1-p) * log(1-p)
        }

        x = seq(-50,50,by = 0.5)
        prob = psycho_ACC(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse)/2)

        bin_pred = rbinom(length(prob),1,prob)

        rt_mu = params$rt_int + params$rt_slope * entropy(prob) + params$rt_stim * x

        # apply inverse logit and scale
        RT_pred = rlnorm(length(rt_mu),rt_mu,exp(params$rt_prec))+constants$rt_ndt

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
                       phi = exp(params$conf_prec),
                       cutpoints = c(constants$c0, constants$c0 + exp(constants$c11)))

          dat = data.frame(pred = q, ACC = conf_dat$ACC[i])

          confidence_pred = rbind(confidence_pred, dat)
        }

        confidence_pred =  confidence_pred %>% mutate(ACC = ifelse(ACC == 1, "Correct","Incorrect")) %>%
          pivot_wider(names_from = "ACC", values_from = pred,values_fn = list) %>% unnest(cols = c(Correct, Incorrect))

        predictions = confidence_pred %>% mutate(bin_pred = bin_pred, RT_pred = RT_pred,
                                                 X = x, draw = d, prob = prob)

        return(predictions)
    },future.seed=TRUE)

    # flatten nested list and create a tidy long dataframe
    return(predictions = map_dfr(pred_list, bind_rows))
  }
  if(model == "full"){
    df$subject = as.numeric(as.factor(df$subject))

    workers = 15
    memory = 10000 * 1024^2

    parameters = c("alpha","beta","lapse","rt_int","rt_slope","rt_prec","rt_stim",
                   "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC","conf_prec",
                   "meta_un","rt_un")

    df_param = as_draws_df(fit$draws("gm")) %>%
      select(-contains(".")) %>%
      rename_with(~parameters) %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(-draw, names_to = "variable")



    constants = as_draws_df(fit$draws(c("rt_ndt","c0","c11"))) %>%
      select(-contains(".")) %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(-draw) %>%
      extract(name, into = c("variable", "subject"),
              regex = "([a-zA-Z0-9_]+)\\[(\\d+)\\]", convert = TRUE) %>% group_by(variable) %>%
      summarize(mean = mean(value)) %>% pivot_wider(names_from = variable,values_from = mean)


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
    pred_list <- future_lapply(draws, function(d) {

        # extract parameter vectors (mu, phi, c0 and c1) for that draw and that subject
        params <- dfq %>%
          filter(draw == d) %>%
          select(variable, value) %>%
          pivot_wider(names_from = "variable", values_from = "value")

        psycho_ACC = function(x,alpha,beta,lapse){
          lapse + (1 - 2 * lapse)*(brms::inv_logit_scaled(beta*(x-alpha)))
        }
        entropy = function(p){
          -p * log(p) - (1-p) * log(1-p)
        }

        x = seq(-50,50,by = 0.5)
        prob = psycho_ACC(x, params$alpha, exp(params$beta), brms::inv_logit_scaled(params$lapse)/2)

        bin_pred = rbinom(length(prob),1,prob)

        rt_mu = params$rt_int + params$rt_slope * entropy(psycho_ACC(x, params$alpha, exp(params$beta + params$rt_un), brms::inv_logit_scaled(params$lapse)/2)) + params$rt_stim * x

        # apply inverse logit and scale
        RT_pred = rlnorm(length(rt_mu),rt_mu,exp(params$rt_prec))+constants$rt_ndt

        conf_prob = entropy(psycho_ACC(x, params$alpha, exp(params$beta + params$meta_un), brms::inv_logit_scaled(params$lapse)/2))

        conf_dat = data.frame(prob = conf_prob) %>%
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
                       phi = exp(params$conf_prec),
                       cutpoints = c(constants$c0, constants$c0 + exp(constants$c11)))

          dat = data.frame(pred = q, ACC = conf_dat$ACC[i])

          confidence_pred = rbind(confidence_pred, dat)
        }

        confidence_pred =  confidence_pred %>% mutate(ACC = ifelse(ACC == 1, "Correct","Incorrect")) %>%
          pivot_wider(names_from = "ACC", values_from = pred,values_fn = list) %>% unnest(cols = c(Correct, Incorrect))

        predictions = confidence_pred %>% mutate(bin_pred = bin_pred, RT_pred = RT_pred,
                                                 X = x, draw = d, prob = prob)

        return(predictions)
    },future.seed=TRUE)

    # flatten nested list and create a tidy long dataframe
    return(predictions = map_dfr(pred_list, bind_rows))

  }

}


group_predictive = function(predictions,df){

    dataq = bind_rows(
      df%>% mutate(Correct = ifelse(Correct == 1, "Correct","Incorrect")) %>%
        group_by(X) %>%
        summarize(name = "Type-1",
                  mean = mean(Y),
                  q5 = mean-2* sqrt(mean(Y)*(1-mean(Y)) / sqrt(n())),
                  q95 = mean+2* sqrt(mean(Y)*(1-mean(Y)) / sqrt(n()))),

      df%>% mutate(Correct = ifelse(Correct == 1, "Correct","Incorrect")) %>%
        group_by(X) %>%
        summarize(name = "RT",
                  mean = mean(RT),
                  q5 = mean(RT) - 2 * (sd(RT) / sqrt(n())),
                  q95 = mean(RT) + 2 * (sd(RT) / sqrt(n()))),

      df%>% mutate(Correct = ifelse(Correct == 1, "Correct","Incorrect")) %>%
        group_by(X, Correct) %>%   # <-- group by correctness only here
        summarize(name = "Confidence",
                  mean = mean(Confidence),
                  q5 = mean(Confidence) - 2 * (sd(Confidence) / sqrt(n())),
                  q95 = mean(Confidence) + 2 * (sd(Confidence) / sqrt(n()))),
    )%>%
      filter(abs(X) < 25)



    predictionsq = bind_rows(
      predictions %>%
        pivot_longer(cols = c("Incorrect","Correct"), values_to = "Confidence",names_to = "Correct") %>%
        filter(Correct == "Correct") %>%
        group_by(X) %>%
        summarize(name = "Type-1",
                  mean = mean(prob),
                  q5 = quantile(prob,0.05),
                  q10 = quantile(prob,0.1),
                  q20 = quantile(prob,0.2),
                  q95 = quantile(prob,0.95),
                  q90 = quantile(prob,0.90),
                  q80 = quantile(prob,0.80)),

      predictions %>%
        pivot_longer(cols = c("Incorrect","Correct"), values_to = "Confidence",names_to = "Correct") %>%
        group_by(X) %>%
        summarize(name = "RT",
                  mean = mean(RT_pred),
                  q5 = quantile(RT_pred,0.05),
                  q10 = quantile(RT_pred,0.1),
                  q20 = quantile(RT_pred,0.2),
                  q95 = quantile(RT_pred,0.95),
                  q90 = quantile(RT_pred,0.90),
                  q80 = quantile(RT_pred,0.80)),

      predictions %>%
        pivot_longer(cols = c("Incorrect","Correct"), values_to = "Confidence",names_to = "Correct") %>%
        group_by(X, Correct) %>%   # <-- group by correctness only here
        summarize(name = "Confidence",
                  mean = mean(Confidence),
                  q5 = quantile(Confidence,0.05),
                  q10 = quantile(Confidence,0.1),
                  q20 = quantile(Confidence,0.2),
                  q95 = quantile(Confidence,0.95),
                  q90 = quantile(Confidence,0.90),
                  q80 = quantile(Confidence,0.80))
    ) %>%
      filter(abs(X) < 25)



   plot =  predictionsq %>%
      ggplot() +
        geom_ribbon(aes(x = X, y = mean,ymin = q5,ymax = q95, fill = Correct),alpha = 0.1)+
        geom_ribbon(aes(x = X, y = mean,ymin = q10,ymax = q90, fill = Correct),alpha = 0.3)+
        geom_ribbon(aes(x = X, y = mean,ymin = q20,ymax = q80, fill = Correct),alpha = 0.5)+
        geom_pointrange(data = dataq , aes(x = X, y = mean, ymin = q5, ymax = q95,color = as.factor(Correct))) +   # only affects responseConf
        geom_line(aes(x = X, y = mean, color = Correct))+
      facet_wrap(~name, scales = "free", ncol = 1) +
      theme_classic(base_size = 16) +
      labs(color = "Correct")+
      geom_vline(xintercept = 0, linetype = 2)+
      theme(legend.position = "top")

  return(plot)


}
