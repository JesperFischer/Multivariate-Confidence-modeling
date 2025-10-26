remover = function(df){

  avg_acc = df %>% group_by(subject) %>%
    summarize(mean = mean(Correct),
              se = (mean(Correct) * (1-mean(Correct))) /sqrt(n())) %>%
    mutate(measure = "A")

  avg_rt = df %>%
    group_by(subject) %>%
    summarize(mean = mean(RT),
              se = sd(RT)/sqrt(n())) %>%
    mutate(measure = "RT")

  avg_conf = df %>% group_by(subject) %>%
    summarize(mean = mean(Confidence),
              se = sd(Confidence)/sqrt(n())) %>%
    mutate(measure = "C")

  avg_confrt = df %>% mutate(ConfidenceRT = as.numeric(ConfidenceRT)) %>%
    group_by(subject) %>%
    summarize(mean = mean(ConfidenceRT, na.rm = T),
              se = sd(ConfidenceRT)/sqrt(n())) %>%
    mutate(measure = "cRT")



  outliers = rbind(avg_acc,avg_rt,avg_conf,avg_confrt) %>%
    group_by(measure) %>% summarize(meann = mean(mean),
                                    se = sd(mean, na.rm = T)) %>%
    mutate(low_int = meann - 2*se,
           high_int = meann + 2*se) %>%
    mutate(se = NULL)


  dd = inner_join(rbind(avg_acc,
                        avg_rt,
                        avg_conf,
                        avg_confrt),outliers, by = "measure") %>%
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0))


  outlier_labels = dd %>%
    filter(outlier == 1) %>%
    group_by(subject) %>%
    summarize(flag = paste0("(", paste(measure, collapse = ","), ")"))



  plotallsubj = inner_join(rbind(avg_acc,avg_rt,avg_conf,avg_confrt),outliers, by = "measure") %>%
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0)) %>%
    ggplot(aes(x = subject, y = mean, ymin = mean-se, ymax = mean+se, col = as.factor(outlier)))+
    facet_wrap(~measure, scales = "free")+theme_minimal()+geom_pointrange()

  badids = dd %>% filter(outlier == 1) %>% .$subject

  plotdf = inner_join(outlier_labels, df %>% filter(subject %in% badids)) %>%
    mutate(subject = paste0(subject," ", flag))

  if(nrow(plotdf) != 0){
    Group_plot(plotdf)

    n_subj <- length(unique(plotdf %>% .$subject))

    subject_chunks <- split(unique(plotdf$subject), ceiling(seq_along(1:n_subj) / 5))

    badsubject_plots <- lapply(subject_chunks, function(chunk) plot_subjects(plotdf, chunk, bin = 7))
    badsubject_plots

    return(list(plotallsubj,badsubject_plots, badids))
  }
  return(list(plotallsubj))


}





check_modeltype = function(ACC,modeltype,conf){

  if(modeltype == "pure" & ACC == T & conf == "ord_beta"){
    mod = cmdstan_model(here::here("Stanmodels","Confidence_ACC.stan"))

  }else if(modeltype == "meta_un" & ACC == T & conf == "ord_beta"){
    mod = NA

  }else if(modeltype == "meta_un_rt_un" & ACC == T & conf == "ord_beta"){
    mod = NA

  }else if(modeltype == "pure" & ACC == F & conf == "ord_beta"){
    mod = cmdstan_model(here::here("Stanmodels","Resp_Bin_RT_Contin_Conf.stan"))

  }else if(modeltype == "meta_un" & ACC == F & conf == "ord_beta"){
    mod = cmdstan_model(here::here("Stanmodels","Resp_Bin_RT_Contin_Conf_metaun.stan"))

  }else if(modeltype == "meta_un_rt_un" & ACC == F & conf == "ord_beta"){
    mod = cmdstan_model(here::here("Stanmodels","Resp_Bin_RT_Contin_Conf_metaun_rt_un.stan"))
  }


  if(modeltype == "pure" & ACC == T & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","ACC_Bin_RT_Discrete_conf.stan"))

  }else if(modeltype == "meta_un" & ACC == T & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","ACC_Bin_RT_Discrete_conf_metaun.stan"))

  }else if(modeltype == "meta_un_rt_un" & ACC == T & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","ACC_Bin_RT_Discrete_conf_metaun_rt_un.stan"))

  }else if(modeltype == "pure" & ACC == F & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","Resp_Bin_RT_Discrete_conf.stan"))

  }else if(modeltype == "meta_un" & ACC == F & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","Resp_Bin_RT_Discrete_conf_metaun.stan"))

  }else if(modeltype == "meta_un_rt_un" & ACC == F & conf == "discrete_conf"){
    mod = cmdstan_model(here::here("Stanmodels","Discrete Confidence","Resp_Bin_RT_Discrete_conf_metaun_rt_un.stan"))
  }



  return(mod)

}



fit_data_copula_rt = function(df,ACC, outputname,modeltype,conf){


  df$subject = as.numeric(as.factor(df$subject))

  t_p_s = df %>% group_by(subject) %>% summarize(n = n())


  ends <- cumsum(t_p_s$n)

  # Calculate the start points
  starts <- c(1, head(ends, -1) + 1)

  mod = check_modeltype(ACC,modeltype,conf)

  datastan = list(N = nrow(df),
                  S = length(unique(df$subject)),
                  starts = starts,
                  minRT = df %>% group_by(subject) %>% summarize(minRT = min(RT)) %>% .$minRT,
                  ends = ends,
                  t_p_s = t_p_s$n,
                  X = df$X,
                  S_id = df$subject,
                  RT = df$RT,
                  ACC = df$Correct,
                  K = length(unique(df$Confidence)),
                  Conf = df$Confidence,
                  binom_y = df$Y)


  cor <-mod$sample(
    data = datastan,
    refresh = 10,
    iter_sampling = 10,
    iter_warmup = 10,
    adapt_delta = 0.95,
    max_treedepth = 12,
    init  = 0,
    parallel_chains = 4)

  name = paste0("N_",outputname,"ACC_",ACC,"modeltype_",modeltype,"conf_",conf)
  cor$save_object(here::here("Saved models",paste0(name,".rds")))

  return(cor)

}
