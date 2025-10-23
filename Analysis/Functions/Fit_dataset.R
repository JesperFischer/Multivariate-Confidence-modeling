remover = function(df){
  avg_acc = df %>% group_by(participant_id) %>% summarize(mean = mean(ResponseCorrect)) %>% mutate(se = NA, measure = "correct")
  avg_rt = df %>% group_by(participant_id) %>% summarize(mean = mean(DecisionRT), se = sd(DecisionRT)/sqrt(n())) %>% mutate(measure = "rt")
  avg_conf = df %>% group_by(participant_id) %>% summarize(mean = mean(Confidence), se = sd(Confidence)/sqrt(n())) %>% mutate(measure = "confidence")
  avg_confrt = df %>% mutate(ConfidenceRT = as.numeric(ConfidenceRT)) %>%
    group_by(participant_id) %>% summarize(mean = mean(ConfidenceRT, na.rm = T), se = sd(ConfidenceRT)/sqrt(n())) %>% mutate(measure = "confidencert")


  outliers = rbind(avg_acc,avg_rt,avg_conf,avg_confrt) %>% group_by(measure) %>% summarize(meann = mean(mean), se = sd(mean, na.rm = T)) %>%
    mutate(low_int = meann - 2*se,high_int = meann + 2*se) %>% mutate(se = NULL)

  dd = inner_join(rbind(avg_acc,avg_rt,avg_conf,avg_confrt),outliers, by = "measure") %>%
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0))


  inner_join(rbind(avg_acc,avg_rt,avg_conf,avg_confrt),outliers, by = "measure") %>%
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0)) %>%
    ggplot(aes(x = participant_id, y = mean, ymin = mean-se, ymax = mean+se, col = as.factor(outlier)))+
    facet_wrap(~measure, scales = "free")+theme_minimal()+geom_pointrange()

  badids = dd %>% filter(measure == "correct" & outlier == 1) %>% .$participant_id

  df%>% filter(cohort == "vmp1") %>%
    filter(subject %in% remover(df%>% filter(cohort == "vmp1")))

  return()
}









fit_data_copula_rt = function(df,ACC, outputname){


  df$subject = as.numeric(as.factor(df$subject))

  t_p_s = df %>% group_by(subject) %>% summarize(n = n())


  ends <- cumsum(t_p_s$n)

  # Calculate the start points
  starts <- c(1, head(ends, -1) + 1)

  if(ACC){
  mod = cmdstan_model(here::here("Stanmodels","Confidence_ACC.stan"))
  }else{
  mod = cmdstan_model(here::here("Stanmodels","Confidence_Standard.stan"))
  }

  datastan = list(N = nrow(df),
                  S = length(unique(df$subject)),
                  starts = starts,
                  minRT = df %>% group_by(subject) %>% summarize(minRT = min(RT)) %>% .$minRT,
                  ends = ends,
                  t_p_s = t_p_s$n,
                  X = df$X,
                  S_id = df$subject,
                  RT = df$RT,
                  binom_y = df$Y)


  cor <-mod$sample(
    data = datastan,
    refresh = 10,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    max_treedepth = 12,
    init  = 0,
    parallel_chains = 4)


  cor$save_object(here::here("Saved models",paste0(outputname,".rds")))

  return(cor)

}
