Plot_bin = function(predictions,df){

  df$subject = as.numeric(as.factor(df$subject))

  sum = df %>% group_by(subject,X) %>% summarize(mean = mean(Y),
                                                 q5 = mean *(Y*(1-Y)) / sqrt(n()),
                                                 q95 = mean *(Y*(1-Y)) / sqrt(n())
  ) %>% mutate(draw = NA)

  psychometrics = predictions %>%
    ggplot() +
    geom_line(aes(x = x, y = prob, group = draw), col = "red") +
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


  accuarcy = predictions %>% group_by(draw,subject) %>%
    mutate(correct = ifelse(bin_pred == 1 & x > 0,1, ifelse(bin_pred == 0 & x < 0,1,0))) %>%
    summarize(p_correct = sum(correct) / n()) %>%
    ggplot() +
    geom_col(data = df %>% mutate(draw = NA) %>% group_by(subject) %>% summarize(p_correct = sum(Correct) / n()),
             aes(x = subject, y = p_correct))+
    geom_point(aes(x = subject,y = p_correct, group = draw), position = position_dodge(width = 0.1), alpha = 0.5)+
    theme_classic(base_size = 16)


  resp = predictions %>% group_by(draw,subject) %>%
    mutate(correct = ifelse(bin_pred == 1 & x > 0,1, ifelse(bin_pred == 0 & x < 0,1,0))) %>%
    summarize(p_correct = sum(bin_pred) / n()) %>%
    ggplot() +
    geom_col(data = df %>% mutate(draw = NA) %>% group_by(subject) %>% summarize(p_correct = sum(Y) / n()),
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


  conditional_X = predictions %>% group_by(subject,x) %>%
    summarize(mean = mean(RT_pred),
              q5 = quantile(RT_pred,0.05),
              q10 = quantile(RT_pred,0.1),
              q20 = quantile(RT_pred,0.2),
              q95 = quantile(RT_pred,0.95),
              q90 = quantile(RT_pred,0.90),
              q80 = quantile(RT_pred,0.80),
    ) %>%
    ggplot() +
    geom_line(aes(x = x, y = mean), col = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
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


  aggregate_conditional_x = predictions %>% group_by(subject,x) %>%
    summarize(mean = mean(RT_pred),
              q5 = quantile(RT_pred,0.05),
              q10 = quantile(RT_pred,0.1),
              q20 = quantile(RT_pred,0.2),
              q95 = quantile(RT_pred,0.95),
              q90 = quantile(RT_pred,0.90),
              q80 = quantile(RT_pred,0.80),
    ) %>%
    ggplot() +
    geom_line(aes(x = x, y = mean), col = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q5,ymax = q95),alpha = 0.1, fill = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q10,ymax = q90),alpha = 0.3, fill = "red")+
    geom_ribbon(aes(x = x, y = mean,ymin = q20,ymax = q80),alpha = 0.5, fill = "red")+
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





Plot_conf = function(predictions,df){

  df$subject = as.numeric(as.factor(df$subject))


  conditional_X = predictions %>%
    mutate(correct = ifelse(bin_pred == 1 & x > 0,1, ifelse(bin_pred == 0 & x < 0,1,0))) %>%
    pivot_longer(cols = c("Incorrect","Correct"), names_to = "Correct",values_to = "Confidence")%>%
    group_by(subject,Correct,x) %>%
    summarize(mean = mean(Confidence),
              q5 = quantile(Confidence,0.05),
              q10 = quantile(Confidence,0.1),
              q20 = quantile(Confidence,0.2),
              q95 = quantile(Confidence,0.95),
              q90 = quantile(Confidence,0.90),
              q80 = quantile(Confidence,0.80)
    ) %>%
    ggplot() +
    geom_line(aes(x = x, y = mean, col = as.factor(Correct)))+
    geom_ribbon(aes(x = x, y = mean,ymin = q5,ymax = q95, fill = as.factor(Correct)),alpha = 0.1)+
    geom_ribbon(aes(x = x, y = mean,ymin = q10,ymax = q90, fill = as.factor(Correct)),alpha = 0.3)+
    geom_ribbon(aes(x = x, y = mean,ymin = q20,ymax = q80, fill = as.factor(Correct)),alpha = 0.5)+
    geom_point(data = df, aes(x = X, y = Confidence / 100, fill = as.factor(Correct)), col = "black",shape = 21)+
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
