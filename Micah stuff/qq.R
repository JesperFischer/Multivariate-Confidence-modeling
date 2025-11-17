
pacman::p_load(tidyverse,posterior,tidybayes,cmdstanr, furrr,ordbetareg, future.apply, patchwork,purrr, bayesplot)


invisible(
  lapply(
    list.files(here::here("Analysis", "Functions"), pattern = "\\.[Rr]$", full.names = TRUE),
    source
  )
)

df = read.csv(here::here("Data","VMP_extero","raw_hrd.csv")) %>%
  mutate(DecisionRT = as.numeric(DecisionRT),
         ConfidenceRT = as.numeric(ConfidenceRT),
         Confidence = as.numeric(Confidence),
         ResponseCorrect = as.numeric(ifelse(ResponseCorrect == "1",1,
                                             ifelse(ResponseCorrect == "1.0",1,
                                                    ifelse(ResponseCorrect == "True",1,
                                                           ifelse(ResponseCorrect == "0",0,
                                                                  ifelse(ResponseCorrect == "0.0",0,
                                                                         ifelse(ResponseCorrect == "False",0,NA))))))),
         Decision = ifelse(Decision == "More",1,ifelse(Decision == "Less",0,NA))) %>%
  mutate(participant_id = as.numeric(as.factor(participant_id))) %>%
  filter(Modality == "Intero",
         session == 1) %>%
  filter(DecisionRT < 8 & DecisionRT > 0.1)%>%
  select(ResponseCorrect,participant_id,DecisionRT,Decision, Confidence, ConfidenceRT, Alpha,cohort)  %>% drop_na() %>%
  rename(
    Confidence = Confidence,
    X  = Alpha,
    Y = Decision,
    Correct = ResponseCorrect,
    RT = DecisionRT,
    subject = participant_id
  ) %>%
  mutate(Confidence = Confidence/100,
         resp = Y)


removers_vmp1 = c("16","167")


df1 = df %>% filter(!subject %in% c("16", "167")) %>%
  filter(cohort == "vmp1") %>%
  mutate(X_scaled = (X - mean(X)) / sd(X))

mod = cmdstanr::cmdstan_model(here::here("Micah stuff","psychos.stan"))

df1$subject = as.numeric(as.factor(df1$subject))


datastan = list(N = nrow(df1),
                S = length(unique(df1$subject)),

                X = df1$X,

                S_id = df1$subject,
                binom_y = df1$Y)


mod_alpha <-mod$sample(
  data = datastan,
  refresh = 100,
  iter_sampling = 500,
  iter_warmup = 500,
  adapt_delta = 0.99,
  max_treedepth = 13,
  init  = 0,
  parallel_chains = 4)



# mod tone:



df = read.csv(here::here("Data","VMP_extero","raw_hrd.csv")) %>%
  mutate(DecisionRT = as.numeric(DecisionRT),
         ConfidenceRT = as.numeric(ConfidenceRT),
         Confidence = as.numeric(Confidence),
         ResponseCorrect = as.numeric(ifelse(ResponseCorrect == "1",1,
                                             ifelse(ResponseCorrect == "1.0",1,
                                                    ifelse(ResponseCorrect == "True",1,
                                                           ifelse(ResponseCorrect == "0",0,
                                                                  ifelse(ResponseCorrect == "0.0",0,
                                                                         ifelse(ResponseCorrect == "False",0,NA))))))),
         Decision = ifelse(Decision == "More",1,ifelse(Decision == "Less",0,NA))) %>%
  mutate(participant_id = as.numeric(as.factor(participant_id))) %>%
  filter(Modality == "Intero",
         session == 1) %>%
  filter(DecisionRT < 8 & DecisionRT > 0.1)%>%
  select(ResponseCorrect,participant_id,DecisionRT,Decision, Confidence, ConfidenceRT, Alpha,cohort,listenBPM)  %>% drop_na() %>%
  rename(
    Confidence = Confidence,
    Y = Decision,
    Correct = ResponseCorrect,
    RT = DecisionRT,
    subject = participant_id
  ) %>%
  mutate(Confidence = Confidence/100,
         resp = Y,
         X  = Alpha + listenBPM)


removers_vmp1 = c("16","167")


df1 = df %>% filter(!subject %in% c("16", "167")) %>%
  filter(cohort == "vmp1") %>%
  mutate(X_scaled = (X - mean(X)) / sd(X))

mod = cmdstanr::cmdstan_model(here::here("Micah stuff","Tone.stan"))

df1$subject = as.numeric(as.factor(df1$subject))


datastan = list(N = nrow(df1),
                S = length(unique(df1$subject)),

                X = df1$X,

                S_id = df1$subject,
                binom_y = df1$Y)


mod_tone <-mod$sample(
  data = datastan,
  refresh = 100,
  iter_sampling = 500,
  iter_warmup = 500,
  adapt_delta = 0.99,
  max_treedepth = 13,
  init  = 0,
  parallel_chains = 4)



toneloo = mod_tone$loo()
alphaloo = mod_alpha$loo()



loo::loo_compare(list(toneloo, alphaloo))

