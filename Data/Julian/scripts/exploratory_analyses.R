## EXPLORATORY ANALYSES ----
# For Matthews et al. (Communications Psychology, 2025).

# Clear environment
rm(list = ls())

# Which analyses?
do_frequentist_analysis <- FALSE
do_accuracy_analysis <- FALSE
do_confidence_analysis <- FALSE
do_confid_database_analysis <- FALSE
do_hmetad_analysis <- FALSE
do_regression_analysis <- FALSE

# Analysis precision
num_iterations <- 1e5 #1e5 is default

# PRELIMINARIES ----

# Load packages
library(tidyverse)
library(see)
library(BayesFactor)
library(bayestestR)

# Load custom functions
source("supporting_functions/my_functions.R")
source("supporting_functions/plot_functions.R")

# Locate folders
data_location <- "../../data/csvs/"
fig_location <- ""

# Colours
baseline_col <- "#fdf8e2"
time_col <- "#93cfc8"
action_col <- "#fbab88"
context_col <- c(baseline_col, time_col, action_col)

# Plot dimensions
regular_width <- 8
regular_height <- 6

# PREPARE DATA ----

the_experiments <- c("experiment1","experiment2")

for(task_version in the_experiments){
  
  # Load data
  g_data <- read.csv(paste0(data_location,"block_data_",task_version,".csv"), na.strings = "NaN")
  t_data <- read.csv(paste0(data_location,"trial_data_",task_version,".csv"), na.strings = "NaN")
  d_data <- read.csv(paste0(data_location,"demographics_",task_version,".csv"), na.strings = "NaN")
  
  # Compute total number of participants before exclusions
  pp_number <- length(d_data$subjID)
  
  # Excluded participants
  excluded_participants <- switch(task_version,
                                  "experiment1" = 
                                    c("008_RM", # This person had to quit after Day #2
                                      "105_KK"), # This person quit after Day #1
                                  "experiment2" = 
                                    c("402_MS") # This person quit after Day #1
  )
  
  # Add any uncompleted participant IDs to excluded list
  exclusions <- unique(g_data[g_data$completed==FALSE,"subjID"])
  excluded_participants <- unique(c(excluded_participants, exclusions))
  
  # Filter out excluded blocks (assign as)
  g_data$include_block <- as.logical(g_data$include_block)
  g_data <- g_data %>%
    mutate(
      efficiency = if_else(include_block, efficiency, NA_real_),
      mean_accuracy = if_else(include_block, mean_accuracy, NA_real_),
      mean_confidence = if_else(include_block, mean_confidence, NA_real_),
      median_confidence = if_else(include_block, median_confidence, NA_real_),
      mean_total_RT = if_else(include_block, mean_total_RT, NA_real_),
      mean_log_RT = if_else(include_block, mean_log_RT, NA_real_),
      sensitivity = if_else(include_block, sensitivity, NA_real_),
      criterion = if_else(include_block, criterion, NA_real_),
      metad = if_else(include_block, metad, NA_real_),
      m_difference = if_else(include_block, m_difference, NA_real_)
    )
  
  g_data <- g_data |> filter(!subjID %in% excluded_participants)
  t_data <- t_data |> filter(!subjID %in% excluded_participants)
  d_data <- d_data |> filter(!subjID %in% excluded_participants)
  
  # Class coding & assignments
  g_data <- assign_group_variables(g_data)
  t_data <- assign_trial_variables(t_data)
  d_data <- assign_demographic_variables(d_data)
  
  # Assign to appropriate experiment
  switch(task_version,
         "experiment1" = {
           group_data.exp1 <- g_data
           trial_data.exp1 <- t_data
           demographics.exp1 <- d_data
           how_many_participants.exp1 <- pp_number
         },
         "experiment2" = {
           group_data.exp2 <- g_data
           trial_data.exp2 <- t_data
           demographics.exp2 <- d_data
           how_many_participants.exp2 <- pp_number
         }
  )
}

# Concatenate data
g_data <- rbind(group_data.exp1, group_data.exp2)
t_data <- rbind(trial_data.exp1, trial_data.exp2)
d_data <- rbind(demographics.exp1, demographics.exp2)

# FREQUENTIST ANALYSES ----

if (do_frequentist_analysis) {
  cat("Running frequentist analyses...\n")
  
  library(lmerTest)
  library(emmeans)
  library(effectsize)
  
  ## HYPOTHESIS 1: ----
  # EXP1: F(1,339.96)=129.21, p<.001, ηp2=.28
  # EXP2: F(1,222.17)=91.92, p<.001, ηp2=.29
  
  # Perception
  test_data <- group_data.exp1 |> 
    filter(include_block == TRUE) |> 
    filter(motor_time_type %in% c("baseline", "time limited"))
  hypothesis1.exp1.mdl <- lmer(mean_total_RT ~ motor_time_type + (1|subjID),
                               data = test_data)
  hypothesis1.exp1.stat <- anova(hypothesis1.exp1.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis1.exp1.eta <- eta_squared(hypothesis1.exp1.stat)
  
  # Memory
  test_data <- group_data.exp2 |> 
    filter(include_block == TRUE) |>
    filter(motor_time_type %in% c("baseline", "time limited"))
  hypothesis1.exp2.mdl <- lmer(mean_total_RT ~ motor_time_type + (1|subjID),
                               data = test_data)
  hypothesis1.exp2.stat <- anova(hypothesis1.exp2.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis1.exp2.eta <- eta_squared(hypothesis1.exp2.stat)
  
  ## HYPOTHESIS 2: ----
  # EXP1 ORDER: F(2,332.46)=1.38, p=.253, ηp2<.01
  # EXP1 ORDER*TIME: F(2,331.05)=2.20, p=.112, ηp2=.01
  
  # EXP2 ORDER: F(2,220.92)=1.23, p=.295, ηp2=.01
  # EXP2 ORDER*TIME: F(1,218.41)=1.10, p=.335, ηp2<.01
  
  test_data <- group_data.exp1 |> 
    filter(include_block == TRUE) |> 
    filter(motor_time_type %in% c("baseline", "time limited"))
  hypothesis234.exp1.mdl <- lmer(
    efficiency ~ response_type * motor_time_type + (1|subjID),
    data = test_data)
  hypothesis234.exp1.stat <- anova(hypothesis234.exp1.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis234.exp1.eta <- eta_squared(hypothesis234.exp1.stat)
  
  test_data <- group_data.exp2 |> 
    filter(include_block == TRUE) |> 
    filter(motor_time_type %in% c("baseline", "time limited"))
  hypothesis234.exp2.mdl <- lmer(
    efficiency ~ response_type * motor_time_type + (1|subjID),
    data = test_data)
  hypothesis234.exp2.stat <- anova(hypothesis234.exp2.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis234.exp2.eta <- eta_squared(hypothesis234.exp2.stat)
  
  ## HYPOTHESIS 3: ----
  # EXP1 TIME: F(1,330.49)=0.16, p=.689, ηp2<.01
  # EXP2 TIME: F(1,218.93)=0.32, p=.571, ηp2<.01
  
  ## HYPOTHESIS 4: ----
  # EXP1 BASELINE D>C | D+C: t(338)=2.372, ptukey=0.048 *
  # EXP1 BASELINE D>C | C>D: t(337)=2.225, ptukey=0.069
  # EXP1 BASELINE D+C | C>D: t(339)=-0.150, ptukey=0.988
  
  # EXP1 TIME D>C | D+C: t(339)=-0.339, ptukey=0.939
  # EXP1 TIME D>C | D+C: t(338)=-0.179, ptukey=0.982
  # EXP1 TIME D>C | D+C: t(339)=-0.150, ptukey=0.986
  
  # EXP2 BASELINE D>C | D+C: t(221)=1.413, ptukey=0.336
  # EXP2 BASELINE D>C | C>D: t(222)=2.028, ptukey=0.108
  # EXP2 BASELINE D+C | C>D: t(220)=0.611, ptukey=0.814
  
  # EXP2 TIME D>C | D+C: t(223)=0.493, ptukey=0.875
  # EXP2 TIME D>C | D+C: t(223)=-0.072, ptukey=0.997
  # EXP2 TIME D>C | D+C: t(221)=-0.566, ptukey=0.838
  
  # Pairwise comparisons of context main effect
  h4.exp1 <- emmeans(hypothesis234.exp1.mdl, ~  response_type | motor_time_type)
  hypothesis234.exp1.pairs <- pairs(h4.exp1)
  
  h4.exp2 <- emmeans(hypothesis234.exp2.mdl, ~  response_type | motor_time_type)
  hypothesis234.exp2.pairs <- pairs(h4.exp2)
  
  ## HYPOTHESIS 5: ----
  # EXP1 ORDER*MOTOR: F(2,334.02)=0.80, p=.449, ηp2<.01
  # EXP2 ORDER*MOTOR: F(2,222.99)=0.43, p=.654, ηp2<.01
  
  # Perception
  test_data <- group_data.exp1 |> 
    filter(include_block == TRUE) |> 
    filter(motor_time_type %in% c("baseline", "motor fixed"))
  hypothesis5.exp1.mdl <- lmer(
    efficiency ~ response_type * motor_time_type + (1|subjID),
    data = test_data)
  hypothesis5.exp1.stat <- anova(hypothesis5.exp1.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis5.exp1.eta <- eta_squared(hypothesis5.exp1.stat)
  
  # Memory
  test_data <- group_data.exp2 |> 
    filter(include_block == TRUE) |> 
    filter(motor_time_type %in% c("baseline", "motor fixed"))
  hypothesis5.exp2.mdl <- lmer(
    efficiency ~ response_type * motor_time_type + (1|subjID),
    data = test_data)
  hypothesis5.exp2.stat <- anova(hypothesis5.exp2.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis5.exp2.eta <- eta_squared(hypothesis5.exp2.stat)
  
  # HYPOTHESIS 6: ----
  # EXP1 ORDER: F(2,134.31)=0.45, p=.639, ηp2<.01
  # EXP1 DIRECTIONAL D>C | C>D: t(63)=0.605, p=0.274, 95%CI[-0.101, INF]
  
  # EXP1 MOTOR D>C | D+C: t(134)=0.652, ptukey=0.792
  # EXP1 MOTOR D>C | C>D: t(132)=0.924, ptukey=0.626
  # EXP1 MOTOR D+C | C>D: t(132)=0.267, ptukey=0.961
  
  # EXP2 ORDER: F(2,87.18)=0.22, p=.806, ηp2<.01
  # EXP2 DIRECTIONAL D>C | C>D: t(45)=0.655, p=0.258, 95%CI[-0.106, INF]
  
  # EXP2 MOTOR D>C | D+C: t(92.4)=0.501, ptukey=0.871
  # EXP2 MOTOR D>C | C>D: t(89.7)=0.618, ptukey=0.811
  # EXP2 MOTOR D+C | C>D: t(91)=0.104, ptukey=0.994
  
  # Perception
  test_data <- group_data.exp1 |> 
    filter(include_block==TRUE) |> 
    filter(motor_time_type %in% c("motor fixed"))
  hypothesis6.exp1.mdl <- lmer(
    efficiency ~ response_type + (1|subjID),
    data = test_data)
  hypothesis6.exp1.stat <- anova(hypothesis6.exp1.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis6.exp1.eta <- eta_squared(hypothesis6.exp1.stat)
  h6.exp1 <- emmeans(hypothesis6.exp1.mdl, ~  response_type)
  hypothesis6.exp1.pairs <- pairs(h6.exp1)
  filtered_data <- test_data |> 
    filter(response_type %in% c("D→C", "C→D")) |> 
    tidyr::pivot_wider(id_cols = subjID, names_from = response_type, values_from = efficiency)
  hypothesis6.exp1.directional <- t.test(
    x = filtered_data[["D→C"]],
    y = filtered_data[["C→D"]],
    paired = TRUE, alternative = "greater")
  
  # Memory
  test_data <- group_data.exp2 |> 
    filter(include_block==TRUE) |> 
    filter(motor_time_type %in% c("motor fixed"))
  hypothesis6.exp2.mdl <- lmer(
    efficiency ~ response_type + (1|subjID),
    data = test_data)
  hypothesis6.exp2.stat <- anova(hypothesis6.exp2.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis6.exp2.eta <- eta_squared(hypothesis6.exp2.stat)
  h6.exp2 <- emmeans(hypothesis6.exp2.mdl, ~  response_type)
  hypothesis6.exp2.pairs <- pairs(h6.exp2)
  filtered_data <- test_data |> 
    filter(response_type %in% c("D→C", "C→D")) |> 
    tidyr::pivot_wider(id_cols = subjID, names_from = response_type, values_from = efficiency)
  hypothesis6.exp2.directional <- t.test(
    x = filtered_data[["D→C"]],
    y = filtered_data[["C→D"]],
    paired = TRUE, alternative = "greater")
  
  # HYPOTHESIS 7: ----
  
  # DOMAIN:   F(1,114.30)=47.66, p<.001, ηp2=0.29 ***
  # ORDER:    F(2,891.23)=3.00, p=.050, ηp2<.01 .
  # CONTEXT:  F(2,888.39)=0.42, p=.659, ηp2<.01
  # DOMAIN*ORDER: F(2,891.23)=0.003, p=.997, ηp2<.01
  # DOMAIN*CONTEXT: F(2,888.39)=0.23, p=.795, ηp2<.01
  # DOMAIN*ORDER*CONTEXT: F(4,887.82)=0.24, p=.913, ηp2<.01
  
  hypothesis7.mdl <- lmer(
    efficiency ~ domain_type * response_type * motor_time_type + (1|subjID),
    data = g_data)
  hypothesis7.stat <- anova(hypothesis7.mdl, ddf = "Satterthwaite", type = 3)
  hypothesis7.eta <- eta_squared(hypothesis7.stat)
  
  # Following-up the ORDER factor with a targeted model that disregards DOMAIN & CONTEXT
  order.mdl <- lmer(
    efficiency ~ response_type + (1|subjID),
    data = g_data)
  order.stat <- anova(order.mdl, ddf = "Satterthwaite")
  eta_squared(order.stat)
  h7 <- emmeans(order.mdl, ~  response_type)
  pairs(h7)
  
  # ORDER: F(2,907)=3.30, p=.037, ηp2<.01 *
  # D>C | D+C: t(910)=2.164, ptukey=0.078
  # D>C | C>D: t(905)=2.277, ptukey=0.060
  # D+C | C>D: t(906)=0.098, ptukey=0.995
  
}

# TASK PERFORMANCE ----

if (do_accuracy_analysis) {
  cat("Running task performance analysis...\n")
  
  # Experiment 1
  test_data <- group_data.exp1 |>
    filter(include_block == TRUE)
  
  accuracy.bf <- anovaBF(mean_accuracy ~ response_type + motor_time_type + subjID,
                         whichRandom = "subjID",
                         data = test_data,
                         progress = FALSE,
                         rscaleFixed = "wide",
                         rscaleRandom = "nuisance",
                         whichModels = "withmain",
                         iterations = num_iterations)
  accuracy.exp1 <- bayesfactor_inclusion(accuracy.bf, match_models = TRUE)
  
  # Experiment 2
  test_data <- group_data.exp2 |>
    filter(include_block == TRUE)
  
  accuracy.bf <- anovaBF(mean_accuracy ~ response_type + motor_time_type + subjID,
                         whichRandom = "subjID",
                         data = test_data,
                         progress = FALSE,
                         rscaleFixed = "wide",
                         rscaleRandom = "nuisance",
                         whichModels = "withmain",
                         iterations = num_iterations)
  accuracy.exp2 <- bayesfactor_inclusion(accuracy.bf, match_models = TRUE)
  
  test_data <- g_data |> 
    mutate(version = fct_recode(version, perception = "experiment1", memory = "experiment2")) |> 
    filter(include_block == TRUE)
  
  accuracy.bf <- anovaBF(mean_accuracy ~ domain_type + subjID,
                         whichRandom = "subjID",
                         data = test_data,
                         progress = FALSE,
                         rscaleFixed = "wide",
                         rscaleRandom = "nuisance",
                         whichModels = "withmain",
                         iterations = num_iterations*10)
  accuracy.group <- bayesfactor_inclusion(accuracy.bf, match_models = TRUE)
  
} else {
  cat("Skipping task performance analysis...\n")
}

# METACOGNITIVE BIAS ----

if (do_confidence_analysis) {
  cat("Running confidence analysis...\n")
  
  # Experiment 1
  test_data <- group_data.exp1 |>
    filter(include_block == TRUE)
  
  confidence.bf <- anovaBF(mean_confidence ~ response_type + motor_time_type + subjID,
                           whichRandom = "subjID",
                           data = test_data,
                           progress = FALSE,
                           rscaleFixed = "wide",
                           rscaleRandom = "nuisance",
                           whichModels = "withmain",
                           iterations = num_iterations)
  confidence.exp1 <- bayesfactor_inclusion(confidence.bf, match_models = TRUE)
  
  # Experiment 2
  test_data <- group_data.exp2 |>
    filter(include_block == TRUE)
  
  confidence.bf <- anovaBF(mean_confidence ~ response_type + motor_time_type + subjID,
                           whichRandom = "subjID",
                           data = test_data,
                           progress = FALSE,
                           rscaleFixed = "wide",
                           rscaleRandom = "nuisance",
                           whichModels = "withmain",
                           iterations = num_iterations)
  confidence.exp2 <- bayesfactor_inclusion(confidence.bf, match_models = TRUE)
  
  test_data <- g_data |> 
    mutate(version = fct_recode(version, perception = "experiment1", memory = "experiment2")) |> 
    filter(include_block == TRUE)
  
  confidence.bf <- anovaBF(mean_confidence ~ domain_type + subjID,
                           whichRandom = "subjID",
                           data = test_data,
                           progress = FALSE,
                           rscaleFixed = "wide",
                           rscaleRandom = "nuisance",
                           whichModels = "withmain",
                           iterations = num_iterations*10)
  confidence.group <- bayesfactor_inclusion(confidence.bf, match_models = TRUE)
  
} else {
  cat("Skipping confidence analysis...\n")
}

# ORDER BIAS IN PAST STUDIES ----

if (do_confid_database_analysis) {
  
  library(lmerTest)
  library(emmeans)
  
  all_trials <- read.csv("../../data/past studies/all_past_trials.csv")
  
  all_trials$study <- factor(all_trials$study,
                           levels = c(
                             "Wierzchon_2014",
                             "Siedlecka_2016",
                             "Siedlecka_2018_bioRxiv",
                             "Siedlecka_2019_Exp1",
                             "Siedlecka_2019_Exp2"
                           ),
                           labels = c(
                             "2014_W",
                             "2016_S",
                             "2018_S",
                             "2019_S exp1",
                             "2019_S exp2"
                           ))
  all_trials$subject_idx <- factor(all_trials$subject_idx)
  all_trials$condition <- factor(all_trials$condition,
                               levels = c("D_then_C","C_then_D"),
                               labels = c("D→C","C→D"))
  
  # Ensure numeric
  all_trials$confidence <- as.numeric(all_trials$confidence)
  
  full_mdl <- glmer(data = all_trials,
                    formula = accuracy ~ confidence * condition + 
                      (confidence|subject_idx),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  
  null_mdl <- glmer(data = all_trials,
                    formula = accuracy ~ confidence + condition + 
                      (confidence|subject_idx),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  anova(full_mdl, null_mdl, test = "Chisq")
  # ALL STUDIES: chisq(1)=17.813, p<.001 #Full model wins
  
  these_trials <- all_trials |> 
    filter(study != "2019_S exp1")
  
  full_mdl <- glmer(data = these_trials,
                    formula = accuracy ~ confidence * condition + 
                      (confidence|subject_idx),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  
  null_mdl <- glmer(data = these_trials,
                    formula = accuracy ~ confidence + condition + 
                      (confidence|subject_idx),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  anova(full_mdl, null_mdl, test = "Chisq")
  # EXCLUDE 2019_S exp1: chisq(1)=6.097, p=.014 #Full model wins
  
  # Load relevant data from the confidence database
  all_data <- read.csv("../../data/past studies/all_past_data.csv")
  
  all_data$study <- factor(all_data$study,
                           levels = c(
                             "Wierzchon_2014",
                             "Siedlecka_2016",
                             "Siedlecka_2018_bioRxiv",
                             "Siedlecka_2019_Exp1",
                             "Siedlecka_2019_Exp2"
                           ),
                           labels = c(
                             "2014_W",
                             "2016_S",
                             "2018_S",
                             "2019_S exp1",
                             "2019_S exp2"
                           ))
  all_data$subject_idx <- factor(all_data$subject_idx)
  all_data$condition <- factor(all_data$condition,
                                levels = c("D_then_C","C_then_D"),
                                labels = c("D→C","C→D"))
  all_data$valid_subj <- factor(all_data$valid_subj)
  
  # "2014_W", "2016_S", "2018_S", "2019_S exp1", "2019_S exp2"
  past_data <- all_data |> 
    filter(study == enter_study_label_here)
  
  D_C <- past_data |> filter(condition == "D→C")
  C_D <- past_data |> filter(condition == "C→D")
  
  # Compute m-ratio
  t.test(D_C$m_ratio,C_D$m_ratio, paired = FALSE)
  ttestBF(D_C$m_ratio,C_D$m_ratio, paired = FALSE, rscale = "medium")
  
  # 2014_W: t(18.171)=-0.949, p=.355, BF10=0.470
  # 2016_S: t(60.736)=0.883, p=.381, BF10=0.355
  # 2018_S: t(57.839)=0.948, p=.347, BF10=0.315
  # 2019_S exp1: t(30.684)=10.514, p<.001, BF10=4691839889
  # 2019_S exp2: t(52.069)=0.470, p=.641, BF10=0.298
  
  past_data <- all_data |> 
    filter(study != "2019_S exp1")
  
  exclusion_num <- nlevels(past_data$subject_idx) - length(unique(past_data$subject_idx))
  
  D_C <- past_data |> filter(condition == "D→C")
  C_D <- past_data |> filter(condition == "C→D")
  
  # Compute m-ratio
  the_test <- t.test(D_C$m_ratio,C_D$m_ratio, paired = FALSE)
  the_bf <- ttestBF(D_C$m_ratio,C_D$m_ratio, paired = FALSE, rscale = "medium")
  
  }

# HMETA-D' ANALYSIS ----

if (do_hmetad_analysis) {
  
  setwd("./metad/")
  
  test_data <- trial_data.exp1 |> 
    filter(include_trial == 1)
  
  # Build counts ----
  # 3x3 repeated measures design:
  # Context (Baseline, Time, Motor) x Order (D→C, Simultaneous, C→D)
  
  # Context: (B,B,B,T,T,T,M,M,M) - Order: (D→,S,C→, D→,S,C→, D→,S,C→)
  nsubj <- nlevels(test_data$subjID)
  nRatings <- length(unique(test_data$confidence))
  ncontexts <- nlevels(test_data$motor_time_type)
  nresps <- nlevels(test_data$response_type)
  
  # Build blank lists within lists (required to input counts)
  nR_S1 <- rep(list(vector(mode = "list", length = ncontexts*nresps)),nsubj)
  nR_S2 <- rep(list(vector(mode = "list", length = ncontexts*nresps)),nsubj)
  
  for (n in 1:(nsubj)) {
    for (c in 1:(ncontexts)) {
      for (o in 1:(nresps)) {
        
        # Select appropriate conditions
        this_subj <- levels(test_data$subjID)[n]
        this_cont <- levels(test_data$motor_time_type)[c]
        this_resp <- levels(test_data$response_type)[o]
        
        # Subset to appropriate data
        this_data <- test_data |>  
          filter(subjID==this_subj & 
                   motor_time_type==this_cont &
                   response_type==this_resp)
        
        this_data$correct_side <- factor(this_data$correct_side,
                                         levels = c("left","right"),
                                         labels = c(0,1))
        this_data$response_side <- factor(this_data$response_side,
                                          levels = c("left","right"),
                                          labels = c(0,1))
        
        # stimID =    list(0, 1, 0, 0, 1, 1, 1, 1)
        stimID <- as.list(as.numeric(this_data$correct_side)-1)
        
        # response =  list(0, 1, 1, 1, 0, 0, 1, 1)
        response <- as.list(as.numeric(this_data$response_side)-1)
        
        # rating =    list(1, 2, 3, 4, 4, 3, 2, 1)
        rating <- as.list(this_data$confidence)
        
        source("trials2counts.R")
        newlist <- trials2counts(stimID,response,rating,nRatings)
        
        nR_S1[[n]][((c-1)*3)+o] <- newlist[1]
        nR_S2[[n]][((c-1)*3)+o] <- newlist[2]
        
      }
    }
  }
  
  # Simply load the outputs of a previous analysis
  # output <- readRDS("output_results_exp1.rds")
  
  # Run a fresh computational analysis (takes a while)
  source("fit_metad_2wayANOVA.R")
  output <- metad_2wayANOVA(nR_S1 = nR_S1, nR_S2 = nR_S2)
  saveRDS(output, file = "output_results_exp1.rds") # Save to an RDS
  
  setwd("../")
  
  # Mean values 
  Value <- summary(output)
  stat <- data.frame(mean = Value[["statistics"]][, "Mean"])
  stat %<>%
    rownames_to_column(var = "name")
  
  # Rhat 
  Value <- gelman.diag(output, confidence = 0.95)
  Rhat <- data.frame(conv = Value$psrf)
  
  # HDI 
  HDI <- data.frame(HPDinterval(output, prob = 0.95))
  HDI %<>%
    rownames_to_column(var = "name")
  
  # All values in the same dataframe
  Fit <- stat %>%
    cbind(lower = HDI$lower,
          upper = HDI$upper)
  
  fixed_stats <- stat
  
  fixed_stats$isMratio <- 0
  
  # Build useful data.frame
  for (n in 1:(nsubj)) {
    for (c in 1:(ncontexts)) {
      for (o in 1:(nresps)) {
        
        this_subj <- levels(test_data$subjID)[n]
        this_cont <- levels(test_data$motor_time_type)[c]
        this_resp <- levels(test_data$response_type)[o]
        
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"subjID"] <- this_subj
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"motor_time_type"] <- this_cont
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"response_type"] <- this_resp
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"isMratio"] <- 1
      }
    }
  }
  
  fixed_stats <- fixed_stats %>% filter(isMratio==1)
  
  fixed_stats$motor_time_type <- factor(fixed_stats$motor_time_type,
                                        levels = motor_time_labels)
  
  fixed_stats$response_type <- factor(fixed_stats$response_type,
                                      levels = response_labels)
  
  # mcmc values in df for plot posterior distributions
  mcmc.sample <- ggs(output)
  
  # Compute mean for plot #1
  the_mean <- stat$mean[stat$name == "muBd_Factor1"]
  mean_label <- sprintf("mu==%1.2f",the_mean)
  
  # Plot posterior distribution for rho value
  PlotCondition1 <- mcmc.sample %>%
    filter(Parameter == "muBd_Factor1") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "slategray2", alpha = 0.7, bins = 100) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid", color='grey') +
    geom_segment(aes(x = HDI$lower[HDI$name == "muBd_Factor1"], 
                     y = 50, xend = HDI$upper[HDI$name == "muBd_Factor1"], yend = 50), 
                 colour = "black", size = 2) +
    annotate("text",label=mean_label, y=500, x = the_mean, parse=TRUE) +
    annotate("text",label="95% HDI", y=150, x = the_mean, colour="black") +
    ylim(0,2000) +
    xlim(-0.25,0.25) +
    ylab("sample count") +
    ggtitle("motor & time context") +
    xlab(expression(paste(log,(beta)))) +
    theme_julian() + theme(axis.title = element_text(face = "plain"))
  
  # Compute mean for plot #2
  the_mean <- stat$mean[stat$name == "muBd_Factor2"]
  mean_label <- sprintf("mu==%1.2f",the_mean)
  
  PlotCondition2 <- mcmc.sample %>%
    filter(Parameter == "muBd_Factor2") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "tan1", alpha = 0.3, bins = 100) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid", color='grey') +
    geom_segment(aes(x = HDI$lower[HDI$name == "muBd_Factor2"], 
                     y = 50, xend = HDI$upper[HDI$name == "muBd_Factor2"], yend = 50), 
                 colour = "black", size = 2) +
    annotate("text",label=mean_label, y=500, x = the_mean, parse=TRUE) +
    annotate("text",label="95% HDI", y=150, x = the_mean, colour="black") +
    ylim(0,2000) +
    xlim(-0.25,0.25) +
    ylab("sample count") +
    ggtitle("report order") +
    xlab(expression(paste(log,(beta)))) +
    theme_julian() + theme(axis.title = element_text(face = "plain"))
  
  library(cowplot)
  plot_grid(PlotCondition1,PlotCondition2 + ylab(""), nrow = 1)
  
  ggsave(paste0(fig_location,"sampled_Hmetad_exp1.png"),width = 8, height = 6)
  
    # Create plot
  sumstats <- Rmisc::summarySE(data = fixed_stats,
                               measurevar = "mean",
                               groupvars = c("motor_time_type","response_type"))
  library(ggbeeswarm)
  mestimates <- ggplot(data = fixed_stats, aes(x = response_type, y = mean, colour = response_type)) +
    facet_grid(. ~ motor_time_type) +
    geom_quasirandom(shape = 20, size = 0.4, alpha = 0.3,
                     dodge.width = 0.3,
                     show.legend = FALSE) +
    geom_line(data = sumstats,
              aes(group = motor_time_type),
              linewidth = 1,
              colour = "black") +
    geom_errorbar(data = sumstats,
                  linewidth = 2,
                  aes(ymin = mean-se,ymax=mean+se),
                  width=0,
                  position = position_dodge(width=0.3),
                  show.legend = FALSE) +
    geom_point(data = sumstats, size = 1, shape = 21,
               fill = "white", stroke = 0.5,
               position = position_dodge(width=0.3),
               show.legend = FALSE) +
    # scale_fill_manual(values = context_col) +
    # scale_colour_manual(values = context_col) +
    ylim(0,1.5) +
    xlab(NULL) + ylab(expression("meta-d'"/"d'")) + labs(colour=NULL) +
    ggtitle("Hmeta-d metacognitive efficiency") +
    theme_julian()
    # theme(legend.position = c(0.1,0.85),
    #       legend.key.size = unit(4,"mm"),
    #       legend.background = element_blank())
  mestimates
  
  ggsave(plot = mestimates, paste0(fig_location,"mratio_Hmetad_exp1.png"),width=6, height=4)
  
  full_plot <- plot_grid(mestimates,PlotCondition1,PlotCondition2 + ylab(""),
                         nrow=1,rel_widths = c(6/12,3/12,3/12),
                         labels = c("a)","b)",""))
  ggsave(plot = full_plot,paste0(fig_location,"full_Hmetad_exp1.png"),width=10, height=4)
  
  # EXP1 RESULTS
  
  rhat_range <- range(Rhat$conv.Point.est)
  rhat_results_exp1 <- sprintf("Rhat: [%0.3g, %0.3g]",
                          rhat_range[1],
                          rhat_range[2])
  
  message(rhat_results_exp1) 
  # Rhat range: [0.99998, 1.30495]
  
  exp1_results <- sprintf("Context: mu=%0.3f, %%95HDI[%0.3f, %0.3f], Order: mu=%0.3f, %%95HDI[%0.3f, %0.3f]",
                          stat$mean[stat$name == "muBd_Factor1"],
                          HDI$lower[HDI$name == "muBd_Factor1"],
                          HDI$upper[HDI$name == "muBd_Factor1"],
                          stat$mean[stat$name == "muBd_Factor2"],
                          HDI$lower[HDI$name == "muBd_Factor2"],
                          HDI$upper[HDI$name == "muBd_Factor2"])
  message(exp1_results)
  # Context: mu=-0.007, %95HDI[-0.082, 0.073], Order: mu=0.014, %95HDI[-0.067, 0.084]
                          
  
  setwd("./metad/")
  
  test_data <- trial_data.exp2 |> 
    filter(include_trial == 1)
  
  # Build counts ----
  # 3x3 repeated measures design:
  # Context (Baseline, Time, Motor) x Order (D→C, Simultaneous, C→D)
  
  # Context: (B,B,B,T,T,T,M,M,M) - Order: (D→,S,C→, D→,S,C→, D→,S,C→)
  nsubj <- nlevels(test_data$subjID)
  nRatings <- length(unique(test_data$confidence))
  ncontexts <- nlevels(test_data$motor_time_type)
  nresps <- nlevels(test_data$response_type)
  
  # Build blank lists within lists (required to input counts)
  nR_S1 <- rep(list(vector(mode = "list", length = ncontexts*nresps)),nsubj)
  nR_S2 <- rep(list(vector(mode = "list", length = ncontexts*nresps)),nsubj)
  
  for (n in 1:(nsubj)) {
    for (c in 1:(ncontexts)) {
      for (o in 1:(nresps)) {
        
        # Select appropriate conditions
        this_subj <- levels(test_data$subjID)[n]
        this_cont <- levels(test_data$motor_time_type)[c]
        this_resp <- levels(test_data$response_type)[o]
        
        # Subset to appropriate data
        this_data <- test_data |>  
          filter(subjID==this_subj & 
                   motor_time_type==this_cont &
                   response_type==this_resp)
        
        this_data$correct_side <- factor(this_data$correct_side,
                                         levels = c("left","right"),
                                         labels = c(0,1))
        this_data$response_side <- factor(this_data$response_side,
                                          levels = c("left","right"),
                                          labels = c(0,1))
        
        # stimID =    list(0, 1, 0, 0, 1, 1, 1, 1)
        stimID <- as.list(as.numeric(this_data$correct_side)-1)
        
        # response =  list(0, 1, 1, 1, 0, 0, 1, 1)
        response <- as.list(as.numeric(this_data$response_side)-1)
        
        # rating =    list(1, 2, 3, 4, 4, 3, 2, 1)
        rating <- as.list(this_data$confidence)
        
        source("trials2counts.R")
        newlist <- trials2counts(stimID,response,rating,nRatings)
        
        nR_S1[[n]][((c-1)*3)+o] <- newlist[1]
        nR_S2[[n]][((c-1)*3)+o] <- newlist[2]
        
      }
    }
  }
  
  # Simply load the outputs of a previous analysis
  # output <- readRDS("output_results_exp2.rds")
  
  # Run a fresh computational analysis (takes a while)
  source("fit_metad_2wayANOVA.R")
  output <- metad_2wayANOVA(nR_S1 = nR_S1, nR_S2 = nR_S2)
  saveRDS(output, file = "output_results_exp2.rds") # Save to an RDS
  
  setwd("../")
  
  # Mean values 
  Value <- summary(output)
  stat <- data.frame(mean = Value[["statistics"]][, "Mean"])
  stat %<>%
    rownames_to_column(var = "name")
  
  # Rhat 
  Value <- gelman.diag(output, confidence = 0.95)
  Rhat <- data.frame(conv = Value$psrf)
  
  # HDI 
  HDI <- data.frame(HPDinterval(output, prob = 0.95))
  HDI %<>%
    rownames_to_column(var = "name")
  
  # All values in the same dataframe
  Fit <- stat %>%
    cbind(lower = HDI$lower,
          upper = HDI$upper)
  
  fixed_stats <- stat
  
  fixed_stats$isMratio <- 0
  
  # Build useful data.frame
  for (n in 1:(nsubj)) {
    for (c in 1:(ncontexts)) {
      for (o in 1:(nresps)) {
        
        this_subj <- levels(test_data$subjID)[n]
        this_cont <- levels(test_data$motor_time_type)[c]
        this_resp <- levels(test_data$response_type)[o]
        
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"subjID"] <- this_subj
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"motor_time_type"] <- this_cont
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"response_type"] <- this_resp
        fixed_stats[fixed_stats$name==paste0("Mratio[",n,",",((c-1)*3)+o,"]"),"isMratio"] <- 1
      }
    }
  }
  
  fixed_stats <- fixed_stats %>% filter(isMratio==1)
  
  fixed_stats$motor_time_type <- factor(fixed_stats$motor_time_type,
                                        levels = motor_time_labels)
  
  fixed_stats$response_type <- factor(fixed_stats$response_type,
                                      levels = response_labels)
  
  # mcmc values in df for plot posterior distributions
  mcmc.sample <- ggs(output)
  
  # Compute mean for plot #1
  the_mean <- stat$mean[stat$name == "muBd_Factor1"]
  mean_label <- sprintf("mu==%1.2f",the_mean)
  
  # Plot posterior distribution for rho value
  PlotCondition1 <- mcmc.sample %>%
    filter(Parameter == "muBd_Factor1") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "slategray2", alpha = 0.7, bins = 100) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid", color='grey') +
    geom_segment(aes(x = HDI$lower[HDI$name == "muBd_Factor1"], 
                     y = 50, xend = HDI$upper[HDI$name == "muBd_Factor1"], yend = 50), 
                 colour = "black", size = 2) +
    annotate("text",label=mean_label, y=500, x = the_mean, parse=TRUE) +
    annotate("text",label="95% HDI", y=150, x = the_mean, colour="black") +
    ylim(0,2500) +
    xlim(-0.25,0.25) +
    ylab("sample count") +
    ggtitle("motor & time context") +
    xlab(expression(paste(log,(beta)))) +
    theme_julian() + theme(axis.title = element_text(face = "plain"))
  
  # Compute mean for plot #2
  the_mean <- stat$mean[stat$name == "muBd_Factor2"]
  mean_label <- sprintf("mu==%1.2f",the_mean)
  
  PlotCondition2 <- mcmc.sample %>%
    filter(Parameter == "muBd_Factor2") %>% 
    ggplot(aes(value)) +
    geom_histogram(fill = "tan1", alpha = 0.3, bins = 100) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid", color='grey') +
    geom_segment(aes(x = HDI$lower[HDI$name == "muBd_Factor2"], 
                     y = 50, xend = HDI$upper[HDI$name == "muBd_Factor2"], yend = 50), 
                 colour = "black", size = 2) +
    annotate("text",label=mean_label, y=500, x = the_mean, parse=TRUE) +
    annotate("text",label="95% HDI", y=150, x = the_mean, colour="black") +
    ylim(0,2500) +
    xlim(-0.25,0.25) +
    ylab("sample count") +
    ggtitle("report order") +
    xlab(expression(paste(log,(beta)))) +
    theme_julian() + theme(axis.title = element_text(face = "plain"))
  
  plot_grid(PlotCondition1,PlotCondition2 + ylab(""), nrow = 1)
  
  ggsave(paste0(fig_location,"sampled_Hmetad_exp2.png"),width = 8, height = 6)
  
  # Create plot
  sumstats <- Rmisc::summarySE(data = fixed_stats,
                               measurevar = "mean",
                               groupvars = c("motor_time_type","response_type"))
  
  mestimates <- ggplot(data = fixed_stats, aes(x = response_type, y = mean, colour = response_type)) +
    facet_grid(. ~ motor_time_type) +
    geom_quasirandom(shape = 20, size = 0.4, alpha = 0.3,
                     dodge.width = 0.3,
                     show.legend = FALSE) +
    geom_line(data = sumstats,
              aes(group = motor_time_type),
              linewidth = 1,
              colour = "black") +
    geom_errorbar(data = sumstats,
                  linewidth = 2,
                  aes(ymin = mean-se,ymax=mean+se),
                  width=0,
                  position = position_dodge(width=0.3),
                  show.legend = FALSE) +
    geom_point(data = sumstats, size = 1, shape = 21,
               fill = "white", stroke = 0.5,
               position = position_dodge(width=0.3),
               show.legend = FALSE) +
    # scale_fill_manual(values = context_col) +
    # scale_colour_manual(values = context_col) +
    ylim(0,1.5) +
    xlab(NULL) + ylab(expression("meta-d'"/"d'")) + labs(colour=NULL) +
    ggtitle("Hmeta-d metacognitive efficiency") +
    theme_julian()
  # theme(legend.position = c(0.1,0.85),
  #       legend.key.size = unit(4,"mm"),
  #       legend.background = element_blank())
  mestimates
  
  ggsave(plot = mestimates, paste0(fig_location,"mratio_Hmetad_exp2.png"),width=6, height=4)
  
  full_plot <- plot_grid(mestimates,PlotCondition1,PlotCondition2 + ylab(""),
                         nrow=1,rel_widths = c(6/12,3/12,3/12),
                         labels = c("a)","b)",""))
  ggsave(plot = full_plot,paste0(fig_location,"full_Hmetad_exp2.png"),width=10, height=4)
  
  # EXP2 RESULTS
  
  rhat_range <- range(Rhat$conv.Point.est)
  rhat_results_exp2 <- sprintf("Rhat: [%0.3g, %0.3g]",
                          rhat_range[1],
                          rhat_range[2])
  
  message(rhat_results_exp2) 
  # Rhat range: [0.99998, 1.30923]
  
  exp2_results <- sprintf("Context: mu=%0.3f, %%95HDI[%0.3f, %0.3f], Order: mu=%0.3f, %%95HDI[%0.3f, %0.3f]",
                          stat$mean[stat$name == "muBd_Factor1"],
                          HDI$lower[HDI$name == "muBd_Factor1"],
                          HDI$upper[HDI$name == "muBd_Factor1"],
                          stat$mean[stat$name == "muBd_Factor2"],
                          HDI$lower[HDI$name == "muBd_Factor2"],
                          HDI$upper[HDI$name == "muBd_Factor2"])
  message(exp2_results)
  # Context: mu=0.013, %95HDI[-0.042, 0.070], Order: mu=0.032, %95HDI[-0.027, 0.094]
  
}

# REGRESSION ANALYSIS OF METACOGNITIVE ACCURACY ----

if (do_regression_analysis) {
  
  # Relevant packages
  library(lmerTest)
  library(emmeans)
  
  # EXPERIMENT 1: PERCEPTION ----
  
  test_data <- trial_data.exp1 |> 
    filter(include_trial == 1)
  
  test_data$confidence <- as.numeric(test_data$confidence)
  
  ## Build models ----
  full_mdl <- glmer(data = test_data,
                    formula = accuracy ~ confidence * response_type * motor_time_type + 
                      (confidence|subjID),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  
  add_mdl <- glmer(data = test_data,
                   formula = accuracy ~ confidence * response_type + confidence * motor_time_type +
                     (confidence|subjID),
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                          optCtrl = list(method = "nlminb",
                                                         starttests = FALSE, kkt = FALSE)))
  rm_context <- glmer(data = test_data,
                      formula = accuracy ~ confidence * response_type + motor_time_type +
                        motor_time_type:response_type +
                        (confidence|subjID),
                      family = binomial(link = "logit"),
                      control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                             optCtrl = list(method = "nlminb",
                                                            starttests = FALSE, kkt = FALSE)))
  rm_response <- glmer(data = test_data,
                       formula = accuracy ~ confidence * motor_time_type + response_type +
                         motor_time_type:response_type +
                         (confidence|subjID),
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                              optCtrl = list(method = "nlminb",
                                                             starttests = FALSE, kkt = FALSE)))
  
  ## Model tests ----
  
  anova(full_mdl, add_mdl, test = "Chisq")
  # EXP1: chisq(8)=16.595, p=.03461 #Full model wins
  
  anova(full_mdl, rm_response, test = "Chisq")
  # EXP1: chisq(6)=19.872, p=.002918 #Full model wins
  
  anova(full_mdl, rm_context, test = "Chisq")
  # EXP1: chisq(6)=17.043, p=.009125 #Full model wins
  
  ## Additional analysis ----
  
  # Determine prediction reference grid
  mdl.refgrid <- ref_grid(full_mdl, at = list(confidence = c(1,2,3,4)), type = "response")
  
  # ref_grid for glmer assumes a fixed SD of 1
  # we can form a more precise estimate of the SD by examining the actual SD
  update_sd <- lme4::VarCorr(full_mdl) # true SD = 0.24885
  update_sd <- as.data.frame(update_sd) # convert to data.frame
  
  # Update SD using true SD and pythagorean theorem
  mdl.refgrid <- update(mdl.refgrid,
                        sigma = sqrt(mean(abs(update_sd$sdcor),na.rm = TRUE)^2))
  
  mdl.predictions <- predict(mdl.refgrid, interval = "prediction")
  
  ## Plot ----
  
  # Convert accuracy to numeric data for plotting
  test_data$accuracy <- as.numeric(test_data$accuracy)
  
  # Compute mean accuracy and counts of each confidence rating for each subject
  subjstats <- test_data  |>  
    group_by(subjID,motor_time_type,response_type,confidence) |> 
    dplyr::summarise(subjAcc = mean(accuracy, na.rm = TRUE),
                     count = n(),
                     .groups = "keep") |> 
    group_by(subjID, motor_time_type, response_type) |> 
    mutate(full_count = sum(count)) |> 
    group_by(subjID, motor_time_type, response_type, confidence) |> 
    mutate(freq = count / full_count)
  
  # Compute mean accuracy and frequency of each confidence rating for each condition
  groupstats <- subjstats |> 
    group_by(motor_time_type,response_type,confidence) |> 
    dplyr::summarise(mAcc = mean(subjAcc),
                     frequency = sum(count)/(100*nlevels(subjstats$subjID)),
                     mFreq = mean(freq),
                     se = sd(subjAcc)/sqrt(nlevels(subjstats$subjID)),
                     .groups = "drop")
  
  # Build polygons to represent confidence intervals
  positions <- mdl.predictions |> 
    group_by(motor_time_type,response_type) |> 
    reframe(x1 = confidence,
              x2 = confidence,
              y1 = prediction+SE,
              y2 = prediction-SE)
  
  polygrid <- data.frame(
    motor_time_type = c(positions$motor_time_type,rev(positions$motor_time_type)),
    response_type = c(positions$response_type, rev(positions$response_type)),
    x = c(positions$x1,rev(positions$x2)),
    y = c(positions$y1,rev(positions$y2))
  )
  
  # Construct plot
  ggplot(data = groupstats, aes(x = confidence, y = mAcc)) +
    facet_grid(motor_time_type~response_type) +
    geom_polygon(data = polygrid,
                 aes(x = x, y = y,
                     group = interaction(motor_time_type,response_type)),
                 fill = "grey92",inherit.aes = FALSE) +
    geom_line(data = mdl.predictions,
              aes(x = confidence, y = prediction),
              size = 1, alpha = 1, colour = "grey") +
    geom_point(aes(size = frequency), shape=21, fill = "black", stroke=0) +
    scale_size(range = c(2,6)) +
    guides(fill = NULL, group = NULL) +
    ylab("mean accuracy") +
    labs(subtitle = "Behaviour + Full Model predictions") +
    theme_julian() +
    theme(legend.key.size = unit(5,"mm"))
  
  ggsave(paste0(fig_location,"EXPLORATORY_logistic_plots_exp1.png"),width=6, height=6)
  
  ## Follow-ups ----
  
  # Subset data to each context
  data_baseline <- test_data |> filter(motor_time_type == "baseline")
  data_time <- test_data |> filter(motor_time_type == "time limited")
  data_motor <- test_data |> filter(motor_time_type == "motor fixed")
  
  # Generate model for each context
  
  exp1.baseline.mdl <- 
    glmer(data = data_baseline,
          formula = accuracy ~ confidence * response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  exp1.baseline.null <- 
    glmer(data = data_baseline,
          formula = accuracy ~ confidence + response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  anova(exp1.baseline.mdl, exp1.baseline.null, test = "Chisq")
  # EXP1 BASELINE: chisq(2)=13.215, p=.00135 #Order model wins
  
  # Perform pairwise contrasts:
  exp1.model.pairs <- emmeans(exp1.baseline.mdl, ~ response_type | confidence, type = "response")
  pairs(exp1.model.pairs)
  
  # Tukey adjusted pairwise contrasts
  # D→C vs. D+C: z=2.360, p=.0479 *
  # D→C vs. C→D: z=0.908, p=.6353
  # D+C vs. C→D: z=-1.456, p=.3125
  
  exp1.time.mdl <- 
    glmer(data = data_time,
          formula = accuracy ~ confidence * response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  exp1.time.null <- 
    glmer(data = data_time,
          formula = accuracy ~ confidence + response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  anova(exp1.time.mdl, exp1.time.null, test = "Chisq")
  # EXP1 TIME LIMITED: chisq(2)=3.1149, p=.2107 #Null model wins
  
  exp1.motor.mdl <- 
    glmer(data = data_motor,
          formula = accuracy ~ confidence * response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  exp1.motor.null <- 
    glmer(data = data_motor,
          formula = accuracy ~ confidence + response_type + 
            (confidence|subjID),
          family = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "optimx", 
            calc.derivs = FALSE,
            optCtrl = list(method = "nlminb",starttests = FALSE, kkt = FALSE)))
  
  anova(exp1.motor.mdl, exp1.motor.null, test = "Chisq")
  # EXP1 MOTOR FIXED: chisq(2)=1.7333, p=.4204 #Null model wins
  
  # EXPERIMENT 2: MEMORY ----
  
  test_data <- trial_data.exp2 |> 
    filter(include_trial == 1)
  
  test_data$confidence <- as.numeric(test_data$confidence)
  
  ## Build models ----
  full_mdl <- glmer(data = test_data,
                    formula = accuracy ~ confidence * response_type * motor_time_type + 
                      (confidence|subjID),
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb",
                                                          starttests = FALSE, kkt = FALSE)))
  
  add_mdl <- glmer(data = test_data,
                   formula = accuracy ~ confidence * response_type + confidence * motor_time_type +
                     (confidence|subjID),
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                          optCtrl = list(method = "nlminb",
                                                         starttests = FALSE, kkt = FALSE)))
  rm_context <- glmer(data = test_data,
                      formula = accuracy ~ confidence * response_type + motor_time_type +
                        (confidence|subjID),
                      family = binomial(link = "logit"),
                      control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                             optCtrl = list(method = "nlminb",
                                                            starttests = FALSE, kkt = FALSE)))
  rm_response <- glmer(data = test_data,
                       formula = accuracy ~ confidence * motor_time_type + response_type +
                         (confidence|subjID),
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                              optCtrl = list(method = "nlminb",
                                                             starttests = FALSE, kkt = FALSE)))
  
  ## Model tests ----
  
  anova(full_mdl, add_mdl, test = "Chisq")
  # EXP2: chisq(8)=9.7524, p=0.2828 #Additive wins
  
  # Against (winning) additive model
  
  anova(add_mdl, rm_response, test = "Chisq")
  # EXP2: chisq(2)=10.349, p=.005659 #Additive model wins
  
  anova(add_mdl, rm_context, test = "Chisq")
  # EXP2: chisq(2)=4.0445, p=.1324 #Context-free model wins
  
  ## Additional analysis ----
  
  # Determine prediction reference grid
  mdl.refgrid <- ref_grid(rm_context, at = list(confidence = c(1,2,3,4)), type = "response")
  
  # ref_grid for glmer assumes a fixed SD of 1
  # we can form a more precise estimate of the SD by examining the actual SD
  update_sd <- lme4::VarCorr(rm_context) # true SD = 0.20182
  update_sd <- as.data.frame(update_sd) # convert to data.frame
  
  # Update SD using true SD and pythagorean theorem
  mdl.refgrid <- update(mdl.refgrid,
                        sigma = sqrt(mean(abs(update_sd$sdcor),na.rm = TRUE)^2))
  
  mdl.predictions <- predict(mdl.refgrid, interval = "prediction")
  
  ## Plot ----
  
  # Convert accuracy to numeric data for plotting
  test_data$accuracy <- as.numeric(test_data$accuracy)
  
  # Compute mean accuracy and counts of each confidence rating for each subject
  subjstats <- test_data  |>  
    group_by(subjID,motor_time_type,response_type,confidence) |> 
    dplyr::summarise(subjAcc = mean(accuracy, na.rm = TRUE),
                     count = n(),
                     .groups = "keep") |> 
    group_by(subjID, motor_time_type, response_type) |> 
    mutate(full_count = sum(count)) |> 
    group_by(subjID, motor_time_type, response_type, confidence) |> 
    mutate(freq = count / full_count)
  
  # Compute mean accuracy and frequency of each confidence rating for each condition
  groupstats <- subjstats |> 
    group_by(motor_time_type,response_type,confidence) |> 
    dplyr::summarise(mAcc = mean(subjAcc),
                     frequency = sum(count)/(100*nlevels(subjstats$subjID)),
                     mFreq = mean(freq),
                     se = sd(subjAcc)/sqrt(nlevels(subjstats$subjID)),
                     .groups = "drop")
  
  # Build polygons to represent confidence intervals
  positions <- mdl.predictions |> 
    group_by(motor_time_type,response_type) |> 
    reframe(x1 = confidence,
            x2 = confidence,
            y1 = prediction+SE,
            y2 = prediction-SE)
  
  polygrid <- data.frame(
    motor_time_type = c(positions$motor_time_type,rev(positions$motor_time_type)),
    response_type = c(positions$response_type, rev(positions$response_type)),
    x = c(positions$x1,rev(positions$x2)),
    y = c(positions$y1,rev(positions$y2))
  )
  
  # Construct plot
  ggplot(data = groupstats, aes(x = confidence, y = mAcc)) +
    facet_grid(motor_time_type~response_type) +
    geom_polygon(data = polygrid,
                 aes(x = x, y = y,
                     group = interaction(motor_time_type,response_type)),
                 fill = "grey92",inherit.aes = FALSE) +
    geom_line(data = mdl.predictions,
              aes(x = confidence, y = prediction),
              size = 1, alpha = 1, colour = "grey") +
    geom_point(aes(size = frequency), shape=21, fill = "black", stroke=0) +
    scale_size(range = c(2,6)) +
    guides(fill = NULL, group = NULL) +
    ylab("mean accuracy") +
    labs(subtitle = "Behaviour + Order Model predictions") +
    theme_julian() +
    theme(legend.key.size = unit(5,"mm"))
  
  ggsave(paste0(fig_location,"EXPLORATORY_logistic_plots_exp2.png"),width=6, height=6)
  
  ## Follow-ups ----
  
  # Perform pairwise contrasts:
  exp2.model.pairs <- emmeans(rm_context, ~ response_type | confidence, type = "response")
  pairs(exp2.model.pairs)
  
  # Tukey adjusted pairwise contrasts
  # D→C vs. D+C: z=3.302, p=.0028 **
  # D→C vs. C→D: z=0.796, p=.7053
  # D+C vs. C→D: z=-2.515, p=.0319 *
  
}

# REPORT RESULTS ----

## Task performance
if (do_accuracy_analysis) {
  message("\n\nPerformance in Experiment 1") 
  print(accuracy.exp1)
  message("\n\nPerformance in Experiment 2") 
  print(accuracy.exp2)
  message("\n\nPerformance between experiments") 
  print(accuracy.group)
}

## Metacognitive bias
if (do_confidence_analysis) {
  message("\n\nConfidence in Experiment 1") 
  print(confidence.exp1)
  message("\n\nConfidence in Experiment 2") 
  print(confidence.exp2)
  message("\n\nConfidence between experiments") 
  print(confidence.group)
}


