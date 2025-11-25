# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Desender_2022_Exp1.csv",
  "data_Desender_2022_Exp1S.csv",
  "data_Desender_2022_Exp2B.csv"
  # "data_Desender_2021_Cognit.csv"
  # "data_Rouault_2018_Expt1.csv",
  # "data_Reyes_2015.csv",
  # "data_Reyes_unpub.csv",
  # "data_Schmidt_2019_memory.csv",
  # "data_Schmidt_2019_perception.csv",
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Function to load and clean a dataset
load_and_clean_dataset <- function(filename) {

  # filename = datasets_wrong_right[4]

  filepath <- file.path(dataPath, filename)

  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filename))
    return(NULL)
  }

  # Read data
  data <- read.csv(filepath, fileEncoding = "UTF-8-BOM")

  # Add dataset identifier
  data$Dataset <- gsub(".csv", "", filename)

  # Ensure required columns exist
  required_cols <- c("Subj_idx", "Stimulus", "Response", "Confidence")
  if (!all(required_cols %in% names(data))) {
    warning(paste("Missing required columns in:", filename))
    return(NULL)
  }

  # Add Correct column if not present
  if (!"Correct" %in% names(data)) {
    data$Correct <- as.numeric(data$Stimulus == data$Response)
  }

  # Clean data
  data <- data %>%
    filter(!is.na(Confidence), !is.na(Correct)) %>%
    mutate(
      Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"),
      Subj_idx = as.factor(Subj_idx)
    )

  return(data)
}

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec","Response")) %>%
  group_by(Condition, Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")


all_data %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec")) %>%
  group_by(Condition, Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")





dataset = 3

all_data %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  # mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Difficulty, -Difficulty)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec","Response")) %>%
  group_by(Condition, Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")


all_data %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  # mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Difficulty, -Difficulty)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec")) %>%
  group_by(Condition, Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")









################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Desender_2021_Cognit.csv"
  # "data_Rouault_2018_Expt1.csv",
  # "data_Reyes_2015.csv",
  # "data_Reyes_unpub.csv",
  # "data_Schmidt_2019_memory.csv",
  # "data_Schmidt_2019_perception.csv",
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Function to load and clean a dataset
load_and_clean_dataset <- function(filename) {

  # filename = datasets_wrong_right[4]

  filepath <- file.path(dataPath, filename)

  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filename))
    return(NULL)
  }

  # Read data
  data <- read.csv(filepath, fileEncoding = "UTF-8-BOM")

  # Add dataset identifier
  data$Dataset <- gsub(".csv", "", filename)

  # Ensure required columns exist
  required_cols <- c("Subj_idx", "Stimulus", "Response", "Confidence")
  if (!all(required_cols %in% names(data))) {
    warning(paste("Missing required columns in:", filename))
    return(NULL)
  }

  # Add Correct column if not present
  if (!"Correct" %in% names(data)) {
    data$Correct <- as.numeric(data$Stimulus == data$Response)
  }

  # Clean data
  data <- data %>%
    filter(!is.na(Confidence), !is.na(Correct)) %>%
    mutate(
      Correct_label = ifelse(Correct == 1, "Correct", "Incorrect"),
      Subj_idx = as.factor(Subj_idx)
    )

  return(data)
}

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec","Response")) %>%
  group_by(Condition, Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")


all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec")) %>%
  group_by(Condition, Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Condition, scales = "free")



all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec")) %>%
  filter(name == "Confidence") %>%
  group_by(Subj_idx, Condition, Difficulty,name,Correct) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  geom_smooth(se = F)+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")


all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) > 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence","RT_conf","RT_dec")) %>%
  filter(name == "Confidence") %>%
  group_by(Subj_idx, Condition, Difficulty,name,Correct) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  geom_smooth(se = F)+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")



all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Response")) %>%
  filter(name == "Response") %>%
  group_by(Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  geom_smooth()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")




all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) > 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("RT_dec")) %>%
  filter(name == "RT_dec") %>%
  group_by(Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_smooth()+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")





all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("RT_dec")) %>%
  filter(name == "RT_dec") %>%
  group_by(Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  geom_smooth()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")




all_data %>%
  filter(RT_dec < 10 & Confidence > -1 & Condition == 1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Response")) %>%
  filter(name == "Response") %>%
  group_by(Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_smooth()+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")


all_data %>% filter(Condition == 2) %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Response")) %>%
  filter(name == "Response") %>%
  group_by(Volatility,Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Volatility)))+
  geom_smooth()+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")



all_data %>% filter(Condition == 2) %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("RT_dec")) %>%
  filter(name == "RT_dec") %>%
  group_by(Volatility,Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Volatility)))+
  geom_smooth()+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(~Subj_idx, scales = "free")

all_data %>% filter(Condition == 2) %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(RT_dec < 10) %>%
  filter(Training != 1 & as.numeric(Subj_idx) < 15) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  pivot_longer(cols = c("Confidence")) %>%
  filter(name == "Confidence") %>%
  group_by(Correct,Volatility,Subj_idx, Condition, Difficulty,name) %>%
  summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_smooth(se = F)+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(Volatility~Subj_idx, scales = "free")




df1 = all_data %>% filter(Condition == 2) %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(RT_dec < 10) %>%
filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, Coherence, -Coherence)) %>%
  rename(
    Confidence = Confidence,
    X  = Difficulty,
    Y = Response,
    Correct = Correct,
    RT = RT_dec,
    subject = Subj_idx
  ) %>%
  mutate(resp = Y) %>%
  filter(Training != 1 & as.numeric(subject)  < 15)



df1$subject = as.numeric(as.factor(df1$subject))

t_p_s = df1 %>% group_by(subject) %>% summarize(n = n())


ends <- cumsum(t_p_s$n)

# Calculate the start points
starts <- c(1, head(ends, -1) + 1)

# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","ACC_Bin_hier_Correct.stan"))
mod = cmdstanr::cmdstan_model(here::here("Stanmodels","Discrete Confidence","discrete_acc_bin_hier.stan"))
# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","Discrete Confidence","discrete_acc_bin_single_subject.stan"))

# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","Pure","psycho.stan"))

datastan = list(N = nrow(df1),
                S = length(unique(df1$subject)),
                starts = starts,
                minRT = df1 %>% group_by(subject) %>% summarize(minRT = min(RT)) %>% .$minRT,
                ends = ends,
                t_p_s = t_p_s$n,
                X = df1$X,
                X_scaled = df1$X,
                S_id = df1$subject,
                RT = df1$RT,
                ACC = df1$Correct,
                K = length(unique(df1$Confidence)),
                Conf = df1$Confidence,
                binom_y = df1$Correct)



cor <-mod$sample(
  data = datastan,
  refresh = 10,
  iter_sampling = 1000,
  iter_warmup = 500,
  adapt_delta = 0.99,
  max_treedepth = 12,
  # init  = 0,
  parallel_chains = 4)


cor$save_object(here::here("Saved models","Desender_Cognit.rds"))





################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Rouault_2018_Expt1.csv"
  # "data_Reyes_2015.csv",
  # "data_Reyes_unpub.csv",
  # "data_Schmidt_2019_memory.csv",
  # "data_Schmidt_2019_perception.csv",
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, DotDiff , -DotDiff )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response")) %>%
  group_by(Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")


all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 1, DotDiff, -DotDiff)) %>%
  pivot_longer(cols = c("Confidence","RT_dec")) %>%
  group_by(Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")









################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Reyes_2015.csv"
  # "data_Reyes_unpub.csv",
  # "data_Schmidt_2019_memory.csv",
  # "data_Schmidt_2019_perception.csv",
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")


all_data %>%
  filter(RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  group_by(Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")








################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Reyes_unpub.csv"
  # "data_Schmidt_2019_memory.csv",
  # "data_Schmidt_2019_perception.csv",
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Condition ,Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(Condition ~name, scales = "free")


all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  group_by(Condition, Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(Condition~name, scales = "free")








################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Schmidt_2019_perception.csv"
  # "data_Pereira_2018.csv",
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Condition ,Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(Condition ~name, scales = "free")


all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast, -Contrast)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  group_by(Condition, Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(Condition~name, scales = "free")






################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Pereira_2018.csv"
  # "data_Prieto_unpub.csv",
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease            )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap( ~name, scales = "free")



all_data %>% filter(as.numeric(Subj_idx) < 14) %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease            )) %>%
  mutate(Difficulty = cut(Difficulty,10)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  filter(name != "RT_dec") %>%
  group_by(Subj_idx,Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name ~Subj_idx, scales = "free")


all_data %>% filter(as.numeric(Subj_idx) > 14) %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease            )) %>%
  mutate(Difficulty = cut(Difficulty,5)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Subj_idx,Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  filter(name != "RT_dec") %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name ~Subj_idx, scales = "free")






all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease           )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  group_by(Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")





all_data %>% filter(as.numeric(Subj_idx) < 10) %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease           )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  filter(name != "RT_conf") %>%
  group_by(Subj_idx,Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~Subj_idx, scales = "free")



all_data %>% filter(as.numeric(Subj_idx) < 14) %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease            )) %>%
  group_by(Subj_idx) %>%
  mutate(Difficulty = cut(Difficulty,10)) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  filter(name == "Response") %>%
  group_by(Subj_idx,Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name ~Subj_idx, scales = "free")


all_data %>% filter(as.numeric(Subj_idx) < 12) %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease           , -Ease            ), Response = Response-1) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  filter(name == "Response") %>%
  ggplot(aes(x = Difficulty, y = value))+
  geom_jitter(height = 0.05, width = 0.3)+
  theme_minimal(base_size = 16)+
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE) +
  facet_grid(name ~Subj_idx, scales = "free")




all_data %>% group_by(Subj_idx) %>% summarize(mean = mean(Correct)) %>% ggplot(aes(x = Subj_idx, y = mean))+geom_point()


badid = "18"

summary(all_data)
hist(all_data$RT_dec)


df1 = all_data %>% filter(Subj_idx != badid & RT_dec > 0 & RT_dec < 500) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Ease, -Ease),
         Response = Response-1) %>%
  rename(
    Confidence = Confidence,
    X  = Difficulty,
    Y = Response,
    Correct = Correct,
    RT = RT_dec,
    subject = Subj_idx
  ) %>%
  mutate(resp = Y)
  # filter(as.numeric(subject) %in% c(1,2,7,10,11))


df1 %>%
  mutate(X = cut(X,7)) %>%
  pivot_longer(cols = c("Confidence","resp","RT")) %>%
  group_by(X,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = X, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap( ~name, scales = "free")

df1 %>%
  mutate(X = cut(X,7)) %>%
  pivot_longer(cols = c("Confidence","resp","RT")) %>%
  group_by(subject,X,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = X, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(subject~name, scales = "free")


df1 %>%
  pivot_longer(cols = c("Confidence","RT")) %>%
  group_by(subject,X,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = X, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(name~subject, scales = "free")


df1 %>%
  pivot_longer(cols = c("Confidence","RT")) %>%
  group_by(X,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = X, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_grid(~name, scales = "free")




df1$subject = as.numeric(as.factor(df1$subject))

t_p_s = df1 %>% group_by(subject) %>% summarize(n = n())


ends <- cumsum(t_p_s$n)

# Calculate the start points
starts <- c(1, head(ends, -1) + 1)

# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","ACC_Bin_hier_Correct.stan"))
mod = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","ACC_Bin_hier_Correct.stan"))

# mod = cmdstanr::cmdstan_model(here::here("Stanmodels","ss","Pure","psycho.stan"))

datastan = list(N = nrow(df1),
                S = length(unique(df1$subject)),
                starts = starts,
                minRT = df1 %>% group_by(subject) %>% summarize(minRT = min(RT)) %>% .$minRT,
                ends = ends,
                t_p_s = t_p_s$n,
                X = df1$X,
                X_scaled = df1$X,
                S_id = df1$subject,
                RT = df1$RT,
                ACC = df1$Correct,
                K = length(unique(df1$Confidence)),
                Conf = df1$Confidence,
                binom_y = df1$Correct)



cor <-mod$sample(
  data = datastan,
  refresh = 10,
  iter_sampling = 1000,
  iter_warmup = 1000,
  adapt_delta = 0.90,
  max_treedepth = 10,
  init  = 0,
  parallel_chains = 4)


# cor$save_object(here::here("Saved models","pereira_all_subs.rds"))

source("~/Multivariate-Confidence-modeling/Analysis/Functions/Correct_models/plot_ACC_hier_simp.R")

preds = Get_predictive_ACC_hier(pereira_all_subs,df1,20)

bin_pure = Plot_psychometric_ACC_hier(preds,df1, bin = 10)

bin_pure

Plot_psychometric_resp_hier(preds,df1, bin = 10)





rt = Plot_RT_hier(preds,df1, bin = 10)
rt


conf = Plot_Conf_ACC_hier(preds,df1, bin = 10)
conf


predictions = Get_predictive_group(pereira_all_subs ,df1,2000)

plots = Plot_group_predictive_psycho(predictions, df1,n_bins = 10)

plots


ggsave(here::here("Extero_VMP1_marginal_means.tiff"), plots$plot_combined, dpi = 100, width = 12, height = 8)








################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_Prieto_unpub.csv"
  # "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast            , -Contrast             )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response","RT_conf")) %>%
  group_by(Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap( ~name, scales = "free")


all_data %>%
  filter(RT_conf < 10, RT_dec < 10 & Confidence > -1) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Contrast            , -Contrast            )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","RT_conf")) %>%
  group_by(Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")






################### next

# Load required libraries
pacman::p_load(tidyverse, patchwork, here)

# Set data path
dataPath <- here::here("Data", "Confidence Database", "Confidence Database")

# Define the datasets with "wrong to right" confidence scales
datasets_wrong_right <- c(
  "data_CalderTravis_unpub.csv"
)

# Load all datasets
cat("Loading datasets...\n")
all_data <- map_dfr(datasets_wrong_right, load_and_clean_dataset)

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total datasets loaded:", length(unique(all_data$Dataset)), "\n")
cat("Total subjects:", length(unique(all_data$Subj_idx)), "\n")
cat("Total trials:", nrow(all_data), "\n\n")

# Print dataset-specific summaries
dataset_summary <- all_data %>%
  group_by(Dataset) %>%
  summarise(
    N_subjects = n_distinct(Subj_idx),
    N_trials = n(),
    Mean_Accuracy = mean(Correct, na.rm = TRUE),
    Has_RT_dec = any(!is.na(RT_dec)),
    # Has_RT_conf = any(!is.na(RT_conf)),
    Conf_min = min(Confidence, na.rm = TRUE),
    Conf_max = max(Confidence, na.rm = TRUE)
  )

print(dataset_summary)

################################################################################
# PLOTTING FUNCTION
################################################################################

dataset = 1

all_data %>%
  filter(RT_dec < 10) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Stim_ref                              , -Stim_ref                               )) %>%
  pivot_longer(cols = c("Confidence","RT_dec","Response")) %>%
  group_by(Difficulty,name) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap( ~name, scales = "free")


all_data %>%
  filter(RT_dec < 10) %>%
  filter(Dataset == gsub(".csv","",datasets_wrong_right[dataset])) %>%
  mutate(Difficulty = ifelse(Stimulus == 2, Stim_ref            , -Stim_ref            )) %>%
  pivot_longer(cols = c("Confidence","RT_dec")) %>%
  group_by(Difficulty,name, Correct) %>% summarize(mean = mean(value,na.rm = T)) %>%
  ggplot(aes(x = Difficulty, y = mean, col = as.factor(Correct)))+
  geom_point()+
  theme_minimal(base_size = 16)+
  facet_wrap(~name, scales = "free")

