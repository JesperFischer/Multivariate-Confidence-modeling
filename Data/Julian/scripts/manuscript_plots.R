# Manuscript plots
## Compiled for 2025-05-30

# PRELIMINARIES ----
# Load data and assign classes

rm(list = ls())

# Load packages
library(tidyverse)
library(see)

# Load custom functions
source("supporting_functions/my_functions.R")
source("supporting_functions/plot_functions.R")

# Locate folders
data_location <- "../../data/csvs/"
fig_location <- ""

# PREPARE DATA ----

the_experiments <- c("experiment1","experiment2")

for(task_version in the_experiments){
  
  # Load data
  g_data <- read.csv(paste0(data_location,"block_data_",task_version,".csv"), na.strings = "NaN")
  t_data <- read.csv(paste0(data_location,"trial_data_",task_version,".csv"), na.strings = "NaN")
  
  # Compute total number of participants before exclusions
  pp_number <- length(length(unique(g_data$subjID)))
  
  # Excluded participants
  excluded_participants <- switch(task_version,
                                  "experiment1" = 
                                    c("008_RM", # This person quit after Day #2
                                      "105_KK", # This person quit after Day #1
                                      "111_RY"), # This person quit after Day #1
                                  "experiment2" = 
                                    c("103_HK", # This person quit after Day #1
                                      "402_MS") # This person quit after Day #1
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
  
  # Class coding & assignments
  g_data <- assign_group_variables(g_data)
  t_data <- assign_trial_variables(t_data)
  
  # Assign to appropriate experiment
  switch(task_version,
         "experiment1" = {
           group_data.exp1 <- g_data
           trial_data.exp1 <- t_data
           how_many_participants.exp1 <- pp_number
         },
         "experiment2" = {
           group_data.exp2 <- g_data
           trial_data.exp2 <- t_data
           how_many_participants.exp2 <- pp_number
         }
  )
}

# Concatenate data
g_data <- rbind(group_data.exp1, group_data.exp2)
t_data <- rbind(trial_data.exp1, trial_data.exp2)

# FIGURES ----

## Main effects of report context ----
a <- context_plot(group_data.exp1, "a) Experiment 1: Perception", "mean m-ratio")
b <- context_plot(group_data.exp2, "b) Experiment 2: Memory", "")
a <- a + theme(plot.subtitle = element_text(hjust = 0))
b <- b + theme(plot.subtitle = element_text(hjust = 0))
cowplot::plot_grid(a,b, nrow = 1, align = "h", rel_widths = c(1,1))
ggsave(paste0(fig_location, "figure4_contexts", ".png"),
       width = regular_width*0.9, height = regular_height*0.7)

## Main effects of report order ----
summary_stats <- g_data |> 
  group_by(domain_type,motor_time_type,response_type) |> 
  reframe(metacognition = efficiency)

subject_data <- summary_stats |> 
  group_by(domain_type,motor_time_type,response_type) |> 
  summarise(m_meta = mean(metacognition, na.rm = TRUE),
            se_meta = sd(metacognition, na.rm = TRUE) / 
              sqrt(nlevels(g_data$subjID)),
            .groups = "drop")

subject_data$motor_time_type <- factor(
  subject_data$motor_time_type,
  levels = c("time limited","baseline","motor fixed")
)

b <- ggplot(data = subject_data, 
            aes(x = response_type, y = m_meta, 
                fill = response_type,
                colour = response_type,
                shape = domain_type)) +
  facet_grid(. ~ motor_time_type) +
  geom_hline(yintercept = 1, colour = "grey", linewidth = 1) +
  geom_line(aes(group = interaction(domain_type,motor_time_type)), 
            linewidth = 1, 
            colour = "black",
            show.legend = FALSE) +
  geom_errorbar(aes(y = m_meta, 
                    ymin = m_meta - se_meta,
                    ymax = m_meta + se_meta),
                width = 0, 
                linewidth = 2,
                show.legend = FALSE) +
  geom_point(stroke = 1,
             colour = "black",
             fill = "white",
             size = 2) +
  labs(x = "report order", 
       y = "mean m-ratio",
       subtitle = "b) Metacognitive efficiency across all conditions") +
  scale_colour_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
  scale_shape_manual(values = c(21,24)) +
  theme_julian() + 
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(t = -0.1, unit='cm'),
    strip.text = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 0)
  )

subject_data <- summary_stats |> 
  group_by(domain_type,response_type) |> 
  summarise(m_meta = mean(metacognition, na.rm = TRUE),
            se_meta = sd(metacognition, na.rm = TRUE) / 
              sqrt(nlevels(g_data$subjID)),
            .groups = "drop")
subject_data$domain_type <- factor(
  subject_data$domain_type,
  levels = c("memory","perception")
)

a <- ggplot(subject_data, aes(x = response_type, y = m_meta, 
                              fill = response_type)) +
  facet_grid(domain_type~.) +
  geom_hline(yintercept = 1, colour = "grey", linewidth = 1) +
  geom_col(show.legend = FALSE, width = 0.8) +
  geom_errorbar(aes(ymin = m_meta - se_meta, ymax = m_meta + se_meta),
                width = 0, linewidth = 2, colour = "black") +
  labs(x = "",
       y = "mean m-ratio",
       subtitle = "a) Report order") +
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb")) +
  theme_julian() + 
  theme(strip.text = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 12, hjust = 0))

cowplot::plot_grid(a,b,nrow = 1, rel_widths = c(1,2.2))
ggsave(paste0(fig_location,"figure5_order", ".png"),
       width = regular_width, height = regular_height*0.8)
