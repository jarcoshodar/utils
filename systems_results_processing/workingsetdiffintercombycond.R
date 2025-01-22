library(ggplot2)
library(dplyr)

setwd("~/Documentos")

files <- list.files("intercomall", recursive = TRUE)

file_parts <- strsplit(files, "/")

df <- do.call(rbind, lapply(file_parts, function(parts) {
  study <- parts[1]
  organ <- parts[2]
  filename <- parts[3]
  
  condition <- ifelse(grepl("Het-O", filename), "Het-O",
                      ifelse(grepl("O", filename), "O",
                             ifelse(grepl("IA", filename), "IA", NA)))
  
  return(data.frame(Study = study, Organ = organ, Condition = condition, stringsAsFactors = FALSE))
}))

process_file <- function(filepath, study, organ, treatment) {
  df <- read.csv(filepath)
  df$Interaction <- paste(df$Ligand, df$Receptor, sep = " ")
  interaction_df <- data.frame(Study = study, Organ = organ, Treatment = treatment, Interaction = df$Interaction, stringsAsFactors = FALSE)
  return(interaction_df)
}

interaction_df <- do.call("rbind", lapply(files, function(filepath) {
  parts <- strsplit(filepath, "/")[[1]]
  study <- parts[1]
  organ <- parts[2]
  filename <- parts[3]
  condition <- ifelse(grepl("Het-O", filename), "Het-O",
                      ifelse(grepl("O", filename), "O",
                             ifelse(grepl("IA", filename), "IA", NA)))
  treatment <- ifelse(condition %in% c("O", "IA"), "untreated", "treated")
  process_file(paste("intercomall", filepath, sep = "/"), study, organ, treatment)
}))

interaction_df_grouped <- interaction_df %>%
  group_by(Study, Organ, Treatment, Interaction) %>%
  summarise(n = n(), .groups = 'drop')

treated <- interaction_df_grouped[interaction_df_grouped$Treatment == "treated", ]
untreated <- interaction_df_grouped[interaction_df_grouped$Treatment == "untreated", ]

all_interactions <- unique(interaction_df_grouped$Interaction)

treated_counts <- treated %>% group_by(Interaction) %>% summarise(n = sum(n))
untreated_counts <- untreated %>% group_by(Interaction) %>% summarise(n = sum(n))

treated_counts[!(treated_counts$Interaction %in% all_interactions), "n"] <- 0
untreated_counts[!(untreated_counts$Interaction %in% all_interactions), "n"] <- 0

delta_counts <- merge(treated_counts, untreated_counts, by = "Interaction", suffixes = c("_treated", "_untreated"))
delta_counts$delta <- delta_counts$n_treated - delta_counts$n_untreated

delta_counts_sorted <- delta_counts[order(delta_counts$delta), ]

# Top 10 interactions increased in treated condition
top_10 <- tail(delta_counts_sorted, 10)

# Top 10 interactions increased in untreated condition
bottom_10 <- head(delta_counts_sorted, 10)

# Compute breakdowns for top 10
top_10_breakdown <- interaction_df_grouped[interaction_df_grouped$Interaction %in% top_10$Interaction, ]

# Compute breakdowns for bottom 10
bottom_10_breakdown <- interaction_df_grouped[interaction_df_grouped$Interaction %in% bottom_10$Interaction, ]

library(dplyr)

filter_by_interaction <- function(df, query) {
  filtered_df <- df %>% filter(Interaction == query)
  return(filtered_df)
}

