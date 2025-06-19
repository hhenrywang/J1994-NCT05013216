# This script calculate the correlation between peripheral and tumor frequency
# Last edited by Henry Wang on 06.18.25

##### Libraries to load
library(ggplot2)
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

##### A. Data Import & Reformat #####
# Sets the working directory.
setwd("C:/Users/henry/Dropbox/Research/Jaffee Lab/J1994/Sequencing/Bulk")

# Setting the data input directory, using Tumor subdirectory for tissue sequencing
Dirwd <- getwd()
DirOutputAll <- paste0(Dirwd, "/Output/")
DirMeta <- paste0(Dirwd, "/Input/Bulk_Putative/")
DirTE <- paste0(Dirwd, "/Input/Bulk_Putative/Treatment-enriched/")
DirTumor <- paste0(Dirwd, "/Input/Bulk_Putative/Tumor-infiltration/")

#Gets the metadata file in the Data_Bulk folder
meta = read_excel(paste0(DirMeta, "metadata.xlsx"))

meta$Patient <- meta$Patient %>%
  str_replace("^J1994\\.", "") %>%
  str_replace("^0+", "")

#### A.1 mKRAS TCRs ####
#Get the merged mKRAS TCR table and merge into a single data frame
TCRs_mKRAS <- read.csv(paste0(DirOutputAll, "All_TCR_Merged.csv"), stringsAsFactors = FALSE)

TCRs_mKRAS <- TCRs_mKRAS %>%
  filter(Patient != "Published") %>%
  left_join(meta %>% select(Patient, Mutation), by = "Patient") %>%
  filter(str_detect(Antigen, fixed(Mutation))) %>%
  pivot_longer(
    cols      = ends_with("freq"),
    names_to  = "freq_type",
    values_to = "Freq"
  ) %>%
  filter(freq_type == paste0(Mutation, "freq")) %>%
  select(CDR3aa, TCRV, TCRJ, Patient, Freq)

#### A.2 Treatment-enriched TCRs ####
#Get treatment-enriched TCRs and merge into a single data frame
file_paths <- list.files(DirTE, pattern = "\\.csv$", full.names = TRUE)

TE_list <- lapply(file_paths, function(file) {
  read.csv(file, header = TRUE)
})

names(TE_list) <- file_paths %>%
  basename() %>%
  substr(1, 2) %>%
  sub("^0+", "", .)

TCRs_TE <- bind_rows(TE_list, .id = "Patient")

TCRs_TE <- TCRs_TE %>%
  filter(mKRAS_tumor_specific == FALSE) %>%
  select(-mKRAS_tumor_specific)

rm(file_paths)
rm(TE_list)

#### A.3 Tumor TCRs ####
#Imports tumor-infiltration data
file_paths <- list.files(DirTumor, pattern = "\\.csv$", full.names = TRUE)

Tumor_list <- lapply(file_paths, function(file) {
  read.csv(file, header = TRUE)
})

names(Tumor_list) <- file_paths %>%
  basename() %>%
  substr(1, 2) %>%
  sub("^0+", "", .)

TCRs_tumor <- map_dfr(
  Tumor_list,
  ~ .x %>%
    select(
      CDR3aa, TCRV, TCRJ,
      Counts_Post, Frequency_Post, Treatment_enriched,
      ends_with("_TCR_all")
    ) %>%
    mutate(TCR = if_any(ends_with("_TCR_all"), identity)) %>%
    select(-ends_with("_TCR_all")),
  .id = "Patient"
) %>%
  filter(Counts_Post != 0) %>%
  select(-Counts_Post)

rm(Tumor_list)

#### A.4 Reformatting ####
# Create Unique columns
create_unique <- function(df) {
  df %>% mutate(Unique = str_c(Patient, CDR3aa, TCRV, TCRJ, sep = "_"))
}

TCRs_tumor <- create_unique(TCRs_tumor)
TCRs_TE    <- create_unique(TCRs_TE)
TCRs_mKRAS <- create_unique(TCRs_mKRAS)

# Merge mKRAS and TE with tumor infiltration
merge_tumor <- function(source, df, name) {
  source %>%
    filter({{name}}) %>%
    full_join(df, by = "Unique") %>%
    select(-Treatment_enriched, -TCR) %>%
    mutate(across(c(Frequency_Post, Freq), ~ replace_na(.x, 0)))
}

TCRs_TE_merged <- merge_tumor(TCRs_tumor, TCRs_TE, Treatment_enriched)
TCRs_mKRAS_merged <- merge_tumor(TCRs_tumor, TCRs_mKRAS, TCR)

#### B. Statistical Testing ####
# Spearman's correlation testing
cor.test(TCRs_mKRAS_merged$Freq, TCRs_mKRAS_merged$Frequency_Post, method = "spearman", exact = FALSE)

cor.test(TCRs_TE_merged$Freq, TCRs_TE_merged$Frequency_Post, method = "spearman", exact = FALSE)

#### C. Plotting ####
# Plotting mKRAS TCRs
TCRs_mKRAS_merged <- TCRs_mKRAS_merged %>%
  mutate(
    Freq_adj = ifelse(Freq == 0, (min(Freq[Freq > 0]) / 2), Freq),
    Frequency_Post_adj = ifelse(Frequency_Post == 0, (min(Frequency_Post[Frequency_Post > 0]) / 2), Frequency_Post)
  )

y_threshold <- 0.7 * min(TCRs_mKRAS_merged$Frequency_Post[TCRs_mKRAS_merged$Frequency_Post > 0])

plot <- ggplot(TCRs_mKRAS_merged, aes(x = Freq_adj, y = Frequency_Post_adj)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = y_threshold, linetype = "dashed", color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Peripheral Frequency after Expansion",
    y = "Tumor Frequency on-treatment",
    title = "mKRAS TCR infiltration"
  ) +
  theme(aspect.ratio = 1)

ggsave(filename = paste0(DirOutputAll, "mKRAS_infiltration.tiff"),
       plot = plot, device = "tiff", width = 6, height = 6)

# Plotting Treatment-enriched TCRs
TCRs_TE_merged <- TCRs_TE_merged %>%
  mutate(
    Freq_adj = ifelse(Freq == 0, (min(Freq[Freq > 0]) / 2), Freq),
    Frequency_Post_adj = ifelse(Frequency_Post == 0, (min(Frequency_Post[Frequency_Post > 0]) / 2), Frequency_Post)
  )
y_threshold <- 0.7 * min(TCRs_TE_merged$Frequency_Post[TCRs_TE_merged$Frequency_Post > 0])

plot <- ggplot(TCRs_TE_merged, aes(x = Freq_adj, y = Frequency_Post_adj)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = y_threshold, linetype = "dashed", color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Peripheral Frequency",
    y = "Tumor Frequency on-treatment",
    title = "Treatment-enriched TCR infiltration"
  ) +
  theme(aspect.ratio = 1)

ggsave(filename = paste0(DirOutputAll, "treatment_enriched_infiltration.tiff"),
       plot = plot, device = "tiff", width = 6, height = 6)