# This script loads Adaptive Bulk TCRSeq Data to identify mKRAS T-cells per patient
# Clonotypes are tracked by their unique TCR Beta CDR3 amino acid sequence
# Last edited by Henry Wang on 06.19.25

# Loading Libraries
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)
library(readxl)

##### 0. Variable Selection ####
# Sets the patient ID for analysis
Patient_ID <- '29'

# Sets FDR for mKRAS TCRs
FDR <- 0.01

# Sets min. fold-exp for mKRAS TCRs
Ratio <- 5

# Sets min. post-exp freq for mKRAS TCRs
exp_freq <- 2.5e-4

# Sets FDR for TE TCRs
FDR_TE <- 0.01

##### A. Loading Data #####
# Set the working directory. Change as needed
setwd("C:/Users/henry/Dropbox/Research/Jaffee Lab/J1994/Sequencing/Bulk")

# Setting input and output directories
Dirwd <- getwd()
DirMeta <- paste0(Dirwd, "/Data_Bulk/") #Location of metadata file
DirDataBulk <- paste0(Dirwd, "/Data_Bulk/", Patient_ID, "/") #Location of input sequencing
DirOutput <- paste0(Dirwd, "/Output/", Patient_ID, "/") #Output folder

# Retrieves file paths for all .tsv files in /Data_Bulk
Dirtsv_file_names <- list.files(path = DirDataBulk, pattern = ".tsv$")
Dirtsv_files <- paste0(DirDataBulk, Dirtsv_file_names)

# Loads metadata
meta_path <- file.path(DirMeta, "metadata.xlsx")
meta <- read_excel(meta_path)
rm(meta_path)

# Identifies tumor-specific KRAS mutation for the patient
mKRAS_tumor_specific <- meta$Tumor_specific_mKRAS[
  meta$Patient == Patient_ID & meta$Sample == 'Baseline']

# Read the .tsv files into a list of data frames - BulkTCR
BulkTCR <- lapply(Dirtsv_files, function(file) {
  read.delim(file, header = TRUE, sep = "\t")
})
names(BulkTCR) <- gsub("\\.tsv$", "", Dirtsv_file_names)
rm(Dirtsv_file_names, Dirtsv_files)

# Renames each file based on metadata
for (i in seq_along(BulkTCR)) {
  current_name <- names(BulkTCR)[i] 
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- matching_row$Sample 
    names(BulkTCR)[i] <- new_name
    rm(matching_row)
  }
  rm(current_name, new_name)
}

#### B. Data Clean-up ####
# Only keep columns that have an in frame CDR3aa (productive TCRs)
BulkTCR <- lapply(BulkTCR, function(df){
  df <- df[df$sequenceStatus == "In",]
  return(df)
})

# Keeping and renaming the indicated columns
columns_to_keep <- c("aminoAcid",
                     "count..templates.reads.",
                     "vMaxResolved",
                     "dMaxResolved",
                     "jMaxResolved")
new_column_names <- c("CDR3aa",
                      "Counts",
                      "vResolved",
                      "dResolved",
                      "jResolved")

# Adjust columns in BulkTCR
BulkTCR <- lapply(BulkTCR, function(df) {
  df <- df %>% select(all_of(columns_to_keep))
  colnames(df) <- new_column_names
  return(df)
})
rm(columns_to_keep, new_column_names)

# Create a unique identifier for each clonotype based on the CDR3aa, TCRV, TCRJ
add_unique_column <- function(df) {
  df$Unique <- paste(df$CDR3aa, df$vResolved, df$jResolved, sep = "_")
  return(df)
}
BulkTCR <- lapply(BulkTCR, add_unique_column)

# Merge rows with same unique identifier and sum the counts
BulkTCR_unique <- lapply(BulkTCR, function(df) {
  df %>%
    group_by(Unique) %>%
    summarise(
      Counts = sum(Counts, na.rm = TRUE),
      Clonotypes = n(),
      across(-Counts, first),
      .groups = "drop"
    )
})

# Create a new summary data frame 'meta_lim'
meta_lim <- data.frame(names(BulkTCR))
colnames(meta_lim) <- 'Sample'

# Renames the data frames in BulkTCR_unique
names(BulkTCR_unique) <- meta_lim$Sample

# Save # productive counts in meta_lim file
for(i in seq_along(BulkTCR_unique)){
  Total <- sum(BulkTCR_unique[[i]]$Counts)
  meta_lim$Total_Counts[[i]] <- Total
  meta_lim$Unique_Clonotypes[[i]] <- nrow(BulkTCR_unique[[i]])
  rm(Total)
}
meta_lim$Total_Counts <- as.numeric(as.character(meta_lim$Total_Counts))
meta_lim$Unique_Clonotypes <- as.numeric(as.character(meta_lim$Unique_Clonotypes))

# For log-log plotting purposes, create a pseudo-frequency count
meta_lim$Unobs_Freq <- (3* meta_lim$Total_Counts)^(-1)

# Calculates clonotype frequencies based on Total_Counts
BulkTCR_unique <- lapply(names(BulkTCR_unique), function(sample_name) {
  df <- BulkTCR_unique[[sample_name]]
  total_counts <- meta_lim$Total_Counts[meta_lim$Sample == sample_name]
  df$Frequency <- df$Counts / total_counts
  df
})
names(BulkTCR_unique) <- meta_lim$Sample
rm(BulkTCR)

#### C. Identifying mKRAS clonotypes ####
#### C.1 Merging Datasets ####
# Merging sequencing for each mKRAS expansion against FLC expansion
# For each row, replace 0 with NA under counts
# Note that Counts.x represents the mKRAS exp count whereas Count.y represents the FLC / baseline / Unstim condition
# BulkTCR_MergedU compares each expansion condition against the Unstim on-treatment sample
# BulkTCR_MergedF compares each expansion condition against the FLC expansion sample
# BulkTCR_MergedP compares the on-treatment Unstim sample vs Baseline pre-treatment sample
peptide_names <- meta_lim$Sample
peptide_names <- peptide_names[!(peptide_names %in% c("Baseline", "Unstim"))]
peptide_names_no_FLC <- setdiff(peptide_names, "FLC")

Baseline <- BulkTCR_unique[["Baseline"]]
Unstim <- BulkTCR_unique[["Unstim"]]
FLC <- BulkTCR_unique[["FLC"]]

BulkTCR_MergedU <- lapply(peptide_names, function(df) {
  merge(BulkTCR_unique[[df]], Unstim, by = c("CDR3aa", "vResolved", "jResolved", "Unique"), all=TRUE)
})
names(BulkTCR_MergedU) <- peptide_names

for(i in seq_along(BulkTCR_MergedU)){
  BulkTCR_MergedU[[i]] <- BulkTCR_MergedU[[i]] %>%
    mutate (Counts.x = replace_na(Counts.x, 0)) %>%
    mutate (Counts.y = replace_na(Counts.y, 0))
}

# For BulkTCR_MergedF, no FLC since comparing against FLC
BulkTCR_MergedF <- lapply(peptide_names, function(df) {
  merge(BulkTCR_unique[[df]], FLC, by = c("CDR3aa", "vResolved", "jResolved", "Unique"), all=TRUE)
})
names(BulkTCR_MergedF) <- peptide_names
BulkTCR_MergedF[["FLC"]] <- NULL

for(i in seq_along(BulkTCR_MergedF)){
  BulkTCR_MergedF[[i]] <- BulkTCR_MergedF[[i]] %>%
    mutate (Counts.x = replace_na(Counts.x, 0)) %>%
    mutate (Counts.y = replace_na(Counts.y, 0))
}

# Finally, create merged data set for non-expansion conditions
BulkTCR_MergedP <- merge(Unstim, Baseline, by = c("CDR3aa", "vResolved", "jResolved", "Unique"), all=TRUE)

BulkTCR_MergedP <- BulkTCR_MergedP %>%
  mutate (Counts.x = replace_na(Counts.x, 0)) %>%
  mutate (Counts.y = replace_na(Counts.y, 0))

# Set unidentified pseudo_freq
Unobs_Freq.y <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Unstim'][1]

# Update each condition where Counts is 0 with pseudo-freq
BulkTCR_MergedU <- lapply(names(BulkTCR_MergedU), function(sample_name) {
  df <- BulkTCR_MergedU[[sample_name]]
  Unobs_Freq.x <- meta_lim$Unobs_Freq[meta_lim$Sample == sample_name][1]
  
  df$Frequency.x <- ifelse(
    df$Counts.x == 0,
    Unobs_Freq.x,
    df$Frequency.x
  )
  df$Frequency.y <- ifelse(
    df$Counts.y == 0,
    Unobs_Freq.y,
    df$Frequency.y
  )
  return(df)
})
names(BulkTCR_MergedU) <- peptide_names

# Next, for BulkTCR_MergedP
Unobs_Freq.x <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Unstim']
Unobs_Freq.y <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Baseline']

BulkTCR_MergedP$Frequency.x <- ifelse(
  BulkTCR_MergedP$Counts.x == 0,
  Unobs_Freq.x,
  BulkTCR_MergedP$Frequency.x
)
BulkTCR_MergedP$Frequency.y <- ifelse(
  BulkTCR_MergedP$Counts.y == 0,
  Unobs_Freq.y,
  BulkTCR_MergedP$Frequency.y
)

# Finally, for BulkTCR_MergedF
Unobs_Freq.y <- meta_lim$Unobs_Freq[meta_lim$Sample == 'FLC']

BulkTCR_MergedF <- lapply(names(BulkTCR_MergedF), function(sample_name) {
  df <- BulkTCR_MergedF[[sample_name]]
  Unobs_Freq.x <- meta_lim$Unobs_Freq[meta_lim$Sample == sample_name][1]

  df$Frequency.x <- ifelse(
    df$Counts.x == 0,
    Unobs_Freq.x,
    df$Frequency.x
  )
  df$Frequency.y <- ifelse(
    df$Counts.y == 0,
    Unobs_Freq.y,
    df$Frequency.y
  )
  return(df)
})
names(BulkTCR_MergedF) <- peptide_names_no_FLC

rm(Unobs_Freq.x, Unobs_Freq.y, FLC)

#### C.2 Stat. Significance ####
# For each mKRAS peptide exp, conduct Fisher's exact test against the FLC control expansion
# Then, calculate p adjusted using Benjamini-Yekutieli's correction for FDR
total_Base <- meta_lim$Total_Counts[meta_lim$Sample == 'Baseline']
total_Unstim <- meta_lim$Total_Counts[meta_lim$Sample == 'Unstim']
total_FLC  <- meta_lim$Total_Counts[meta_lim$Sample == 'FLC']

# First, BulkTCR_MergedF
for (j in seq_along(BulkTCR_MergedF)) {
  
  df <- BulkTCR_MergedF[[j]]
  total_temp <- meta_lim$Total_Counts[meta_lim$Sample == names(BulkTCR_MergedF)[j]]
  
  contig_tables <- mapply(
    function(count_x, count_y) {
      matrix(c(
        count_x, total_temp - count_x,
        count_y, total_FLC - count_y
      ), nrow = 2)
    },
    df$Counts.x, df$Counts.y, SIMPLIFY = FALSE
  )
  p_values <- sapply(contig_tables, function(tbl) fisher.test(tbl, alternative="greater")$p.value)
  df$P_value <- p_values
  df$P_value_BY <- p.adjust(p_values, method = "BY")
  BulkTCR_MergedF[[j]] <- df
}
rm(df)

# Now, repeating the same but for unstim on-treatment vs baseline by testing BulkTCR_MergedP
contig_tables <- mapply(
  function(count_x, count_y) {
    matrix(c(
      count_x, total_Unstim - count_x,
      count_y, total_Base - count_y
    ), nrow = 2)
  },
  BulkTCR_MergedP$Counts.x, BulkTCR_MergedP$Counts.y, SIMPLIFY = FALSE
)
p_values <- sapply(contig_tables, function(tbl) fisher.test(tbl, alternative="greater")$p.value)
BulkTCR_MergedP$P_value <- p_values
BulkTCR_MergedP$P_value_BY <- p.adjust(p_values, method = "BY")

rm(contig_tables,
   p_values,
   total_Base,
   total_Unstim,
   total_FLC,
   total_temp)

#### C.3 FDR & Ratio Filter ####
# Separate out the FLC expansion to a separate data frame and indicate if Sig
BulkTCR_Merged_FvUnstim <- BulkTCR_MergedU$FLC
BulkTCR_MergedU$FLC <- NULL

# Now filter on the list of data frames
BulkTCR_MergedF_sig <- list()
for( i in seq_along(BulkTCR_MergedF)){
  df <- BulkTCR_MergedF[[i]] %>% filter(P_value_BY < FDR)
  df <- df %>% filter(Frequency.x / Frequency.y > Ratio)
  BulkTCR_MergedF_sig[[i]] <- df
  rm(df)
}
names(BulkTCR_MergedF_sig) <- peptide_names_no_FLC

# For the Treatment-enriched clonotypes, now only include clonotypes that were not present at baseline in the blood
BulkTCR_MergedP_sig <- BulkTCR_MergedP %>%
  filter(P_value_BY < FDR_TE) %>%
  filter(!(Unique %in% Baseline$Unique))

#### C.4 Min Freq & De Novo Filter ####
# First, filter by post-exp frequency
BulkTCR_MergedF_final <- list()
for( i in seq_along(BulkTCR_MergedF_sig)){
  BulkTCR_MergedF_final[[i]] <- BulkTCR_MergedF_sig[[i]] %>% filter(Frequency.x >= exp_freq)
}
names(BulkTCR_MergedF_final) <- peptide_names_no_FLC

# Remove clonotypes that are present in Baseline PBMC
BulkTCR_MergedF_sig <- lapply(BulkTCR_MergedF_sig, function(df) {
  df %>% filter(!(Unique %in% Baseline$Unique))
})
BulkTCR_MergedF_final <- lapply(BulkTCR_MergedF_final, function(df) {
  df %>% filter(!(Unique %in% Baseline$Unique))
})

# Now annotate the original data frames for plotting purposes
add_significance <- function(df_list, input_list, col_name) {
  Map(function(df, input_list) {
    df_filt <- input_list$Unique
    df[[col_name]] <- df$Unique %in% df_filt
    df
  }, df_list, input_list)
}

BulkTCR_MergedU <- add_significance(BulkTCR_MergedU, BulkTCR_MergedF_sig, "Sig")
BulkTCR_MergedU <- add_significance(BulkTCR_MergedU, BulkTCR_MergedF_final, "Sig_final")

BulkTCR_MergedF <- add_significance(BulkTCR_MergedF, BulkTCR_MergedF_sig, "Sig")
BulkTCR_MergedF <- add_significance(BulkTCR_MergedF, BulkTCR_MergedF_final, "Sig_final")

# Now identify mKRAS clonotypes in the on-treatment vs baseline samples
BulkTCR_MergedP$Sig <- BulkTCR_MergedP$Unique %in% BulkTCR_MergedP_sig$Unique

temp_df <- BulkTCR_MergedF_final[[mKRAS_tumor_specific]]
BulkTCR_MergedP <- BulkTCR_MergedP %>%
  mutate(
    mKRAS_tumor_specific = Unique %in% temp_df$Unique
  )
rm(temp_df)

#### D. Pairwise Plots####
Baseline_plots <- list()
Unstim_plots <- list()

#### D.1 Unstim vs Baseline ####
# Extract Unobs_Freq values for Baseline and Unstim
Unobs_Freq_Baseline <- meta_lim$Unobs_Freq[meta_lim$Sample == "Baseline"]
Unobs_Freq_Unstim <- meta_lim$Unobs_Freq[meta_lim$Sample == "Unstim"]

# Plotting treatment-enriched TCRs
BulkTCR_MergedP <- BulkTCR_MergedP %>%
  mutate(clonotype_group = ifelse(Sig, "deepskyblue", "black"))

plot_data <- BulkTCR_MergedP %>%
  group_by(Frequency.x, Frequency.y, clonotype_group) %>%
  summarise(n = n(), .groups = "drop")

Periph_plot_sig <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
  geom_point(data = filter(plot_data, clonotype_group == "black"),
             color = "black", alpha = 0.3) +
  
  geom_point(data = filter(plot_data, clonotype_group == "deepskyblue"),
             color = "deepskyblue", alpha = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 2 * Unobs_Freq_Baseline, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2 * Unobs_Freq_Unstim,   color = "black", linetype = "dashed", linewidth = 1) +
  scale_size_continuous(trans = "log2",
                        name = "# Clonotypes",
                        breaks = c(1, 10, 100, 1000),
                        labels = c("1", "10", "100", "1000"),
                        range = c(2,6)) +
  labs(x = "Baseline Peripheral Frequency",
    y = "On-treatment Peripheral Frequency",
    title = "Peripheral TCRs - Baseline vs On-treatment (Treatment-enriched)") +
  theme(legend.position = "right", 
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  coord_equal()

ggsave(filename = paste0(DirOutput, Patient_ID, "_Periph_treatment_induced.tiff"),
  plot = Periph_plot_sig,
  device = "tiff", width = 8, height = 6)

# Plotting tumor-specific mKRAS
BulkTCR_MergedP <- BulkTCR_MergedP %>%
  mutate(clonotype_group = ifelse(mKRAS_tumor_specific, "orange", "black"))

plot_data <- BulkTCR_MergedP %>%
  group_by(Frequency.x, Frequency.y, clonotype_group) %>%
  summarise(n = n(), .groups = "drop")

Periph_plot_mKRAS <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
  geom_point(data = filter(plot_data, clonotype_group == "black"),
             color = "black", alpha = 0.3) +
  
  geom_point(data = filter(plot_data, clonotype_group == "orange"),
             color = "orange", alpha = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 2 * Unobs_Freq_Baseline, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2 * Unobs_Freq_Unstim, color = "black", linetype = "dashed", linewidth = 1) +
  scale_size_continuous(
    trans = "log2",
    name = "# Clonotypes",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000"),
    range = c(2,6)) +
  labs(x = "Baseline Peripheral Frequency",
    y = "On-treatment Peripheral Frequency",
    title = paste("Peripheral TCRs - Baseline vs On-treatment -", mKRAS_tumor_specific, "TCRs")) +
  theme(legend.position = "right", 
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  coord_equal()

ggsave(filename = paste0(DirOutput, Patient_ID, "_Periph_Tumor_Specific_mKRAS_TCRs.tiff"),
  plot = Periph_plot_mKRAS,
  device = "tiff", width = 8, height = 6)

#### D.2 mKRAS vs FLC ####
Unobs_Freq_FLC <- meta_lim$Unobs_Freq[meta_lim$Sample == "FLC"]
Baseline_plots <- vector("list", length(BulkTCR_MergedF))
names(Baseline_plots) <- peptide_names_no_FLC

for (i in seq_along(BulkTCR_MergedF)) {
  peptide_name <- names(BulkTCR_MergedF)[i]
  df <- BulkTCR_MergedF[[i]]
  Unobs_Freq <- meta_lim$Unobs_Freq[meta_lim$Sample == peptide_name]
  
  df <- df %>%
    mutate(clonotype_group = ifelse(Sig_final, "orange", "black"))
  
  plot_data <- df %>%
    group_by(Frequency.x, Frequency.y, clonotype_group) %>%
    summarise(n = n(), .groups = "drop")
  
  plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
    geom_point(data = filter(plot_data, clonotype_group == "black"),
               color = "black", alpha = 0.3) +
    geom_point(data = filter(plot_data, clonotype_group == "orange"),
               color = "orange", alpha = 0.5, show.legend = FALSE) +
    geom_vline(xintercept = 2 * Unobs_Freq_FLC, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 2 * Unobs_Freq,     color = "black", linetype = "dashed", linewidth = 1) +
    scale_size_continuous(
      trans = "log2",
      name = "# Clonotypes",
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000"),
      range = c(2,6)) +
    labs(x = "Frequency with FLC Expansion",
      y = paste("Frequency with", peptide_name, "Expansion"),
      title = paste("TCR Frequencies -", peptide_name, "vs FLC")) +
    theme(legend.position = "right", 
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    coord_equal()
  Baseline_plots[[i]] <- plot
}

for (i in seq_along(Baseline_plots)) {
  ggsave(filename = paste0(DirOutput, Patient_ID, "_", names(Baseline_plots)[i], "_vs_FLC.tiff"),
    plot = Baseline_plots[[i]],
    device = "tiff", width = 8, height = 6)
}

#### D.3 FLC vs On-treatment ####
BulkTCR_Merged_FvUnstim <- BulkTCR_Merged_FvUnstim %>%
  mutate(clonotype_group = "black")

plot_data <- BulkTCR_Merged_FvUnstim %>%
  group_by(Frequency.x, Frequency.y, clonotype_group) %>%
  summarise(n = n(), .groups = "drop")

Periph_plot_sig <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
  geom_point(color = "black", alpha = 0.3) +
  geom_vline(xintercept = 2 * Unobs_Freq_Unstim, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2 * Unobs_Freq_FLC, color = "black", linetype = "dashed", linewidth = 1) +
  scale_size_continuous(
    trans = "log2",
    name = "# Clonotypes",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000"),
    range = c(2,6)) +
  labs(x = "Frequency on Treatment",
    y = "Frequency with FLC Stim",
    title = "TCR Frequencies - On Treatment vs FLC Stim") +
  theme(legend.position = "right", 
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  coord_equal()

ggsave(filename = paste0(DirOutput, Patient_ID, "_FLC.tiff"),
  plot = Periph_plot_sig,
  device = "tiff", width = 8, height = 6)

rm(Baseline_plots,
   Unstim_plots,
   Periph_plot_sig,
   Unobs_Freq,
   Unobs_Freq_Baseline,
   Unobs_Freq_Unstim,
   Unobs_Freq_FLC,
   plot)

#### E. Additional Stats ####
# All additional statistics are stored in meta_lim table for export

#### E.1 Total Frequencies ####
# Tracks total mKRAS clonotype frequencies across conditions
# First, extract the total frequencies of expanded clonotypes in Baseline, Unstim, FLC conditions
Frequency <- matrix(nrow = length(peptide_names_no_FLC), ncol = 3)
rownames(Frequency) <- peptide_names_no_FLC
colnames(Frequency) <- c('Unstim Freq', 'FLC Stim Freq', 'Stim Freq')

# Now in BulkTCR_MergedF_final, include the Unstim and FLC frequencies for each clonotype
for (i in seq_along(BulkTCR_MergedF_final)) {
  df_U <- BulkTCR_MergedF_final[[i]]
  df_B <- BulkTCR_MergedP %>%
    filter(Counts.x != 0) %>%
    select(Unique, Frequency.x) %>%
    rename(Frequency.unstim = Frequency.x)
  df_F <- BulkTCR_MergedF[[i]] %>%
    filter(Counts.y != 0) %>%
    select(Unique, Frequency.y) %>%
    rename(Frequency.FLC = Frequency.y)
  df_U <- merge(df_U, df_B, by = "Unique", all.x = TRUE)
  df_U <- merge(df_U, df_F, by = "Unique", all.x = TRUE)
  BulkTCR_MergedF_final[[i]] <- df_U
}
rm(df_B, df_U, df_F)

# Now calculate cumulative clonotype frequencies and store in meta_lim
for (i in seq_along(BulkTCR_MergedF_final)) {
  df <- BulkTCR_MergedF_final[[i]]
  sample_name <- names(BulkTCR_MergedF_final)[i]
  
  Frequency[sample_name, 'Unstim Freq'] <- sum(df$Frequency.unstim, na.rm = TRUE)
  Frequency[sample_name, 'FLC Stim Freq'] <- sum(df$Frequency.FLC, na.rm = TRUE)
  Frequency[sample_name, 'Stim Freq'] <- sum(df$Frequency.x, na.rm = TRUE)
}

Frequency <- as.data.frame(Frequency)
Frequency$Sample <- rownames(Frequency)
meta_lim <- merge(meta_lim, Frequency, by = "Sample", all.x = TRUE)

# Identify frequencies for the treatment-enriched clonotypes
meta_lim <- bind_rows(
  meta_lim,
  data.frame(Sample = "Treatment-enriched", stringsAsFactors = FALSE)
  )

meta_lim$`Unstim Freq`[meta_lim$Sample == "Treatment-enriched"] <-
  sum(BulkTCR_MergedP_sig$Frequency.x, na.rm = TRUE)

#### E.2 Number of Clonotypes ####
# Adding additional columns
meta_lim <- meta_lim %>%
  mutate(`# Clonotypes` = NA_integer_) %>%
  relocate(`# Clonotypes`, .after = Unobs_Freq)

# Now count the number of mKRAS clonotypes
for (i in seq_along(BulkTCR_MergedF_final)) {
  df <- BulkTCR_MergedF_final[[i]]
  sample_name <- names(BulkTCR_MergedF_final)[i]
  meta_lim$`# Clonotypes`[meta_lim$Sample == sample_name] <- nrow(df)
}
rm(df, sample_name)

# Now count the number of treatment-enriched clonotypes
meta_lim$`# Clonotypes`[meta_lim$Sample == "Treatment-enriched"] <- nrow(BulkTCR_MergedP_sig)

#### E.3 Clonality ####
# Create a function for calculating clonality
calculate_clonality <- function(df) {
  freq <- df$Counts / sum(df$Counts)
  entropy <- -sum(freq * log(freq))
  clonality <- 1 - (entropy / log(sum(df$Counts)))
  return(clonality)
}

# Calculate clonality for all conditions and store in meta_lim
clonality_calc <- tibble(
  Sample    = names(BulkTCR_unique),
  Clonality = map_dbl(BulkTCR_unique, calculate_clonality)
)

meta_lim <- meta_lim %>%
  left_join(clonality_calc, by = "Sample") %>%
  relocate(Clonality, .after = Unobs_Freq)
rm(clonality_calc)

#### F. Exporting Data ####
# First, reformat
TCRSummary <- list()

for (i in seq_along(BulkTCR_MergedF_final)) {
  sample_name <- names(BulkTCR_MergedF_final)[[i]]
  df_name <- paste0(sample_name)
  
  df <- BulkTCR_MergedF_final[[i]] %>%
    select(CDR3aa, vResolved, jResolved, Frequency.x, Counts.x) %>%
    rename(
      TCRV = vResolved,
      TCRJ = jResolved,
      !!paste0(df_name, "freq") := Frequency.x,
      !!paste0(df_name, "count") := Counts.x
    )
  
  TCRSummary[[sample_name]] <- df 
}
names(TCRSummary) <- peptide_names_no_FLC

# Create a Summary table with the CDR3, VDJ, and frequencies under each condition
merge_function <- function(x, y) {
  merge(x, y, by = c("CDR3aa", "TCRV", "TCRJ"), all = TRUE)
}
Summary_table <- Reduce(merge_function, TCRSummary)
rm(TCRSummary)

# Now create a column displaying which antigens each clonotype is reactive to as well as monoreactivity
freq_columns <- names(Summary_table)[grepl("freq$", names(Summary_table))]

Summary_table <- Summary_table %>%
  mutate(
    Antigen = apply(select(., all_of(freq_columns)), 1, function(row) {
      paste(names(row)[!is.na(row)], collapse = ", ")
    }),
    Antigen = gsub("freq", "", Antigen),
    Monoreactive = !grepl(",", Antigen)
  )
rm(freq_columns)
 
# Now set all NA values to 0
Summary_table[is.na(Summary_table)] <- 0

# Relocating columns
Summary_table <- Summary_table %>% relocate('Antigen', .after = 'TCRJ')
Summary_table <- Summary_table %>% relocate('Monoreactive', .after = 'TCRJ')

# This create a treatment-enriched TCR table
Treatment_enriched_table <- BulkTCR_MergedP %>% filter(Sig)

Treatment_enriched_table <- Treatment_enriched_table %>%
  select(CDR3aa, TCRV = vResolved, TCRJ = jResolved, Freq = Frequency.x, mKRAS_tumor_specific)

#### Final Export ####
write.csv(meta_lim, paste0(DirOutput, Patient_ID, "_Bulk_TCR_Numbers.csv"), row.names=FALSE)
write.csv(Summary_table, paste0(DirOutput, Patient_ID, "_Bulk_TCR_Summary.csv"), row.names=FALSE)
write.csv(Treatment_enriched_table, paste0(DirOutput, Patient_ID, "_Bulk_Treatment_Enriched_Summary.csv"), row.names=FALSE)