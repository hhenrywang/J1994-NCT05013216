# This script loads previously identified mKRAS & TE TCRs
# It identifies tumor infiltrating TCRs and additional Public TCRs
# Last edited by Henry Wang on 06.19.25

##### Libraries to load
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)
library(readxl)
library(stringr)

##### 0. Variable Selection ####
# Sets the patient ID that analysis is performed on
Patient_ID <- '29'

# Sets whether outlier TCR is removed for J1994.29
Outlier_TCR <- TRUE

##### A. Loading Data #####
# Set the working directory. Change as needed
setwd("C:/Users/henry/Dropbox/Research/Jaffee Lab/J1994/Sequencing/Bulk")

# Setting input and output directories
Dirwd <- getwd()
DirMeta <- paste0(Dirwd, "/Data_Bulk/")
DirDataBulkTumor <- paste0(Dirwd, "/Data_Bulk/", Patient_ID, "/Tumor/")
DirDataBulk <- paste0(Dirwd, "/Data_Bulk/", Patient_ID, "/")
DirOutput <- paste0(Dirwd, "/Output/", Patient_ID, "/")
DirOutputAll <- paste0(Dirwd, "/Output/")

# Retrieves tumor .tsv file paths
Dirtsv_file_names <- list.files(path = DirDataBulkTumor, pattern = ".tsv$")
Dirtsv_files <- paste0(DirDataBulkTumor, Dirtsv_file_names)

# Retrieves PBMC .tsv file paths
Dirtsv_file_names_periph <- list.files(path = DirDataBulk, pattern = ".tsv$")
Dirtsv_files_periph <- paste0(DirDataBulk, Dirtsv_file_names_periph)

#Gets the metadata file
meta = read_excel(paste0(DirMeta, "metadata.xlsx"))

#Identifies tumor-specific KRAS mutation from metadata file
mKRAS_tumorspecific <- meta$Tumor_specific_mKRAS[
  meta$Patient == Patient_ID & meta$Sample == 'Baseline']

# Now read all of the tumor .tsv files into a list of data frames called BulkTCR
BulkTCR <- lapply(Dirtsv_files, function(file) {
  read.delim(file, header = TRUE, sep = "\t")
})
names(BulkTCR) <- gsub("\\.tsv$", "", Dirtsv_file_names)

# Now read the PBMC .tsv files into a list of data frames called BulkTCRPeriph
BulkTCRPeriph <- lapply(Dirtsv_files_periph, function(file) {
  read.delim(file, header = TRUE, sep = "\t")
})
names(BulkTCRPeriph) <- gsub("\\.tsv$", "", Dirtsv_file_names_periph)

# Renames each tumor sample based on the metadata file
for (i in seq_along(BulkTCR)) {
  current_name <- names(BulkTCR)[i] 
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- matching_row$Sample 
    names(BulkTCR)[i] <- new_name
  }
  rm(current_name, matching_row, new_name)
}

# Renames each peripheral sample based on the metadata file,
for (i in seq_along(BulkTCRPeriph)) {
  current_name <- names(BulkTCRPeriph)[i] 
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- matching_row$Sample 
    names(BulkTCRPeriph)[i] <- new_name
  }
  rm(current_name, matching_row, new_name)
}
rm(meta)

# Keep only "Baseline" and "Unstim" PBMC sequencing
BulkTCRPeriph <- BulkTCRPeriph[names(BulkTCRPeriph) %in% c("Baseline", "Unstim")]

# Checks whether there is an extra Tumor sample. Only relevant for J1994.21
Tumor_extra <- "Tumor_Extra" %in% names(BulkTCR)

# Merging peripheral sequencing into BulkTCR
BulkTCR$Baseline <- BulkTCRPeriph$Baseline
BulkTCR$Unstim <- BulkTCRPeriph$Unstim
rm(BulkTCRPeriph)

#### B. Data Clean-Up ####
# Indicating which columns to keep and rename
columns_to_keep <- c("aminoAcid",
                     "count..templates.reads.",
                     "vMaxResolved",
                     "jMaxResolved")
new_column_names <- c("CDR3aa",
                      "Counts",
                      "TCRV",
                      "TCRJ")

# Only keep columns that have in frame CDR3aa
BulkTCR <- lapply(BulkTCR, function(df){
  df <- df[df$sequenceStatus == "In",]
  return(df)
})

# Keep indicated columns and rename them
BulkTCR <- lapply(BulkTCR, function(df) {
  df <- df %>% select(all_of(columns_to_keep))
  colnames(df) <- new_column_names
  return(df)
})
rm(columns_to_keep, new_column_names)

# Create unique identifiers
add_unique_column <- function(df) {
  df$Unique <- paste(df$CDR3aa, df$TCRV, df$TCRJ, sep = "_")
  return(df)
}
BulkTCR <- lapply(BulkTCR, add_unique_column)

# Merge clonotypes with same unique identifier, but sum the counts
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

# Create meta_lim to store summary information
meta_lim <- data.frame(names(BulkTCR))
colnames(meta_lim) <- 'Sample'
rm(BulkTCR)

# Drop outlier TCR for J1994.29
if (Patient_ID == "29" & Outlier_TCR == TRUE) {
BulkTCR_unique <- lapply(BulkTCR_unique, function(df) {
  df[df$Unique != "CASSVAGGGQETQYF_TCRBV09-01*01_TCRBJ02-05*01", ]
})}

# Save the total number of counts across clonotypes
for(i in seq_along(BulkTCR_unique)){
  Total <- sum(BulkTCR_unique[[i]]$Counts)
  meta_lim$Total_Counts[[i]] <- Total
  meta_lim$Unique_Clonotypes[[i]] <- nrow(BulkTCR_unique[[i]])
  rm(Total)
}
meta_lim$Total_Counts <- as.numeric(as.character(meta_lim$Total_Counts))
meta_lim$Unique_Clonotypes <- as.numeric(as.character(meta_lim$Unique_Clonotypes))

# Create pseudo-frequency if count is 0 for log-log plots
meta_lim$Unobs_Freq <- (3* meta_lim$Total_Counts)^(-1)

# Calculates new clonotype frequencies based on Total_Counts
BulkTCR_unique <- lapply(names(BulkTCR_unique), function(sample_name) {
  df <- BulkTCR_unique[[sample_name]]
  total_counts <- meta_lim$Total_Counts[meta_lim$Sample == sample_name]
  df$Frequency <- df$Counts / total_counts
  df
})
names(BulkTCR_unique) <- meta_lim$Sample

#### C. Tumor Infiltration ####
# Creating separate data frames for each condition
Baseline <- BulkTCR_unique[["Baseline"]]
Unstim <- BulkTCR_unique[["Unstim"]]
Tumor_Pre <- BulkTCR_unique[["Tumor_Pre"]]
Tumor_Post <- BulkTCR_unique[["Tumor_Post"]]

if (Tumor_extra){
  Tumor_Extra<- BulkTCR_unique[["Tumor_Extra"]]
}
rm(BulkTCR_unique)

# Function to merge tumor with PBMC sequencing
# Only keep TCR clones that are found in the Tumor
merge_and_fill <- function(df1, df2) {
  merge(df1, df2, by = c("CDR3aa", "TCRV", "TCRJ", "Unique"), all = TRUE) %>%
    mutate(across(starts_with("Counts"), ~ replace_na(.x, 0)))
}

# Create merged data franes
BulkTCR_MergedB <- merge_and_fill(Tumor_Pre, Baseline)
BulkTCR_MergedU <- merge_and_fill(Tumor_Post, Unstim)
BulkTCR_MergedT <- merge_and_fill(Tumor_Post, Tumor_Pre)

if (Tumor_extra){
  BulkTCR_MergedTadd_Pre <- merge_and_fill(Tumor_Extra, Tumor_Pre)
  BulkTCR_MergedTadd_Post <- merge_and_fill(Tumor_Extra, Tumor_Post)
  BulkTCR_MergedTadd_U <- merge_and_fill(Tumor_Extra, Unstim)
}

# Now setting unidentified frequencies based on Unobs_Freq found in meta_lim
Unob_Freq_Baseline <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Baseline']
Unob_Freq_Unstim <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Unstim']
Unob_Freq_Tumor_Pre <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Tumor_Pre']
Unob_Freq_Tumor_Post <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Tumor_Post']
if (Tumor_extra){
  Unob_Freq_Tumor_Extra <- meta_lim$Unobs_Freq[meta_lim$Sample == 'Tumor_Extra']
  }

# Function to update frequency values based on unidentified frequencies
update_frequency <- function(df, counts_col, freq_col, unobs_value) {
  df[[freq_col]] <- ifelse(df[[counts_col]] == 0, unobs_value, df[[freq_col]])
  return(df)
}

# Now update frequency values for when Counts are 0
BulkTCR_MergedB <- BulkTCR_MergedB %>%
  update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Pre) %>%
  update_frequency("Counts.y", "Frequency.y", Unob_Freq_Baseline)

BulkTCR_MergedU <- BulkTCR_MergedU %>%
  update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Post) %>%
  update_frequency("Counts.y", "Frequency.y", Unob_Freq_Unstim)

BulkTCR_MergedT <- BulkTCR_MergedT %>%
  update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Post) %>%
  update_frequency("Counts.y", "Frequency.y", Unob_Freq_Tumor_Pre)

if (Tumor_extra){
  BulkTCR_MergedTadd_U <- BulkTCR_MergedTadd_U %>%
  update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Extra) %>%
  update_frequency("Counts.y", "Frequency.y", Unob_Freq_Unstim)
  
  BulkTCR_MergedTadd_Pre <- BulkTCR_MergedTadd_Pre %>%
    update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Extra) %>%
    update_frequency("Counts.y", "Frequency.y", Unob_Freq_Tumor_Pre)
  
  BulkTCR_MergedTadd_Post <- BulkTCR_MergedTadd_Post %>%
    update_frequency("Counts.x", "Frequency.x", Unob_Freq_Tumor_Extra) %>%
    update_frequency("Counts.y", "Frequency.y", Unob_Freq_Tumor_Post)
}

# Identifying tumor clonotypes present in baseline blood
Tumor_Pre$Periph_baseline <- Tumor_Pre$Unique %in% Baseline$Unique
Tumor_Post$Periph_unstim <- Tumor_Post$Unique %in% Unstim$Unique

if (Tumor_extra){
  Tumor_Extra$Periph_unstim <- Tumor_Extra$Unique %in% Unstim$Unique
}

#### D. mKRAS TCR clonotypes ####
#### D.1 Identifying TCRs ####
# Importing the .csv file for annotated TCRs for the patient
# Also imports the summary of all mKRAS TCRs across patients
Summary_table <- read.csv(paste0(DirOutput, Patient_ID, "_Bulk_TCR_Summary.csv"), stringsAsFactors = FALSE)
Treatment_enriched_table <- read.csv(paste0(DirOutput, Patient_ID, "_Bulk_Treatment_Enriched_Summary.csv"), stringsAsFactors = FALSE)
Summary_All_Merged_Table <- read.csv(paste0(DirOutputAll, "All_TCR_Merged.csv"), stringsAsFactors = FALSE)

# Remove treatment-enriched TCRs that are also mKRAS TCRs
Treatment_enriched_table <- Treatment_enriched_table %>%
  filter(mKRAS_tumor_specific == FALSE)

# Creating an unique TCR identifier column for the imported data frames
Summary_table <- add_unique_column(Summary_table)
Treatment_enriched_table <- add_unique_column(Treatment_enriched_table)
Summary_All_Merged_Table <- add_unique_column(Summary_All_Merged_Table)

# Create a data frame of mKRAS TCRs from other patients that are reactive to tumor-specific mKRAS
Patient_ID_pattern <- paste0("(^|, )", as.numeric(Patient_ID), "($|, )")

Summary_All_other_patients <- Summary_All_Merged_Table %>%
  filter(!str_detect(Patient, Patient_ID_pattern))

Summary_All_tumor_specific <- Summary_All_Merged_Table %>%
  filter(!str_detect(Patient, Patient_ID_pattern)) %>%
  filter(str_detect(Antigen, mKRAS_tumorspecific))
rm(Patient_ID_pattern)

# Creating a data frame of tumor-specific mKRAS TCRs
Summary_tumor_specific_TCRs <- Summary_table %>%
  filter(str_detect(Antigen, mKRAS_tumorspecific))

# Function to check if `Unique` values exist in a reference list
# If there is no TCRV, then check based on CDR3aa only
tcr_present <- function(df, ref_df, new_col) {
  unique_refs <- ref_df %>%
    filter(TCRV != "") %>%
    pull(Unique)
  
  cdr3_refs <- ref_df %>%
    filter(TCRV == "") %>%
    pull(CDR3aa)
  
  df[[new_col]] <- df$Unique %in% unique_refs | df$CDR3aa %in% cdr3_refs
  return(df)
}

# Define the patient's tumor-specific mKRAS column
mKRAS_TS <- paste0(mKRAS_tumorspecific,"_TCR")

# Now, check for the indicated clonotypes
BulkTCR_MergedT <- BulkTCR_MergedT %>%
  tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
  tcr_present(Summary_All_other_patients, "Public_TCR") %>%
  tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
  tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)

BulkTCR_MergedB <- BulkTCR_MergedB %>%
  tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
  tcr_present(Summary_All_other_patients, "Public_TCR") %>%
  tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
  tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)

BulkTCR_MergedU <- BulkTCR_MergedU %>%
  tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
  tcr_present(Summary_All_other_patients, "Public_TCR") %>%
  tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
  tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)

Tumor_Pre <- Tumor_Pre %>%
  tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
  tcr_present(Summary_All_other_patients, "Public_TCR") %>%
  tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
  tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)

Tumor_Post <- Tumor_Post %>%
  tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
  tcr_present(Summary_All_other_patients, "Public_TCR") %>%
  tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
  tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)

if (Tumor_extra){
  BulkTCR_MergedTadd_Pre <- BulkTCR_MergedTadd_Pre %>%
    tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
    tcr_present(Summary_All_other_patients, "Public_TCR") %>%
    tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
    tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)
  
  BulkTCR_MergedTadd_Post <- BulkTCR_MergedTadd_Post %>%
    tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
    tcr_present(Summary_All_other_patients, "Public_TCR") %>%
    tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
    tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)
  
  BulkTCR_MergedTadd_U <- BulkTCR_MergedTadd_U %>%
    tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
    tcr_present(Summary_All_other_patients, "Public_TCR") %>%
    tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
    tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)
  
  Tumor_Extra <- Tumor_Extra %>%
    tcr_present(Treatment_enriched_table, "Treatment_enriched_TCR") %>%
    tcr_present(Summary_All_other_patients, "Public_TCR") %>%
    tcr_present(Summary_All_tumor_specific, "Public_TCR_tumor_specific") %>%
    tcr_present(Summary_tumor_specific_TCRs, mKRAS_TS)
}

#### D.2 New Public TCRs ####
# Function to identify new Public TCRs
find_public_tcr <- function(df_tumor, df_ref, Patient_ID) {
  df_tumor %>%
    filter(Public_TCR == TRUE) %>%
    select(CDR3aa, TCRV, TCRJ, Unique, Counts.x, Frequency.x, Counts.y, Frequency.y) %>%
    rename(
      Counts_Post = Counts.x,
      Frequency_Post = Frequency.x,
      Counts_Pre = Counts.y,
      Frequency_Pre = Frequency.y
    ) %>%
    merge(df_ref, by = c("CDR3aa", "TCRV", "TCRJ", "Unique")) %>%
    select(CDR3aa,
           TCRV,
           TCRJ,
           Counts_Post,
           Frequency_Post,
           Counts_Pre,
           Frequency_Pre,
           Patient,
           Antigen)
}

# Apply the function
Summary_tumor_public <- find_public_tcr(BulkTCR_MergedT, Summary_All_other_patients, Patient_ID)

# If there additional tumor sequencing
if (Tumor_extra){
  Summary_tumor_add_public <-
    find_public_tcr(BulkTCR_MergedTadd_Pre, Summary_All_other_patients, Patient_ID)
  Summary_tumor_public <- bind_rows(Summary_tumor_public, Summary_tumor_add_public)
}
rm(Summary_tumor_specific_TCRs,
   Summary_All_other_patients,
   Summary_All_tumor_specific)

#### E. Update tumor TCR flags ####
# Now, create a new column that indicates whether a tumor TCR is:
# tumor-specific mKRAS-reactive for expansions or as a Public TCR

add_tcr_all_flag <- function(df, mKRAS_tumorspecific) {
  col_TCR             <- paste0(mKRAS_tumorspecific, "_TCR")
  col_TCR_all         <- paste0(mKRAS_tumorspecific, "_TCR_all")
  
  df <- df %>%
    mutate(
      !!col_TCR_all := coalesce(.data[[col_TCR]], FALSE) |
        coalesce(Public_TCR_tumor_specific, FALSE)
    )
  return(df)
}

Tumor_Pre <- add_tcr_all_flag(Tumor_Pre, mKRAS_tumorspecific)
Tumor_Post <- add_tcr_all_flag(Tumor_Post, mKRAS_tumorspecific)
BulkTCR_MergedT <- add_tcr_all_flag(BulkTCR_MergedT, mKRAS_tumorspecific)
BulkTCR_MergedU <- add_tcr_all_flag(BulkTCR_MergedU, mKRAS_tumorspecific)

if (Tumor_extra){
  Tumor_Extra <- add_tcr_all_flag(Tumor_Extra, mKRAS_tumorspecific)
  BulkTCR_MergedTadd_Pre <- add_tcr_all_flag(BulkTCR_MergedTadd_Pre, mKRAS_tumorspecific)
  BulkTCR_MergedTadd_Post <- add_tcr_all_flag(BulkTCR_MergedTadd_Post, mKRAS_tumorspecific)
  BulkTCR_MergedTadd_U <- add_tcr_all_flag(BulkTCR_MergedTadd_U, mKRAS_tumorspecific)
}

#### F. Statistics####
#### F.1 # Clonotypes & Freq ####
# Calculating total number and frequency of clonotypes for each condition
mKRAS_TS_all <- paste0(mKRAS_TS, "_all")

conditions <- c("Treatment_enriched_TCR",
                "Public_TCR",
                "Public_TCR_tumor_specific",
                mKRAS_TS,
                mKRAS_TS_all)
column_names <- c("Number_Treatment_enriched_TCR",
                "Number_Public_TCR",
                "Number_Public_TCR_tumor_specific",
                paste0("Number_", mKRAS_TS),
                paste0("Number_", mKRAS_TS_all))
freq_column_names <- c("Freq_Treatment_enriched_TCR",
                  "Freq_Public_TCR",
                  "Freq_Public_TCR_tumor_specific",
                  paste0("Freq_", mKRAS_TS),
                  paste0("Freq_", mKRAS_TS_all))

for (i in seq_along(conditions)) {
  # Count total number of clonotypes
  meta_lim[[column_names[i]]][meta_lim$Sample == "Tumor_Pre"] <- sum(Tumor_Pre[[conditions[i]]], na.rm = TRUE)
  meta_lim[[column_names[i]]][meta_lim$Sample == "Tumor_Post"] <- sum(Tumor_Post[[conditions[i]]], na.rm = TRUE)
  
  # Compute total frequency
  meta_lim[[freq_column_names[i]]][meta_lim$Sample == "Tumor_Pre"] <- sum(Tumor_Pre$Counts[Tumor_Pre[[conditions[i]]]], na.rm = TRUE) / meta_lim$Total_Counts[meta_lim$Sample == 'Tumor_Pre']
  meta_lim[[freq_column_names[i]]][meta_lim$Sample == "Tumor_Post"] <- sum(Tumor_Post$Counts[Tumor_Post[[conditions[i]]]], na.rm = TRUE) / meta_lim$Total_Counts[meta_lim$Sample == 'Tumor_Post']
}

# If there is additional tumor sequencing data
if (Tumor_extra){
    for (i in seq_along(conditions)) {
      meta_lim[[column_names[i]]][meta_lim$Sample == "Tumor_Extra"] <- sum(Tumor_Extra[[conditions[i]]], na.rm = TRUE)
      
      meta_lim[[freq_column_names[i]]][meta_lim$Sample == "Tumor_Extra"] <- sum(Tumor_Extra$Counts[Tumor_Extra[[conditions[i]]]], na.rm = TRUE) / meta_lim$Total_Counts[meta_lim$Sample == 'Tumor_Extra']
    }
}

#### F.2 Clonality ####
# Define function to calculate clonality
calculate_clonality <- function(df) {
  freq <- df$Counts / sum(df$Counts)
  entropy <- -sum(freq * log(freq))
  clonality <- 1 - (entropy / log(sum(df$Counts)))
  return(clonality)
}

# Calculate clonality for entire set of tumor TCRs
meta_lim[meta_lim$Sample == "Tumor_Pre", "Clonality"] <- calculate_clonality(Tumor_Pre)
meta_lim[meta_lim$Sample == "Tumor_Post", "Clonality"] <- calculate_clonality(Tumor_Post)

if(Tumor_extra){
  meta_lim[meta_lim$Sample == "Tumor_Extra", "Clonality"] <- calculate_clonality(Tumor_Extra)
}

#### H. Pairwise Plots ####
#### H.1 Tumor - Pre vs Post ####
# Create a column for group assignment
BulkTCR_MergedT <- BulkTCR_MergedT %>%
  mutate(clonotype_group = case_when(
      !!as.name(mKRAS_TS_all) ~ "orange",
      Treatment_enriched_TCR == TRUE ~ "deepskyblue",
      TRUE ~ "black")
  )

# Aggregate by Frequency.x and Frequency.y for each group
plot_data <- BulkTCR_MergedT %>%
  group_by(Frequency.x, Frequency.y, clonotype_group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
  geom_point(
    data = filter(plot_data, clonotype_group == "black"),
    color = "black", alpha = 0.3
  ) +
  geom_point(
    data = filter(plot_data, clonotype_group == "orange"),
    color = "orange", alpha = 0.8, show.legend = FALSE
  ) +
  geom_point(
    data = filter(plot_data, clonotype_group == "deepskyblue"),
    color = "deepskyblue", alpha = 0.8, show.legend = FALSE
  ) +
  geom_vline(xintercept = 2 * Unob_Freq_Tumor_Pre, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2 * Unob_Freq_Tumor_Post, color = "black", linetype = "dashed", linewidth = 1) +
  scale_size_continuous(trans = "log2",
                        name = "# Clonotypes",
                        breaks = c(1, 10, 100, 1000),
                        labels = c("1", "10", "100", "1000"),
                        range = c(2,8))+
  labs(x = "Baseline Tumor Frequency",
    y = "On-treatment Tumor Frequency",
    title = paste0("Tumor TCR - Baseline vs On-treatment (", mKRAS_tumorspecific, ")")) +
  theme(legend.position = "right",
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  coord_equal()

ggsave(filename = paste0(DirOutput, Patient_ID, "_Tumor_", mKRAS_tumorspecific, ".tiff"),
       plot = plot, device = "tiff", width = 8, height = 6)
rm(plot, plot_data)

if (Tumor_extra) {
  # Plot pre vs extra first
  plot_df <- BulkTCR_MergedTadd_Pre
  
  plot_df <- plot_df %>%
    mutate(clonotype_group = case_when(
        !!as.name(mKRAS_TS_all) ~ "orange",
        Treatment_enriched_TCR == TRUE ~ "deepskyblue",
        TRUE ~ "black"
      )
    )
  
  plot_data <- plot_df %>%
    group_by(Frequency.x, Frequency.y, clonotype_group) %>%
    summarise(n = n(), .groups = "drop")
  
  plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
    geom_point(data = filter(plot_data, clonotype_group == "black"),
               color = "black", alpha = 0.3) +
    geom_point(data = filter(plot_data, clonotype_group == "orange"),
               color = "orange", alpha = 0.8, show.legend = FALSE) +
    geom_point(data = filter(plot_data, clonotype_group == "deepskyblue"),
               color = "deepskyblue", alpha = 0.8, show.legend = FALSE) +
    geom_vline(xintercept = 2 * Unob_Freq_Tumor_Pre, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 2 * Unob_Freq_Tumor_Extra, color = "black", linetype = "dashed", linewidth = 1) +
    scale_size_continuous(trans = "log2",
                          name = "# Clonotypes",
                          breaks = c(1, 10, 100, 1000),
                          labels = c("1", "10", "100", "1000"),
                          range = c(2,8))+
    labs(x = "Baseline Tumor Frequency",
      y = "Extra Tumor Frequency",
      title = paste0("Tumor TCR - Extra vs Pre (", mKRAS_tumorspecific, ")")) +
    theme(legend.position = "right",
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    coord_equal()
  
  ggsave(filename = paste0(DirOutput, Patient_ID, "_Tumor_Extra_vs_Pre", mKRAS_tumorspecific, ".tiff"),
         plot = plot, device = "tiff", width = 8, height = 6)

  # Next, plot post vs extra
  plot_df <- BulkTCR_MergedTadd_Post
  
  plot_df <- plot_df %>%
    mutate(clonotype_group = case_when(
        !!as.name(mKRAS_TS_all) ~ "orange",
        Treatment_enriched_TCR == TRUE ~ "deepskyblue",
        TRUE ~ "black"
      )
    )
  
  plot_data <- plot_df %>%
    group_by(Frequency.x, Frequency.y, clonotype_group) %>%
    summarise(n = n(), .groups = "drop")
  
  plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
    geom_point(data = filter(plot_data, clonotype_group == "black"),
               color = "black", alpha = 0.3) +
    geom_point(data = filter(plot_data, clonotype_group == "orange"),
               color = "orange", alpha = 0.8, show.legend = FALSE) +
    geom_point(data = filter(plot_data, clonotype_group == "deepskyblue"),
               color = "deepskyblue", alpha = 0.8, show.legend = FALSE) +
    geom_vline(xintercept = 2 * Unob_Freq_Tumor_Post, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 2 * Unob_Freq_Tumor_Extra, color = "black", linetype = "dashed", linewidth = 1) +
    scale_size_continuous(trans = "log2",
                          name = "# Clonotypes",
                          breaks = c(1, 10, 100, 1000),
                          labels = c("1", "10", "100", "1000"),
                          range = c(2,8))+
    labs(x = "On-treatment Tumor Frequency",
      y = "Extra Tumor Frequency",
      title = paste0("Tumor TCR - Extra vs Post (", mKRAS_tumorspecific, ")")
    ) +
    theme(legend.position = "right",
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    coord_equal()
  
  ggsave(filename = paste0(DirOutput, Patient_ID, "_Tumor_Extra_vs_Post", mKRAS_tumorspecific, ".tiff"),
         plot = plot, device = "tiff", width = 8, height = 6)
  rm(plot, plot_data)
}

#### H.2 Periph vs Tumor ####
# Plot on treatment periph vs tumor
BulkTCR_MergedU <- BulkTCR_MergedU %>%
  mutate(clonotype_group = case_when(
      !!as.name(mKRAS_TS_all) ~ "orange",
      Treatment_enriched_TCR == TRUE ~ "deepskyblue",
      TRUE ~ "black"
    )
  )

# Aggregate by Frequency.x and Frequency.y for each group
plot_data <- BulkTCR_MergedU %>%
  group_by(Frequency.x, Frequency.y, clonotype_group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
  geom_point(
    data = filter(plot_data, clonotype_group == "black"),
    color = "black", alpha = 0.3
  ) +
  geom_point(
    data = filter(plot_data, clonotype_group == "orange"),
    color = "orange", alpha = 0.8, show.legend = FALSE
  ) +
  geom_point(
    data = filter(plot_data, clonotype_group == "deepskyblue"),
    color = "deepskyblue", alpha = 0.8, show.legend = FALSE
  ) +
  geom_vline(xintercept = 2 * Unob_Freq_Unstim, color = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2 * Unob_Freq_Tumor_Post, color = "black", linetype = "dashed", linewidth = 1) +
  scale_size_continuous(trans = "log2",
                        name = "# Clonotypes",
                        breaks = c(1, 10, 100, 1000),
                        labels = c("1", "10", "100", "1000"),
                        range = c(2,8))+
  labs(x = "On-treatment Peripheral Blood Frequency",
    y = "On-treatment Tumor Frequency",
    title = paste0("Tumor TCR - On-treatment Peripheral vs Tumor (", mKRAS_tumorspecific, ")")) +
  theme(legend.position = "right",
        legend.key = element_rect(fill = "white"),
        aspect.ratio = 1) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  coord_equal()

ggsave(filename = paste0(DirOutput, Patient_ID, "_Tumor_Periph_", mKRAS_tumorspecific, ".tiff"),
       plot = plot, device = "tiff", width = 8, height = 6)
rm(plot, plot_data)

if (Tumor_extra) {
  plot_df <- BulkTCR_MergedTadd_U
  
  plot_df <- plot_df %>%
    mutate(clonotype_group = case_when(
        !!as.name(mKRAS_TS_all) ~ "orange",
        Treatment_enriched_TCR == TRUE ~ "deepskyblue",
        TRUE ~ "black"
      )
    )
  
  plot_data <- plot_df %>%
    group_by(Frequency.x, Frequency.y, clonotype_group) %>%
    summarise(n = n(), .groups = "drop")
  
  plot <- ggplot(plot_data, aes(x = Frequency.y, y = Frequency.x, size = n)) +
    geom_point(data = filter(plot_data, clonotype_group == "black"),
               color = "black", alpha = 0.3) +
    geom_point(data = filter(plot_data, clonotype_group == "orange"),
               color = "orange", alpha = 0.8, show.legend = FALSE) +
    geom_point(data = filter(plot_data, clonotype_group == "deepskyblue"),
               color = "deepskyblue", alpha = 0.8, show.legend = FALSE) +
    geom_vline(xintercept = 2 * Unob_Freq_Unstim, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 2 * Unob_Freq_Tumor_Extra, color = "black", linetype = "dashed", linewidth = 1) +
    scale_size_continuous(trans = "log2",
                          name = "# Clonotypes",
                          breaks = c(1, 10, 100, 1000),
                          labels = c("1", "10", "100", "1000"),
                          range = c(2,8))+
    labs(x = "Baseline Peripheral Blood Frequency",
      y = "Extra Tumor Frequency",
      title = paste0("Tumor TCR - On-treatment Peripheral vs Extra (", mKRAS_tumorspecific, ")")) +
    theme(legend.position = "right",
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    coord_equal()
  
  ggsave(filename = paste0(DirOutput, Patient_ID, "_Tumor_Periph_vs_Extra", mKRAS_tumorspecific, ".tiff"),
         plot = plot, device = "tiff", width = 8, height = 6)
  
  rm(plot, plot_data)
}

#### I. Data Export ####
# Reformatting
Summary_table_tumor <- BulkTCR_MergedT

# Removing excess columns
tumor_columns_to_keep <- c("CDR3aa",
                           "TCRV",
                           "TCRJ",
                           "Counts.x",
                           "Frequency.x",
                           "Counts.y",
                           "Frequency.y",
                           "Treatment_enriched_TCR",
                           "Public_TCR",
                           mKRAS_TS,
                           mKRAS_TS_all)
tumor_new_column_names <- c("CDR3aa",
                            "TCRV",
                            "TCRJ",
                            "Counts_Post",
                            "Frequency_Post",
                            "Counts_Pre",
                            "Frequency_Pre",
                            "Treatment_enriched",
                            "Public_TCR",
                            mKRAS_TS,
                            mKRAS_TS_all)

# Removes extraneous columns and then renames the remaining columns
Summary_table_tumor <- Summary_table_tumor[, tumor_columns_to_keep, drop = FALSE]
colnames(Summary_table_tumor) <- tumor_new_column_names

if (Tumor_extra){
  Summary_table_tumor_extra <- BulkTCR_MergedTadd_Pre
  Summary_table_tumor_extra <- Summary_table_tumor_extra[, tumor_columns_to_keep, drop = FALSE]
  colnames(Summary_table_tumor_extra) <- tumor_new_column_names
}

# Filtering for tumor TCRs where at least one TRUE exists across selected columns
logical_cols <- sapply(Summary_table_tumor, is.logical)
Summary_table_tumor <- Summary_table_tumor[rowSums(Summary_table_tumor[, logical_cols], na.rm = TRUE) > 0, ]

if (Tumor_extra){
  logical_cols <- sapply(Summary_table_tumor_extra, is.logical)
  Summary_table_tumor_extra <- Summary_table_tumor_extra[rowSums(Summary_table_tumor_extra[, logical_cols], na.rm = TRUE) > 0, ]
  Summary_table_tumor_extra <- Summary_table_tumor_extra %>%
    rename(
      Counts_Extra = Counts_Post,
      Frequency_Extra = Frequency_Post
    )
  
  Summary_table_tumor <- full_join(Summary_table_tumor, Summary_table_tumor_extra)
  
  Summary_table_tumor <- Summary_table_tumor %>%
    relocate(Counts_Extra, .after = TCRJ) %>%
    relocate(Frequency_Extra, .after = Counts_Extra)
  
  Summary_table_tumor <- Summary_table_tumor %>%
    mutate(across(
      c(Counts_Extra, Frequency_Extra, Counts_Post, Frequency_Post),
      ~ replace_na(.x, 0)
    ))
}

#### Final Export ####
write.csv(meta_lim, paste0(DirOutput, Patient_ID, "_Tumor_TCR_Numbers.csv"), row.names=FALSE)
write.csv(Summary_table_tumor, paste0(DirOutput, Patient_ID, "_Tumor_TCR_Summary.csv"), row.names=FALSE)
write.csv(Summary_tumor_public, paste0(DirOutput, Patient_ID, "_Tumor_New_Public_TCR.csv"), row.names=FALSE)