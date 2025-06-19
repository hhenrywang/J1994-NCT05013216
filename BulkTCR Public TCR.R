# This imports the summary of each patient's mKRAS clonotypes
# to create a master summary table and identify public clonotypes
# Last edited by Henry Wang on 06.19.25

##### Libraries to load
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(writexl)
library(readxl)
library(readr)

##### A. Data Import #####
# Set the working directory
setwd("C:/Users/henry/Dropbox/Research/Jaffee Lab/J1994/Sequencing/Bulk")

# Set data directories
Dirwd <- getwd()
DirInput <- paste0(Dirwd, "/Input/Bulk_Putative/")
DirOutput <- paste0(Dirwd, "/Output/")
DirTumor <- paste0(Dirwd, "/Data_Bulk/")

# Getting the full path of .csv files in the folder
Dirtsv_file_names <- list.files(path = DirInput, pattern = ".csv$")
Dirtsv_files <- paste0(DirInput, Dirtsv_file_names)

# Gets the metadata file in the Data_Bulk folder
meta = read_excel(paste0(DirInput, "metadata.xlsx"))

# Inputs all of the .csv files into a list of data frames called BulkTCR
BulkTCR <- lapply(Dirtsv_files, function(file) {
  read.csv(file, header = TRUE)
})
names(BulkTCR) <- gsub("\\.csv$", "", Dirtsv_file_names)

# Renames each sample based on the metadata file
for (i in seq_along(BulkTCR)) {
  current_name <- names(BulkTCR)[i] 
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- matching_row$Patient 
    names(BulkTCR)[i] <- new_name
  }
  rm(current_name)
}
rm(matching_row)

#### B. Data Cleanup ####
# Create a new column called Patient with the Patient number
for (i in seq_along(BulkTCR)) {
  BulkTCR[[i]]$Patient <- names(BulkTCR)[i]
}

# Make sure all data frames have the same columns
# Get the column names from the summary for J1994.08
reference_cols <- colnames(BulkTCR[["J1994.08"]])

# Add in missing columns
for (i in seq_along(BulkTCR)) {
  missing_cols <- setdiff(reference_cols, colnames(BulkTCR[[i]]))
  for (col in missing_cols) {
    BulkTCR[[i]][[col]] <- 0
  }
  BulkTCR[[i]] <- BulkTCR[[i]][, reference_cols, drop = FALSE]
}
rm(col)

# Setting logical columns for published TCRs
BulkTCR$Published$Monoreactive <- TRUE

#### C. Data Merging ####
# Merge separate data frames, clean-up patient IDs & TCRV/J
# Create unique clonotype identifiers
TCR_Merged <- bind_rows(BulkTCR) %>%
  relocate(TCRJ, Patient, .after = TCRV) %>%
  mutate(patient_clean = str_remove(Patient, "^J1994\\.")) %>%
  mutate(Patient = case_when(
      patient_clean == "Published" ~ "Published",
      TRUE~ as.character(as.numeric(patient_clean)))) %>%
  select(-patient_clean) %>%
  mutate(TCRV_clean = TCRV %>%
      str_replace("TCRBV(\\d)(?!\\d)", "TCRBV0\\1") %>%
      str_remove("-.*"),
    TCRJ_clean = TCRJ %>%
      str_replace("TCRBJ(\\d)(?!\\d)", "TCRBJ0\\1") %>%
      str_remove("-.*"),
    Unique = str_c(CDR3aa, TCRV_clean, TCRJ_clean, sep = "_")
  ) %>%
  select(-TCRV_clean, -TCRJ_clean)
rm(BulkTCR)

# Split into two objects with and without TCRV/TCRJ
TCR_with_V    <- TCR_Merged %>% filter(TCRV != "")
TCR_missing_V <- TCR_Merged %>% filter(TCRV == "")

#### D. Public Clonotypes ####
# Identify public clonotypes via identical Unique
Public <- TCR_with_V %>%
  group_by(Unique) %>%
  filter(n_distinct(Patient) > 1) %>%
  ungroup()

# Now, for clonotypes without TCRV/TCRJ, match by CDR3aa only
miss_public <- TCR_missing_V %>%
  semi_join(TCR_with_V, by = "CDR3aa")

with_public <- TCR_with_V %>%
  semi_join(TCR_missing_V, by = "CDR3aa") %>%
  anti_join(TCR_missing_V, by = c("CDR3aa", "Patient"))

Public_noV <- bind_rows(miss_public, with_public)

# Now merging both sets
TCR_Public <- bind_rows(Public, Public_noV)
rm(Public, Public_noV, miss_public, with_public, TCR_with_V, TCR_missing_V)

# Annotate TCR_Merged with Public_ID
annotated <- TCR_Merged %>%
  left_join(
    TCR_Public %>%
      group_by(CDR3aa) %>%
      mutate(Public_ID = cur_group_id()) %>%
      ungroup() %>%
      select(Unique, Patient, Public_ID),
    by = c("Unique", "Patient")
  ) %>%
  mutate(Public = !is.na(Public_ID))

non_public <- annotated %>%
  filter(!Public)

# Collapse the public clonotypes into one row
public_collapsed <- annotated %>%
  filter(Public, !is.na(TCRV), TCRV != "") %>%
  group_by(Public_ID) %>%
  summarize(
    CDR3aa  = first(CDR3aa),
    TCRV    = first(TCRV),
    TCRJ    = first(TCRJ),
    Unique  = first(Unique),
    Patient = toString(unique(Patient)),
    Antigen = toString(sort(unique(unlist(str_split(Antigen, ",\\s*"))))),
    Public  = TRUE,
    across(ends_with("freq"),  ~ max(.x, na.rm = TRUE)),
    across(ends_with("count"), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Reassemble TCR_merged
TCR_Merged <- bind_rows(non_public, public_collapsed) %>%
  select(-Public_ID) %>%
  relocate(Public, .after = TCRJ)

rm(annotated, non_public, public_collapsed)

#### E. Merge Tumor Public TCR ####
# Now, adds on additional Public TCRs found within the tumor
# Imports the files first
DirInputTumor <- paste0(Dirwd, "/Input/Bulk_Putative/Tumor/")

Tumor_file_list <- list.files(path = DirInputTumor, pattern = "\\.csv$", full.names = TRUE)

Tumor_file_names <- basename(Tumor_file_list)
Tumor_patient_ids <- substr(Tumor_file_names, 1, 2)

Tumor_Public_TCR <- mapply(function(file, id) {
  df <- read.csv(file)
  df$Patient_tumor <- as.integer(id)
  return(df)
}, Tumor_file_list, Tumor_patient_ids, SIMPLIFY = FALSE)

TCR_Public_Tumor <- bind_rows(Tumor_Public_TCR)
rm(Tumor_Public_TCR)

# Clean TCRV/J and build Unique clonotype ID
TCR_Public_Tumor <- TCR_Public_Tumor %>%
  mutate(
    TCRV_clean = TCRV %>%
      str_replace("TCRBV(\\d)(?!\\d)", "TCRBV0\\1") %>%
      str_replace("-.*", ""),
    
    TCRJ_clean = TCRJ %>%
      str_replace("TCRBJ(\\d)(?!\\d)", "TCRBJ0\\1") %>%
      str_replace("-.*", ""),                            
    
    Unique = paste(CDR3aa, TCRV_clean, TCRJ_clean, sep = "_")
  ) %>%
  select(-TCRV_clean, -TCRJ_clean, -Patient, -Antigen)

# Collapse Tumor TCR rows
tumor_collapsed <- TCR_Public_Tumor %>%
  group_by(Patient_tumor, Unique) %>%
  summarise(
    across(everything(), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Patient = Patient_tumor) %>%
  mutate(Patient = as.character(Patient))

# Identify matched TCRs
matched_master <- TCR_Merged %>%
  semi_join(tumor_collapsed, by = "Unique")

# Merge paired rows
TCR_Public_Tumor <- bind_rows(
  TCR_Public,
  tumor_collapsed,
  matched_master) %>%
  distinct() %>%
  group_by(Unique) %>%
  mutate(Public_ID = cur_group_id()) %>%
  ungroup()

rm(tumor_collapsed, matched_master)

#### F. Exporting Data ####
write.csv(TCR_Public, paste0(DirOutput, "All_TCR_Public.csv"), row.names=FALSE)
write.csv(TCR_Public_Tumor, paste0(DirOutput, "All_TCR_Public_Tumor.csv"), row.names=FALSE)
write.csv(TCR_Merged, paste0(DirOutput, "All_TCR_Merged.csv"), row.names=FALSE)