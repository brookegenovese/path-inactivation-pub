# 12-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * SCRIPT FOR PEPTIDE INTERSECTION ANALYSIS (SHARED AND UNIQUE) *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(ggplot2)
library(ggview)
library(readr)
library(readxl) 


##============================================================================##
# Repeating the process for peptide analysis
##============================================================================##

# --- import msdap data --- #
file.pep <- "./results/msdap_results/2024-12-05_16-30-09/peptide_abundance__filter by group independently.tsv" # peptide-level data
pep <- read.delim(file.pep, check.names = FALSE)

pdata.file.use <- "./results/msdap_results/2024-12-05_16-30-09/samples.xlsx" # metadata 
pdata <- as.data.frame(read_excel(pdata.file.use))
rownames(pdata) <- pdata$sample_id

# --- rename columns for ease --- #
cols_to_rename <- colnames(pep)[grepl("^FL", colnames(pep))]

for (col in cols_to_rename) {
  match_row <- pdata[pdata$sample_id == col, ]
  
  if (nrow(match_row) > 0) {
    new_name <- paste0(match_row$shortname, "_", match_row$group)
    
    colnames(pep)[colnames(pep) == col] <- new_name
  }
}
colnames(pep)
# --- group col by condition --- #

pcontrol_cols <- grep("_CONTROL$", names(pep))
pgamma_cols <- grep("_GAMMA$", names(pep))
pheat56_cols <- grep("_HEAT56$", names(pep))
pheat95_cols <- grep("_HEAT95$", names(pep))
ptrizol_cols <- grep("_TRIZOL$", names(pep))

pep_present_in_control <- rowSums(!is.na(pep[, pcontrol_cols])) > 0
pep_present_in_gamma <- rowSums(!is.na(pep[, pgamma_cols])) > 0
pep_present_in_heat56 <- rowSums(!is.na(pep[, pheat56_cols])) > 0
pep_present_in_heat95 <- rowSums(!is.na(pep[, pheat95_cols])) > 0
pep_present_in_trizol <- rowSums(!is.na(pep[, ptrizol_cols])) > 0


# --- create sets of peptide IDs per tx --- #

peptides_control <- pep$peptide_id[pep_present_in_control] 
peptides_gamma <- pep$peptide_id[pep_present_in_gamma]
peptides_heat56 <- pep$peptide_id[pep_present_in_heat56]
peptides_heat95 <- pep$peptide_id[pep_present_in_heat95]
peptides_trizol <- pep$peptide_id[pep_present_in_trizol] 

peptide_groups <- list(
  control = peptides_control,
  gamma = peptides_gamma,
  heat56 = peptides_heat56,
  heat95 = peptides_heat95,
  trizol = peptides_trizol
) 

peptide_lengths <- sapply(peptide_groups, length)
print(peptide_lengths)

# --- peptide stats table compared to control --- #
control_peptides <- peptide_groups$control
calculate_peptide_stats <- function(group, control) {
  total_shared <- sum(group %in% control)
  percent_shared <- 100 * total_shared / length(control)
  unique_to_method <- sum(!group %in% control)
  
  list(
    total_shared = total_shared,
    percent_shared = percent_shared,
    unique_to_method = unique_to_method
  )
}
peptide_stats <- lapply(names(peptide_groups)[-1], function(method) {
  stats <- calculate_peptide_stats(peptide_groups[[method]], control_peptides)
  c(method = method, stats)
})

peptide_stats_df <- do.call(rbind, lapply(peptide_stats, as.data.frame))
print(peptide_stats_df)


# --- sort peptides by presence in control & not other tx --- #

