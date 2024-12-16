# 12-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * SCRIPT FOR PROTEIN INTERSECTION ANALYSIS (SHARED AND UNIQUE) *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(ggplot2)
library(ggview)
library(readr)
library(readxl) 

##============================================================================##
# Protein-level analyses of control-tx comparisons
##============================================================================##

# --- import msdap data and format --- #
file.prot <- "./results/msdap_results/2024-12-05_16-30-09/protein_abundance__filter by group independently.tsv" # protein-level data
pdata.file.use <- "./results/msdap_results/2024-12-05_16-30-09/samples.xlsx" # metadata 

dat0 <- read.delim(file.prot, check.names = FALSE)
rownames(dat0) <- dat0$protein_id
dat <- select(dat0, starts_with("FL"))
anno <- select(dat0, -starts_with("FL"))
pdata <- as.data.frame(read_excel(pdata.file.use))
rownames(pdata) <- pdata$sample_id

# How many proteins are shared btw control and each tx?
##============================================================================##

# --- rename columns for ease --- #
dat0 <- dat0 %>%
  mutate(protein_id = stringr::str_extract(protein_id, "^[^;]+"))
rename_map <- pdata %>%
  select(sample_id, shortname, group) %>%
  mutate(new_name = paste0(shortname, "_", group)) %>%
  select(sample_id, new_name) %>%
  deframe()
colnames(dat0) <- colnames(dat0) %>%
  recode(!!!rename_map)
colnames(dat0)

# --- group col by condition --- #

control_cols <- grep("_CONTROL$", names(dat0))
gamma_cols <- grep("_GAMMA$", names(dat0))
heat56_cols <- grep("_HEAT56$", names(dat0))
heat95_cols <- grep("_HEAT95$", names(dat0))
trizol_cols <- grep("_TRIZOL$", names(dat0))

is_present_in_control <- rowSums(!is.na(dat0[, control_cols])) > 0
is_present_in_gamma <- rowSums(!is.na(dat0[, gamma_cols])) > 0
is_present_in_heat56 <- rowSums(!is.na(dat0[, heat56_cols])) > 0
is_present_in_heat95 <- rowSums(!is.na(dat0[, heat95_cols])) > 0
is_present_in_trizol <- rowSums(!is.na(dat0[, trizol_cols])) > 0

# --- what proteins are present in control & not other tx --- #
only_in_control <- is_present_in_control & 
  !is_present_in_gamma & 
  !is_present_in_heat56 & 
  !is_present_in_heat95 & 
  !is_present_in_trizol
proteins_only_in_control <- dat0$protein_id[only_in_control]
print(proteins_only_in_control) # [1]  "A0A1D5QWN5" "A0A5F7ZA98" "A0A1D5QEN9" "F7FWN4"  
length(proteins_only_in_control) # [1] 4

# --- create sets of protein IDs per tx --- #

proteins_control <- dat0$protein_id[is_present_in_control] 
proteins_gamma <- dat0$protein_id[is_present_in_gamma]
proteins_heat56 <- dat0$protein_id[is_present_in_heat56]
proteins_heat95 <- dat0$protein_id[is_present_in_heat95]
proteins_trizol <- dat0$protein_id[is_present_in_trizol] 

protein_groups <- list(
  control = proteins_control,
  gamma = proteins_gamma,
  heat56 = proteins_heat56,
  heat95 = proteins_heat95,
  trizol = proteins_trizol
) 

protein_lengths <- sapply(protein_groups, length)
print(protein_lengths)

# --- protein stats table compared to control --- #
# X% of the tx group's proteins are also found in the ctrl group
control_proteins <- protein_groups$control
calculate_protein_stats <- function(group, control) {
  total_shared <- sum(group %in% control)
  percent_shared <- 100 * total_shared / length(control)
  unique_to_method <- sum(!group %in% control)
  
  list(
    total_shared = total_shared,
    percent_shared = percent_shared,
    unique_to_method = unique_to_method
  )
}
protein_stats <- lapply(names(protein_groups)[-1], function(method) {
  stats <- calculate_protein_stats(protein_groups[[method]], control_proteins)
  c(method = method, stats)
})

protein_stats_df <- do.call(rbind, lapply(protein_stats, as.data.frame))
print(protein_stats_df)

# What proteins are unique to the control group?
##============================================================================##

treatment_proteins <- unique(c(
  protein_groups$gamma,
  protein_groups$heat56,
  protein_groups$heat95,
  protein_groups$trizol
))
unique_to_control <- setdiff(protein_groups$control, treatment_proteins) # 4 prot

print(unique_to_control) # [1] "A0A1D5QWN5" "A0A5F7ZA98" "A0A1D5QEN9" "F7FWN4"  

# What proteins are shared by gamma and heat56?
##============================================================================##

# --- identify unique proteins for gamma and heat56 --- #

unique_to_gamma <- setdiff(protein_groups$gamma, control_proteins)
unique_to_heat56 <- setdiff(protein_groups$heat56, control_proteins)
shared_unique_proteins <- intersect(unique_to_gamma, unique_to_heat56) # 29
length(unique(shared_unique_proteins))
