# 05-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * CLEANING RAW SPECTRONAUT OUTPUT *
##============================================================================##

# --- load libraries --- #
library(dplyr)

# --- import spn datasets --- #
# two experimental reports (called 'exp1 and exp2' for ease)
# these can be found in corresponding PRIDE repo
exp1_raw <- read.delim("raw_proteomics_tsv/20240108_inactivation_set1_full.tsv")
exp2_raw <- read.delim("raw_proteomics_tsv/20240131_inactivation_set2_full.tsv")
spectronaut_raw <- rbind(exp1_raw, exp2_raw) 
head(spectronaut_raw) # make sure you have one col called "PG.ProteinGroups"

# --- filter out contaminants --- #
# the spn report marks all contaminants with "Cont_" in PG.ProteinGroups  
spectronaut_filtered <- spectronaut_raw[!grepl("^Cont_", 
                                               spectronaut_raw$PG.ProteinGroups), ]

# quick check 
cat("Unique values after filtering:\n")
cat(unique(spectronaut_filtered$PG.ProteinGroup), "\n")

# --- obtain day 0 samples only --- #

day0  <- spectronaut_filtered %>%
filter(grepl("([aAfFkK]$)|([aAfFkK].$)", R.FileName))

write_tsv(day0, "./data/clean/22oct24_day0.tsv")
