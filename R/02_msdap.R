# 05-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * MS-DAP ALGO SET UP AND ANALYSIS *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(msdap)
library(readr)
library(readxl) 


# Section 1: Data importation and setting up MS-DAP analysis parameters
##============================================================================##

# --- load output from upstream raw data processor and fasta files --- #

dataset =  msdap::import_dataset_spectronaut(filename = "./data/clean/22oct24_day0.tsv") 
# info: 9761/14685 precursors remain after selecting the 'best' precursor for each modified sequence

# - fasta files deposited in corresponding PRIDE repo
dataset = import_fasta(
  dataset,
  files = c("~/Desktop/Projects/2023-dissertation/2023-inactivationexperiments/data/msdap/gg_rhesusMacaque_tUP06718.fasta",
            "~/Desktop/Projects/2023-dissertation/2023-inactivationexperiments/data/msdap/Universal Contaminant Protein FASTA.fasta")
)

# info: 2421/2421 protein accessions and 1222/1222 protein groups were mapped to provided fasta file(s)

# --- create temp file for metadata --- #
#write_template_for_sample_metadata(dataset, "day0_meta_22oct.xlsx") 

# --- edit template in Excel --- #
# - describe the sample group of each sample in the "group" column
# - add additional columns with any metadata that varies between samples
# - further documentation is in the "instructions" tab within the Excel file

# --- load metadata from file --- #
import_dataset = import_sample_metadata(dataset, filename = "./data/metadata/day0_meta_22oct.xlsx")

# --- statistical contrasts --- #
# - if you don't want to do stats, simply remove or comment this line (e.g. just look at QC report, or maybe your dataset has 1 experimental group only).
contrast_dataset = setup_contrasts(import_dataset, 
                                   contrast_list = list(  c("GAMMA","CONTROL"), 
                                                          c("HEAT56","CONTROL"), 
                                                          c("HEAT95","CONTROL"), 
                                                          c("TRIZOL","CONTROL") 
                                                          ))


# Section 2: Running main MS-DAP algorithm (analysis_quickstart)
##============================================================================##

# --- main functions --- #
model_dataset = analysis_quickstart(
  contrast_dataset,
  filter_min_detect = 2,            
  filter_min_quant = 3,             
  filter_fraction_detect = 0.75,   
  filter_fraction_quant = 0.75,     
  filter_min_peptide_per_prot = 1,  
  filter_by_contrast = TRUE,        
  norm_algorithm = c("vsn", "modebetween_protein"), 
  dea_algorithm = c("deqms"),       
  dea_qvalue_threshold = 0.01,                      
  dea_log2foldchange_threshold = NA,                
  output_qc_report = TRUE,                          
  output_abundance_tables = TRUE,                   
  output_dir = "./results/msdap_results",                    
  output_within_timestamped_subdirectory = TRUE
) 

# print a short summary of results at the end
print_dataset_summary(model_dataset)

# optional plot to examine variance 
plot_variance_explained(model_dataset)

