# 15-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * SCRIPT FOR EXPLORING ABUNDANCE CHANGES FOR 5 PROTEINS *
##============================================================================##

# --- libraries --- #
library(ggplot2)
library(ggview)
library(tidyr)

# --- import msdap data for proteins and format --- #
file.prot <- "./results/msdap_results/2024-12-18_pub/protein_abundance__filter by group independently.tsv" # protein-level data
pdata.file.use <- "./results/msdap_results/2024-12-18_pub/samples.xlsx" # metadata 

dat0 <- read.delim(file.prot, check.names = FALSE)
rownames(dat0) <- dat0$protein_id
pdata <- as.data.frame(read_excel(pdata.file.use))
rownames(pdata) <- pdata$sample_id

# --- rename columns for ease --- #
dat0 <- dat0 %>%
  mutate(protein_id = stringr::str_extract(protein_id, "^[^;]+"))
cols_to_rename <- colnames(dat0)[grepl("^FL", colnames(dat0))]

for (col in cols_to_rename) {
  match_row <- pdata[pdata$sample_id == col, ]
  

  if (nrow(match_row) > 0) {
    new_name <- match_row$group
    colnames(dat0)[colnames(dat0) == col] <- new_name
  }
}
colnames(dat0)

# --- boxplot colors --- #

group_colors <- c(
  "TRIZOL" = "#a00f83",
  "HEAT95" = "#dd4124",
  "HEAT56" = "#F88500",
  "GAMMA" = "#BCCF21",
  "CONTROL" = "#0f85a0"
)

# Albumin: F7HCH0;Q28522
##============================================================================##

albumin <- dat0[dat0$protein_id %in% c("F7HCH0"), ] ## correct to use, manual check
control_alb <- grep("CONTROL$", colnames(albumin), value = TRUE)
gamma_alb   <- grep("GAMMA$", colnames(albumin), value = TRUE)
heat56_alb  <- grep("HEAT56$", colnames(albumin), value = TRUE)
heat95_alb  <- grep("HEAT95$", colnames(albumin), value = TRUE)
trizol_alb  <- grep("TRIZOL$", colnames(albumin), value = TRUE)
all_method_alb <- c(control_alb, gamma_alb, heat56_alb, heat95_alb, trizol_alb)

albumin_long <- albumin %>%
  pivot_longer(
    cols = all_of(all_method_alb),  
    names_to = "method",              
    values_to = "intensity"           
  ) %>%
  filter(!is.na(intensity))           # Remove any rows where intensity is NA


alb <- ggplot(albumin_long, aes(x = method, y = intensity, fill = method)) +  
  geom_boxplot(color = "black") +  
  labs(x = "Inactivation Method", y = "Intensity", title = "Albumin (ALB)") + 
  theme_minimal() +
  theme(
    text = element_text(size = 12), 
    axis.title.x = element_text(size = 10, color = "#3b3b3b", hjust = 0.5, margin = margin(t = 11)),  
    axis.text.x = element_text(hjust = 0.5),  
    axis.title.y = element_text(size = 10, color = "#3b3b3b", margin = margin(r = 10)),  
    legend.position = "none",                
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = "white", colour = NA), 
    plot.background = element_rect(fill = "white"),  
    panel.grid.major = element_line(color = "#D7DBDD"), 
    panel.grid.minor = element_line(color = "#D7DBDD")  
  ) +
  scale_fill_manual(values = group_colors)  # Apply custom fill colors  

print(alb)
plot_alb <- alb  + canvas(1200, 1600, units = "px")
save_ggplot(plot_alb, "./figs/alb_d0_boxplot.png")

# SEROTRANSFERRIN; F7DHR8
##============================================================================##

tf <- dat0[dat0$protein_id %in% c("F7DHR8"), ] ## correct to use, manual check
control_tf <- grep("CONTROL$", colnames(tf), value = TRUE)
gamma_tf   <- grep("GAMMA$", colnames(tf), value = TRUE)
heat56_tf  <- grep("HEAT56$", colnames(tf), value = TRUE)
heat95_tf  <- grep("HEAT95$", colnames(tf), value = TRUE)
trizol_tf  <- grep("TRIZOL$", colnames(tf), value = TRUE)
all_method_tf <- c(control_tf, gamma_tf, heat56_tf, heat95_tf, trizol_tf)

tf_long <- tf %>%
  pivot_longer(
    cols = all_of(all_method_tf),  
    names_to = "method",              
    values_to = "intensity"           
  ) %>%
  filter(!is.na(intensity))           

tf <- ggplot(tf_long, aes(x = method, y = intensity, fill = method)) +  
  geom_boxplot(color = "black") +  
  labs(x = "Inactivation Method", y = "Intensity", title = "Serotransferrin (TF)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12), 
    axis.title.x = element_text(size = 10, color = "#3b3b3b", hjust = 0.5, margin = margin(t = 11)),  
    axis.text.x = element_text(hjust = 0.5),  
    axis.title.y = element_text(size = 10, color = "#3b3b3b", margin = margin(r = 10)),  
    legend.position = "none",                
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = "white", colour = NA), 
    plot.background = element_rect(fill = "white"),  
    panel.grid.major = element_line(color = "#D7DBDD"), 
    panel.grid.minor = element_line(color = "#D7DBDD")  
  ) +
  scale_fill_manual(values = group_colors)  

print(tf)
plot_tf <- tf  + canvas(1200, 1600, units = "px")
save_ggplot(plot_tf, "./figs/tf_d0_boxplot.png")

#  APOLIPOPROTEIN A-1; F7FQM7
##============================================================================##

a1 <- dat0[dat0$protein_id %in% c("F7FQM7"), ] ## correct to use, manual check
control_a1 <- grep("CONTROL$", colnames(a1), value = TRUE)
gamma_a1   <- grep("GAMMA$", colnames(a1), value = TRUE)
heat56_a1  <- grep("HEAT56$", colnames(a1), value = TRUE)
heat95_a1  <- grep("HEAT95$", colnames(a1), value = TRUE)
trizol_a1  <- grep("TRIZOL$", colnames(a1), value = TRUE)
all_method_a1 <- c(control_a1, gamma_a1, heat56_a1, heat95_a1, trizol_a1)

a1_long <- a1 %>%
  pivot_longer(
    cols = all_of(all_method_a1),  
    names_to = "method",              
    values_to = "intensity"           
  ) %>%
  filter(!is.na(intensity))           

a1 <- ggplot(a1_long, aes(x = method, y = intensity, fill = method)) +  
  geom_boxplot(color = "black") +  
  labs(x = "Inactivation Method", y = "Intensity", title = "Apolipoprotein A-1 (APOA1)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12), 
    axis.title.x = element_text(size = 10, color = "#3b3b3b", hjust = 0.5, margin = margin(t = 11)),  
    axis.text.x = element_text(hjust = 0.5),  
    axis.title.y = element_text(size = 10, color = "#3b3b3b", margin = margin(r = 10)),  
    legend.position = "none",                
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = "white", colour = NA), 
    plot.background = element_rect(fill = "white"),  
    panel.grid.major = element_line(color = "#D7DBDD"), 
    panel.grid.minor = element_line(color = "#D7DBDD")  
  ) +
  scale_fill_manual(values = group_colors)  

print(a1)
plot_a1 <- a1  + canvas(1200, 1600, units = "px")
save_ggplot(plot_a1, "./figs/a1_d0_boxplot.png")

#  COMPLEMENT C3; F7EV32
##============================================================================##

c3 <- dat0[dat0$protein_id %in% c("F7EV32"), ] ## correct to use, manual check
control_c3 <- grep("CONTROL$", colnames(c3), value = TRUE)
gamma_c3   <- grep("GAMMA$", colnames(c3), value = TRUE)
heat56_c3  <- grep("HEAT56$", colnames(c3), value = TRUE)
heat95_c3  <- grep("HEAT95$", colnames(c3), value = TRUE)
trizol_c3  <- grep("TRIZOL$", colnames(c3), value = TRUE)
all_method_c3 <- c(control_c3, gamma_c3, heat56_c3, heat95_c3, trizol_c3)

c3_long <- c3 %>%
  pivot_longer(
    cols = all_of(all_method_c3),  
    names_to = "method",              
    values_to = "intensity"           
  ) %>%
  filter(!is.na(intensity))           

c3 <- ggplot(c3_long, aes(x = method, y = intensity, fill = method)) +  
  geom_boxplot(color = "black") +  
  labs(x = "Inactivation Method", y = "Intensity", title = "Complement C3 (C3)") +
  theme_minimal() +
  theme(
      text = element_text(size = 12), 
      axis.title.x = element_text(size = 10, color = "#3b3b3b", hjust = 0.5, margin = margin(t = 11)),  
      axis.text.x = element_text(hjust = 0.5),  
      axis.title.y = element_text(size = 10, color = "#3b3b3b", margin = margin(r = 10)),  
      legend.position = "none",                
      plot.title = element_text(hjust = 0.5), 
      panel.background = element_rect(fill = "white", colour = NA), 
      plot.background = element_rect(fill = "white"),  
      panel.grid.major = element_line(color = "#D7DBDD"), 
      panel.grid.minor = element_line(color = "#D7DBDD")  
  ) +
  scale_fill_manual(values = group_colors)  

print(c3)
plot_c3 <- c3  + canvas(1200, 1600, units = "px")
save_ggplot(plot_c3, "./figs/c3_d0_boxplot.png")

#  PENTRAXIN FAMILY MEMBER / C-REACTIVE PROTEIN; F7DHQ1
##============================================================================##

crp <- dat0[dat0$protein_id %in% c("F7DHQ1"), ] ## correct to use, manual check
control_crp <- grep("CONTROL$", colnames(crp), value = TRUE)
gamma_crp   <- grep("GAMMA$", colnames(crp), value = TRUE)
heat56_crp  <- grep("HEAT56$", colnames(crp), value = TRUE)
heat95_crp  <- grep("HEAT95$", colnames(crp), value = TRUE)
trizol_crp  <- grep("TRIZOL$", colnames(crp), value = TRUE)
all_method_crp <- c(control_crp, gamma_crp, heat56_crp, heat95_crp, trizol_crp)

crp_long <- crp %>%
  pivot_longer(
    cols = all_of(all_method_crp),  
    names_to = "method",              
    values_to = "intensity"           
  ) %>%
  filter(!is.na(intensity))      

# Joining to keep all methods, fill NAs for missing intensities
# Handle NAs by setting intensity to a finite value (for plotting purposes)
all_methods_df <- data.frame(method = unique(all_method_crp))
crp_long <- all_methods_df %>%
  left_join(crp_long, by = "method")
crp_long <- crp_long %>%
  mutate(intensity = ifelse(is.na(intensity), NA, intensity))  
crp_long$method <- factor(crp_long$method, levels = unique(all_method_crp)) 

c <- ggplot(crp_long, aes(x = method, y = intensity, fill = method)) +
  geom_boxplot(color = "black", na.rm = TRUE) + 
  labs(x = "Inactivation Method", y = "Intensity", title = "Pentraxin family member (CRP)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12), 
    axis.title.x = element_text(size = 10, color = "#3b3b3b", hjust = 0.5, margin = margin(t = 11)),  
    axis.text.x = element_text(hjust = 0.5),  
    axis.title.y = element_text(size = 10, color = "#3b3b3b", margin = margin(r = 10)),  
    legend.position = "none",                
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = "white", colour = NA), 
    plot.background = element_rect(fill = "white"),  
    panel.grid.major = element_line(color = "#D7DBDD"), 
    panel.grid.minor = element_line(color = "#D7DBDD")  
  ) +
  scale_fill_manual(values = group_colors)   

print(c) 
plot_crp <- c  + canvas(1200, 1600, units = "px")
save_ggplot(plot_crp, "./figs/crp_d0_boxplot.png")