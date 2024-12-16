# 09-Dec-2024 ==== author: Brooke Genovese
##============================================================================##
# * SCRIPT FOR MDS PLOT  *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(ggplot2)
library(ggview)
library(readr)
library(readxl) 

# --- import msdap data --- #
file.prot <- "./results/msdap_results/2024-12-05_16-30-09/protein_abundance__filter by group independently.tsv" # protein-level data
pdata.file.use <- "./results/msdap_results/2024-12-05_16-30-09/samples.xlsx" # metadata 


# Generate MDS plots to look at underlying data similarity
##============================================================================##

dat0 <- read.delim(file.prot, check.names = FALSE)
rownames(dat0) <- dat0$protein_id
dat <- select(dat0, starts_with("FL"))
anno <- select(dat0, -starts_with("FL"))
pdata <- as.data.frame(read_excel(pdata.file.use))
rownames(pdata) <- pdata$sample_id

# multidimensional scaling plot (MDS) 
coords <- limma::plotMDS(dat, plot = FALSE)
plotdat <- pdata
plotdat$`MDS Dim 1` <- coords$x
plotdat$`MDS Dim 2` <- coords$y
group_colors <- c(
  "TRIZOL" = "#a00f83",   
  "HEAT95" = "#dd4124",   
  "HEAT56" = "#F88500",
  "GAMMA" = "#BCCF21",  
  "CONTROL" = "#0f85a0"    
  
) 
# Create the base MDS plot with layered points
group_MDS <- ggplot(plotdat, aes(x = `MDS Dim 1`, y = `MDS Dim 2`)) + 
  # Plot all points except "gamma" and "control"
  geom_point(data = subset(plotdat, !(group %in% c("gamma", "control"))), aes(color = group), size = 4, alpha = 0.7) + 
  # Plot "control" group points next
  geom_point(data = subset(plotdat, group == "control"), aes(color = group), size = 4, alpha = 0.7) + 
  # Plot "gamma" group points on top
  geom_point(data = subset(plotdat, group == "gamma"), aes(color = group), size = 4, alpha = 0.7) + 
  # Apply custom colors to each group
  scale_color_manual(values = group_colors) +  
  labs(title = "MDS Plot by Inactivation Method",  # Update title
       x = "MDS Dimension 1",
       y = "MDS Dimension 2",
       color = "Method") +  # Update legend title
  theme_minimal(base_size = 14) +  # Adjust font size
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        legend.position = "right")  # Position the legend to the right

print(group_MDS)
plot_group_MDS  <- group_MDS  + canvas(1800, 1400, units = "px")
save_ggplot(plot_group_MDS, "./figs/group_MDS.png") # save the plot as an image

