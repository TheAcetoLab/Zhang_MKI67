# Load libraries
library(RColorBrewer)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(reshape2)
library(gplots)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)
library(doParallel)
library(ggbeeswarm)
library(gridExtra)
library(viridis)
library(grid)
library(scran)
library(venn)
library(EnhancedVolcano)
library(waterfalls)
#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

setwd("file")
analysis <- sub(".R", "", stringr::str_split(rstudioapi::getSourceEditorContext()$path, "/", simplify = T)[,10])
analysisid <- unlist(lapply(analysis, function(x){strsplit(as.character(x), '\\_')[[1]][1]}))
an.desc <- unlist(lapply(analysis, function(x){strsplit(as.character(x), '\\_')[[1]][2]}))
Exp_ID <- stringr::str_split(rstudioapi::getSourceEditorContext()$path, "/", simplify = T)[,8]

dir.create(paste0("plots"))
dir.create(paste0("files"))
dir.create(paste0("plots/", analysisid,"_",an.desc))
dir.create(paste0("files/", analysisid,"_",an.desc))
dir.figs <- paste0("plots/", analysisid, "_", an.desc, "/", analysisid, "_")
dir.files <- paste0("files/", analysisid, "_", an.desc, "/", analysisid, "_")

#=================================================================================================================#
# Supplementary Fig2B
#=================================================================================================================#
selected_df <- read_csv("CRISPR_screen/scrren_raw_count.csv")
library_sizes <- colSums(selected_df)
# Calculate CPM
cpm <- t(t(selected_df) / library_sizes * 1e6)
# Optionally, convert CPM to a data frame and inspect it
cpm_dfx <- as.data.frame(cpm)

# Reorder the columns in the data frame
cpm_dfx$gene <- rownames(cpm_dfx)
mixed_rows <- melt(cpm_dfx, id.vars = "gene")
# Replace variable if gene contains 'Targeting'
mixed_rows$variable <- as.character(mixed_rows$variable)
# Update the variable column where condition is met
mixed_rows$variable[mixed_rows$variable == "CTC_mixed" & grepl("Targeting", mixed_rows$gene)] <- "CTC_non targeting control"
# Optionally, convert back to factor if needed
mixed_rows$variable <- factor(mixed_rows$variable)

den.zc <- ggplot(mixed_rows, aes(x = log2(value+1), color = variable, fill = variable)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = c("#9400D3", "#8071A8","#D30DCA",  "#B6A9CC")) +
  scale_fill_manual(values = c("#9400D3", "#8071A8","#D30DCA",  "#B6A9CC")) +
  labs(title = "", 
       x = "log2(CPM+1)", 
       y = "Density") +
  theme_custom
ggsave(paste0(dir.figs, 'density.plot.all.log2.non.tar', '.pdf'), den.zc, height = 7, width = 7)

#=================================================================================================================#
# Supplementary Fig2D
#=================================================================================================================#
# Load data
guide_summary <- read.table(paste0("CRISPR_screen/ctc_vs_tumor.sgrna_summary.txt"), header = TRUE, fill = TRUE)

# Genes of interest
genes <- c("MKI67")

# Loop over each gene
for (gene in genes) {
  gene_subset <- subset(guide_summary, Gene == gene)
  
  gene_long <- melt(gene_subset, 
                    id.vars = "sgrna", 
                    measure.vars = c("control_count", "treatment_count"),
                    variable.name = "Condition", 
                    value.name = "Count")
  
  gene_long$Count <- as.numeric(gene_long$Count)
  
  p <- ggplot(gene_long, aes(x = Condition, y = Count, color = sgrna, group = sgrna)) +
    geom_line() +
    geom_point(size = 3) +
    scale_y_continuous(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 10000)) +
    scale_x_discrete(labels = c("control_count" = "Tumor", "treatment_count" = "CTC")) +
    labs(title = gene,
         x = "",
         y = "log CPM") +
    scale_color_manual(values = rep("#756bb1", length(unique(gene_long$sgrna)))) +
    theme_custom +
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Save each plot
  ggsave(paste0(dir.figs, "line.plot.", gene, ".same.color.pdf"), p, height = 4, width = 7)
}

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to file/sessioninfo_detailed.txt

