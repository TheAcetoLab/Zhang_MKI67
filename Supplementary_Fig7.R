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
# FigS7.G top
#=================================================================================================================#
#read chip-seq data
mark <- "H3K27ac"
file_name <- paste0("file/peak_df.", mark, ".csv")
peak_df <- read.csv(file_name, header = TRUE)

rnaseq_path <- "file/an0035_Guide4 VS Ctrl/Guide4 VS Ctrl.all.csv"
rnaseq_df <- read.csv(rnaseq_path, header = TRUE)

# Connect to Ensembl mouse
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")

# Get ensembl_gene_id for your mouse gene symbols (safe query)
mouse_ids <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "mgi_symbol",
  values = unique(peak_df$geneSymbol),
  mart = ensembl_mouse
)

ortholog_map <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = mouse_ids$ensembl_gene_id,
  mart = ensembl_mouse
)
ortholog_map <- merge(mouse_ids, ortholog_map, by = "ensembl_gene_id")
colnames(ortholog_map) <- c("mouse_ensembl", "mouse_symbol", "human_ensembl", "human_symbol")

# Remove entries with no human match
ortholog_map <- ortholog_map[ortholog_map$human_symbol != "", ]
chipseq_df <- merge(peak_df, ortholog_map, by.x = "geneSymbol", by.y = "mouse_symbol")
rnaseq_mapped <- merge(rnaseq_df, ortholog_map,
                       by.x = "SYMBOL", by.y = "human_symbol")

# Now merge with ChIP-seq results using mouse gene names
merged_df <- merge(peak_df, rnaseq_mapped,
                   by.x = "geneSymbol", by.y = "mouse_symbol")

# Add a significance column to the full dataset
merged_df$Significant <- with(merged_df,
                              ifelse(abs(Fold) > 1 & p.value < 0.05 & abs(log2FoldChange) >= 0.5 & padj < 0.05,
                                     "Significant", "Not Significant")
)

# Subset label genes

# Split data into significant and non-significant for controlled plotting order
nonsig_df <- merged_df %>% filter(Significant == "Not Significant")
sig_df <- merged_df %>% filter(Significant == "Significant")

label_genes <- sig_df %>%
  filter(geneSymbol %in% c(
    "Cd47", "Klf4", "Unc13d",
    "Adam8"
  ))

# Plot
scatter.lm.g4 <- ggplot() +
  # Plot non-significant points first (grey, background)
  geom_point(data = nonsig_df, aes(x = Fold, y = log2FoldChange),
             color = "grey80", alpha = 0.5, size = 2,shape = 16) +
  
  # Then plot significant points (purple, foreground)
  geom_point(data = sig_df, aes(x = Fold, y = log2FoldChange),
             color = "#986DB2", alpha = 0.5, size = 2,shape = 16) +
  
  # Dashed lines for 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(ylim = c(-5, NA)) +
  
  # Highlight Cd47 and Klf4 with border
  geom_point(data = label_genes,
             aes(x = Fold, y = log2FoldChange),
             shape = 21, size = 2.5, stroke = 0.5,
             fill = "#80287D", color = "black") +
  
  geom_text_repel(data = label_genes,
                  aes(x = Fold, y = log2FoldChange, label = geneSymbol),
                  size = 4, color = "black", max.overlaps = 100) +
  
  labs(
    title = "",
    x = "log2 Fold Change (H3K27ac)",
    y = "log2 Fold Change (LM2)"
  ) +
  theme_custom +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# Save plot
ggsave(paste0(dir.figs, "scatter.LM2.G4", ".pdf"), scatter.lm.g4, height = 4, width = 4)


#=================================================================================================================#
# FigS7.G bottom
#=================================================================================================================#
rnaseq_path <- "file/Guide1 VS Ctrl.all.csv"
rnaseq_df <- read.csv(rnaseq_path, header = TRUE)

# Now merge with ChIP-seq results using mouse gene names
merged_df <- merge(peak_df, rnaseq_df,
                   by.x = "geneSymbol", by.y = "SYMBOL")

# Add a significance column to the full dataset
merged_df$Significant <- with(merged_df,
                              ifelse(abs(Fold) > 1 & p.value < 0.05 & abs(log2FoldChange) >= 0.5 & padj < 0.05,
                                     "Significant", "Not Significant")
)

# Subset label genes
label_genes <- merged_df %>%
  filter(geneSymbol %in% c("Cd47", "Klf4"))

# Split data into significant and non-significant for controlled plotting order
nonsig_df <- merged_df %>% filter(Significant == "Not Significant")
sig_df <- merged_df %>% filter(Significant == "Significant")
label_genes <- sig_df %>%
  filter(geneSymbol %in% c(
    "Cd47", "Klf4", "Unc13d",
    "Adam8"
  ))

# Plot
scatter.lm.g1 <- ggplot() +
  # Plot non-significant points first (grey, background)
  geom_point(data = nonsig_df, aes(x = Fold, y = log2FoldChange),
             color = "grey80", alpha = 0.5, size = 2,shape = 16) +
  
  # Then plot significant points (purple, foreground)
  geom_point(data = sig_df, aes(x = Fold, y = log2FoldChange),
             color = "#986DB2", alpha = 0.5, size = 2,shape = 16) +
  
  # Dashed lines for 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  #coord_cartesian(ylim = c(-5, NA)) +
  
  # Highlight Cd47 and Klf4 with border
  geom_point(data = label_genes,
             aes(x = Fold, y = log2FoldChange),
             shape = 21, size = 2.5, stroke = 0.5,
             fill = "#80287D", color = "black") +
  
  geom_text_repel(data = label_genes,
                  aes(x = Fold, y = log2FoldChange, label = geneSymbol),
                  size = 4, color = "black", max.overlaps = 100) +
  
  labs(
    title = "",
    x = "log2 Fold Change (H3K27ac)",
    y = "log2 Fold Change (MVT1)"
  ) +
  theme_custom +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# Save plot
ggsave(paste0(dir.figs, "scatter.MVT1.G1", ".pdf"), scatter.lm.g1, height = 4, width = 4)

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to files/sessioninfo_detailed.txt

