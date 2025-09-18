# Load libraries
library(RUVSeq)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(xlsx)
#library(calibrate)
library(biomaRt)
library(ggrepel)
#library(devtools)
#library(ggpubr)
library(reshape2)
library(gplots)
library(org.Hs.eg.db)
#library(topGO)
library(clusterProfiler)
library(tidyr)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(tidyverse)
library(doParallel)
library(ggbeeswarm)
library(pheatmap)
library(gridExtra)
library(viridis)
library(grid)
library(scran)
library(venn)
library(ggVennDiagram)

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
# Supplementary  Fig6A
#=================================================================================================================#
# Define your comparison groups
groups <- list(
  LM2.featurecount  = c("Guide4", "Guide6"),
  MVT1.featurecount = c("Guide1", "Guide2")
)
# Base path to input files
# Initialize storage
dge_results <- list()

# Load differential expression results
for (group_name in names(groups)) {
  guides <- groups[[group_name]]
  for (guide in guides) {
    file_path <- file.path(paste0("file/an0035_", guide, " VS Ctrl/", guide, " VS Ctrl.DGE.csv"))
    key_name <- paste(group_name, guide, sep = "_")
    
    if (file.exists(file_path)) {
      dge_data <- read.csv(file_path)
      dge_results[[key_name]] <- dge_data
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}

# Step 1: Create combined up/down gene list
gene_list_combined <- list()
for (set_name in names(dge_results)) {
  dge <- dge_results[[set_name]]
  gene_list_combined[[paste0(set_name, "_up")]] <- unique(dge$SYMBOL[dge$diff == "up"])
  gene_list_combined[[paste0(set_name, "_down")]] <- unique(dge$SYMBOL[dge$diff == "down"])
}

# Step 2: Clean invalid gene symbols
clean_genes <- function(genes) {
  genes[!grepl("^ENSMUSG00", genes) & genes != "" & !is.na(genes)]
}
gene_list_combined <- lapply(gene_list_combined, clean_genes)

# Step 3: Mouse to human conversion for LM2 and MVT1 sets
mouse_sets <- grep("^LM2|^MVT1", names(gene_list_combined), value = TRUE)
mouse_genes <- unique(unlist(gene_list_combined[mouse_sets]))

mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia")
conversion_table <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = mouse_genes,
  mart = mouse_mart
)
mouse_to_human <- setNames(conversion_table$hsapiens_homolog_associated_gene_name,
                           conversion_table$external_gene_name)

# Replace mouse symbols with human homologs
for (set_name in mouse_sets) {
  original_genes <- gene_list_combined[[set_name]]
  mapped_genes <- mouse_to_human[original_genes]
  mapped_genes[is.na(mapped_genes)] <- original_genes[is.na(mapped_genes)]
  gene_list_combined[[set_name]] <- unique(mapped_genes)
}

# Step 4: Create binary presence/absence matrix
all_genes <- unique(unlist(gene_list_combined))
binary_df <- data.frame(Gene = all_genes)
for (set_name in names(gene_list_combined)) {
  binary_df[[set_name]] <- binary_df$Gene %in% gene_list_combined[[set_name]]
}

up.dot <- ComplexUpset::upset(
  binary_df,
  intersect = names(gene_list_combined),
  min_size = 40,  # <- filters out intersections with <40 genes
  name = "Gene Intersections",
  base_annotations = list(
    'Intersection size' = intersection_size()
  ))+theme(panel.grid = element_blank(),
           panel.background = element_rect(colour = NA, fill = NA))

ggsave(paste0(dir.figs, 'upset.up.down.filtered.dot', '.pdf'), up.dot, height = 5, width =18 )


#=================================================================================================================#
# Supplementary  Fig6B
#=================================================================================================================#

ego_df <- read_csv("file/Common_Symbols.down.GO_enrichment.LM2.csv")
  # Convert GeneRatio to numeric
ego_df <- ego_df %>%
  mutate(
    GeneRatioNum = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    LogP = -log10(p.adjust)
  )

  # Keep top 15 by GeneRatioNum (or all if fewer)
  num_terms_to_show <- min(15, nrow(ego_df))
  ego_top <- ego_df %>%
    arrange(desc(GeneRatioNum)) %>%
    slice(1:num_terms_to_show)
  # Color palette
    colo <- colorRampPalette(c("#7DB9DE", "#113285"))(100)

  # Plot
  dotplot <- ggplot(ego_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum))) +
    geom_point(aes(size = Count, color = LogP)) +
    scale_color_gradientn(colours = colo, name = expression(-log[10](padj))) +
    scale_size_continuous(name = "Gene Count") +
    labs(y = NULL, x = "Gene Ratio", title = NULL) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1, 4, 1, 1), "lines"),
      strip.background = element_blank()
    )
ggsave(paste0(dir.figs, 'dotplot.LM2.pdf'), dotplot, height = 5, width = 9)

#=================================================================================================================#
# Supplementary  Fig6C
#=================================================================================================================#
ego_df <- read_csv("file/Common_Symbols.down.GO_enrichment.MVT1.csv")
# Convert GeneRatio to numeric
ego_df <- ego_df %>%
  mutate(
    GeneRatioNum = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    LogP = -log10(p.adjust)
  )

# Keep top 15 by GeneRatioNum (or all if fewer)
num_terms_to_show <- min(15, nrow(ego_df))
ego_top <- ego_df %>%
  arrange(desc(GeneRatioNum)) %>%
  slice(1:num_terms_to_show)
# Color palette
colo <- colorRampPalette(c("#7DB9DE", "#113285"))(100)

# Plot
dotplot <- ggplot(ego_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum))) +
  geom_point(aes(size = Count, color = LogP)) +
  scale_color_gradientn(colours = colo, name = expression(-log[10](padj))) +
  scale_size_continuous(name = "Gene Count") +
  labs(y = NULL, x = "Gene Ratio", title = NULL) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 4, 1, 1), "lines"),
    strip.background = element_blank()
  )
ggsave(paste0(dir.figs, 'dotplot.MVT1.pdf'), dotplot, height = 5, width = 9)

#=================================================================================================================#
# Supplementary  Fig6D
#=================================================================================================================#
ego_df_lm2.g4 <- read.csv(file = "files/Guide4 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_lm2.g6 <- read.csv(file = "files/Guide6 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g1 <- read.csv(file = "files/Guide1 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g2 <- read.csv(file = "filesGuide2 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)

# Apply transformation to both dataframes
for (df_name in c("ego_df_lm2.g4", "ego_df_lm2.g6", "ego_df_mvt1.g1","ego_df_mvt1.g2")) {
  df <- get(df_name)
  # Handle zero p.adjust values to avoid -Inf
  df$p.adjust[df$p.adjust == 0] <- 1e-300
  
  # Compute LogP
  df$LogP <- -log10(df$p.adjust)
  
  # Assign the modified dataframe back
  assign(df_name, df)
  
}
#-----------------------------------------------------------------------------------------------------------------#
#subset and combine
# List of dataframes and their labels
df_list <- list(
  "Lm2.KO1"  = ego_df_lm2.g4,
  "lm2.KO2" = ego_df_lm2.g6,
  "mvt1.KO1" = ego_df_mvt1.g1,
  "mvt1.KO2" = ego_df_mvt1.g2
)
filtered_dfs <- list()

# Loop through each dataframe
for (label in names(df_list)) {
  df <- df_list[[label]]
  
  # Some datasets use "GO" instead of "ID"
  if (!"ID" %in% colnames(df)) {
    df$ID <- df$GO
  }
  
  # Filter based on keywords: include "adhesion"/"adhesive", exclude "spreading"/"host"
  df_filtered <- df %>%
    dplyr::select(ID, LogP, GeneRatio, Description) %>%
    dplyr::filter(
      grepl("adhesion|adhesive", Description, ignore.case = TRUE),
      !grepl("spreading|host", Description, ignore.case = TRUE)
    ) %>%
    dplyr::mutate(Source = label)
  # Add to list
  filtered_dfs[[label]] <- df_filtered
}

# Combine all filtered dataframes into one
combined_df <- do.call(rbind, filtered_dfs)
combined_df$GeneRatio <- sapply(combined_df$GeneRatio, function(x) {
  if (grepl("/", x)) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  } else {
    as.numeric(x)
  }
})

#-----------------------------------------------------------------------------------------------------------------#
#plot
source_levels <- levels(factor(combined_df$Source))
vline_positions <- seq(1.5, length(source_levels) - 0.5, by = 1)
combined_df$LogP_capped <- pmin(combined_df$LogP, 10)
# go term summary plot with label
go.dot.label <- ggplot(combined_df, aes(x = Source, y = Description, size = GeneRatio, color = -LogP_capped)) +
  geom_point() +
  scale_color_gradient(high = "#ecd0ff", low = "#6A0572", limits = c(-10, 0)) +
  labs(title = "",
       x = "",
       size = "Gene Ratio",
       color = "-Log P.adjusted") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_blank())
ggsave(paste0(dir.figs, 'dot.go.term.cell.type.order.with.label', '.pdf'), go.dot.label, height = 5.3, width = 6.7)

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to file/sessioninfo_detailed.txt
