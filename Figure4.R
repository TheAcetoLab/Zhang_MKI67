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
# Fig4C left
#=================================================================================================================#
mat <- read_csv("Bulk_RNA_seq/pca.mat.lm2.csv")
condition <- unlist(lapply(colnames(mat), function(x){strsplit(as.character(x), '\\_')[[1]][1]}))
metadata_df <- data.frame(condition)
rownames(metadata_df) <- colnames(mat)
condition_colors <- c("Guide4" = "#9E9AC7", "Ctrl" = "black","Guide6" = "#552E8D")
p <- PCAtools::pca(mat, metadata = metadata_df)

bi.500 <- PCAtools::biplot(
  p, 
  colby = "condition", 
  colkey = condition_colors,
  colLegendTitle = '',
  #shape = "condition",
  #shapekey = condition_shape,
  #shapeLegendTitle = '',
  lab = NULL, 
  encircle = TRUE,
  encircleAlpha = 1/3,
  legendPosition = 'right',
  pointSize = 3     
) + theme_custom + #xlim(-70,160) 
  #xlim(-45,55) + 
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.border = element_blank(),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.line = element_line(size=1,color = "black"))

ggsave(paste0(dir.figs, 'PCA.all.condition.top500', '.pdf'), bi.500, height = 5, width = 5)

#=================================================================================================================#
# Fig4C right
#=================================================================================================================#
mat <- read_csv("Bulk_RNA_seq/pca.mat.MVT1.csv")
colnames(mat) <- sapply(strsplit(colnames(mat), "_"), function(x) paste(x[3], x[4], sep = "_"))
condition <- unlist(lapply(colnames(mat), function(x){strsplit(as.character(x), '\\_')[[1]][1]}))
metadata_df <- data.frame(condition)
rownames(metadata_df) <- colnames(mat)
condition_colors <- c("G1" = "#D6A1C9", "Ctrl" = "black","G2" = "#80287F")
p <- PCAtools::pca(mat, metadata = metadata_df)

bi.500 <- PCAtools::biplot(
  p, 
  colby = "condition", 
  colkey = condition_colors,
  colLegendTitle = '',
  #shape = "condition",
  #shapekey = condition_shape,
  #shapeLegendTitle = '',
  lab = NULL, 
  encircle = TRUE,
  encircleAlpha = 1/3,
  legendPosition = 'right',
  pointSize = 3     
) + theme_custom + #xlim(-70,160) 
  #xlim(-45,55) + 
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.border = element_blank(),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.line = element_line(size=1,color = "black"))

ggsave(paste0(dir.figs, 'PCA.all.condition.top500', '.pdf'), bi.500, height = 5, width = 5)


#=================================================================================================================#
# Fig4D
#=================================================================================================================#
#generating data needed for the plot, the plot itself is done by Prism 10 to make sure the style of pie charts remain consistant
# Define groups and corresponding file identifiers
groups <- list(
  LM2.featurecount  = c("Guide4", "Guide6"),
  #LM2         = c("Guide4", "Guide6"),
  MVT1.featurecount = c("Guide1", "Guide2")
)

# Initialize storage for results and summary
dge_results <- list()
summary_list <- list()

# Loop through each group and read the corresponding files
for (group_name in names(groups)) {
  guides <- groups[[group_name]]
  for (guide in guides) {
    # Construct file path
    file_path <- file.path(paste0("Bulk_RNA_seq/", guide, " VS Ctrl.DGE.csv"))
    key_name <- paste(group_name, guide, sep = "_")
    
    if (file.exists(file_path)) {
      # Read the CSV
      dge_data <- read.csv(file_path)
      
      # Save raw data
      dge_results[[key_name]] <- dge_data
      
      # Summarize up/down genes
      up_count <- sum(dge_data$diff == "up", na.rm = TRUE)
      down_count <- sum(dge_data$diff == "down", na.rm = TRUE)
      
      # Store summary
      summary_list[[key_name]] <- data.frame(
        Comparison = key_name,
        Up = up_count,
        Down = down_count
      )
    } else {
      warning(paste("File not found:", file_path))
    }
  }
}

# Combine all summaries into one data frame
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL
write.csv(summary_df, file=paste0(dir.files,"gene.count_separate_for_prism.csv"), row.names = FALSE)

#=================================================================================================================#
# Fig4E
#=================================================================================================================#
ego_df_lm2.g4 <- read.csv(file = "Bulk_RNA_seq/Guide4 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_lm2.g6 <- read.csv(file = "Bulk_RNA_seq/Guide6 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g1 <- read.csv(file = "Bulk_RNA_seq/Guide1 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g2 <- read.csv(file = "Bulk_RNA_seq/Guide2 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)

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
go.dot <- ggplot(combined_df, aes(x = Source, y = Description, size = GeneRatio, color = -LogP_capped)) +
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
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_blank())
ggsave(paste0(dir.figs, 'dot.go.term.cell.type.order.new', '.pdf'), go.dot, height = 4.5, width = 3)

#=================================================================================================================#
# Fig4H
#=================================================================================================================#
# select shared downregulated in at least one KO group for LM2 and MVT1 

# Define groups and paths
groups <- list(
  LM2.featurecount  = c("Guide4", "Guide6"),
  MVT1.featurecount = c("Guide1", "Guide2")
)
dge_results <- list()

# Step 1: Load DE results from all files
for (group_name in names(groups)) {
  for (guide in groups[[group_name]]) {
    file_path <- file.path(paste0("Bulk_RNA_seq/", guide, " VS Ctrl.DGE.csv"))
    key_name <- paste(group_name, guide, sep = "_")
    
    if (file.exists(file_path)) {
      df <- read.csv(file_path)
      dge_results[[key_name]] <- df
    }
  }
}

# Split into human and mouse
lm2_results <- dge_results[grep("^LM2", names(dge_results))]
mvt1_results <- dge_results[grep("^MVT1", names(dge_results))]

# Get downregulated genes in each dataset
get_down <- function(df) {
  df %>% filter(!is.na(SYMBOL), padj < 0.05, log2FoldChange < 0) %>% pull(SYMBOL) %>% unique()
}

lm2_down_list <- lapply(lm2_results, get_down)
mvt1_down_list <- lapply(mvt1_results, get_down)

# Convert mouse gene symbols to human
mouse_symbols <- unique(unlist(mvt1_down_list))
mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
conversion_table <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = mouse_symbols,
  mart = mouse_mart
)
conversion_table <- conversion_table[conversion_table$hsapiens_homolog_associated_gene_name != "", ]

# Map all mouse gene symbols to human
symbol_map <- setNames(conversion_table$hsapiens_homolog_associated_gene_name,
                       conversion_table$external_gene_name)

mvt1_down_list_human <- lapply(mvt1_down_list, function(x) {
  mapped <- symbol_map[x]
  mapped[!is.na(mapped)]
})

# Genes that are downregulated in all files (after conversion)
all_down_genes <- Reduce(intersect, c(lm2_down_list, mvt1_down_list_human))

# significant genes 
get_strong_down <- function(df) {
  df %>% filter(!is.na(SYMBOL), padj < 0.05, log2FoldChange <= -0.5) %>% pull(SYMBOL) %>% unique()
}
lm2_strong <- unique(unlist(lapply(lm2_results, get_strong_down)))
mvt1_strong <- unique(unlist(lapply(mvt1_results, get_strong_down)))
mvt1_strong_human <- unique(symbol_map[mvt1_strong])
mvt1_strong_human <- mvt1_strong_human[!is.na(mvt1_strong_human)]

# Final result: genes in all downregulated sets + strong in both species
shared_de_genes <- intersect(all_down_genes, intersect(lm2_strong, mvt1_strong_human))

#=================================================================================================================#
# select adhesion related genes for validation for LM2 and MVT1

ego_df_lm2.g4 <- read.csv(file = "Bulk_RNA_seq/Guide4 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_lm2.g6 <- read.csv(file = "Bulk_RNA_seq/Guide6 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g1 <- read.csv(file = "Bulk_RNA_seq/Guide1 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)
ego_df_mvt1.g2 <- read.csv(file = "Bulk_RNA_seq/Guide2 VS Ctrl.down.GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE)

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
# Initialize list to store filtered results
filtered_dfs <- list()

# Loop through each dataframe
for (label in names(df_list)) {
  df <- df_list[[label]]
  
  # Some of the lm2 datasets use "GO" instead of "ID"
  if (!"ID" %in% colnames(df)) {
    df$ID <- df$GO
  }
  
  # Filter based on keywords
  df_filtered <- df %>%
    dplyr::select(ID, LogP, GeneRatio, Description, geneID, p.adjust) %>%
    dplyr::mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
    dplyr::filter(
      grepl("cell-cell adhesion|cell adhesion", Description, ignore.case = TRUE),
      GeneRatio > 0.05
    ) %>%
    dplyr::mutate(Source = label)
  
  # Add to list
  filtered_dfs[[label]] <- df_filtered
}
# Combine all filtered dataframes into one
combined_df <- do.call(rbind, filtered_dfs)


#-----------------------------------------------------------------------------------------------------------------#
#  Separate LM2 (human) and MVT1 (mouse) rows
combined_df$Species <- ifelse(grepl("^Lm2", combined_df$Source, ignore.case = TRUE), "human", "mouse")

#  Split geneID strings into individual Entrez IDs
combined_df_long <- combined_df %>%
  separate_rows(geneID, sep = "/") %>%
  rename(entrez_id = geneID)

#  Convert Entrez IDs to symbols using biomaRt
# Human
mart_human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
human_genes <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = unique(combined_df_long$entrez_id[combined_df_long$Species == "human"]),
  mart = mart_human
)

# Mouse
mart_mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_genes <- getBM(
  attributes = c("entrezgene_id", "mgi_symbol"),
  filters = "entrezgene_id",
  values = unique(combined_df_long$entrez_id[combined_df_long$Species == "mouse"]),
  mart = mart_mouse
)

#  Combine and merge back
colnames(human_genes) <- c("entrez_id", "gene_symbol")
colnames(mouse_genes) <- c("entrez_id", "gene_symbol")
gene_symbol_map <- bind_rows(human_genes, mouse_genes)
combined_df_long$entrez_id <- as.character(combined_df_long$entrez_id)
gene_symbol_map$entrez_id <- as.character(gene_symbol_map$entrez_id)

combined_df_long <- left_join(combined_df_long, gene_symbol_map, by = "entrez_id")

# Collapse back into geneID string (but now using gene symbols)
combined_df_converted <- combined_df_long %>%
  group_by(across(-c(entrez_id, gene_symbol))) %>%
  summarise(geneID = paste(unique(gene_symbol[!is.na(gene_symbol)]), collapse = "/"), .groups = "drop")

#  Get unique mouse gene symbols
mouse_genes <- unique(combined_df_long$gene_symbol[combined_df_long$Species == "mouse"])

#  Get mouse → human orthologs
mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
conversion_table <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = mouse_genes,
  mart = mouse_mart
)
colnames(conversion_table) <- c("mouse_symbol", "human_symbol")
conversion_table <- conversion_table[conversion_table$human_symbol != "", ]
# Step 3: Convert mouse genes in combined_df_long to human
combined_df_long <- combined_df_long %>%
  mutate(gene_symbol_converted = case_when(
    Species == "mouse" ~ conversion_table$human_symbol[match(gene_symbol, conversion_table$mouse_symbol)],
    TRUE ~ gene_symbol
  )) %>%
  filter(!is.na(gene_symbol_converted))

#  Find shared human genes
genes_human <- unique(combined_df_long$gene_symbol_converted[combined_df_long$Species == "human"])
genes_mouse <- unique(combined_df_long$gene_symbol_converted[combined_df_long$Species == "mouse"])
shared_genes <- intersect(genes_human, genes_mouse)

#  Keep only rows with shared genes
combined_df_long_shared <- combined_df_long %>%
  filter(gene_symbol_converted %in% shared_de_genes)

adhesion.select <- unique(combined_df_long_shared$gene_symbol_converted)

# Read list of human gene symbols
gene.for.plot <- adhesion.select
gene.for.plot.lower <- tolower(gene.for.plot)

# File paths and labels
count_files <- list(
  lm2_vivo = "file/an0035_raw.counts.LM2.csv",
  mvt1_vivo = ".file/an0035_raw.counts.MVT1.csv"
)

# Normalize counts
normalized_counts_list <- list()

for (name in names(count_files)) {
  raw_counts <- read.csv(count_files[[name]], header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  coldata <- data.frame(row.names = colnames(raw_counts),
                        condition = rep("A", ncol(raw_counts)))
  
  dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                colData = coldata,
                                design = ~1)
  
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_log2 <- log2(norm_counts + 1)
  
  normalized_counts_list[[name]] <- norm_counts_log2
}

# Human genes from your list
human_genes <- unique(final.list$V1)

# Use the mouse mart
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get mouse genes and their human orthologs
conversion_table <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  mart = mouse_mart
)

# Filter: keep only mouse genes whose human ortholog matches your list
conversion_table_filtered <- subset(
  conversion_table,
  hsapiens_homolog_associated_gene_name %in% human_genes
)
# Remove unwanted mouse genes
conversion_table_filtered <- subset(
  conversion_table_filtered,
  !external_gene_name %in% c("H2-DMb2", "Cd59b","Tnfrsf10b")
)
# Lowercase for matching
mouse_genes_converted <- tolower(conversion_table_filtered$external_gene_name)

# Subset each normalized matrix
normalized_subset_list <- list()

for (name in names(normalized_counts_list)) {
  mat <- normalized_counts_list[[name]]
  gene_symbols <- sub(".*_", "", rownames(mat))
  gene_symbols_lower <- tolower(gene_symbols)
  
  if (name == "mvt1_vivo") {
    matched_rows <- gene_symbols_lower %in% mouse_genes_converted
  } else {
    matched_rows <- gene_symbols_lower %in% gene.for.plot.lower
  }
  
  normalized_subset_list[[name]] <- mat[matched_rows, , drop = FALSE]
}

# Final matrices
mat1 <- normalized_subset_list$lm2_vivo
mat2 <- normalized_subset_list$mvt1_vivo

# Build a named vector for mapping: human → mouse
human_to_mouse <- setNames(
  conversion_table_filtered$external_gene_name,
  conversion_table_filtered$hsapiens_homolog_associated_gene_name
)

# Now: apply mapping to mat1 (human gene matrix)
rownames(mat1) <- sub(".*_", "", rownames(mat1))  # Strip any prefix from gene IDs

human_genes <- rownames(mat1)

# Create unified names like "SLK/Slk"
unified_names <- ifelse(human_genes %in% names(human_to_mouse),
                        paste0(human_genes, "/", human_to_mouse[human_genes]),
                        human_genes)

# Update rownames
rownames(mat1) <- unified_names


#  Create mouse → unified mapping
# Remove any NA values in human_to_mouse[human_genes]
valid_mappings <- !is.na(human_to_mouse[human_genes])
human_genes_valid <- human_genes[valid_mappings]
unified_names_valid <- unified_names[valid_mappings]

# Now build mouse → unified mapping cleanly
mouse_to_unified <- setNames(unified_names_valid, human_to_mouse[human_genes_valid])

#  Subset and rename mat2 rows
gene_symbols_mat2 <- sub(".*_", "", rownames(mat2))  # Clean up mat2 gene names
rownames(mat2) <- gene_symbols_mat2

matched_mouse_genes <- intersect(rownames(mat2), names(mouse_to_unified))
mat2 <- mat2[matched_mouse_genes, , drop = FALSE]
rownames(mat2) <- mouse_to_unified[rownames(mat2)]

#  Ensure row order matches mat1
mat2 <- mat2[match(rownames(mat1), rownames(mat2)), , drop = FALSE]

# Sanity check
stopifnot(all(rownames(mat1) == rownames(mat2)))

#  Z-score normalization (row-wise)
mat1_z <- t(scale(t(mat1)))
mat2_z <- t(scale(t(mat2)))

#  Combine matrices
combined_mat_z <- cbind(mat1_z, mat2_z)

# Column split
column_split <- factor(c(rep("LM2 (Human)", ncol(mat1)),
                         rep("MVT1 (Mouse)", ncol(mat2))),
                       levels = c("LM2 (Human)", "MVT1 (Mouse)"))

# Define colors
#heatmap_colors <- colorRamp2(c(min(combined_mat_z), 0, max(combined_mat_z)),
#                            c("dodgerblue3", "white", "firebrick3"))
heatmap_colors <- colorRamp2(c(min(combined_mat_z), 0, max(combined_mat_z)),
                             c("#2B6C93", "white", "#C83521"))

# Define groupings for annotation
column_groups <- c(
  rep("Ctrl", 3),
  rep("Guide4", 3),
  rep("Guide6", 2),
  rep("Ctrl", 3),
  rep("Guide1", 4),
  rep("Guide2", 3)
)

# Create a data frame for annotation
column_annotation_df <- data.frame(Group = column_groups)
rownames(column_annotation_df) <- colnames(combined_mat_z)

# Define colors for groups
group_colors <- c(
  "Ctrl" = "#483C32",
  "Guide4" = "#9E9AC6",
  "Guide6" = "#552E8C",
  "Guide1" = "#D6A1C9",
  "Guide2" = "#80287E"
)

# Create annotation
top_annotation <- HeatmapAnnotation(
  Group = column_annotation_df$Group,
  col = list(Group = group_colors),
  show_annotation_name = TRUE
)

# Build heatmap
ht <- Heatmap(
  combined_mat_z,
  name = "Z-score",
  col = heatmap_colors,
  #column_title = "Human vs Mouse",
  column_split = column_split,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,   # ⛔️ Remove column names
  gap = unit(10, "mm"),
  column_gap = unit(2.5, "mm"),  # <--- This sets the gap between column_split panels
  bottom_annotation = top_annotation
)

pdf(paste0(dir.figs, 'heatmap.adhesion.genes.color2', '.pdf'), width = 6, height = 4.5)  # adjust size as needed
# Draw with borders
draw(ht)
dev.off()

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to file/sessioninfo_detailed.txt

