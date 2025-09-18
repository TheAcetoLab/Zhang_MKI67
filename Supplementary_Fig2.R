#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

library(Seurat)
library(ggplot2)
library(ggpubr)
library(pagoda2)
library(ggrepel)
library(clustree)
library(DESeq2)
library(SeuratWrappers)
library(dplyr)
library(stringr)
library(sctransform)
library(plyr)
library(ggcorrplot)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scCustomize)
library(viridis)
library(RColorBrewer)
library(scplotter)
library(EnhancedVolcano)
library(patchwork)

#=================================================================================================================#
# Prepare custom theme for plot
#=================================================================================================================#
theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = "black"),
                      axis.text.y = element_text(size = 12.8, color = "black"), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,15,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = "black"), 
                      axis.ticks.y = element_line(size = 0.4, colour = "black"),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour="black"),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = "black", size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = "bold", margin = margin(0,0,5,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      #plot.background = element_rect(colour = NA, fill = NA), 
                      aspect.ratio = 1)

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

#-----------------------------------------------------------------------------------------------------------------#
# read variable needed
#-----------------------------------------------------------------------------------------------------------------#
lm2.mg <- readRDS("files/an0013.filtered.merge.sample2.3.rds")
mv.mg <- readRDS("files/an0026.filtered.merge.sample2.3.rds")
Br.mg <- readRDS("files/an0042.filtered.merge.batch2.rds")

# create metadata without sample number
c.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][1]}))
hx.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][3]}))
h.sme <- paste(c.sme,hx.sme,sep = "-")
lm2.mg$cell.type <- h.sme

#=================================================================================================================#
# Supplementary Fig2A upper
#=================================================================================================================#
df.mg <- data.frame(#tx = mv.mg@reductions$tsne@cell.embeddings[,1],
  ux = mv.mg@reductions$umap@cell.embeddings[,1],
  uy = mv.mg@reductions$umap@cell.embeddings[,2],
  cells = mv.mg$cell.type,
  cluster = mv.mg@active.ident)
df.mg <- df.mg %>% 
  arrange(desc(cells == "MVT1-Primary tumor"),desc(cells == "MVT1-CTC Single"),desc(cells == "MVT1-Lung Metastasis"),
          desc(cells == "MVT1-CTC cluster")
  )
df.mg$alpha_val <- ifelse(df.mg$cells == "MVT1-CTC cluster", 1, 0.6)
ts.all.bd <- ggplot(df.mg, aes(x = ux, y = uy,color = cells)) + 
  geom_point(shape = 16, size = 0.8,alpha=0.5) +
  scale_color_manual(values = c("MVT1-CTC Single" = "#7CAE00","MVT1-CTC cluster" = "#E16B8C", "MVT1-Lung Metastasis" = "#6FA0C6", "MVT1-Primary tumor" = "gray70")) +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme_custom  +theme(legend.title = element_blank())+guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove horizontal gridlines
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))
ggsave(paste0(dir.figs, 'merge_sample_together_umap_update', '.pdf'), ts.all.bd, height = 5, width = 8)

#=================================================================================================================#
# Supplementary Fig2A lower
#=================================================================================================================#

df.mg <- data.frame(#tx = Br.mg@reductions$tsne@cell.embeddings[,1],
  ux = Br.mg@reductions$umap@cell.embeddings[,1],
  uy = Br.mg@reductions$umap@cell.embeddings[,2],
  cells = Br.mg$orig.ident,
  cluster = Br.mg@active.ident)
df.mg <- df.mg %>% 
  arrange(desc(cells == "BR16-2-Primary tumor"),desc(cells == "BR16-2-CTC Single"),desc(cells == "BR16-2-Lung Metastasis"),
          desc(cells == "BR16-2-CTC cluster")
  )

df.mg$alpha_val <- ifelse(df.mg$cells == "BR16-2-CTC cluster", 1, 0.6)

ts.all.bd <- ggplot(df.mg, aes(x = ux, y = uy, color = cells, alpha = alpha_val)) + 
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = c(
    "BR16-2-CTC Single" = "#7CAE00",
    "BR16-2-CTC cluster" = "#E16B8C",
    "BR16-2-Lung Metastasis" = "#6FA0C6",
    "BR16-2-Primary tumor" = "gray70")) +
  scale_alpha_continuous(range = c(0.6, 1), guide = "none") +  # Hide legend for alpha
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme_custom +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black")
  )
ggsave(paste0(dir.figs, 'merge_sample_together_umap_update', '.pdf'), ts.all.bd, height = 5, width = 8)

#=================================================================================================================#
# Supplementary Fig2B
#=================================================================================================================#
av.exp <- AverageExpression(lm2.mg,group.by="cell.type")$SCT
# Ensure the data is numeric and in matrix form before calculating correlation
av.exp <- as.matrix(av.exp)  # Convert to matrix
# Calculate correlation matrix
cor.exp <- cor(av.exp)
# Convert to data frame
cor.exp <- as.data.frame(cor.exp)
h.dist <- pheatmap(cor.exp,
                   col = hcl.colors(50, "Purples", rev = TRUE),  # Purple color scale
                   fontsize_row = 14,  # Set font size for row labels to 14
                   fontsize_col = 14
)
ggsave(paste0(dir.figs, 'heatmap.correlation.all.samples.purple.LM2', '.pdf'), plot=h.dist$gtable, height = 7, width = 7)

#=================================================================================================================#
av.exp <- AverageExpression(mv.mg,group.by="cell.type")$SCT
av.exp <- as.matrix(av.exp)  # Convert to matrix
cor.exp <- cor(av.exp)
cor.exp <- as.data.frame(cor.exp)
h.dist <- pheatmap(cor.exp,
                   #col = hcl.colors(50),
                   fontsize_row = 14,  # Set font size for row labels to 14
)
ggsave(paste0(dir.figs, 'heatmap.correlation.all.samples.MVT1', '.pdf'), plot=h.dist$gtable, height = 7, width = 7)

#=================================================================================================================#
av.exp <- AggregateExpression(Br.mg,group.by="orig.ident")$SCT
# Ensure the data is numeric and in matrix form before calculating correlation
av.exp <- as.matrix(av.exp)  # Convert to matrix
cor.exp <- cor(av.exp)
cor.exp <- as.data.frame(cor.exp)
ordered_samples <- c(
  "BR16-2-Primary tumor",
  "BR16-2-Lung Metastasis",
  "BR16-2-CTC cluster",
  "BR16-2-CTC Single"
)
# Reorder rows and columns of the correlation matrix
cor.exp <- cor.exp[ordered_samples, ordered_samples]
h.dist <- pheatmap(cor.exp,
                   col = hcl.colors(50, "Purples", rev = TRUE),  # Purple color scale
                   fontsize_row = 14,  # Set font size for row labels to 14
                   fontsize_col = 14,
                   #annotation_col = ann_df, 
                   #annotation_colors = ann_colors,
                   #main = 'Cell type correlation using all genes detected.'
)
ggsave(paste0(dir.figs, 'heatmap.correlation.all.samples.Br16', '.pdf'), plot=h.dist$gtable, height = 7, width = 7)

#=================================================================================================================#
# Supplementary Fig2C
#=================================================================================================================#
df.bd <- data.frame(
  cells = lm2.mg$cell.type,
  cluster = lm2.mg@active.ident)

cell.number <- df.bd %>% group_by(cells,cluster) %>% 
  dplyr::summarize(count = n()) %>% group_by(cluster) %>% 
  dplyr::mutate(freq = 100*(count / sum(count)))  %>% # find percent of total
  dplyr::mutate(freq.whole = 100*(count / 1623))  %>% # find percent of total
  dplyr::mutate(sum=sum(count))

cell.number <- ddply(cell.number, "cluster",
                     transform, label_ypos=(109-cumsum(freq)))
# Create an ordering based on LM2-CTC cluster frequency
cluster_order <- cell.number %>%
  filter(cells == "LM2-CTC cluster") %>%
  arrange(desc(freq)) %>%
  pull(cluster)

# Apply the ordering to the factor
cell.number$cluster <- factor(cell.number$cluster, levels = cluster_order)
#plot bar plot
freq.dist <- ggplot(cell.number, aes(x = cluster, y = freq)) + 
  geom_bar(aes(fill = cells), stat = "identity", width = 0.7, alpha = 0.85,show.legend = FALSE) + 
  scale_fill_manual(values = c("LM2-CTC Single" = "#7CAE00","LM2-CTC cluster" = "#E16B8C", "LM2-Lung Metastasis" = "#6FA0C6", "LM2-Primary tumor" = "gray70")) +
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)+theme_custom+ theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  ) + labs(x = NULL, y = NULL)
ggsave(paste0(dir.figs, 'bar.percent.cell.population.update', '.pdf'), freq.dist, height = 8, width = 8) 
#=================================================================================================================#
# Supplementary Fig2D
#=================================================================================================================#
# Separate up and down-regulated genes
cluster9.15.markers <- cluster9.15.markers %>%
  dplyr::group_by(sign(avg_log2FC)) %>%    # Group by sign of y (negative/positive)
  dplyr::arrange(abs(p_val_adj)) %>%    # Arrange by absolute value of y in descending order
  dplyr::arrange(desc(abs(avg_log2FC)))  # Arrange by absolute value of y in descending order

gene_sets <- list(
  Up = cluster9.15.markers %>% filter(avg_log2FC > 2, p_val_adj < 0.001) %>% pull(genes),
  Down = cluster9.15.markers %>% filter(avg_log2FC < -2, p_val_adj < 0.001) %>% pull(genes)
)
#=================================================================================================================#
# uploade above files to https://metascape.org/ for GO term analysis 
down_file <- paste0("file/an0013_", "9.15_","down","/Enrichment_GO/_FINAL_GO.csv")

preprocess_data <- function(file_path, regulation_type) {
  data <- read.csv(file_path, sep=',', header=TRUE) %>%
    dplyr::mutate(regulation = regulation_type) %>%
    dplyr::arrange(LogP) %>%
    head(15) %>%
    dplyr:: mutate(Description = sub("(.)", "\\U\\1", Description, perl=TRUE),
                   Count=X.GeneInGOAndHitList,
                   p.adjust= 10^(-LogP),
                   Gene.ratio = as.numeric(`X.InGO`) / 100,
                   GeneRatio = paste0(X.GeneInGOAndHitList, "/", X.GeneInGO))
  # Handling duplicated Description values
  duplicated_desc <- duplicated(data$Description)
  data$Description[duplicated_desc] <- paste0(data$Description[duplicated_desc], "-2")
  
  # Sort by -LogP (p-value) in descending order for ranking in the plot
  data <- data[order(data$Gene.ratio, decreasing = T),]
  
  # Convert Description to factor with levels based on -LogP ordering
  data$Description <- str_pad(data$Description, 65, side = "left")
  data$Description <- factor(data$Description, levels = data$Description[order(data$Gene.ratio, decreasing = F)])
  data
}
enrich.down <- preprocess_data(down_file, "Down")

enrich.down <- enrich.down %>%
  dplyr::mutate(Group = case_when(
    grepl("matrix|structure|extracellular|encapsulating|NABA", Description, ignore.case = TRUE) ~ "Extracellular Matrix & Structure",
    grepl("vessel|vasculature|morphogenesis|migration|tube", Description, ignore.case = TRUE) ~ "Vasculature & Morphogenesis",
    grepl("development|differentiation|response", Description, ignore.case = TRUE) ~ "Development & Differentiation",
    TRUE ~ "Uncategorized"  # fallback in case any are missed
  ))
go.dot.d <- ggplot(enrich.down, aes(Gene.ratio, Description, color = -LogP, size = X.GeneInGOAndHitList)) +
  geom_point() +
  #scale_color_viridis_c(option = "D", direction = -1, trans = 'log10') +
  #scale_color_gradient(high = "#ecd0ff", low = "#6A0572", trans = 'log10')+
  scale_color_gradient(low = "#7DB9DE", high = "#113285", trans = 'log10')+
  labs(y = '', x = 'GeneRatio',size = "Count") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() + scale_x_continuous(limits = c(0.02, 0.08), expand = c(0.02, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_blank())
#ggsave(paste0(dir.figs, 'dot.go.term.cell.type', '.pdf'), go.dot, height = 5, width = 5.5)
ggsave(paste0(dir.figs, 'meta.go_term.down.9.15.dotplot.arranged', '.pdf'), go.dot.d, height = 5.4, width = 4.5)

#=================================================================================================================#
# Supplementary Fig2E
#=================================================================================================================#
# uploade above files to https://metascape.org/ for GO term analysis 
up_file <- paste0("file/an0013_", "9.15_","up","/Enrichment_GO/_FINAL_GO.csv")
enrich.up <- preprocess_data(up_file, "Up")
# group in to 3 groups for up regulated
enrich.up <- enrich.up %>%
  dplyr::mutate(Group = case_when(
    grepl("spindle|segregation|metaphase|anaphase|chromosome|mitosis", Description, ignore.case = TRUE) ~ "Mitosis & Chromosome Segregation",
    grepl("cell cycle|cell division|cycle process", Description, ignore.case = TRUE) ~ "Cell Cycle Progression",
    grepl("nuclear|fission", Description, ignore.case = TRUE) ~ "Nuclear Division ",
    TRUE ~ "Other"
  ))

go.dot <- ggplot(enrich.up, aes(Gene.ratio, Description, color = -LogP, size = X.GeneInGOAndHitList)) +
  geom_point() +
  scale_color_gradient(low = "#DC9FB4", high = "#D0104C") +
  labs(y = 'GO Term', x = 'Gene ratio', size = "Count") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 0.17), expand = c(0.05, 0)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(),               # Show y-axis title
    axis.text.y = element_text(),                # Show y-axis text
    axis.ticks.y = element_line(),               # Optional: show ticks
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    strip.text.y.left = element_text(angle = 0),
    plot.title = element_blank()
  )

ggsave(paste0(dir.figs, 'meta.go_term.up.9.15.dotplot.with.termName', '.pdf'), go.dot, height = 5.5, width = 8.4)

#=================================================================================================================#
# Supplementary Fig2F
#=================================================================================================================#
genes_to_fetch <- c("PBK", "NMU", "ARL6IP1", "MKI67", "CENPF", "CCNB1", "CEP55", "PTTG1", "ASPM", "CKS2", 
                    "TPX2", "HMGB2", "PRC1", "UBE2S", "BIRC5", "ANLN", "CDKN3", "AURKA", "DLGAP5", "CKAP2")

# Add additional metadata to be fetched
metadata_to_fetch <- c('cell.type', 'integrated_snn_res.1.2')

# Combine genes and metadata into a single vector
vars_to_fetch <- c(genes_to_fetch, metadata_to_fetch)

# Fetch the data
d.top20 <- FetchData(lm2.mg, vars = vars_to_fetch)
# Create a new column based on integrated_snn_res.1.2 values
d.top20$cluster_category <- ifelse(d.top30$integrated_snn_res.1.2 %in% c(9, 15), 'Combined Clusters', 'Other Clusters')

# List of genes
genes <- genes_to_fetch

# Function to calculate average expression and proportion of cells expressing the gene
calculate_metrics <- function(data, gene, group_var) {
  data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      avg_expression = mean(!!sym(gene), na.rm = TRUE),
      proportion_expressing = sum(!!sym(gene) > 0, na.rm = TRUE) / n()
    ) %>%
    mutate(gene = gene, group = !!sym(group_var)) %>%
    dplyr::select(group, gene, avg_expression, proportion_expressing)
}

# Calculate metrics for each gene and grouping variable (cell.type and cluster_category)
results_list <- lapply(genes, function(gene) {
  list(
    calculate_metrics(d.top20, gene, "cell.type"),
    calculate_metrics(d.top20, gene, "cluster_category")
  )
})

# Flatten the list of results and combine into a single data frame
bubble_data <- do.call(rbind, lapply(results_list, function(x) do.call(rbind, x)))
bubble_data$group <- gsub("^LM2-", "", bubble_data$group)  # Remove "LM2-"
bubble_data$group <- tools::toTitleCase(bubble_data$group)  # Capitalize the first letter of each word

# Specify the order of the groups, with "LM2-CTC cluster" first
desired_order <- c("Primary Tumor", "Lung Metastasis", "CTC Single", "CTC Cluster", 'Other Clusters', 'Combined Clusters')
# Convert the 'group' column to a factor with levels in the desired order
bubble_data$group <- factor(bubble_data$group, levels = desired_order)
genes_to_fetch.r <- rev(genes_to_fetch)
bubble_data$gene <- factor(bubble_data$gene, levels = genes_to_fetch.r)
# Create the ggplot
bb_plot <- ggplot(bubble_data, aes(x = group, y = gene, size = proportion_expressing, color = avg_expression)) +
  geom_point(alpha = 0.7, shape = 16) +
  scale_size_continuous(name = "Proportion Expressing") +
  scale_color_gradient(high = "#6A0572", low = "#ecd0ff", name = "Avg Expression") +
  theme_minimal() +
  labs(title = NULL,  # Remove plot title
       x = NULL,      # Remove x-axis label
       y = NULL) +    # Remove y-axis label
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),    # Remove major grid lines
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +    # Add border
  geom_vline(xintercept = 4.5, color = "black", linetype = "dashed")    # Add vertical line after 4th group

ggsave(paste0(dir.figs,'bubble.plot.top20.diff.order', '.pdf'), bb_plot, height = 5.2, width = 4.5)

#=================================================================================================================#
# Supplementary Fig2G
#=================================================================================================================#
# Initialize an empty list to store the markers for each cell type
all_markers <- list()
all_markers[["LM2-CTC cluster"]] <- read.csv('files/an0013_lm2.sample2.3.merge.differentiation1/an0013_an0013_de.marker_LM2-CTC cluster.csv')
all_markers[["LM2-CTC Single"]] <- read.csv('files/an0013_lm2.sample2.3.merge.differentiation1/an0013_an0013_de.marker_LM2-CTC Single.csv')
all_markers[["LM2-Lung Metastasis"]] <- read.csv('files/an0013_lm2.sample2.3.merge.differentiation1/an0013_an0013_de.marker_LM2-Lung Metastasis.csv')
all_markers[["LM2-Primary tumor"]] <-  read.csv('files/an0013_lm2.sample2.3.merge.differentiation1/an0013_an0013_de.marker_LM2-Primary tumor.csv')
# Initialize an empty list to store GO enrichment results
go_results <- list()
go_results_df <- list()
reactome_results <- list()
reactome_results_df <- list()
bind_results_df <- list()
all_markers <- lapply(all_markers, function(df) {
  df %>% filter(avg_log2FC > 2) %>% filter(p_val_adj < 0.001)
})

# Loop through each entry in all_markers
for (marker_name in names(all_markers)) {
  # Extract gene list (assuming gene symbols are in the first column named "X")
  gene_list <- all_markers[[marker_name]]$X
  # Convert gene symbols to Entrez IDs
  entrez_genes <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene = entrez_genes$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "BP",  # Biological Process
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  reactome_result <- enrichPathway(gene = entrez_genes$ENTREZID,
                                   organism = "human", 
                                   pvalueCutoff = 0.05, 
                                   readable = TRUE)
  # Store the results in the go_results list
  go_results[[marker_name]] <- go_enrichment
  # Convert enrichResult to data frame and store in go_results_df list
  go_results_df[[marker_name]] <- as.data.frame(go_enrichment)
  reactome_results[[marker_name]] <- reactome_result
  # Convert enrichResult to data frame and store in go_results_df list
  reactome_results_df[[marker_name]] <- as.data.frame(reactome_result)
  bind_results_df[[marker_name]] <- rbind(go_results_df[[marker_name]][, c("ID", "Description", "GeneRatio", "p.adjust","Count")],reactome_results_df[[marker_name]][, c("ID", "Description", "GeneRatio", "p.adjust","Count")])
}

# Combine all data frames into one
merged_df.initial <- do.call(rbind, bind_results_df)
merged_df.initial <- merged_df.initial %>%
  mutate(cluster = str_extract(rownames(.), "^[^.]+"))
# Extract the top 20 rows from the LM2-CTC cluster list
top15 <- head(bind_results_df[["LM2-CTC cluster"]][order(bind_results_df[["LM2-CTC cluster"]]$p.adjust), ], 15)
top15_ids <- top15$ID

# selected data 
merged_df <- merged_df.initial %>%
  filter(ID %in% top15_ids)
# make the df for plot 
# Copy first 15 rows
template_rows <- merged_df[1:15, ]
# Create version for LM2-Lung Metastasis
lung_metastasis <- template_rows
lung_metastasis$GeneRatio <- NA
lung_metastasis$p.adjust <- NA
lung_metastasis$Count <- NA
lung_metastasis$cluster <- "LM2-Lung Metastasis"
# Create version for LM2-Primary tumor
primary_tumor <- template_rows
primary_tumor$GeneRatio <- NA
primary_tumor$p.adjust <- NA
primary_tumor$Count <- NA
primary_tumor$cluster <- "LM2-Primary tumor"

# Combine all into an expanded data frame
merged_df <- rbind(merged_df, lung_metastasis, primary_tumor)
## for code for loop, return back the R file differential1

# Split the GeneRatio column and convert it to numeric
merged_df$p.adjust <- as.numeric(merged_df$p.adjust)
merged_df$GeneRatio <- as.numeric(sapply(strsplit(as.character(merged_df$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
#merged_df$term <- merged_df$Description

####
# Create a new column 'term' to summarize the Description
merged_df$term <- dplyr::case_when(
  grepl("chromosome|chromatid|segregation|alignment", merged_df$Description, ignore.case = TRUE) ~ "Chromosome segregation",
  grepl("regulation|transition|phase|nuclear division", merged_df$Description, ignore.case = TRUE) ~ "Mitotic cell cycle",
  grepl("spindle|kinetochore|organelle", merged_df$Description, ignore.case = TRUE) ~ "Spindle Organization",
  TRUE ~ NA_character_
)

order_of_descriptions <- rev(unique(top15$Description))  # Reverse the order
# Set Description factor levels in merged_df based on the reversed order from LM2-CTC cluster
merged_df$Description <- factor(merged_df$Description, levels = order_of_descriptions)
# Reorder the levels in one line
merged_df$cluster <- factor(merged_df$cluster, levels = c("LM2-Primary tumor", "LM2-Lung Metastasis", "LM2-CTC Single", "LM2-CTC cluster"))
# Create the dot plot with viridis color scale, reversing the color direction
go.dot <- ggplot(merged_df, aes(x = cluster, y = Description, size = GeneRatio, color = -log(p.adjust))) +
  geom_point() +
  #scale_color_viridis_c(option = "D", direction = -1, trans = 'log10') +
  scale_color_gradient(low = "#ecd0ff", high = "#6A0572", trans = 'log10')+
  labs(title = "",
       x = "",
       size = "Gene Ratio",
       color = "Adjusted p-value") +
  facet_grid(rows = vars(term), scales = "free_y", space = "free_y", switch = "y") +
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
ggsave(paste0(dir.figs, 'dot.go.term.cell.type.order', '.pdf'), go.dot, height = 5.6, width = 4.7)

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to files/sessioninfo_detailed.txt
