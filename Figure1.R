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
lm2.mg <- readRDS("10X_single_cell/merged.sample.LM2.rds")
mv.mg <- readRDS("10X_single_cell/merged.sample.MVT1.rds")
Br.mg <- readRDS("10X_single_cell/merged.sample.MVT1.BR16.rds")

# create metadata without sample number
c.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][1]}))
hx.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][3]}))
h.sme <- paste(c.sme,hx.sme,sep = "-")
lm2.mg$cell.type <- h.sme

#=================================================================================================================#
# Fig1B
#=================================================================================================================#
df.mg <- data.frame(
  ux = lm2.mg@reductions$umap@cell.embeddings[,1],
  uy = lm2.mg@reductions$umap@cell.embeddings[,2],
  cells = lm2.mg$cell.type,
  cluster = lm2.mg@active.ident)
df.mg <- df.mg %>% 
  arrange(desc(cells == "LM2-Primary tumor"),desc(cells == "LM2-CTC Single"),desc(cells == "LM2-Lung Metastasis"),
          desc(cells == "LM2-CTC cluster")
  )
ts.all.bd <- ggplot(df.mg, aes(x = ux, y = uy,color = cells)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.6) +
  scale_color_manual(
    values = c("LM2-CTC Single" = "#7CAE00","LM2-CTC cluster" = "#e3242b", "LM2-Lung Metastasis" = "#6FA0C6", "LM2-Primary tumor" = "gray70")) +
  coord_cartesian(xlim = c(min(df.mg$ux) - 5, max(df.mg$ux) + 1), ylim = c(min(df.mg$uy) - 1, max(df.mg$uy) + 1)) +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme_custom  +theme(legend.title = element_blank())+guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove horizontal gridlines
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))
ggsave(paste0(dir.figs, 'merge_sample_together_umap_update', '.pdf'), ts.all.bd, height = 5, width = 8)

#=================================================================================================================#
# Fig1B left
#=================================================================================================================#
df.mg <- data.frame(
  ux = lm2.mg@reductions$umap@cell.embeddings[,1],
  uy = lm2.mg@reductions$umap@cell.embeddings[,2],
  cells = lm2.mg$cell.type,
  cluster = lm2.mg@active.ident)
df.mg <- df.mg %>% 
  arrange(desc(cells == "LM2-Primary tumor"),desc(cells == "LM2-CTC Single"),desc(cells == "LM2-Lung Metastasis"),
          desc(cells == "LM2-CTC cluster")
  )
ts.all.bd <- ggplot(df.mg, aes(x = ux, y = uy,color = cells)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.6) +
  scale_color_manual(
    values = c("LM2-CTC Single" = "#7CAE00","LM2-CTC cluster" = "#E16B8C", "LM2-Lung Metastasis" = "#6FA0C6", "LM2-Primary tumor" = "gray70")) +
  coord_cartesian(xlim = c(min(df.mg$ux) - 5, max(df.mg$ux) + 1), ylim = c(min(df.mg$uy) - 1, max(df.mg$uy) + 1)) +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme_custom  +theme(legend.title = element_blank())+guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),  # Remove horizontal gridlines
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"))
    
ggsave(paste0(dir.figs, 'merge_sample_together_umap_update', '.pdf'), ts.all.bd, height = 5, width = 8)
    

#=================================================================================================================#
# Fig1B Right
#=================================================================================================================#
umap_plot <- CellDimPlot(lm2.mg, group_by = "integrated_snn_res.1.2", reduction = "umap",
                         label = F,theme = "theme_blank", legend.position = "none",pt_alpha = 0.6)
ggsave(paste0(dir.figs, 'umap.by.cluster', '.pdf'), plot = umap_plot, width = 6, height = 6)

#=================================================================================================================#
# Fig1C
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
subset_data <- cell.number[cell.number$cells == "LM2-CTC cluster", ]

# Create a bubble plot
bubble.pct <- ggplot(subset_data, aes(x = sum, y = freq, size = sum,color=freq)) +
  geom_point(alpha = 1) + # Adjust transparency with alpha
  scale_size_continuous(range = c(1, 10)) + # Adjust bubble sizes
  scale_color_gradient(high = "#6A0572", low = "#ecd0ff", name = "% CTC cluster") +
  geom_vline(xintercept = 900, linetype = "dashed", color = "black") + # Vertical line at x=1000
  geom_hline(yintercept = 15, linetype = "dashed", color = "black") + # Horizontal line at y=15
  geom_text(data = subset_data[subset_data$cluster %in% c(9, 15), ],
            aes(label = cluster), nudge_x = 0.05, nudge_y = 0.05, color = "black") +
  labs(title = "",
       x = "Cell count",
       y = "% CTC cluster",
       size = "Cell count") +
  theme_custom + theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.border = element_blank(), # Remove all panel borders
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.line = element_line(color = "black"), # Add axis line only on bottom and left
    axis.ticks = element_line(color = "black"), # Ensure ticks are visible
    axis.ticks.length = unit(0.2, "cm"), # Adjust tick length
    axis.line.x = element_line(color = "black", size = 0.5), # Bottom axis line
    axis.line.y = element_line(color = "black", size = 0.5)  # Left axis line
  )
ggsave(paste0(dir.figs,'bubble.pct.cluster.all', '.pdf'), bubble.pct, height = 5, width = 5)

#=================================================================================================================#
# Fig1D
#=================================================================================================================#
lm2.mg$combined_cluster <- ifelse(
  lm2.mg$integrated_snn_res.1.2 %in% c(9, 15),
  "combined cluster",
  "other clusters"
)
pseudo_ifnb <- AggregateExpression(lm2.mg, assays = "SCT", return.seurat = T, group.by = c("combined_cluster"))
print(pseudo_ifnb@assays)
gene_aggre <- as.data.frame(GetAssayData(pseudo_ifnb, assay = "SCT", slot = "data"))
gene_avg <- rowMeans(gene_aggre[, c("combined cluster", "other clusters")])
avg_df <- data.frame(average = gene_avg)

cluster9.15.markers <- FindMarkers(lm2.mg, ident.1 = c(9,15),logfc.threshold = 0,group.by = 'integrated_snn_res.1.2')
cluster9.15.markers <- cluster9.15.markers[grep("\\.", rownames(cluster9.15.markers), invert = TRUE), ]

# Merge the average expression into the markers data frame
cluster9.15.markers$average_expression <- avg_df[rownames(cluster9.15.markers), "average"]

colo1 <- ifelse(
  cluster9.15.markers$p_val_adj < 0.05 & cluster9.15.markers$avg_log2FC > 0.5, "firebrick3",
  ifelse(cluster9.15.markers$p_val_adj < 0.05 & cluster9.15.markers$avg_log2FC < -0.5, "dodgerblue3", "grey")
)

# Assign color category
cluster9.15.markers$color_group <- ifelse(
  cluster9.15.markers$p_val_adj < 0.05 & cluster9.15.markers$avg_log2FC > 0.5, "firebrick3",
  ifelse(cluster9.15.markers$p_val_adj < 0.05 & cluster9.15.markers$avg_log2FC < -0.5, "dodgerblue3", "grey")
)

# Split by color
grey_df <- cluster9.15.markers[cluster9.15.markers$color_group == "grey", ]
red_df  <- cluster9.15.markers[cluster9.15.markers$color_group == "firebrick3", ]
blue_df <- cluster9.15.markers[cluster9.15.markers$color_group == "dodgerblue3", ]

# Plot: grey first (beneath), then blue/red
MA.de <- ggplot() +
  geom_point(data = grey_df, aes(x = average_expression, y = avg_log2FC), color = "grey", alpha = 0.4,shape=16,size=3) +
  geom_point(data = blue_df, aes(x = average_expression, y = avg_log2FC), color = "dodgerblue3", alpha = 0.4, shape=16,size=3) +
  geom_point(data = red_df, aes(x = average_expression, y = avg_log2FC), color = "firebrick3", alpha = 0.4, shape=16,size=3) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  labs(
    x = "Average of SCT Normalized Counts",
    y = "Log2 Fold Change",
    title = ""
  ) +
  scale_x_continuous(limits = c(0, 5.3), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# Save
ggsave(paste0(dir.figs, 'MA.deg.update', '.pdf'), MA.de, height = 4, width = 6)

#=================================================================================================================#
# Fig1E
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
up_file <- paste0("10X_single_cell/metascape_", "9.15_","up","/Enrichment_GO/_FINAL_GO.csv")
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
  scale_color_gradient(low = "#DC9FB4", high = "#D0104C")+ #, trans = 'log10'
  labs(y = '', x = 'Gene ratio',size = "Count") +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() + scale_x_continuous(limits = c(0, 0.17), expand = c(0.05, 0)) +
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
ggsave(paste0(dir.figs, 'meta.go_term.up.9.15.dotplot.arranged', '.pdf'), go.dot, height = 5.5, width = 4.8)

#=================================================================================================================#
# Fig1F
#=================================================================================================================#
my_colors <- c("G1" = "#D4D9E9",# Dark purple
               "S" = "#8C90B5",      # Light purple 696D93"#A7A9C8"
               "G2M" = "#81318A"
) # Purple
umap_plot <- CellDimPlot(lm2.mg, reduction = "umap", group_by = "Phase",
                         theme = "theme_blank") +
  scale_color_manual(values = my_colors)+ 
  scale_fill_manual(values = my_colors)
# save as PDF
ggsave(paste0(dir.figs, 'umap.cell.cycle', '.pdf'), plot = umap_plot, width = 6, height = 6)

#=================================================================================================================#
# Fig1G
#=================================================================================================================#
# Generating data needed for the plot, the plot itself is done by Prism 10 to make sure the style of pie charts remain consistant
# Create a new dataframe with Phase and cell.type
#=================================================================================================================#
#For LM2-NSG
cc.com <- lm2.mg@meta.data[, c('Phase', 'cell.type')]
# Calculate the percentage of each phase within each cell type
percentage_df <- cc.com %>%
  group_by(cell.type, Phase) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
write.csv(percentage_df, file = paste0(dir.files, analysisid, "_", "cell.cycle.phase.LM2", '.csv'), quote=F, row.names=T)
#=================================================================================================================#
#For MVT1-FVB
# Create a new dataframe with Phase and cell.type
cc.com <- mv.mg@meta.data[, c('Phase', 'cell.type')]
# Calculate the percentage of each phase within each cell type
percentage_df <- cc.com %>%
  dplyr::group_by(cell.type, Phase) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
write.csv(percentage_df, file = paste0(dir.files, analysisid, "_", "cell.cycle.phase.MVT1", '.csv'), quote=F, row.names=T)

#=================================================================================================================#
#For Br16-NSG
cc.com <- Br.mg@meta.data[, c('Phase', 'orig.ident')]
# Calculate the percentage of each phase within each cell type
percentage_df <- cc.com %>%
  dplyr::group_by(orig.ident, Phase) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
write.csv(percentage_df, file = paste0(dir.files, analysisid, "_", "cell.cycle.phase.Br16", '.csv'), quote=F, row.names=T)

#=================================================================================================================#
# Fig1H
#=================================================================================================================#
# Plot generated with Prism 10 using data above

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to files/sessioninfo_detailed.txt
