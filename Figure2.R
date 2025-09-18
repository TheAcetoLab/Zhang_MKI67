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
# Fig2B
#=================================================================================================================#
selected_df <- read_csv("files/scrren_raw_count.csv")
selected_df <- selected_df[!grepl("Targeting", rownames(selected_df)), ]

baseline_count <- sum(selected_df$Injection_control_1 > 5)

non_zero_counts      <- sapply(selected_df, function(x) sum(x > 5))
non_zero_percentage  <- (non_zero_counts / baseline_count) * 100

non_zero_percentage_df <- data.frame(
  Condition  = names(non_zero_percentage),
  Percentage = non_zero_percentage
) |>
  dplyr::mutate(Drop = Percentage - 100) |>
  dplyr::mutate(Drop = dplyr::case_when(                # set controls to 0
    Condition %in% c("Injection_control_1", "Injection_control_2") ~ 0,
    TRUE ~ Drop
  ))
non_zero_percentage_df$Drop[2] <- non_zero_percentage_df$Drop[2]-non_zero_percentage_df$Drop[1]
# non targeting control has 34% drop
## ── split CTC_mixed into a 66 % + 34 % pair ────────────────
ctc_row <- dplyr::filter(non_zero_percentage_df, Condition == "CTC_mixed")

ctc_66  <- ctc_row
ctc_34  <- ctc_row
ctc_66$Condition <- "CTC_mixed_66"
ctc_34$Condition <- "CTC_mixed_34"

non_df <- count.new[, c("Tumor_mixed", "CTC_mixed",
                        "Injection_control_1", "Injection_control_2")]
non_df <- non_df[grepl("Targeting", rownames(non_df)), ]

baseline_count <- sum(non_df$Injection_control_1 > 5)

non_zero_counts      <- sapply(non_df, function(x) sum(x > 5))
non_targeting_percentage  <- 100-((non_zero_counts / baseline_count) * 100)

ctc_66$Drop <- -non_targeting_percentage[2]    # 66 % of the original drop
ctc_34$Drop <- non_zero_percentage_df$Drop[2]+non_targeting_percentage[2]    # 34 % of the original drop

## rebuild the data frame (remove old CTC row, add the two new ones)
wf_df <- dplyr::filter(non_zero_percentage_df, Condition != "CTC_mixed") |>
  dplyr::bind_rows(ctc_66, ctc_34)

## ── 🔶 enforce the desired order *in the rows themselves* ──
desired_order <- c("Injection_control_1", "Injection_control_2",
                   "Tumor_mixed", "CTC_mixed_66", "CTC_mixed_34")

wf_df <- wf_df |>
  dplyr::mutate(Condition = factor(Condition, levels = desired_order)) |>
  dplyr::arrange(Condition)                      # 🔶 this does the trick
# ------------------------------------------------------------

## pick colours (same length as desired_order)
fill_vec <- c("#efbbff",  # Injection_control_1
              "#d896ff",  # Injection_control_2
              "#efbbff",  # Tumor_mixed
              "#BF99C7",  # first 66 % of CTC bar
              "#A56FB2")  # remaining 34 % of CTC bar


## plot  (two-decimal labels on bars **and** y-axis)
wf.drop2 <- waterfall(
  wf_df,
  values            = wf_df$Drop,
  labels            = wf_df$Condition,
  rect_text_labels  = sprintf("%.2f", wf_df$Drop),   # <-- bar labels, 2 dp
  fill_by_sign      = FALSE,
  fill_colours      = fill_vec
) +
  theme_minimal() +
  labs(title = "", x = "", y = "Retention rate") +
  scale_y_continuous(
    limits = c(-50, 0)
  ) +
  theme(
    panel.grid         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks         = element_line(colour = "black"),
    axis.line          = element_line(colour = "black")
  )
ggsave(paste0(dir.figs, 'waterfall.plot.no.lung.non.targeting', '.pdf'), wf.drop2, height = 4, width = 4)
#=================================================================================================================#
# Fig2C
#=================================================================================================================#
# Define a function to create and save Enhanced Volcano plots
create_and_save_volcano <- function(input_file, plot_title, output_file) {
  gene_summary <- read.table(paste0(dir.filer, input_file), header = TRUE)
  
  volcano_plot <- EnhancedVolcano(
    gene_summary,
    lab = gene_summary$id,  # Labels for the points
    x = 'lfc',          # Column for log fold change
    y = 'p.value',      # Column for p-value
    title = plot_title,
    shape = 16,
    arrowheads = F,
    subtitle = NULL,
    caption = NULL,
    xlab = 'log2FoldChange',
    ylab = '-Log10(p.value)',
    pCutoff = 0.05,         # Threshold for p-value
    FCcutoff = 1,           # Threshold for fold change
    pointSize = 2.0,        # Size of the points
    labSize = 4.0, 
    colAlpha = 0.8,
    col = c('grey', 'grey', 'grey', 'firebrick3'), # Custom colors
    legendPosition = 'none', # Remove legend
    drawConnectors = TRUE,   # Draw connectors for labels
    widthConnectors = 0.6, 
    typeConnectors = "open",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    ylim = c(0, 2),        # Set y-axis limits
    colConnectors = 'black'  # Color of the connectors
  )
  
  ggsave(paste0(dir.figs, output_file, '.pdf'), volcano_plot, height = 5, width = 5)
}

# Use the function to create and save the plots
create_and_save_volcano("files/ctc_vs_tumor.gene_summary.combine.txt", "CTC vs Primary Tumor", "Volcano.genes.ctc.vs.tumor")

#=================================================================================================================#
# Fig2D
#=================================================================================================================#

features <- c("MKI67","UBE2C","CKS1B","CIT","KIF5B","PCNT","EVC2","TCTE1","TAGLN")
c.order <- c("LM2-Primary tumor","LM2-Lung Metastasis","LM2-CTC Single","LM2-CTC cluster")
lm2.mg$cell.type <- factor(lm2.mg$cell.type, levels = c.order)
dot_plot <- FeatureStatPlot(lm2.mg, features = features, ident = "cell.type",
                            plot_type = "dot",cluster_columns=F,palette = "Purples")
ggsave(paste0(dir.figs, 'dot.screen', '.pdf'), plot = dot_plot, width = 4, height = 5)
#=================================================================================================================#
# Fig2E
#=================================================================================================================#
# Generating data needed for the plot, the plot itself is done by Prism 10 to make sure the style of pie charts remain consistant
# Create a new dataframe with Phase and cell.type
#=================================================================================================================#
#For LM2-NSG
# get scaled data for MKI67
FP.MKI67 <- FeaturePlot(lm2.mg, features = c("MKI67"))
FP.MKI67.df <- as.data.frame(FP.MKI67$data)
FP.MKI67.df$type <- unlist(lapply(rownames(FP.MKI67.df), function(x){strsplit(as.character(x), '\\_')[[1]][2]}))
# Calculate the maximum expression of the gene
FP.MKI67.df <- FP.MKI67.df %>%
  mutate(Expression = case_when(
    MKI67 >= 2 ~ ">=2",
    MKI67 <= 0 ~ "<=0",
    between(MKI67, 0, 2) ~ "0 to 2"
  ))
FP.MKI67.dfx <- FP.MKI67.df %>% 
  dplyr::group_by(type) %>%
  dplyr::count(Expression) %>%
  dplyr::summarize(count=n,across()) %>%  # count records by celltype,across() will keep all the columns
  dplyr::mutate(pct = count/sum(count)*100) %>% # find percent of total
  dplyr::ungroup() %>%
  dplyr::arrange(match(type, c("Primary", "Metastasis", "single","cluster")), desc(Expression))
write.csv(FP.MKI67.dfx, file = paste0(dir.files, analysisid, "_", "KI67.percentage.LM2", '.csv'), sep=',', quote=F, col.names=T, row.names=T)
#=================================================================================================================#
#For MVT1-FVB
# get scaled data for MKI67
FP.MKI67 <- FeaturePlot(mv.mg, features = c("Mki67"))
FP.MKI67.df <- as.data.frame(FP.MKI67$data)
FP.MKI67.df$type <- unlist(lapply(rownames(FP.MKI67.df), function(x){strsplit(as.character(x), '\\_')[[1]][2]}))
# Calculate the maximum expression of the gene
FP.MKI67.df <- FP.MKI67.df %>%
  mutate(Expression = case_when(
    Mki67 >= 2 ~ ">=2",
    Mki67 <= 0 ~ "<=0",
    between(Mki67, 0, 2) ~ "0 to 2"))
FP.MKI67.dfx <- FP.MKI67.df %>% 
  dplyr::group_by(type) %>%
  dplyr::count(Expression) %>%
  dplyr::summarize(count=n,across()) %>%  # count records by celltype,across() will keep all the columns
  dplyr::mutate(pct = count/sum(count)*100) %>% # find percent of total
  dplyr::ungroup() %>%
  dplyr::arrange(match(type, c("Primary", "Metastasis", "single","cluster")), desc(Expression))
write.csv(FP.MKI67.dfx, file = paste0(dir.files, analysisid, "_", "KI67.percentage.MVT1", '.csv'), sep=',', quote=F, col.names=T, row.names=T)

#=================================================================================================================#
#For Br16-NSG
FP.MKI67 <- FeaturePlot(Br.mg, features = c("MKI67"))
FP.MKI67.df <- as.data.frame(FP.MKI67$data)
FP.MKI67.df$type <- unlist(lapply(rownames(FP.MKI67.df), function(x){strsplit(as.character(x), '\\_')[[1]][1]}))
FP.MKI67.df <- FP.MKI67.df %>%
  mutate(Expression = case_when(
    #MKI67 >= 1.5 ~ "High",
    #MKI67 < 1.5 ~ "Low"
    MKI67 >= 2 ~ ">=2",
    MKI67 <= 0 ~ "<=0",
    between(MKI67, 0, 2) ~ "0 to 2"
  ))
FP.MKI67.dfx <- FP.MKI67.df %>% 
  dplyr::group_by(type) %>%
  dplyr::count(Expression) %>%
  dplyr::summarize(count=n,across()) %>%  # count records by celltype,across() will keep all the columns
  dplyr::mutate(pct = count/sum(count)*100) %>% # find percent of total
  dplyr::ungroup() %>%
  dplyr::arrange(match(type, c("Primary", "Metastasis", "Single","Cluster")), desc(Expression))
write.csv(FP.MKI67.dfx, file = paste0(dir.files, analysisid, "_", "KI67.percentage.Br16", '.csv'), sep=',', quote=F, col.names=T, row.names=T)

#=================================================================================================================#
# Fig2F
#=================================================================================================================#
# Plot generated with Prism 10 using data above

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to files/sessioninfo_detailed.txt

