# 04_sample_distance_analysis.R
# Purpose: Sample-to-sample distance heatmap (QC / clustering)

# Prerequisites:
# - source("scripts/00_setup_and_import.R")
# - source("scripts/03_pca.R")  (or re-create rld here)
# Objects expected:
# - dds (DESeqDataSet)
# - tp_cols (named vector for timepoint colours)

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

## Transform (if not already created elsewhere) --------------------------
# If you already create `rld` in 03_pca.R, you can remove this and instead `source("scripts/03_pca.R")`.
rld <- rlog(dds)

## Sample distances -------------------------------------------------------
# Note: method is "euclidean" (spelling matters)
sample_distance <- dist(t(assay(rld)), method = "euclidean")
sample_distance_matrix <- as.matrix(sample_distance)

## Annotation (Timepoint labels) -----------------------------------------
group_to_timepoint <- function(group) {
  dplyr::case_when(
    group == "Allo2h"  ~ "2h post-reperfusion",
    group == "Allo24h" ~ "24h post-reperfusion",
    TRUE               ~ "Naive (control)"
  )
}

timepoint_levels <- c("Naive (control)", "2h post-reperfusion", "24h post-reperfusion")

heatmap_annotation <- tibble(
  Group = colData(dds)$Group
) %>%
  mutate(Timepoint = factor(group_to_timepoint(Group), levels = timepoint_levels)) %>%
  select(Timepoint) %>%
  as.data.frame()

rownames(heatmap_annotation) <- colnames(sample_distance_matrix)

## Colour scale ----------------------------------------------------------
blue_cols <- colorRampPalette(brewer.pal(9, "Blues"))(100)
blue_cols_rev <- rev(blue_cols)

## Plot heatmap ----------------------------------------------------------
pheatmap(
  sample_distance_matrix,
  color = blue_cols_rev,
  border_color = NA,
  clustering_distance_rows = sample_distance,
  clustering_distance_cols = sample_distance,
  annotation_col = heatmap_annotation,
  annotation_row = heatmap_annotation,
  annotation_colors = list(Timepoint = tp_cols),
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_col = FALSE,
  annotation_names_row = FALSE,
  fontsize = 14,
  main = "Euclidean distance heatmap of RNA-seq samples"
)

## Save figure (optional) ------------------------------------------------
# png("figures/fig_sample_distance_heatmap.png", width = 1800, height = 1600, res = 300)
# pheatmap(
#   sample_distance_matrix,
#   color = blue_cols_rev,
#   border_color = NA,
#   clustering_distance_rows = sample_distance,
#   clustering_distance_cols = sample_distance,
#   annotation_col = heatmap_annotation,
#   annotation_row = heatmap_annotation,
#   annotation_colors = list(Timepoint = tp_cols),
#   show_rownames = FALSE,
#   show_colnames = FALSE,
#   annotation_names_col = FALSE,
#   annotation_names_row = FALSE,
#   fontsize = 14,
#   main = "Euclidean distance heatmap of RNA-seq samples"
# )
# dev.off()
