# 03_principal_component_analysis.R
# Purpose: PCA of RNA-seq samples using DESeq2 rlog-transformed counts

# Prerequisites:
# - source("scripts/00_setup_and_import.R")
# Objects expected:
# - dds (DESeqDataSet)
# - tp_cols (named vector for timepoint colours, matching labels below)

library(DESeq2)
library(tidyverse)

## rlog transform --------------------------------------------------------
# rlog is suitable for small sample sizes; vst is an alternative for larger datasets
rld <- rlog(dds)

## Get PCA coordinates + variance explained ------------------------------
pca_data <- plotPCA(rld, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

## Add readable timepoint labels ----------------------------------------
group_to_timepoint <- function(group) {
  dplyr::case_when(
    group == "Allo2h"  ~ "2h post-reperfusion",
    group == "Allo24h" ~ "24h post-reperfusion",
    TRUE               ~ "Naive (control)"
  )
}

timepoint_levels <- c("Naive (control)", "2h post-reperfusion", "24h post-reperfusion")

pca_data <- pca_data %>%
  mutate(
    Timepoint = factor(group_to_timepoint(group), levels = timepoint_levels)
  )

## Plot PCA --------------------------------------------------------------
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Timepoint, shape = Timepoint)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = tp_cols) +
  labs(
    x = bquote(bold("PC1:") ~ .(percentVar[1]) * "% variance"),
    y = bquote(bold("PC2:") ~ .(percentVar[2]) * "% variance"),
    title = "Principal component analysis of RNA-seq samples",
    colour = "Timepoint",
    shape = "Timepoint"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 12)),
    legend.position = c(0.87, 0.12),
    legend.background = element_rect(fill = "grey93", color = "white")
  )

## Save figure (optional) -----------------------------------------------
# ggsave("figures/fig_pca.png", p_pca, width = 6.5, height = 5.5, dpi = 300)

## Display ---------------------------------------------------------------
p_pca
