# 05_differential_expression.R
# Purpose: Differential expression (2h vs control; 24h vs control) + side-by-side volcano plots

# Prerequisites:
# - source("scripts/00_setup_and_import.R")
# Objects expected:
# - dds (DESeqDataSet) that has been run through DESeq()
#
# Notes:
# - This script uses ggrepel for labels and patchwork to combine plots.
# - If you want gene symbols on labels, provide `annot_results` with:
#     ensembl_gene_id, external_gene_name
#   Otherwise the script will label ENSEMBL IDs.

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(patchwork)

## Optional: annotate ENSEMBL -> gene symbol -----------------------------
# If you already have `annot_results` created elsewhere, this will use it.
# If not, volcano labels will fall back to ENSEMBL IDs.
have_annot <- exists("annot_results") &&
  all(c("ensembl_gene_id", "external_gene_name") %in% colnames(annot_results))

## Helper: build volcano-ready tibble -----------------------------------
prep_volcano <- function(dds, group_num, group_den = "Naive") {
  res <- results(dds, contrast = c("Group", group_num, group_den))
  res_tbl <- as_tibble(res, rownames = "ensembl_gene_id") %>%
    filter(complete.cases(.)) %>%
    mutate(
      logPVal = -log10(padj),
      Expression = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE                              ~ "Not significant"
      ),
      Expression = factor(Expression, levels = c("Not significant", "Upregulated", "Downregulated"))
    )
  list(res = res, tbl = res_tbl)
}

## Helper: choose top labels --------------------------------------------
get_top_labels <- function(volcano_tbl, n = 10) {
  top <- volcano_tbl %>%
    filter(Expression != "Not significant") %>%
    arrange(padj) %>%
    slice_head(n = n)

  if (have_annot) {
    top <- top %>%
      left_join(
        annot_results %>% select(ensembl_gene_id, external_gene_name),
        by = "ensembl_gene_id"
      ) %>%
      mutate(label = ifelse(is.na(external_gene_name) | external_gene_name == "",
                            ensembl_gene_id,
                            external_gene_name))
  } else {
    top <- top %>% mutate(label = ensembl_gene_id)
  }

  top
}

## Helper: volcano plot --------------------------------------------------
make_volcano_plot <- function(volcano_tbl, top_labels, plot_title) {
  ggplot(volcano_tbl, aes(x = log2FoldChange, y = logPVal, colour = Expression)) +
    geom_point(alpha = 0.8, size = 1.5) +
    geom_text_repel(
      data = top_labels,
      aes(label = label),
      size = 3,
      max.overlaps = Inf
    ) +
    scale_colour_manual(values = c(
      "Not significant" = "grey70",
      "Upregulated"     = "#D55E00",
      "Downregulated"   = "#0072B2"
    )) +
    labs(
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "adjusted p-value")
    ) +
    ggtitle(plot_title) +
    theme_minimal() +
    theme(
      legend.position = c(0.88, 0.88),
      legend.background = element_rect(fill = "grey93", color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10))
    )
}

## 1) 2h vs control ------------------------------------------------------
de_2h <- prep_volcano(dds, group_num = "Allo2h", group_den = "Naive")
message("DE summary: 2h vs control")
print(summary(de_2h$res))

top10_2h <- get_top_labels(de_2h$tbl, n = 10)

p_vol_2h <- make_volcano_plot(
  volcano_tbl = de_2h$tbl,
  top_labels  = top10_2h,
  plot_title  = "Differential gene expression: 2h post-reperfusion vs control"
)

## 2) 24h vs control -----------------------------------------------------
de_24h <- prep_volcano(dds, group_num = "Allo24h", group_den = "Naive")
message("DE summary: 24h vs control")
print(summary(de_24h$res))

top10_24h <- get_top_labels(de_24h$tbl, n = 10)

p_vol_24h <- make_volcano_plot(
  volcano_tbl = de_24h$tbl,
  top_labels  = top10_24h,
  plot_title  = "Differential gene expression: 24h post-reperfusion vs control"
)

## Combine (2h left, 24h right) -----------------------------------------
p_volcano_panel <- p_vol_2h + p_vol_24h + plot_layout(guides = "collect") &
  theme(legend.position = "right")

## Save (optional) -------------------------------------------------------
# ggsave("figures/fig_volcano_panel_2h_and_24h.png", p_volcano_panel,
#        width = 12, height = 5, dpi = 300)

## Display ---------------------------------------------------------------
p_volcano_panel
