# 01_library_size_and_count_distribution.R
# Purpose: QC plots for sequencing depth (library size) and gene count distributions

# Prerequisites:
# - source("scripts/00_setup_and_import.R")
# Objects expected:
# - dds (DESeqDataSet)
# - txi (tximport output)
# - sample_table (must contain Run and Group)
# - tp_cols (named vector for timepoint colours)

## Helper: map Group -> Timepoint label ---------------------------------
group_to_timepoint <- function(group) {
  dplyr::case_when(
    group == "Allo2h"  ~ "2h post-reperfusion",
    group == "Allo24h" ~ "24h post-reperfusion",
    TRUE               ~ "Naive (control)"
  )
}

timepoint_levels <- c("Naive (control)", "2h post-reperfusion", "24h post-reperfusion")

## Sanity check: DESeq2 counts vs tximport counts ------------------------
check_ok <- isTRUE(all.equal(colSums(counts(dds)), colSums(txi$counts)))
message("Check: DESeq2 column sums equal tximport column sums? -> ", check_ok)

## 1) Sequencing depth (library size) -----------------------------------
lib_sizes <- colSums(counts(dds))

lib_df <- tibble(
  Run          = names(lib_sizes),
  library_size = as.numeric(lib_sizes)
) %>%
  left_join(sample_table %>% select(Run, Group), by = "Run") %>%
  mutate(Timepoint = factor(group_to_timepoint(Group), levels = timepoint_levels))

# Summaries (console)
print(range(lib_df$library_size))
print(summary(lib_df$library_size))
print(aggregate(library_size ~ Timepoint, data = lib_df, FUN = summary))

p_libsize <- ggplot(lib_df, aes(x = Timepoint, y = library_size, fill = Timepoint)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
  scale_fill_manual(values = tp_cols) +
  labs(
    x = "Timepoint",
    y = "Library size (total counts)",
    title = "Distribution of RNA-seq library sizes by timepoint"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 12))
  )

## 2) Gene-level count distribution (log10(counts + 1)) ------------------
# Convert counts matrix to long format: one row per (gene, sample)
count_long <- as.data.frame(counts(dds)) %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Run",
    values_to = "count"
  ) %>%
  mutate(Value = log10(count + 1)) %>%
  left_join(sample_table %>% select(Run, Group), by = "Run") %>%
  mutate(Timepoint = factor(group_to_timepoint(Group), levels = timepoint_levels))

# Order samples by timepoint, then Run (so plots are consistent)
sample_order <- count_long %>%
  distinct(Run, Timepoint) %>%
  arrange(Timepoint, Run) %>%
  pull(Run)

count_long <- count_long %>%
  mutate(Run = factor(Run, levels = sample_order))

p_countdist <- ggplot(count_long, aes(x = Run, y = Value, fill = Timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(values = tp_cols) +
  labs(
    x = "Sample",
    y = expression("Gene counts (" * log[10] * ")"),
    title = "Distribution of RNA-seq gene counts across samples by timepoint",
    fill = "Timepoint"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.background = element_rect(fill = "grey93", color = "white"),
    legend.position = c(0.12, 0.88),
    legend.justification = c(0, 1)
  )

## Combine plots ---------------------------------------------------------
p_combined <- p_libsize + p_countdist

## Save figures (optional) ----------------------------------------------
# ggsave("figures/fig_library_size_by_timepoint.png", p_libsize, width = 6, height = 4.5, dpi = 300)
# ggsave("figures/fig_count_distribution_by_sample.png", p_countdist, width = 7, height = 4.5, dpi = 300)
# ggsave("figures/fig_qc_libsize_and_countdist.png", p_combined, width = 12, height = 4.5, dpi = 300)

## Display ----------------------------------------------------------------
p_combined
