# 02_dispersion_estimates.R
# Purpose: Visualise DESeq2 dispersion estimates (model fit / QC)

# Prerequisites:
# - source("scripts/00_setup_and_import.R")
# Objects expected:
# - dds (DESeqDataSet) that has been run through DESeq()

library(DESeq2)

# Plot dispersion estimates
plotDispEsts(dds)

## Save figure (optional) -----------------------------------------------
# png("figures/fig_dispersion_estimates.png", width = 1600, height = 1200, res = 300)
# plotDispEsts(dds)
# dev.off()
