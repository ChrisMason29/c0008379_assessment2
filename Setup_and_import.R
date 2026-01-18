## Packages -------------------------------------------------------------
# Install (run once, then comment out)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("tximport", "DESeq2", "biomaRt"))
# install.packages(c("pheatmap", "tidyverse", "ggrepel", "patchwork"))

library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(ggrepel)
library(patchwork)



## Import data ----------------------------------------------------------

sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')

files = pull(sample_table, Run)

files = paste0('counts/', files, '/quant.sf')

names(files) = pull(sample_table, Run)

gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')

txi = tximport(files, 
                 type='salmon',
                 tx2gene=gene_map,
                 ignoreTxVersion=TRUE)
