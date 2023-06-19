# script to look at mapping statistics, RNA content,
# and community composition of metatranscriptomes.

# load libraries
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(viridis)

# percentages of mapped read types

# read in data
cds_counts <- unlist(read_tsv("./metatranscriptomics_2023_04/data/rna_align/counttable_bb.txt.summary")[1,-1])
rrna_counts <- unlist(read_tsv("./metatranscriptomics_2023_04/data/rna_align/rRNA_counttable_bb.txt.summary")[1,-1])
trna_counts <- unlist(read_tsv("./metatranscriptomics_2023_04/data/rna_align/tRNA_counttable_bb.txt.summary")[1,-1])

# Now I create a dataframe and plot a stacked barplot of the data:
# create dataframe
df_rnatypes <-  tibble(
  rna_type = c(rep("rRNA", length(rrna_counts)),
               rep("tRNA", length(trna_counts)),
               rep("mRNA", length(cds_counts))),
  counts = c(rrna_counts,  trna_counts, cds_counts),
  sample = rep(names(rrna_counts), 3)
)

# rename samples:
df_rnatypes$sample <- gsub(".*rna_align/(.+)\\.bam", "\\1", df_rnatypes$sample)
df_rnatypes$sample <- gsub("^\\d_", "", df_rnatypes$sample)

# protocol
df_rnatypes$protocol <- gsub(".*((NEB)|(RZ)).*", "\\1", df_rnatypes$sample)

# plot
g_rnatype <- ggplot(df_rnatypes) +
  aes(x = sample, y = counts, fill = rna_type) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)",
                     breaks = seq(0,20000000, by = 1000000)) +
  scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))+
  facet_grid(. ~ protocol, scales = "free_x")

g_rnatype

# save as pdf:
ggsave("./metatranscriptomics_2023_04/analysis/fig_rnatypes.pdf", g_rnatype, width = 10, height = 5)




