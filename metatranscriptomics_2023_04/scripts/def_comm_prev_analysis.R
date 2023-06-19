# script to look at mapping statistics, RNA content,
# and community composition of metatranscriptomes.

# load libraries
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(viridis)
library(MetBrewer)



# first, look at mapping statistics
mapping_stats <- read_tsv("./metatranscriptomics_2023_04/data/mappingstats.tab")

# make table longer. keep first column and make all other columns into rows. name new columns "sample" and "value"
mapping_stats_long <- mapping_stats %>%
  pivot_longer(-"read_type", names_to = "sample", values_to = "reads") %>%
  # add column for protocol
  mutate(protocol = gsub(".*((NEB)|(RZ)).*", "\\1", sample)) %>%
  # remove nr before sample names
  mutate(sample = gsub("^\\d_", "", sample)) %>%
  # change read type order and make factor
  mutate(read_type = factor(read_type, levels = c("Input reads", "after trimming", "mapped")))

# now create a stacked barplot with ggplot. It should have one bar per sample and all bars should ordered by the
# protocol (NEB first, then RZ). The bars should be stacked by read type: input reads, after trimming,
# mapped (whole genome). make bars transparent and start on y=0 for all bars, so we have a common baseline and
# overlaps.

# plot
g_mapping <- mapping_stats_long %>%
  ggplot() +
  # create a barplot with x = sample, y = reads, fill = read_type
    aes(x = sample, y = reads, fill = read_type) +
    # make a stacked barplot
    geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "reads (in million)",
                     breaks = seq(0,50000000, by = 5000000)) +
  scale_fill_manual(values = c("midnightblue", "steelblue", "skyblue")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))+
  facet_grid(. ~ protocol, scales = "free_x")
g_mapping

# save as pdf:
ggsave("./metatranscriptomics_2023_04/analysis/fig_mapping.pdf", g_mapping, width = 10, height = 5)

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




