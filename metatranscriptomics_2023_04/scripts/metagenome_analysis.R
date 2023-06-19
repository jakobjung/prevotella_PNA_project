# script to analyse metatranscriptomics data. This script is part of the metatranscriptomics_2023_04 project.
#

# load libraries
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)
library(viridis)
library(MetBrewer)
library(tidybulk)
library(edgeR)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(KEGGREST)
library("EDASeq")
library(xlsx)



# get the count table
count_table_mrna <- read_tsv("./metatranscriptomics_2023_04/data/rna_align/counttable_bb.txt", skip = 1)
count_table_rrna <- read_tsv("./metatranscriptomics_2023_04/data/rna_align/rRNA_counttable_bb.txt", skip = 1)
count_table_trna <- read_tsv("./metatranscriptomics_2023_04/data/rna_align/tRNA_counttable_bb.txt", skip = 1)

# merge by row
count_table <- rbind(count_table_mrna, count_table_rrna, count_table_trna)

# change colnames:
pnapat <- ".*.rna_align/(\\d_)?(.+)\\.bam"
colnames(count_table) <- gsub(pnapat,"\\2", colnames(count_table))

# make table pivot_longer for columns 7 to the last.
count_table_long <- count_table %>%
  pivot_longer(-c("Geneid", "Chr", "Start", "End", "Strand", "Length"), names_to = "sample", values_to = "counts") %>%
  # remove RZ protocol
  filter(grepl("NEB", sample)) %>%
  # add column for organism (part before first "_" in Geneid
  mutate(organism = gsub("^(.*)_.*", "\\1", Geneid)) %>%
  mutate(organism = ifelse(organism == "FE838", "B. theta",
                           ifelse(organism =="Bovatus", "B. ovatus", "P. copri"))) %>%
  mutate(sample = gsub("_NEB_DPLEX", "", sample))

# create stacked barplot with ggplot. It should have one bar per sample and all bars should ordered by the protocol
# (NEB first, then RZ). The bars should be stacked by organism. The y axis should be the sum of all counts per organism.

# plot
g_counts <- count_table_long %>%
  # group by organism and get sum of counts per organism
    group_by(organism, sample) %>% summarise(counts = sum(counts)) %>%
  ggplot() +
  # create a barplot with x = sample, y = counts, fill = organism
  aes(x = sample, y = counts, fill = organism) +
  # make a stacked barplot
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6), name = "counts (in million)",
                     breaks = seq(0,50000000, by = 1000000)) +
  scale_fill_manual(values = c("midnightblue", "steelblue", "skyblue")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))
g_counts

# save plot as pdf
ggsave("./metatranscriptomics_2023_04/analysis/counts_per_organism.pdf", g_counts, width = 6, height = 5)

# make same plot but with relative counts summing up to 100% per sample
g_counts_rel <- count_table_long %>%
  # group by organism and get sum of counts per organism
  group_by(organism, sample) %>% summarise(counts = sum(counts)) %>%
  # group by sample and get sum of counts per sample
  group_by(sample) %>% mutate(counts = counts/sum(counts)) %>%
  ggplot() +
  # create a barplot with x = sample, y = counts, fill = organism
  aes(x = sample, y = counts, fill = organism) +
  # make a stacked barplot
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_y_continuous(labels = scales::percent_format(), name = "relative counts",
                     breaks = seq(0,1, by = 0.1)) +
  scale_fill_manual(values = c("midnightblue", "steelblue", "skyblue")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))
g_counts_rel

# save plot as pdf
ggsave("./metatranscriptomics_2023_04/analysis/counts_per_organism_rel.pdf", g_counts_rel, width = 6, height = 5)


#read in the count_table_mrna to tidybulk

# remove RZ samples
count_table_mrna <- count_table_mrna %>% select(-contains("RZ"))

# rename colnames
colnames(count_table_mrna) <- gsub(pnapat,"\\2", colnames(count_table_mrna))
colnames(count_table_mrna) <- gsub("_NEB_DPLEX","", colnames(count_table_mrna))


# use normal edger for DE analysis:
gwc <- as.data.frame(count_table_mrna[,6:length(count_table_mrna[1,])])
rownames(gwc) <- count_table_mrna$Geneid

gene_lengths <- gwc$Length
raw_counts <- gwc[,-1]
norm_length <- data.frame(sapply(raw_counts, function(x) x / gene_lengths))
tpm_matrix <- data.frame(sapply(norm_length, function(x) x * 1e6 / sum(x)), row.names = rownames(raw_counts))
# pairs(log(tpm_matrix))

# separate dataframe gwc to 3 dataframes for each organism
gwc_btheta <- gwc[grep("FE838", rownames(gwc)),]
gwc_bovatus <- gwc[grep("Bovatus", rownames(gwc)),]
gwc_pcopri <- gwc[grep("LK433", rownames(gwc)),]

# create sample (test variable) for all samples in gwc (no length)
test <- as.factor(gsub("_\\d$", "",colnames(gwc)[-1]))

# create edgeR objects for each organism & remove ordganism_name from rownames
y_btheta <- DGEList(counts = gwc_btheta[-1], group = test, genes =  gwc_btheta[,1,drop=FALSE])
y_bovatus <- DGEList(counts = gwc_bovatus[-1], group = test, genes =  gwc_bovatus[,1,drop=FALSE])
y_pcopri <- DGEList(counts = gwc_pcopri[-1], group = test, genes =  gwc_pcopri[,1,drop=FALSE])

# filter out genes with low counts
keep <- filterByExpr(y_btheta)
y_btheta <- y_btheta[keep,,keep.lib.sizes=FALSE]
keep <- filterByExpr(y_bovatus)
y_bovatus <- y_bovatus[keep,,keep.lib.sizes=FALSE]
keep <- filterByExpr(y_pcopri)
y_pcopri <- y_pcopri[keep,,keep.lib.sizes=FALSE]



# create a design matrix that can be used for all three organisms
design_matrix <- model.matrix(~0+test)
colnames(design_matrix) <- levels(test)
rownames(design_matrix) <- colnames(y_btheta)
design_matrix

# do tmm normalization
y_btheta <- calcNormFactors(y_btheta, method = "TMM")
y_bovatus <- calcNormFactors(y_bovatus, method = "TMM")
y_pcopri <- calcNormFactors(y_pcopri, method = "TMM")

# estimate dispersion
y_btheta <- estimateDisp(y_btheta, design_matrix, robust = TRUE)
y_bovatus <- estimateDisp(y_bovatus, design_matrix, robust = TRUE)
y_pcopri <- estimateDisp(y_pcopri, design_matrix, robust = TRUE)


# define colors:
colors <- c( "navy", "blue", "darkgreen", "green")

# create PCA and RLE plots in a nicer way:

# get logcpms for all organisms
logCPMs <- list( b_theta = cpm(y_btheta, log=TRUE, prior.count=2), b_ovatus = cpm(y_bovatus, log=TRUE, prior.count=2),
                 p_copri = cpm(y_pcopri, log=TRUE, prior.count=2))

# make theme for plots:
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             title=element_text(colour="black", size=20),
             axis.text=element_text(colour="black", size=15),axis.ticks=element_line(colour="black"),
             axis.title=element_text(colour="black", size=20),
             plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none")

# creaate PCA plots and save them as pdf:
lapply(names(logCPMs), function(i) {
  pca <- prcomp(t(logCPMs[[i]]))
  pca_df <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste( colnames(pca_df), "(", paste( as.character(percentage), "%", ")", sep="") )
  pca_df$group <-test
  p<-ggplot(pca_df,aes(x=PC1,y=PC2,group=group,label=rownames(pca_df), colour=group))
  p<-p+geom_point(size=3)+ scale_shape_identity()+
    geom_text_repel(size=8, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 20)+
    theme + xlab(percentage[1]) +
    ggtitle(paste0("PCA after TMM in ", i))+
    ylab(percentage[2])+ scale_color_manual(values = colors)
  p
  # save plot as pdf
  pdf(paste0("./metatranscriptomics_2023_04/analysis/pca_", i, ".pdf"), width = 11, height = 10)
  print(p)
  dev.off()
})

# do DE analysis for all organisms. start by making contrasts betweeen test conditions that
# can be used for all organisms:
contrasts <- makeContrasts( con_2_vs_0 = CON_2 - CON_0,
                            pna_2_vs_0 = PNA_2 - PNA_0,
                            pna_0_vs_con_0 = PNA_0 - CON_0,
                            pna_2_vs_con_2 = PNA_2 - CON_2,
                            pna_2_vs_0_vs_ctrl_2_vs_0 = (PNA_2 - PNA_0) - (CON_2 - CON_0),
                            levels = design_matrix)

# create list of all three organisms
y_list <- list(y_btheta=y_btheta, y_bovatus=y_bovatus, y_pcopri=y_pcopri)

# do DE analysis for all organisms
de_list <- lapply(y_list, function(y) {
  fit <- glmQLFit(y, design_matrix, robust = TRUE)
  res <- list(con_2_vs_0 = glmQLFTest(fit, contrast = contrasts[,1]),
              pna_2_vs_0 = glmQLFTest(fit, contrast = contrasts[,2]),
              pna_0_vs_con_0 = glmQLFTest(fit, contrast = contrasts[,3]),
              pna_2_vs_con_2 = glmQLFTest(fit, contrast = contrasts[,4]),
              pna_2_vs_0_vs_ctrl_2_vs_0 = glmQLFTest(fit, contrast = contrasts[,5]))
  res <- lapply(res, function (r) {
    r$table$FDR <- p.adjust(r$table$PValue, method = "fdr")
    r
  })
  res
})

names(de_list) <- c("b_theta", "b_ovatus", "p_copri")

# save all tables in de_list as csv file:
lapply(names(de_list), function(i) {
  lapply(names(de_list[[i]]), function(j) {
    # sort by FDR:
    de_list[[i]][[j]]$table <- de_list[[i]][[j]]$table[order(de_list[[i]][[j]]$table$FDR),]
    write.csv(de_list[[i]][[j]]$table, file = paste0("./metatranscriptomics_2023_04/data/DE_raw_data/", i, "_", j, ".csv"))
    write.xlsx(de_list[[i]][[j]]$table, file = paste0("./metatranscriptomics_2023_04/data/DE_raw_data/", i, "_", j, ".xlsx"))
  })
})


de_list$p_copri$pna_2_vs_con_2$table[c("LK433_RS02385", "LK433_RS02380"),]
de_list$p_copri$pna_0_vs_con_0$table[c("LK433_RS02385", "LK433_RS02380"),]

# rename locus tags LK433_RS02385, LK433_RS02380 to acpP and fabF, resp. in rownames of all p_copri tables in delist.
# use a for-loop for this.
for (i in 1:length(de_list$p_copri)) {
  rownames(de_list$p_copri[[i]]$table)[rownames(de_list$p_copri[[i]]$table) == "LK433_RS02385"] <- "acpP"
  rownames(de_list$p_copri[[i]]$table)[rownames(de_list$p_copri[[i]]$table) == "LK433_RS02380"] <- "fabF"
}

# same in b_theta
for (i in 1:length(de_list$b_theta)) {
  rownames(de_list$b_theta[[i]]$table)[rownames(de_list$b_theta[[i]]$table) == "FE838_RS06320"] <- "acpP"
  rownames(de_list$b_theta[[i]]$table)[rownames(de_list$b_theta[[i]]$table) == "FE838_RS06315"] <- "fabF"
}

# same in b_ovatus
for (i in 1:length(de_list$b_ovatus)) {
  rownames(de_list$b_ovatus[[i]]$table)[rownames(de_list$b_ovatus[[i]]$table) == "Bovatus_RS17010"] <- "acpP"
  rownames(de_list$b_ovatus[[i]]$table)[rownames(de_list$b_ovatus[[i]]$table) == "Bovatus_RS17005"] <- "fabF"
}

# create volcano plots for all organisms and all contrasts. save them as pdfs in the analysis folder.
# always mark the acpP and fabF genes (locus tags RS02385, RS02380, resp.) in the plots .
# use EnhancedVolcano package for this.

lapply(names(de_list), function(i) {
  names(de_list[[i]]) <- c("CTRL 2h - 0h", "PNA 2h - 0h", "PNA - CTRL 0h", "PNA - CTRL 2h", "PNA - CTRL 2h - 0h")
  lapply(names(de_list[[i]]), function(j) {
    res <- de_list[[i]][[j]]$table
    # save volc as pdf
    g <- EnhancedVolcano(res, lab = NA, x = "logFC", y = "FDR",
                          pCutoff = 0.05, FCcutoff = 1,
                          ylab = bquote(-log[10]~FDR),
                          pointSize = 4,
                          #boxedLabels = TRUE,
                          #parseLabels = TRUE,
                          legendPosition = 'none',
                          subtitle = NULL,
                          xlim = c(-7, 7),
                          ylim = c(0, 20),
                          title = paste0("Volcano plot of ", i, " for ", j),
                          #labSize = 5, labCol = "black",
                         col = c("grey", "gray48", "gray48", "darkred"),
                         # remove caption:
                            caption = NULL,
                          gridlines.major = FALSE, gridlines.minor = FALSE) +
      coord_cartesian(xlim=c(-7, 7)) +
      scale_x_continuous(breaks=seq(-6,6, 2)) +
      geom_point(aes(x = logFC, y = -log10(FDR)), data = res[rownames(res) %in% c("acpP", "fabF"),], size = 4,
                 fill = "darkorange", shape=21) + geom_label_repel(aes(x = logFC, y = -log10(FDR)),
                                                   data = res[rownames(res) %in% c("acpP", "fabF"),], size = 4,
                                                   color = "black", label = rownames(res[rownames(res) %in% c("acpP", "fabF"),]),
                                                    fontface = "italic")

    # save plot as svg
    svg(paste0("./metatranscriptomics_2023_04/analysis/volcano_", i, "_", j, ".svg"), width = 11, height = 10)
    print(g)
    dev.off()

  })
})

# create a ggplot barplot with log2FC on y-axis and sample on x-axis, for all 2h time point PNA vs CTRL for all three
# organisms. split bars by organism. add significance stars.

# start by getting logfc and fdr for all 2h time points for all organisms
df_bar_pna <- as_tibble(de_list$p_copri$pna_2_vs_con_2$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
  mutate(organism = "P. copri")) %>%
    mutate(gene= c("acpP", "fabF")) %>% mutate(sample = "PNA 2h - CTRL 2h") %>%
    rbind(as_tibble(de_list$b_ovatus$pna_2_vs_con_2$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
    mutate(organism = "B. ovatus")) %>% mutate(sample = "PNA 2h - CTRL 2h") %>%
    mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_theta$pna_2_vs_con_2$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
    mutate(organism = "B. theta")) %>% mutate(sample = "PNA 2h - CTRL 2h") %>%
    mutate(gene= c("acpP", "fabF"))) %>%
  # add same but for PNA 2h vs 0h
    rbind(as_tibble(de_list$p_copri$pna_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "P. copri")) %>% mutate(sample = "PNA 2h - PNA 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_ovatus$pna_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. ovatus")) %>% mutate(sample = "PNA 2h - PNA 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_theta$pna_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. theta")) %>% mutate(sample = "PNA 2h - PNA 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
  # add same but for CTRL 2h vs 0h
    rbind(as_tibble(de_list$p_copri$con_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "P. copri")) %>% mutate(sample = "CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_ovatus$con_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. ovatus")) %>% mutate(sample = "CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_theta$con_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. theta")) %>% mutate(sample = "CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
  # add same but for PNA 0h vs CTRL 0h
    rbind(as_tibble(de_list$p_copri$pna_0_vs_con_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "P. copri")) %>% mutate(sample = "PNA 0h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_ovatus$pna_0_vs_con_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. ovatus")) %>% mutate(sample = "PNA 0h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_theta$pna_0_vs_con_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "B. theta")) %>% mutate(sample = "PNA 0h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    # add same but for PNA 2h vs 0h - CTRL 2h vs 0h
    rbind(as_tibble(de_list$p_copri$pna_2_vs_0_vs_ctrl_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
        mutate(organism = "P. copri")) %>% mutate(sample = "PNA 2h - PNA 0h - CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_ovatus$pna_2_vs_0_vs_ctrl_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
      mutate(organism = "B. ovatus")) %>% mutate(sample = "PNA 2h - PNA 0h - CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF"))) %>%
    rbind(as_tibble(de_list$b_theta$pna_2_vs_0_vs_ctrl_2_vs_0$table[c("acpP", "fabF"), c("logFC", "FDR")] %>%
      mutate(organism = "B. theta")) %>% mutate(sample = "PNA 2h - PNA 0h - CTRL 2h - CTRL 0h") %>%
        mutate(gene= c("acpP", "fabF")))


# plot barplots for each sample

# create a function creating above plots for a defined sample of df_bar_pna
plot_bar <- function(df, s, y_lim, y_lab, vj=0.2){
  df %>% filter(sample==s) %>%
    ggplot(aes(x=organism, y=logFC, fill=gene)) +
    geom_hline(yintercept=1, color="black", linetype = "dotted") +
    geom_hline(yintercept=-1, color="black", linetype = "dotted") +
    geom_hline(yintercept=0, color="black" ) +
    geom_bar(stat="identity",width = 0.8, position=position_dodge(), color="black") +
    geom_text(aes(label=ifelse(FDR<0.01, "*", "")), vjust=vj,hjust=0.5 , position=position_dodge(.8, ), size = 10,
              color="black") +
    theme_classic() + ylim(y_lim)+
    scale_fill_manual(values=c("steelblue", "lightblue")) +
    theme(axis.text.x = element_text(face="italic",angle = 0,  vjust=1, size = 16),
          axis.text.y = element_text(size=14),
          axis.title = element_text(size=16),
          title = element_text(size=15),
          #put legend on upper right corner
          legend.position = c(0.9, 0.92),
          legend.title = element_text(size=12),
          legend.text = element_text(face="italic", size = 12)) +
    labs(x="", y=y_lab, fill="target gene")
}

# generate plots for all samples:
g_tc_2_vs_0 <- plot_bar(df_bar_pna, "PNA 2h - PNA 0h", c(-2.5,2.5), expression(log[2]~FC~2~h~PNA~vs.~0~h~PNA))
g_tc_2_vs_0

g_tc_2_vs_ctrl <- plot_bar(df_bar_pna, "PNA 2h - CTRL 2h", c(-2.5,2.5),
                           expression(log[2]~FC~2~h~PNA~vs.~2~h~CTRL), 1.2)
g_tc_2_vs_ctrl

g_tc_ctrl_2h_vs_ctrl_0h <- plot_bar(df_bar_pna, "CTRL 2h - CTRL 0h", c(-2.5,2.5),
                           expression(log[2]~FC~2~h~CTRL~vs.~0~h~CTRL))
g_tc_ctrl_2h_vs_ctrl_0h

g_tc_0_pna_vs_ctrl_0h <- plot_bar(df_bar_pna, "PNA 0h - CTRL 0h", c(-2.5,2.5),
                           expression(log[2]~FC~0~h~PNA~vs.~0~h~CTRL))
g_tc_0_pna_vs_ctrl_0h

g_tc_pna_2_vs_0_vs_ctrl_2_vs_0 <- plot_bar(df_bar_pna, "PNA 2h - PNA 0h - CTRL 2h - CTRL 0h", c(-2.5,2.5),
                           expression(log[2]~FC~PNA~2~h~-~0~h~vs.~CTRL~2~h~-~0~h))
g_tc_pna_2_vs_0_vs_ctrl_2_vs_0

# now use cowplot to arrange plots in a grid. add breaks between and headers:
library(cowplot)
cp <- plot_grid(g_tc_2_vs_0,g_tc_ctrl_2h_vs_ctrl_0h,NULL, g_tc_2_vs_ctrl, g_tc_0_pna_vs_ctrl_0h,
                g_tc_pna_2_vs_0_vs_ctrl_2_vs_0,
                # add space between plots
                scale = 0.8,
                ncol=3, nrow=2, labels = c("A", "B", "","C", "D", "E"), label_size = 20)

# save plot as svg
svg("./metatranscriptomics_2023_04/analysis/TC_pna.svg", width = 18, height = 13)
print(cp)
dev.off()

# Do pathway analysis in p. copri for now. use keggrest package for this.

# check whether
# get link and list of p. copri
link_kegg <- keggLink("pathway", "pcoi")
list_kegg <- keggList("pathway", "pcoi")

kegg_pw_ids <- names(list_kegg)

#rename genes, remove ones which arent in our data:
names(link_kegg) <- gsub("pcoi:", "", names(link_kegg)) #rename genes as locus tags

# get old locus tags:
old_lts <- read.delim("./metatranscriptomics_2023_04/data/old_locus_tags.tsv")
rownames(old_lts) <- old_lts$old_locus_tag
# remove ones without old lt:
link_kegg <- link_kegg[names(link_kegg) %in% old_lts$old_locus_tag]
names(link_kegg) <- old_lts[names(link_kegg),1]


link_kegg <- link_kegg[names(link_kegg) %in% c(rownames(de_list$p_copri$con_2_vs_0))] #remove genes not in data

link_kegg <- gsub("path:","", link_kegg)
idx_kegg <- sapply(kegg_pw_ids, function(x){
  x <- unique(names(link_kegg[link_kegg == x])) # choose all genes, except duplucates
})

# now we have the links of the pathways in p. copri.
# we can run fry gene set analysis on the log2FCs of the genes in the pathways.
l <- length(colnames(contrasts))
kegg_fry <- lapply(1:l, function(x) fry(y_list$y_pcopri,idx_kegg, design_matrix, contrasts[,x]))
names(kegg_fry) <- colnames(contrasts)

# add kegg tems:
for (fryres in names(kegg_fry)) {
  kegg_fry[[fryres]][["TERM"]] <- ifelse(grepl("pcoi",rownames(kegg_fry[[fryres]])),
                                              list_kegg[rownames(kegg_fry[[fryres]])],
                                              rownames(kegg_fry[[fryres]]))
  kegg_fry[[fryres]][["TERM"]] <- gsub("(.*) - Prevotella copri",
                                            "\\1", kegg_fry[[fryres]][["TERM"]])
  write.csv(kegg_fry[[fryres]], paste("./metatranscriptomics_2023_04/analysis/pathway_analysis/", fryres, ".csv", sep = ""))
}


kegg_frysig <- lapply(kegg_fry, function(x) x[x[["FDR"]]<0.1 & x[["NGenes"]]>3,])
kegg_siggos <- c()


for (i in names(kegg_frysig)) {
  print(i)
  print(dim(kegg_frysig[[i]]))
  print(kegg_frysig[[i]][,c(1,2,4,7)])
  kegg_siggos <- c(kegg_siggos, rownames(kegg_frysig[[i]][1:15,]))  # can be modified
}

kegg_siggos <- unique(kegg_siggos[!grepl("NA", kegg_siggos)])

# now generate a heatmap with complexheatmap
idx_kegg_char <- lapply(idx_kegg, as.character)


# I create a dataframe with mean logFC values for each significant GO-term:
hm_kegg <- t(as.data.frame(lapply(idx_kegg_char[kegg_siggos], function(x){
  sapply(names(de_list$p_copri), function(y){
    mean(de_list$p_copri[[y]]$table[x,]$logFC)
  })
})))

hm_kegg <- as.data.frame(hm_kegg)

# rownames(hm_kegg) <- gsub("\\.", "\\:", rownames(hm_kegg))


# make hm
hm_kegg <- hm_kegg[order(hm_kegg[,1], decreasing = T),]

kegg_sizes <- sapply(idx_kegg_char[rownames(hm_kegg)], function(x) length(x))

pvals <- data.frame(sapply(names(kegg_fry),
                           function(x) kegg_fry[[x]][rownames(hm_kegg),"FDR"]),
                    row.names = rownames(hm_kegg))

#select only significant ones:
pvals <-sapply(pvals, function(x) ifelse(x<0.05, x <- "*", x<-"") )

keggpws <- kegg_fry$con_2_vs_0[rownames(hm_kegg),] [["TERM"]]


rownames(hm_kegg) <- ifelse(!is.na(keggpws),keggpws, rownames(hm_kegg) )


# plot hm (save as pdf):
library(circlize)
col_fun = colorRamp2(c(-1.5,0, 1.5), c("darkblue", "white", "darkred"))

g_samples <- factor(c(rep("PNA vs. PNA0", 4), rep("SCR vs. SCR0", 4), rep("PNA vs. CON", 5), rep("PNA vs. SCR", 5)),
                    levels = c("PNA vs. PNA0", "SCR vs. SCR0", "PNA vs. CON", "PNA vs. SCR"))


rownames(hm_kegg) <- gsub(" - Prevotella copri", "", rownames(hm_kegg))
rownames(hm_kegg) <- gsub("; Including.*", "", rownames(hm_kegg))

ht_vert <- Heatmap(hm_kegg, cluster_rows = T, cluster_columns = T,
               name = "GO-analysis", col = col_fun,
               show_heatmap_legend = F,
               #column_split = g_samples,
               row_title_side = "right", row_title_rot = 0,
               border = TRUE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.1s", pvals[i, j]), x, y)
               },
               column_names_gp = gpar(fontsize = 11),
               column_title_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               row_title = NULL,
               width = unit(12, "cm"), height = unit(15, "cm"),

               right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes)))

ht_vert
lgd = Legend(col_fun = col_fun, title = expression("mean log"[2]*" FC"), #direction = "horizontal",
             title_gp = gpar(fontsize = 12), labels = c("-1.5", " 0"," 1.5"), legend_height = unit(6, "cm"),
             at = c(-2, 0, 2), border = "black",
             title_position = "leftcenter-rot")
draw(lgd)

svg("./metatranscriptomics_2023_04/analysis/pathway_analysis/hm_KEGG_mod.svg", width = unit(12, "cm"),  height = unit(10, "cm"))
draw(ht_vert)
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"), just = c("left", "bottom"))
dev.off()

# do pathway analysis for B. theta genes as well:
# first, get the i=link and list from keggrest
link_kegg <- keggLink("pathway", "bth")
list_kegg <- keggList("pathway", "bth")

# for link_kegg, we need to adjust the locus tags:
names(link_kegg) <- gsub("bth:", "", names(link_kegg))

# import orthologs
orthomapper <- read_delim("./metatranscriptomics_2023_04/data/reference_sequences/portho_b_theta.proteinortho.tsv", delim = "\t")
mappings_ortho <- orthomapper %>% select(b_theta_cds.fasta , b_theta_cds_orig.fasta) %>%
  filter(nchar(b_theta_cds.fasta) < 10)

names(link_kegg)  <- mappings_ortho$b_theta_cds_orig.fasta[match(names(link_kegg), mappings_ortho$b_theta_cds.fasta)]
#remove entries with NA name from link_kegg
link_kegg <- link_kegg[!is.na(names(link_kegg))]


kegg_pw_ids <- names(list_kegg)


link_kegg <- link_kegg[names(link_kegg) %in% c(rownames(de_list$b_theta$con_2_vs_0))] #remove genes not in data

link_kegg <- gsub("path:","", link_kegg)
idx_kegg <- sapply(kegg_pw_ids, function(x){
  x <- unique(names(link_kegg[link_kegg == x])) # choose all genes, except duplucates
})

# now we have the links of the pathways in p. copri.
# we can run fry gene set analysis on the log2FCs of the genes in the pathways.
l <- length(colnames(contrasts))
kegg_fry <- lapply(1:l, function(x) fry(y_list$y_btheta,idx_kegg, design_matrix, contrasts[,x]))
names(kegg_fry) <- colnames(contrasts)

# add kegg tems:
for (fryres in names(kegg_fry)) {
  kegg_fry[[fryres]][["TERM"]] <- ifelse(grepl("bth",rownames(kegg_fry[[fryres]])),
                                              list_kegg[rownames(kegg_fry[[fryres]])],
                                              rownames(kegg_fry[[fryres]]))
  kegg_fry[[fryres]][["TERM"]] <- gsub("(.*) - Bacteroides thetaiotaomicron VPI-5482",
                                            "\\1", kegg_fry[[fryres]][["TERM"]])
  write.csv(kegg_fry[[fryres]], paste("./metatranscriptomics_2023_04/analysis/pathway_analysis/btheta/", fryres, ".csv", sep = ""))
}


kegg_frysig <- lapply(kegg_fry, function(x) x[x[["FDR"]]<0.1 & x[["NGenes"]]>3,])
kegg_siggos <- c()


for (i in names(kegg_frysig)) {
  print(i)
  print(dim(kegg_frysig[[i]]))
  print(kegg_frysig[[i]][,c(1,2,4,7)])
  kegg_siggos <- c(kegg_siggos, rownames(kegg_frysig[[i]][1:15,]))  # can be modified
}

kegg_siggos <- unique(kegg_siggos[!grepl("NA", kegg_siggos)])

# now generate a heatmap with complexheatmap
idx_kegg_char <- lapply(idx_kegg, as.character)


# I create a dataframe with mean logFC values for each significant GO-term:
hm_kegg <- t(as.data.frame(lapply(idx_kegg_char[kegg_siggos], function(x){
  sapply(names(de_list$b_theta), function(y){
    mean(de_list$b_theta[[y]]$table[x,]$logFC)
  })
})))

hm_kegg <- as.data.frame(hm_kegg)

# make hm
hm_kegg <- hm_kegg[order(hm_kegg[,1], decreasing = T),]

kegg_sizes <- sapply(idx_kegg_char[rownames(hm_kegg)], function(x) length(x))

pvals <- data.frame(sapply(names(kegg_fry),
                           function(x) kegg_fry[[x]][rownames(hm_kegg),"FDR"]),
                    row.names = rownames(hm_kegg))

#select only significant ones:
pvals <-sapply(pvals, function(x) ifelse(x<0.05, x <- "*", x<-"") )

keggpws <- kegg_fry$con_2_vs_0[rownames(hm_kegg),] [["TERM"]]


rownames(hm_kegg) <- ifelse(!is.na(keggpws),keggpws, rownames(hm_kegg) )


# plot hm (save as pdf):
library(circlize)
col_fun = colorRamp2(c(-1.5,0, 1.5), c("darkblue", "white", "darkred"))

g_samples <- factor(c(rep("PNA vs. PNA0", 4), rep("SCR vs. SCR0", 4), rep("PNA vs. CON", 5), rep("PNA vs. SCR", 5)),
                    levels = c("PNA vs. PNA0", "SCR vs. SCR0", "PNA vs. CON", "PNA vs. SCR"))

ht_b_theta <- Heatmap(hm_kegg, cluster_rows = T, cluster_columns = T,
               name = "GO-analysis", col = col_fun,
               show_heatmap_legend = F,
               #column_split = g_samples,
               row_title_side = "right", row_title_rot = 0,
               border = TRUE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.1s", pvals[i, j]), x, y)
               },
               column_names_gp = gpar(fontsize = 11),
               column_title_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               row_title = NULL,
               width = unit(12, "cm"), height = unit(15, "cm"),

               right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes)))

ht_b_theta
lgd_b_theta = Legend(col_fun = col_fun, title = expression("mean log"[2]*" FC"), #direction = "horizontal",
             title_gp = gpar(fontsize = 12), labels = c("-1.5", " 0"," 1.5"), legend_height = unit(6, "cm"),
             at = c(-2, 0, 2), border = "black",
             title_position = "leftcenter-rot")
draw(lgd_b_theta)

svg("./metatranscriptomics_2023_04/analysis/pathway_analysis/hm_KEGG_b_theta.svg", width = unit(12, "cm"),  height = unit(10, "cm"))
draw(ht_b_theta)
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"), just = c("left", "bottom"))
dev.off()


# same for b. ovatus::
# Do pathway analysis in p. copri for now. use keggrest package for this.

# check whether
# get link and list of p. copri
link_kegg <- keggLink("pathway", "boa")
list_kegg <- keggList("pathway", "boa")

kegg_pw_ids <- names(list_kegg)

#rename genes, remove ones which arent in our data:
names(link_kegg) <- gsub("boa:", "", names(link_kegg)) #rename genes as locus tags

# get old locus tags:
old_lts <- read.delim("./metatranscriptomics_2023_04/data/reference_sequences/b_ovatus_old_lt.tsv", sep = "\t", header = F)
colnames(old_lts) <- c("locus_tag","old_locus_tag")

rownames(old_lts) <- old_lts$old_locus_tag

# remove ones without old lt:
link_kegg <- link_kegg[names(link_kegg) %in% old_lts$old_locus_tag]
names(link_kegg) <- old_lts[names(link_kegg),1]


link_kegg <- link_kegg[names(link_kegg) %in% c(rownames(de_list$b_ovatus$con_2_vs_0))] #remove genes not in data

link_kegg <- gsub("path:","", link_kegg)
idx_kegg <- sapply(kegg_pw_ids, function(x){
  x <- unique(names(link_kegg[link_kegg == x])) # choose all genes, except duplucates
})

# now we have the links of the pathways in p. copri.
# we can run fry gene set analysis on the log2FCs of the genes in the pathways.
l <- length(colnames(contrasts))
kegg_fry <- lapply(1:l, function(x) fry(y_list$y_bovatus,idx_kegg, design_matrix, contrasts[,x]))
names(kegg_fry) <- colnames(contrasts)

# add kegg tems:
for (fryres in names(kegg_fry)) {
  kegg_fry[[fryres]][["TERM"]] <- ifelse(grepl("oa",rownames(kegg_fry[[fryres]])),
                                              list_kegg[rownames(kegg_fry[[fryres]])],
                                              rownames(kegg_fry[[fryres]]))
  kegg_fry[[fryres]][["TERM"]] <- gsub("(.*) - Bacteroides ovatus",
                                            "\\1", kegg_fry[[fryres]][["TERM"]])
  write.csv(kegg_fry[[fryres]], paste("./metatranscriptomics_2023_04/analysis/pathway_analysis/bovatus/", fryres, ".csv", sep = ""))
}


kegg_frysig <- lapply(kegg_fry, function(x) x[x[["FDR"]]<0.1 & x[["NGenes"]]>3,])
kegg_siggos <- c()


for (i in names(kegg_frysig)) {
  print(i)
  print(dim(kegg_frysig[[i]]))
  print(kegg_frysig[[i]][,c(1,2,4,7)])
  kegg_siggos <- c(kegg_siggos, rownames(kegg_frysig[[i]][1:15,]))  # can be modified
}

kegg_siggos <- unique(kegg_siggos[!grepl("NA", kegg_siggos)])

# now generate a heatmap with complexheatmap
idx_kegg_char <- lapply(idx_kegg, as.character)


# I create a dataframe with mean logFC values for each significant GO-term:
hm_kegg <- t(as.data.frame(lapply(idx_kegg_char[kegg_siggos], function(x){
  sapply(names(de_list$b_ovatus), function(y){
    mean(de_list$b_ovatus[[y]]$table[x,]$logFC)
  })
})))

hm_kegg <- as.data.frame(hm_kegg)

# rownames(hm_kegg) <- gsub("\\.", "\\:", rownames(hm_kegg))


# make hm
hm_kegg <- hm_kegg[order(hm_kegg[,1], decreasing = T),]

kegg_sizes <- sapply(idx_kegg_char[rownames(hm_kegg)], function(x) length(x))

pvals <- data.frame(sapply(names(kegg_fry),
                           function(x) kegg_fry[[x]][rownames(hm_kegg),"FDR"]),
                    row.names = rownames(hm_kegg))

#select only significant ones:
pvals <-sapply(pvals, function(x) ifelse(x<0.05, x <- "*", x<-"") )

keggpws <- kegg_fry$con_2_vs_0[rownames(hm_kegg),] [["TERM"]]


rownames(hm_kegg) <- ifelse(!is.na(keggpws),keggpws, rownames(hm_kegg) )


# plot hm (save as pdf):
col_fun = colorRamp2(c(-1.5,0, 1.5), c("darkblue", "white", "darkred"))

ht_bovatus <- Heatmap(hm_kegg, cluster_rows = T, cluster_columns = T,
               name = "GO-analysis", col = col_fun,
               show_heatmap_legend = F,
               #column_split = g_samples,
               row_title_side = "right", row_title_rot = 0,
               border = TRUE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.1s", pvals[i, j]), x, y)
               },
               column_names_gp = gpar(fontsize = 11),
               column_title_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               row_title = NULL,
               width = unit(12, "cm"), height = unit(15, "cm"),

               right_annotation = rowAnnotation(genes = anno_barplot(kegg_sizes)))

ht_bovatus
lgd = Legend(col_fun = col_fun, title = expression("mean log"[2]*" FC"), #direction = "horizontal",
             title_gp = gpar(fontsize = 12), labels = c("-1.5", " 0"," 1.5"), legend_height = unit(6, "cm"),
             at = c(-2, 0, 2), border = "black",
             title_position = "leftcenter-rot")
draw(lgd)

svg("./metatranscriptomics_2023_04/analysis/pathway_analysis/hm_KEGG_bovatus.svg", width = unit(12, "cm"),  height = unit(10, "cm"))
draw(ht_bovatus)
draw(lgd, x = unit(2, "cm"), y = unit(10, "cm"), just = c("left", "bottom"))
dev.off()

ht_vert
ht_b_theta
ht_bovatus
draw(lgd, x = unit(0.1, "cm"), y = unit(10, "cm"), just = c("left", "bottom"))


ht_list



