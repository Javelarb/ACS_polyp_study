#library(DESeq2)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(rstatix)
library(vegan)
library(IHW)
library(ggpubr)
library(lemon)
library(rfPermute)
library(pROC)
library(matrixStats)
library(svglite)

setwd("/media/julio/Storage/CRC/github/")
options(scipen = 100000)
source("/media/julio/Storage/Software/ANCOM-master/programs/ancom.R")

#register(MulticoreParam(4)) #Use four cores

metadata <- read.delim("Shotgun_metadata.tsv", row.names=1, comment.char="#", check.names = F)
contig_table <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/contig_table.txt", row.names=1, check.names = F)
anots <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/eggnog.emapper.annotations", comment.char = "#", header=FALSE, row.names = 1)
CAZyme_annotations = read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/dbcan/overview.txt", row.names=1, check.names = F) %>% filter(`#ofTools` == 3)

metadata$Patient <- as.factor(metadata$Patient)
contig_table <- contig_table[rowSums(contig_table) >= 10,] #Get rid of genes with low read counts

#gene lengths
gene_lengths = read.delim("/media/julio/Storage/CRC/DataDryad/gene_lengths.txt", header=FALSE, comment.char="#") %>% .[,1:2]
colnames(gene_lengths) = c("ORF", "gene_length")
gene_lengths$gene_length = gene_lengths$gene_length/1000 #in Kb

#From microbe census
genome_equivs = read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/mic_cense_genome_equivs_merged.txt")

#Divide reads hitting to ORF by gene length & genome equivalents.
RPKG = merge(gene_lengths, contig_table, by.x = "ORF", by.y = "row.names") %>% column_to_rownames(var = "ORF")
RPKG = RPKG[,-1]/RPKG$gene_length
RPKG2 = merge(genome_equivs, t(RPKG), by.x = "V1", by = "row.names") %>% column_to_rownames(var = "V1")
RPKG2 = as.data.frame(t(RPKG2[,-1]/RPKG2$V2))

eggnog_table = merge(anots, RPKG2, by = "row.names")
cazyme_table = merge(CAZyme_annotations, RPKG2, by = "row.names") %>% column_to_rownames(var = "Row.names")

#### Cazyme diversity ####
cazy_rich = as.data.frame(specnumber(t(cazyme_table[,6:ncol(cazyme_table)]))) %>% merge(., metadata, by = "row.names")
cazy_rich$Shannon = diversity(t(cazyme_table[,6:ncol(cazyme_table)]))
cazy_rich = cazy_rich %>% filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), SampleType == "aspirate")

ggplot(data = cazy_rich) +
  aes(x = as.factor(MostMalignantPolypType), y = cazy_rich$`specnumber(t(cazyme_table[, 6:ncol(cazyme_table)]))`, fill = as.factor(MostMalignantPolypType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .25) +
  theme_bw() +
  labs(x = NULL, y = "CAZyme richness") +
  scale_fill_manual(values=c("forestgreen", "steelblue", "firebrick3")) +
  theme(legend.position = "none")

cazy_nmds = metaMDS(t(cazyme_table[,6:ncol(cazyme_table)]), trymax = 999, parallel = 32, k = 3)
cazy_beta = as.data.frame(cazy_nmds$points[,1:2]) %>% merge(., metadata, by = "row.names") %>% filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), SampleType == "aspirate")

ggplot(data = cazy_beta) +
  aes(x = MDS1, y = MDS2, color = MostMalignantPolypType) +
  theme_bw() +
  geom_text(label = cazy_beta$Patient, size = 4) +
  stat_ellipse(linetype = 2, aes(group = MostMalignantPolypType), show.legend = F) +
  annotate("text", x = 1.35, y = -1, size = 4, label = bquote("K = 3, Stress ="~.(round(cazy_nmds$stress, digits = 2)))) +
  scale_color_manual(values=c("forestgreen", "steelblue", "firebrick3")) +
  labs(title = "CAZyme-beta diversity", color = "Polyp type")

#Permanova
cazy_merged = merge(metadata, t(cazyme_table[,6:ncol(cazyme_table)]), by = "row.names")
cazy_merged[cazy_merged == "unknown"] <- NA
cazy_merged = cazy_merged %>% filter(complete.cases(.), !MostMalignantPolypType == "adenocarcinoma", SampleType == "aspirate") %>% column_to_rownames(var = "Row.names")

adonis2(formula = cazy_merged[,16:ncol(cazy_merged)] ~ as.numeric(cazy_merged$BMI) + 
         as.numeric(cazy_merged$Age) + cazy_merged$Ethnicity + cazy_merged$Gender + 
         (cazy_merged$MostMalignantPolypType /as.character(cazy_merged$Patient) ), 
       data = cazy_merged, method = "bray", permutations = 999, parallel = 32, strata = as.factor(cazy_merged$Plate))

#### Bar plot ####
#brite <- read.delim("/media/julio/Storage/DBs/brite.txt", header=FALSE)

#lastValue <- function(x) tail(x[!is.na(x)], 1)
#brite_list <- tidyr::separate(anots, col = V14, into = c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10", 
#                                                         "L11","L12","L13", "L14","L15","L16","L17","L18","L19",
#                                                         "L20","L21","L22"), sep = "\\,", remove = T, extra = "drop")  %>% .[,13:(13+21)]
#brite_list <- as.data.frame(apply(brite_list, 1, lastValue))
#brite_list <- subset(brite_list, !(brite_list == "")) 

#colnames(brite_list) <- "V1"
#brite_list$Row.names <- rownames(brite_list)
#brite_list <- merge(brite_list, brite, by = "V1")

#contig_relab = t(RPKG2)/rowSums(t(RPKG2))
#contig_relab <- contig_relab %>% reshape2::melt()
#contig_relab <- contig_relab[!(contig_relab$value == 0),]

#plot_df <- subset(contig_relab, contig_relab$Var2 %in% brite_list$Row.names)
#plot_df2 <- merge(plot_df,brite_list[,c("Row.names","V2")], by.x = "Var2", by.y = "Row.names") %>% merge(., metadata, by.x = "Var1", by.y = "row.names")
#plot_df2$V2 = as.character(plot_df2$V2)

#Take top 10 genes.
#top_genes <- group_by(plot_df2, V2) %>% summarise(., top_genes_tmp = sum(value)) %>% arrange(., desc(top_genes_tmp)) %>% slice(., 1:10)
#high_abundance <- split(top_genes$V2, 1:NROW(top_genes))

#Change non top hits to other.
#plot_df2$V2[plot_df2$V2 %in% high_abundance != "TRUE"] <- "Other"
#plot_df2 <- plot_df2[order(plot_df2$V2),] #Re order
#plot_df2 <- rbind(plot_df2[!(plot_df2$V2 == "Other"),],plot_df2[(plot_df2$V2 == "Other"),]) #Move other to bottom
#plot_df2$V2 <- factor(plot_df2$V2, levels = unique(plot_df2$V2)) #Fix the order

#IF you remove enzyme then you need this next line.
#relab <- aggregate(plot_df2$value, by=list(Var1=plot_df2$Var1), FUN=sum)
#plot_df2 <- merge(plot_df2, relab, by = "Var1")
#plot_df2$relab <- plot_df2$value/plot_df2$x
#plot_df2 = plot_df2 %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))

#sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947", "gray")

#brite_plot = ggplot(data = plot_df2, aes(x = Var1, weight = relab, fill = V2)) +
#  geom_bar(width = 1, color = "black", size = .2) +
#  theme_classic(base_size = 16) +
#  facet_wrap(SampleType~MostMalignantPolypType, scales = "free", labeller = 
#               labeller(SampleType = c(`aspirate` = "Aspirate", `fecal` = "Fecal", `lavage` = "Lavage"), 
#                        MostMalignantPolypType = c(`healthy`="Polyp free", `TA`="TA", `serrated`="Serrated"))) +
#  scale_fill_manual(values = sarah_color) +
#  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 10), 
#        strip.background = element_rect(fill="lightblue")) +
#  labs(x = NULL,
#       y = "Relative abundance", fill = "KEGG Brite pathway", title = "Shotgun - Individuals: 104")
#brite_plot

#ggsave("Brite_plot.png", plot = brite_plot, device = "png", units = "in", dpi = 300, height = 10, width = 14)

#Alpha diversity
eggnog_alpha = as.data.frame(diversity(t(eggnog_table[,23:ncol(eggnog_table)]), index = "shannon")) %>% 
  rename(Shannon = `diversity(t(eggnog_table[, 23:ncol(eggnog_table)]), index = "shannon")`)
eggnog_alpha$Richness = specnumber(t(eggnog_table[,23:ncol(eggnog_table)]))
eggnog_alpha = merge(eggnog_alpha, metadata, by = "row.names") %>% filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE")

fig5b1 = ggplot(data = eggnog_alpha) +
  aes(x = SampleType, y = Shannon, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .1), alpha = .5, size = .1) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), labels = c("Polyp free", "Serrated polyp", "Tubular adenoma")) +
  scale_x_discrete(labels = c("Aspirate", "Fecal", "Lavage")) +
  theme(legend.position = "none") +
  ylim(8, 11.5)

alpha_lm <- NULL
alpha_lm$x <- as.numeric(eggnog_alpha$Shannon)
alpha_lm$r <- as.numeric(eggnog_alpha$Richness)
alpha_lm$y <- as.factor(eggnog_alpha$SampleType)
alpha_lm$z <- as.factor(eggnog_alpha$MostMalignantPolypType)
alpha_lm$p <- as.factor(eggnog_alpha$Plate)
alpha_lm$i <- as.factor(eggnog_alpha$Patient)
alpha_lm <- as.data.frame(alpha_lm)
alpha_lm <- within(alpha_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))
summary(lme(r ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))

fig5b2 = ggplot(data = eggnog_alpha) +
  aes(x = SampleType, y = Richness, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .1), alpha = .5, size = .1) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), labels = c("Polyp free (78)", "Serrated polyp (59)", "Tubular adenoma (75)")) +
  scale_x_discrete(labels = c("Aspirate", "Fecal", "Lavage")) +
  theme(legend.position = "none") +
  ylim(0,150000)

fig5b = plot_grid(fig5b1, fig5b2, rel_widths = c(1,1), nrow = 1)

#PCoA
eggnog_bray = vegdist(t(eggnog_table[,23:ncol(eggnog_table)]), method = "bray")
eggnog_pcoa = cmdscale(eggnog_bray, eig = T, k = nrow(t(eggnog_table[,23:ncol(eggnog_table)]))-1, add = T)
eggnog_eig = eigenvals(eggnog_pcoa)
pcoa_var = eggnog_eig/sum(eggnog_eig)
pcoa_var

eggnog_points = as.data.frame(eggnog_pcoa$points[,1:2]) %>% merge(., metadata, by = "row.names") %>% 
  dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE")

plot5c = ggplot(data = eggnog_points) +
  aes(x = V1, y = V2) +
  theme_bw() +
  geom_point(aes(pch = SampleType, fill = MostMalignantPolypType), size = 2, alpha = 3/4) +
  stat_ellipse(linetype = 2, aes(group = SampleType, color = SampleType), show.legend = F) +
  scale_color_manual(values=c("gold", "darkorange4", "plum2")) +
  scale_shape_manual(values = c(21,22,24), name = "Sample type", labels = c("Aspirate (157)", "Fecal (35)", "Lavage (20)")) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), 
                    labels = c("Polyp free (78)", "Serrated polyp (59)", "Tubular adenoma (75)")) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  labs(color = "Subject type", x = bquote("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) "), y = bquote("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) "))
plot5c

#Permanova
eggnog_merged = merge(metadata, as.matrix(eggnog_bray), by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))
eggnog_merged[eggnog_merged == "unknown"] <- NA
eggnog_merged = eggnog_merged %>% filter(complete.cases(.), !MostMalignantPolypType == "adenocarcinoma")

adonis2(formula = eggnog_merged[,15:ncol(eggnog_merged)] ~ as.numeric(eggnog_merged$BMI) + 
         as.numeric(eggnog_merged$Age) + eggnog_merged$Ethnicity + eggnog_merged$Gender + 
         (eggnog_merged$MostMalignantPolypType /as.character(eggnog_merged$Patient) / eggnog_merged$SampleType), 
       data = eggnog_merged, method = "bray", permutations = 999, parallel = 32, strata = as.factor(eggnog_merged$Plate))

#### Differential abundance (Pre-reviewer comments) ####
#merged_cont_table <- merge(metadata, t(contig_table), by = "row.names")
#merged_cont_table <- merged_cont_table[!(merged_cont_table$MostMalignantPolypType %in% c("adenocarcinoma", "unknown")),]
#merged_cont_table <- merged_cont_table[!(row.names(merged_cont_table) %in% c("0837")),]

#rownames(merged_cont_table) <- merged_cont_table$Row.names

#l_asps <- merged_cont_table[(merged_cont_table$SampleType == "aspirate"),] %>% .[(.[,6] == "left"),]
#r_asps <- merged_cont_table[(merged_cont_table$SampleType == "aspirate"),] %>% .[(.[,6] == "right"),]
#fecal <- merged_cont_table[(merged_cont_table$SampleType == "fecal"),]
#sampletype <- rbind(r_asps[!(r_asps$Patient %in% fecal$Patient),], fecal) 

#l_asps_deseq <- DESeqDataSetFromMatrix(countData = t(l_asps[,(ncol(metadata)+2):ncol(l_asps)]), 
#                                       colData = l_asps[,2:(ncol(metadata)+1)],
#                                       design = ~ MostMalignantPolypType)

#r_asps_deseq <- DESeqDataSetFromMatrix(countData = t(r_asps[,(ncol(metadata)+2):ncol(r_asps)]), 
#                                       colData = r_asps[,2:(ncol(metadata)+1)],
#                                       design = ~ MostMalignantPolypType)

#fecal_deseq <- DESeqDataSetFromMatrix(countData = t(fecal[,(ncol(metadata)+2):ncol(fecal)]), 
#                                      colData = fecal[,2:(ncol(metadata)+1)],
#                                      design = ~ MostMalignantPolypType)

#sample_deseq <- DESeqDataSetFromMatrix(countData = t(sampletype[,(ncol(metadata)+2):ncol(sampletype)]), 
#                                       colData = sampletype[,2:(ncol(metadata)+1)],
#                                       design = ~ SampleType)

#l_asps_deseq_analysis <- DESeq(l_asps_deseq, parallel = T)
#r_asps_deseq_analysis <- DESeq(r_asps_deseq, parallel = T)
#fecal_deseq_analysis <- DESeq(fecal_deseq, parallel = T)
#sample_deseq_analysis <- DESeq(fecal_deseq, parallel = T)

#l_asps_results1 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
#l_asps_results2 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
#l_asps_results3 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)

#r_asps_results1 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
#r_asps_results2 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
#r_asps_results3 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)

#fecal_results1 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
#fecal_results2 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
#fecal_results3 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)

#sample_deseq_results <- results(sample_deseq_analysis, alpha = 0.05) %>% subset(., padj < 0.05)

#summary(l_asps_results1)
#summary(l_asps_results2)
#summary(l_asps_results3)

#summary(r_asps_results1)
#summary(r_asps_results2)
#summary(r_asps_results3)

#summary(fecal_results1)
#summary(fecal_results2)
#summary(fecal_results3)

#summary(sample_deseq_results)

#l_asps_summary1 <- bind_rows(l_asps_results1@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results1@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#l_asps_summary2 <- bind_rows(l_asps_results2@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results2@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#l_asps_summary3 <- bind_rows(l_asps_results3@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results3@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)

#r_asps_summary1 <- bind_rows(r_asps_results1@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results1@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#r_asps_summary2 <- bind_rows(r_asps_results2@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results2@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#r_asps_summary3 <- bind_rows(r_asps_results3@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results3@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("aspirate") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)

#fecal_summary1 <- bind_rows(fecal_results1@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results1@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("fecal") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#fecal_summary2 <- bind_rows(fecal_results2@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results2@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("fecal") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)
#fecal_summary3 <- bind_rows(fecal_results3@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results3@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("fecal") %>%
#  dplyr::rename(comparison = 7, SampleType = 8)

#sample_summary <- bind_rows(sample_deseq_results@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(sample_deseq_results@rownames) %>% 
#  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("aspirate to fecal") %>% rename(comparison = 6)

#master_summary <- bind_rows(l_asps_summary1, l_asps_summary2, l_asps_summary3, r_asps_summary1, r_asps_summary2, r_asps_summary3)
#master_summary <- master_summary[!master_summary$V5 == "Planctomycetes",]

#master_summary2 <- subset(master_summary, master_summary$comparison == "serrated to healthy")
#master_summary2$color <- master_summary2$log2FoldChange > 0
#labels1 <- master_summary2 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]

#ggplot(data = master_summary2) +
#  theme_classic(base_size = 14, base_line_size = 1) +
#  geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), size = 3, pch = 21, alpha = .5) +
#  geom_text_repel(data = master_summary2, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
#  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
#  labs(fill ="Subject type") +
#  scale_fill_manual(values = c("forestgreen", "firebrick3"), labels = c("Healthy", "Tubular adenoma"))
#plot5c1

#master_summary3 <- subset(master_summary, master_summary$comparison == "serrated to healthy")
#master_summary3$color <- master_summary3$log2FoldChange > 0
#labels2 <- master_summary3 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]

#ggplot(data = master_summary3) +
#  theme_classic(base_size = 14, base_line_size = 1) +
#  geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), pch = 21, size = 3, alpha = .5) +
#  geom_text_repel(data = master_summary3, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
#  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
#  labs(fill ="Subject type") +
#  scale_fill_manual(values = c("forestgreen", "steelblue3"), labels =  c("Healthy", "Serrated"))
#  theme(legend.position = "none")

#plot_grid(plot5c1, plot5c2, rel_widths = c(.65,1))
#ggsave("figure_5c.svg", plot = plot5c3, device = "svg", units = "in", dpi = 1000, height = 3.5, width = 5)

#master_summary4 <- subset(master_summary, master_summary$comparison == "TA to serrated")
#master_summary4$color <- master_summary4$log2FoldChange > 0
#labels3 <- master_summary4 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]

#ggplot(data = master_summary4) +
#  theme_classic(base_size = 14, base_line_size = 1) +
#  geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), pch = 21, size = 3, alpha = .5) +
#  geom_text_repel(data = master_summary4, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
#  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
#  labs(fill ="Subject type") +
#  scale_fill_manual(values = c("steelblue3", "firebrick3"), labels = c("Serrated", "Tubular adenoma"))

#Checked gene abundances to eliminate false positives.
#Should be present in at least 10% of subjects.
#k111_458716_2 (tonB), k111_2791185_4 (GH31), k111_1237370_2 (sdaAA)

#for (i in unique(master_summary$Row.names)) {
#  assign(paste0("L_", i), plotCounts(l_asps_deseq_analysis, gene = i, intgroup = "MostMalignantPolypType", pc = 0, returnData = T))
#  assign(paste0("R_", i), plotCounts(r_asps_deseq_analysis, gene = i, intgroup = "MostMalignantPolypType", pc = 0, returnData = T))
#  assign(paste0("Both_", i), rbind(get(paste0("L_", i)), get(paste0("R_", i))) %>% merge(., metadata, by = "row.names"))
  
#  assign(paste0("plot_", i), ggplot(data = get(paste0("Both_", i))) +
#           aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
#           geom_boxplot(outlier.shape = NA) +
#           geom_jitter(alpha = .5) +
#           labs(title = i, x = NULL, y = "DESeq2 normalized counts") +
#           scale_x_discrete(labels = c("Healthy (63)", "Serrated (45)", "TA (49)")) +
#           scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
#           theme_classic() +
#           theme(legend.position = "none"))
#}

#tonB <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_458716_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
#              plotCounts(r_asps_deseq_analysis, gene = "k111_458716_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
#  merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
#tonB <- tonB[!duplicated(tonB),]
#dunn_test(data = tonB, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")

#plot5d <- ggplot(data = Both_k111_458716_2) +
#  aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
#  geom_boxplot(outlier.shape = NA, lwd = 1.5) +
#  geom_jitter(width = .25) +
#  labs(title = "", subtitle = expression(italic("tonB")), x = NULL, y = NULL) +
#  scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
#  scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
#  theme_classic() +
#  theme(legend.position = "none") +
#  scale_y_log10() +
#  theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))

#GH31 <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_2791185_4", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
#              plotCounts(r_asps_deseq_analysis, gene = "k111_2791185_4", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
#  merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
#GH31 <- GH31[!duplicated(GH31),]
#dunn_test(data = GH31, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")

#plot5c <- ggplot(data = Both_k111_2791185_4) +
#  aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
#  geom_boxplot(outlier.shape = NA, lwd = 1.5) +
#  geom_jitter(width = .25) +
#  labs(title = "", subtitle = "GH31", x = NULL, y = NULL) +
#  scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
#  scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
#  theme_classic() +
#  theme(legend.position = "none") +
#  scale_y_log10() +
# theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))

#sdaAA <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_1237370_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
#               plotCounts(r_asps_deseq_analysis, gene = "k111_1237370_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
#  merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
#sdaAA <- sdaAA[!duplicated(sdaAA),]
#dunn_test(data = sdaAA, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")

#plot5b <- ggplot(data = Both_k111_1237370_2) +
#  aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
#  geom_boxplot(outlier.shape = NA, lwd = 1.5) +
#  geom_jitter(width = .25) +
#  labs(title = "Mucosal aspirates only", subtitle = expression(italic("sdaA")), x = NULL, y = "DESeq2 normalized reads") +
#  scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
#  scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
#  theme_classic() +
#  theme(legend.position = "none") +
#  scale_y_log10() +
#  theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))

#plot5b_d <- plot_grid(plot5b,plot5c,plot5d, nrow = 1, rel_widths = c(1,.95,.95))
#plot5b_d

#ggsave(filename = "fig5.svg", plot = plot_grid(plot5a, plot5b_d, ncol = 1, rel_heights = c(.6,.4)), device = "svg", dpi = 300, height = 8, width = 8.5)

#### IHW ####
#Filter OTUs that are not present across 33% of individuals.
#RPKG2 = contig_table[1:50000,] #Delete this when you want to run the full dataset (takes forever)

RPKG2_merged = as.matrix(t(RPKG2)) %>% merge(metadata, ., by = "row.names") 
Asps_only = RPKG2_merged %>% filter(SampleType == "aspirate") %>% column_to_rownames(var = "Row.names") %>% select(!Patient:InFinalAnalysis)

variance_filter = as.matrix(Asps_only) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(Asps_only)/3)) #33% filter

#Remove columns which appear in the variance filter list
Asps_filtered = Asps_only[,!(colnames(Asps_only) %in% variance_filter$Var2)] %>% merge(metadata, ., by = "row.names")

#Differential abundance analysis will be performed using kruskal wallis test with independent hypothesis weighting

#Healthy to tubular adenoma comparison
healthy_v_TA = Asps_filtered %>% filter(!MostMalignantPolypType == "serrated") %>% 
  group_by(as.factor(Patient)) %>% summarise(across((ncol(metadata)+2):ncol(Asps_filtered), mean)) %>% 
  remove_rownames() %>% merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% 
  .[!duplicated(.),]

KW_healthy_v_TA = list()

for (i in 3:ncol(healthy_v_TA)) {
  tmp = kruskal.test(healthy_v_TA[,i]~healthy_v_TA$MostMalignantPolypType)
  tmp$taxa = names(healthy_v_TA)[i]
  KW_healthy_v_TA[[length(KW_healthy_v_TA)+1]] = tmp
}

pvals_KW_healthy_v_TA = as.data.frame(tibble(taxa = map(KW_healthy_v_TA, "taxa"), pval = map(KW_healthy_v_TA, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pvals_KW_healthy_v_TA2 = as.data.frame(colSums(healthy_v_TA[,3:ncol(healthy_v_TA)])) %>%
  merge(., pvals_KW_healthy_v_TA, by.x = "row.names", by.y = "taxa") %>% column_to_rownames(var = "Row.names")
names(pvals_KW_healthy_v_TA2)[1] = "rpkg_sum"
pvals_KW_healthy_v_TA2$pval = as.numeric(pvals_KW_healthy_v_TA2$pval)

IHW_healthy_v_TA = ihw(pvalues = pvals_KW_healthy_v_TA2$pval, covariates = pvals_KW_healthy_v_TA2$rpkg_sum, alpha = 0.05, null_proportion = T, nbins = 5)
IHW_healthy_v_TA@df$taxa = rownames(pvals_KW_healthy_v_TA2)
IHW_healthy_v_TA@df$group1 = "healthy"
IHW_healthy_v_TA@df$group2 = "TA"

#Healthy to serrated comparison
healthy_v_ser = Asps_filtered %>% filter(!MostMalignantPolypType == "TA") %>% group_by(as.factor(Patient)) %>% 
  summarise(across((ncol(metadata)+2):ncol(Asps_filtered), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]

KW_healthy_v_ser = list()

for (i in 3:ncol(healthy_v_ser)) {
  tmp = kruskal.test(healthy_v_ser[,i]~healthy_v_ser$MostMalignantPolypType)
  tmp$taxa = names(healthy_v_ser)[i]
  KW_healthy_v_ser[[length(KW_healthy_v_ser)+1]] = tmp
}

pvals_KW_healthy_v_ser = as.data.frame(tibble(taxa = map(KW_healthy_v_ser, "taxa"), pval = map(KW_healthy_v_ser, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pvals_KW_healthy_v_ser2 = as.data.frame(colSums(healthy_v_ser[,3:ncol(healthy_v_ser)])) %>%
  merge(., pvals_KW_healthy_v_ser, by.x = "row.names", by.y = "taxa") %>% column_to_rownames(var = "Row.names")
names(pvals_KW_healthy_v_ser2)[1] = "rpkg_sum"
pvals_KW_healthy_v_ser2$pval = as.numeric(pvals_KW_healthy_v_ser2$pval)

IHW_healthy_v_ser = ihw(pvalues = pvals_KW_healthy_v_ser2$pval, covariates = pvals_KW_healthy_v_ser2$rpkg_sum, alpha = 0.05, null_proportion = T, nbins = 5)
IHW_healthy_v_ser@df$taxa = rownames(pvals_KW_healthy_v_ser2)
IHW_healthy_v_ser@df$group1 = "healthy"
IHW_healthy_v_ser@df$group2 = "serrated"

#TA to serrated comparison
TA_v_ser = Asps_filtered %>% filter(!MostMalignantPolypType == "healthy") %>% group_by(as.factor(Patient)) %>% 
  summarise(across((ncol(metadata)+2):ncol(Asps_filtered), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]

KW_TA_v_ser = list()

for (i in 3:ncol(TA_v_ser)) {
  tmp = kruskal.test(TA_v_ser[,i]~TA_v_ser$MostMalignantPolypType)
  tmp$taxa = names(TA_v_ser)[i]
  KW_TA_v_ser[[length(KW_TA_v_ser)+1]] = tmp
}

pvals_KW_TA_v_ser = as.data.frame(tibble(taxa = map(KW_TA_v_ser, "taxa"), pval = map(KW_TA_v_ser, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pvals_KW_TA_v_ser2 = as.data.frame(colSums(TA_v_ser[,3:ncol(TA_v_ser)])) %>%
  merge(., pvals_KW_TA_v_ser, by.x = "row.names", by.y = "taxa") %>% column_to_rownames(var = "Row.names")
names(pvals_KW_TA_v_ser2)[1] = "rpkg_sum"
pvals_KW_TA_v_ser2$pval = as.numeric(pvals_KW_TA_v_ser2$pval)

IHW_TA_v_ser = ihw(pvalues = pvals_KW_TA_v_ser2$pval, covariates = pvals_KW_TA_v_ser2$rpkg_sum, alpha = 0.05, null_proportion = T, nbins = 5)
IHW_TA_v_ser@df$taxa = rownames(pvals_KW_TA_v_ser2)
IHW_TA_v_ser@df$group1 = "TA"
IHW_TA_v_ser@df$group2 = "serrated"

#Visualizing
#The p values reported are not FDR adjusted but are for multiple comparisons (Sort of like a Dunn's post-hoc test)
Polyp_type_DAb = bind_rows(IHW_healthy_v_TA@df, IHW_healthy_v_ser@df, IHW_TA_v_ser@df) %>% 
  mutate(final_pval = pvalue*3) #%>% filter(final_pval <= 0.05)

RPKG_means = RPKG2_merged %>% filter(SampleType == "aspirate", !MostMalignantPolypType %in% c("adenocarcinoma", "unknown")) %>% 
  select(!variance_filter$Var2) %>% group_by(MostMalignantPolypType) %>% 
  summarise(across((ncol(metadata)+2):(ncol(.)-1), mean)) %>% column_to_rownames(var = "MostMalignantPolypType")
RPKG_means = as.data.frame(t(RPKG_means)) + 0.01
RPKG_means$names = rownames(RPKG_means)

#Healthy vs TA
FC_healthy_TA = log2(as.data.frame(RPKG_means$TA/RPKG_means$healthy)) %>% dplyr::rename(Log2FoldChange = 1) %>% 
  magrittr::set_rownames(rownames(RPKG_means))
FC_healthy_TA = Polyp_type_DAb %>% filter(group1 == "healthy", group2 == "TA") %>% merge(., FC_healthy_TA, by.x = "taxa", by.y = "row.names") %>% mutate(color = case_when(
  final_pval < 0.05 & Log2FoldChange > 0 ~ "color1",
  final_pval < 0.05 & Log2FoldChange < 0 ~ "color2",
  final_pval > 0.05 ~ "color3"
))

volc1 = ggplot(data = FC_healthy_TA) +
  aes(x = Log2FoldChange, y = -log10(final_pval), color = color) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.5, size = .25) + 
  scale_color_manual(values = c("firebrick3", "forestgreen", "gray")) +
  labs(title = "Polyp free vs. Tubular adenoma", x = expression("log"[2]*" fold change"), y = expression("-log"[10]*"(p-value)")) +
  theme_bw() +
  #geom_text_repel(data = volcano_labels[1:10,], aes(x = log2FoldChange, y = -log10(pvalue), label = V6), color = "black", size = 3) +
  theme(legend.position = "none") +
  annotate("text", x = -4.5, y = 3.25, size = 4, hjust = 0, label = bquote("Total:"~.(sum(FC_healthy_TA$final_pval < 0.05)))) +
  annotate("text", x = -4.5, y = 3.05, size = 4, hjust = 0, label = bquote("Positive:"~.(sum(FC_healthy_TA$final_pval < 0.05 & FC_healthy_TA$Log2FoldChange > 0)))) +
  annotate("text", x = -4.5, y = 2.85, size = 4, hjust = 0, label = bquote("Negative:"~.(sum(FC_healthy_TA$final_pval < 0.05 & FC_healthy_TA$Log2FoldChange < 0))))
volc1

#Healthy vs serrated
FC_healthy_serrated = log2(as.data.frame(RPKG_means$serrated/RPKG_means$healthy)) %>% dplyr::rename(Log2FoldChange = 1) %>% 
  magrittr::set_rownames(rownames(RPKG_means))
FC_healthy_serrated = Polyp_type_DAb %>% filter(group1 == "healthy", group2 == "serrated") %>% merge(., FC_healthy_serrated, by.x = "taxa", by.y = "row.names") %>% mutate(color = case_when(
  final_pval < 0.05 & Log2FoldChange > 0 ~ "color1",
  final_pval < 0.05 & Log2FoldChange < 0 ~ "color2",
  final_pval > 0.05 ~ "color3"
))

volc2 = ggplot(data = FC_healthy_serrated) +
  aes(x = Log2FoldChange, y = -log10(final_pval), color = color) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.5, size = .25) + 
  scale_color_manual(values = c("steelblue3", "forestgreen", "gray")) +
  labs(title = "Polyp free vs. Serrated polyp", x = expression("log"[2]*" fold change"), y = expression("-log"[10]*"(p-value)")) +
  theme_bw() +
  #geom_text_repel(data = volcano_labels[1:10,], aes(x = log2FoldChange, y = -log10(pvalue), label = V6), color = "black", size = 3) +
  theme(legend.position = "none") +
  annotate("text", x = 3.9, y = 4, size = 4, hjust = 1, label = bquote("Total:"~.(sum(FC_healthy_serrated$final_pval < 0.05)))) +
  annotate("text", x = 3.9, y = 3.75, size = 4, hjust = 1, label = bquote("Positive:"~.(sum(FC_healthy_serrated$final_pval < 0.05 & FC_healthy_serrated$Log2FoldChange > 0)))) +
  annotate("text", x = 3.9, y = 3.5, size = 4, hjust = 1, label = bquote("Negative:"~.(sum(FC_healthy_serrated$final_pval < 0.05 & FC_healthy_serrated$Log2FoldChange < 0))))

#TA vs serrated
FC_TA_serrated = log2(as.data.frame(RPKG_means$TA/RPKG_means$serrated)) %>% dplyr::rename(Log2FoldChange = 1) %>% 
  magrittr::set_rownames(rownames(RPKG_means))
FC_TA_serrated = Polyp_type_DAb %>% filter(group1 == "TA", group2 == "serrated") %>% merge(., FC_TA_serrated, by.x = "taxa", by.y = "row.names") %>% mutate(color = case_when(
  final_pval < 0.05 & Log2FoldChange > 0 ~ "color1",
  final_pval < 0.05 & Log2FoldChange < 0 ~ "color2",
  final_pval > 0.05 ~ "color3"
))

volc3 = ggplot(data = FC_TA_serrated) +
  aes(x = Log2FoldChange, y = -log10(final_pval), color = color) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.5, size = .25) + 
  scale_color_manual(values = c("firebrick3", "steelblue3", "gray")) +
  labs(title = "Serrated polyp \n vs. Tubular adenoma", x = expression("log"[2]*" fold change"), y = expression("-log"[10]*"(p-value)")) +
  theme_bw() +
  #geom_text_repel(data = volcano_labels[1:10,], aes(x = log2FoldChange, y = -log10(pvalue), label = V6), color = "black", size = 3) +
  theme(legend.position = "none") +
  annotate("text", x = -4.7, y = 3.75, size = 4, hjust = 0, label = bquote("Total:"~.(sum(FC_TA_serrated$final_pval < 0.05)))) +
  annotate("text", x = -4.7, y = 3.5, size = 4, hjust = 0, label = bquote("Positive:"~.(sum(FC_TA_serrated$final_pval < 0.05 & FC_TA_serrated$Log2FoldChange > 0)))) +
  annotate("text", x = -4.7, y = 3.25, size = 4, hjust = 0, label = bquote("Negative:"~.(sum(FC_TA_serrated$final_pval < 0.05 & FC_TA_serrated$Log2FoldChange < 0))))

DAb_supplemental_figure1 = ggarrange(volc1, volc2, volc3, nrow = 1, labels = c("a.", "b.", "c."))
DAb_supplemental_figure1

# Taxonomy visualizations
# Healthy v tubular adenoma
taxa_ht = FC_healthy_TA %>% filter(final_pval < 0.05) %>% merge(., anots, by.x = "taxa", by.y = "row.names")
taxa_ht2 = taxa_ht %>% group_by(V5) %>% summarise(up = sum(Log2FoldChange > 0), down = sum(Log2FoldChange < 0))
taxa_ht2$sum = rowSums(taxa_ht2[,-1])
taxa_ht2 = taxa_ht2[rev(order(abs(taxa_ht2$sum)))[1:10],] #%>% arrange(sum)

taxa_ht3 = as.data.frame(reshape2::melt(taxa_ht2)) %>% filter(!variable == "sum")
taxa_ht3$V5 = factor(taxa_ht3$V5, levels = unique(taxa_ht3$V5))

fig5d1 = ggplot(data = taxa_ht3, 
       aes(x = ifelse(test = variable == "up", yes = -value, no = value), y = V5, fill = variable)) +
  geom_col(color = "black") +
  geom_text(aes(label = value, x= ifelse(test = variable == "up", yes = -value-3, no = value+3))) +
  labs(x = "Number of differentially \n abundant genes", y = NULL, title = "Polyp free \n vs. Tubular adenoma", fill = NULL) +
  scale_fill_manual(values = c("firebrick3","forestgreen"), labels = c("Tubular adenoma", "Polyp free")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_symmetric(labels = abs) +
  coord_flip()

# Healthy v serrated polyp
taxa_hs = FC_healthy_serrated %>% filter(final_pval < 0.05) %>% merge(., anots, by.x = "taxa", by.y = "row.names")
taxa_hs2 = taxa_hs %>% group_by(V5) %>% summarise(up = sum(Log2FoldChange > 0), down = sum(Log2FoldChange < 0))
taxa_hs2$sum = rowSums(taxa_hs2[,-1])
taxa_hs2 = taxa_hs2[rev(order(abs(taxa_hs2$sum)))[1:10],] #%>% arrange(sum)

taxa_hs3 = as.data.frame(reshape2::melt(taxa_hs2)) %>% filter(!variable == "sum")
taxa_hs3$V5 = factor(taxa_hs3$V5, levels = unique(taxa_hs3$V5))

fig5d2 = ggplot(data = taxa_hs3, 
       aes(x = ifelse(test = variable == "up", yes = -value, no = value), y = V5, fill = variable)) +
  geom_col(color = "black") +
  geom_text(aes(label = value, x= ifelse(test = variable == "up", yes = -value-15, no = value+15))) +
  labs(x = NULL, y = NULL, title = "Polyp free \n vs. Serrated polyp", fill = NULL) +
  scale_fill_manual(values = c("steelblue3","forestgreen"), labels = c("Serrated polyp", "Polyp free")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_symmetric(labels = abs) +
  coord_flip()

# TA v serrateed polyp
taxa_ts = FC_TA_serrated %>% filter(final_pval < 0.05) %>% merge(., anots, by.x = "taxa", by.y = "row.names")
taxa_ts2 = taxa_ts %>% group_by(V5) %>% summarise(up = sum(Log2FoldChange > 0), down = sum(Log2FoldChange < 0))
taxa_ts2$sum = rowSums(taxa_ts2[,-1])
taxa_ts2 = taxa_ts2[rev(order(abs(taxa_ts2$sum)))[1:10],] #%>% arrange(sum)

taxa_ts3 = as.data.frame(reshape2::melt(taxa_ts2)) %>% filter(!variable == "sum")
taxa_ts3$V5 = factor(taxa_ts3$V5, levels = unique(taxa_ts3$V5))

fig5d3 = ggplot(data = taxa_ts3, 
       aes(x = ifelse(test = variable == "up", yes = -value, no = value), y = V5, fill = variable)) +
  geom_col(color = "black") +
  geom_text(aes(label = value, x= ifelse(test = variable == "up", yes = -value-20, no = value+20))) +
  labs(x = NULL, y = NULL, title = "Serrated polyp \n vs. Tubular adenoma", fill = NULL) +
  scale_fill_manual(values = c("firebrick3", "steelblue3"), labels = c("Serrated polyp", "Tubular adenoma")) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_symmetric(labels = abs) +
  coord_flip()

DAb_supplemental_figure2 = ggarrange(fig5d1, fig5d2, fig5d3, nrow = 1, labels = c("d.", "e.", "f."))
DAb_supplemental_figure = ggarrange(DAb_supplemental_figure1, DAb_supplemental_figure2, ncol = 1)

ggsave(filename = "DAB_suppl_plot.png", plot = DAb_supplemental_figure, device = "png", dpi = 600, width = 10, height = 8)

#Supplemental table
suppl_table = bind_rows(FC_healthy_TA, FC_healthy_serrated, FC_TA_serrated) %>% filter(final_pval < 0.05) %>% 
  select(taxa, group1, group2, final_pval, Log2FoldChange) %>% 
  merge(., anots[,c(4,5,21)], by.x = "taxa", by.y = "row.names") %>% 
  rename(ORF = taxa, unadj_pval = final_pval, Taxa = V5, Gene = V6, Function = V22)
write.table(suppl_table, file = "Differentially_abundant_genes.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#### Humann3 ####
pathway_table = read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/RPK_pathabundance.tsv", check.names = F, row.names = 1)
pathway_table = as.data.frame(t(pathway_table[!grepl("\\|", rownames(pathway_table)),] %>% 
                         rename_with(~ gsub('.merged_Abundance', '', .x)))) %>% select(!UNMAPPED:UNINTEGRATED) %>% filter(!rowSums(.) == 0)

pathway_merged = merge(metadata, pathway_table, by = "row.names") %>% 
  filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE")

pathway_merged[pathway_merged == "unknown"] <- NA
row.names(pathway_merged) = pathway_merged$Row.names

#Permanova with all sample types included
perma2_df = pathway_merged[complete.cases(pathway_merged),]

adonis2(formula = perma2_df[,16:ncol(perma2_df)] ~ as.numeric(perma2_df$BMI) + as.numeric(perma2_df$Age) + perma2_df$Ethnicity 
        + perma2_df$Gender + perma2_df$MostMalignantPolypType / as.character(perma2_df$Patient) / perma2_df$SampleType, 
        data = perma2_df, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df$Plate))

#Permanova with aspirates only
perma2_df_asps = perma2_df %>% filter(SampleType == "aspirate")

adonis2(formula = perma2_df_asps[,16:ncol(perma2_df_asps)] ~ as.numeric(perma2_df_asps$BMI) + as.numeric(perma2_df_asps$Age) + perma2_df_asps$Ethnicity 
        + perma2_df_asps$Gender + perma2_df_asps$ColonLocation + as.character(perma2_df_asps$PrepType) + 
          perma2_df_asps$MostMalignantPolypType / as.character(perma2_df_asps$Patient), 
        data = perma2_df_asps, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_asps$Plate))

#Pairwise permanovas with aspirates only
#Healthy vs. TA
perma2_df_hta = perma2_df_asps %>% filter(!MostMalignantPolypType == "serrated")
adonis2(formula = perma2_df_hta[,16:ncol(perma2_df_hta)] ~ as.numeric(perma2_df_hta$BMI) + as.numeric(perma2_df_hta$Age) + perma2_df_hta$Ethnicity 
        + perma2_df_hta$Gender + perma2_df_hta$ColonLocation + as.character(perma2_df_hta$PrepType) + perma2_df_hta$MostMalignantPolypType / as.character(perma2_df_hta$Patient), 
        data = perma2_df_hta, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_hta$Plate))
#Polyp type R2 = 1.0%

#Healthy vs. Serrated
perma2_df_hs = perma2_df_asps %>% filter(!MostMalignantPolypType == "TA")
adonis2(formula = perma2_df_hs[,16:ncol(perma2_df_hs)] ~ as.numeric(perma2_df_hs$BMI) + as.numeric(perma2_df_hs$Age) + perma2_df_hs$Ethnicity 
        + perma2_df_hs$Gender + perma2_df_hs$ColonLocation + as.character(perma2_df_hs$PrepType) + perma2_df_hs$MostMalignantPolypType / as.character(perma2_df_hs$Patient), 
        data = perma2_df_hs, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_hs$Plate))
#Polyp type R2 = 1.8%

#TA vs. Serrated
perma2_df_tas = perma2_df_asps %>% filter(!MostMalignantPolypType == "healthy")
adonis2(formula = perma2_df_tas[,16:ncol(perma2_df_tas)] ~ as.numeric(perma2_df_tas$BMI) + as.numeric(perma2_df_tas$Age) + perma2_df_tas$Ethnicity 
        + perma2_df_tas$Gender + perma2_df_tas$ColonLocation + as.character(perma2_df_tas$PrepType) + perma2_df_tas$MostMalignantPolypType / as.character(perma2_df_tas$Patient), 
        data = perma2_df_tas, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_tas$Plate))
#Polyp type R2 = 1.1%
#Serrated vs healthy has the LARGEST effect size. TA vs healthy has the smallest.

#pcoa of pathways
pathway_bray = vegdist(pathway_table, method = "bray")
pathway_pcoa = cmdscale(pathway_bray, eig = T, k = (nrow(pathway_table)-1), add = T)
pathway_eig = eigenvals(pathway_pcoa)
pathway_var = pathway_eig/sum(pathway_eig)
pathway_var

pathway_points = as.data.frame(pathway_pcoa$points[,1:2]) %>% merge(., metadata, by = "row.names") %>% 
  dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE")

ggplot(data = pathway_points) +
  aes(x = V1, y = V2) +
  theme_bw() +
  geom_point(aes(pch = SampleType, fill = MostMalignantPolypType), size = 2, alpha = 3/4) +
  stat_ellipse(linetype = 2, aes(group = SampleType, color = SampleType), show.legend = F) +
  scale_color_manual(values=c("gold", "darkorange4", "plum2")) +
  scale_shape_manual(values = c(21,22,24), name = "Sample type", labels = c("Aspirate", "Fecal", "Lavage")) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), 
                    labels = c("Polyp free", "Serrated polyp", "Tubular adenoma")) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  labs(color = "Subject type", x = bquote("PCo1 ("~.(round(pathway_var[1]*100, digits = 1))~"%) "), y = bquote("PCo2 ("~.(round(pathway_var[2]*100, digits = 1))~"%) "))

#Differential abundance testing of pathways.
pathway_merged_avg = pathway_merged %>% group_by(as.factor(Patient)) %>% 
  summarise(across((ncol(metadata)+2):ncol(pathway_merged), mean)) %>% 
  remove_rownames() %>% merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% 
  .[!duplicated(.),]

pw_filter_list = pathway_merged_avg %>% dplyr::select(!Patient:MostMalignantPolypType) %>% as.matrix(.) %>% 
  reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(pathway_merged_avg)/3))

pathway_merged_avg = pathway_merged_avg %>% dplyr::select(Patient, MostMalignantPolypType, !pw_filter_list$Var2)

#Healthy to tubular adenoma comparison
pathway_healthy_TA = pathway_merged_avg %>% filter(!MostMalignantPolypType == "serrated")

pathway_KW_healthy_v_TA = list()

for (i in 3:ncol(pathway_healthy_TA)) {
  tmp = kruskal.test(pathway_healthy_TA[,i]~pathway_healthy_TA$MostMalignantPolypType)
  tmp$pathway = names(pathway_healthy_TA)[i]
  pathway_KW_healthy_v_TA[[length(pathway_KW_healthy_v_TA)+1]] = tmp
}

pathway_pval_healthy_v_TA = as.data.frame(tibble(pathway = map(pathway_KW_healthy_v_TA, "pathway"), pval = map(pathway_KW_healthy_v_TA, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pathway_pval_healthy_v_TA2 = as.data.frame(colSums(pathway_healthy_TA[,3:ncol(pathway_healthy_TA)])) %>%
  merge(., pathway_pval_healthy_v_TA, by.x = "row.names", by.y = "pathway") %>% column_to_rownames(var = "Row.names")
names(pathway_pval_healthy_v_TA2)[1] = "rpm_sum"
pathway_pval_healthy_v_TA2$pval = as.numeric(pathway_pval_healthy_v_TA2$pval)

pathway_IHW_healthy_v_TA = ihw(pvalues = pathway_pval_healthy_v_TA2$pval, covariates = pathway_pval_healthy_v_TA2$rpm_sum, alpha = 0.05, 
                               nbins = 5, null_proportion = T)

pathway_IHW_healthy_v_TA@df$pathway = rownames(pathway_pval_healthy_v_TA2)
pathway_IHW_healthy_v_TA@df$group1 = "healthy"
pathway_IHW_healthy_v_TA@df$group2 = "TA"

#Healthy to serrated comparison
pathway_healthy_serrated = pathway_merged_avg %>% filter(!MostMalignantPolypType == "TA")

pathway_KW_healthy_v_serrated = list()

for (i in 3:ncol(pathway_healthy_serrated)) {
  tmp = kruskal.test(pathway_healthy_serrated[,i]~pathway_healthy_serrated$MostMalignantPolypType)
  tmp$pathway = names(pathway_healthy_serrated)[i]
  pathway_KW_healthy_v_serrated[[length(pathway_KW_healthy_v_serrated)+1]] = tmp
}

pathway_pval_healthy_v_serrated = as.data.frame(tibble(pathway = map(pathway_KW_healthy_v_serrated, "pathway"), pval = map(pathway_KW_healthy_v_serrated, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pathway_pval_healthy_v_serrated2 = as.data.frame(colSums(pathway_healthy_serrated[,3:ncol(pathway_healthy_serrated)])) %>%
  merge(., pathway_pval_healthy_v_serrated, by.x = "row.names", by.y = "pathway") %>% column_to_rownames(var = "Row.names")
names(pathway_pval_healthy_v_serrated2)[1] = "rpm_sum"
pathway_pval_healthy_v_serrated2$pval = as.numeric(pathway_pval_healthy_v_serrated2$pval)

pathway_IHW_healthy_v_serrated = ihw(pvalues = pathway_pval_healthy_v_serrated2$pval, covariates = pathway_pval_healthy_v_serrated2$rpm_sum, alpha = 0.05,
                                     nbins = 5, null_proportion = T)
pathway_IHW_healthy_v_serrated@df$pathway = rownames(pathway_pval_healthy_v_serrated2)
pathway_IHW_healthy_v_serrated@df$group1 = "healthy"
pathway_IHW_healthy_v_serrated@df$group2 = "serrated"

#Serrated to tubular adenoma comparison
pathway_serrated_TA = pathway_merged_avg %>% filter(!MostMalignantPolypType == "healthy")

pathway_KW_serrated_v_TA = list()

for (i in 3:ncol(pathway_serrated_TA)) {
  tmp = kruskal.test(pathway_serrated_TA[,i]~pathway_serrated_TA$MostMalignantPolypType)
  tmp$pathway = names(pathway_serrated_TA)[i]
  pathway_KW_serrated_v_TA[[length(pathway_KW_serrated_v_TA)+1]] = tmp
}

pathway_pval_serrated_v_TA = as.data.frame(tibble(pathway = map(pathway_KW_serrated_v_TA, "pathway"), pval = map(pathway_KW_serrated_v_TA, "p.value"))) %>%
  filter(!pval %in% NaN) %>% arrange(as.numeric(pval))
pathway_pval_serrated_v_TA2 = as.data.frame(colSums(pathway_serrated_TA[,3:ncol(pathway_serrated_TA)])) %>%
  merge(., pathway_pval_serrated_v_TA, by.x = "row.names", by.y = "pathway") %>% column_to_rownames(var = "Row.names")
names(pathway_pval_serrated_v_TA2)[1] = "rpm_sum"
pathway_pval_serrated_v_TA2$pval = as.numeric(pathway_pval_serrated_v_TA2$pval)

pathway_IHW_serrated_v_TA = ihw(pvalues = pathway_pval_serrated_v_TA2$pval, covariates = pathway_pval_serrated_v_TA2$rpm_sum, alpha = 0.05,
                                nbins = 5, null_proportion = T)
pathway_IHW_serrated_v_TA@df$pathway = rownames(pathway_pval_serrated_v_TA2)
pathway_IHW_serrated_v_TA@df$group1 = "serrated"
pathway_IHW_serrated_v_TA@df$group2 = "TA"

Pathway_DAb = bind_rows(pathway_IHW_healthy_v_TA@df, pathway_IHW_healthy_v_serrated@df, pathway_IHW_serrated_v_TA@df) %>% 
  mutate(final_pval = pvalue*3) %>% filter(final_pval <= 0.05) #Zero pathways were determined to be differentially abundant, even without FDR correction.

#Heatmap
heatmap_df = merge(metadata, pathway_table, by = "row.names") %>% 
  filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE") %>%
  select(!Patient:InFinalAnalysis) %>% column_to_rownames(var = "Row.names") 

heatmap_filter = as.matrix(heatmap_df) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(heatmap_df)/3))

zscore_df = sapply(heatmap_df, function(heatmap_df) (heatmap_df-mean(heatmap_df))/sd(heatmap_df)) %>% as.data.frame(.)
zscore_df = zscore_df[colSums(!is.na(zscore_df)) > 0]
zscore_df = zscore_df[,!names(zscore_df) %in% heatmap_filter$Var2]
row.names(zscore_df) = row.names(heatmap_df)

heatmap_df_melt = reshape2::melt(as.matrix(zscore_df)) %>% merge(., metadata, by.x = "Var1", by.y = "row.names") %>% 
  group_by(SampleType, MostMalignantPolypType) %>% arrange()

clust_order = hclust(dist(zscore_df, method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt$Var1 = factor(heatmap_df_melt$Var1, levels = row.names(zscore_df)[clust_order$order])

clust_order2 = hclust(dist(t(zscore_df), method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt$Var2 = factor(heatmap_df_melt$Var2, levels = colnames(zscore_df)[clust_order2$order])
heatmap_df_melt = heatmap_df_melt %>% mutate(Var2 = gsub(".*: ", "", Var2))

heat = ggplot(data = heatmap_df_melt) +
  aes(x = Var1, y = Var2, fill = value) +
  geom_tile() +
  theme_classic() +
  facet_grid(~SampleType+MostMalignantPolypType, scales = "free", space = "free") +
  labs(x = NULL, y = NULL, fill = "Z-score", title = "All functional pathways") +
  scale_fill_gradient2(low="steelblue3", high="firebrick3") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 3, margin = margin(l = 30)))
heat

top50 = as.data.frame(sort(colSums(abs(heatmap_df)), decreasing = T)[1:50])
zscore_df2 = zscore_df[, names(zscore_df) %in% row.names(top50)]

heatmap_df_melt2 = reshape2::melt(as.matrix(zscore_df2)) %>% merge(., metadata, by.x = "Var1", by.y = "row.names") %>% 
  group_by(SampleType, MostMalignantPolypType) %>% arrange()

clust_order3 = hclust(dist(zscore_df2, method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt2$Var1 = factor(heatmap_df_melt2$Var1, levels = row.names(zscore_df2)[clust_order3$order])

clust_order4 = hclust(dist(t(zscore_df2), method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt2$Var2 = factor(heatmap_df_melt2$Var2, levels = colnames(zscore_df2)[clust_order4$order])
heatmap_df_melt2 = heatmap_df_melt2 %>% mutate(Var2 = gsub(".*: ", "", Var2))

heat2 = ggplot(data = heatmap_df_melt2) +
  aes(x = Var1, y = Var2, fill = value) +
  geom_tile() +
  theme_classic() +
  facet_grid(~SampleType+MostMalignantPolypType, scales = "free", space = "free") +
  labs(x = NULL, y = NULL, fill = "Z-score", title = "Top 50 most abundant functional pathways") +
  scale_fill_gradient2(low="steelblue3", high="firebrick3") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 8, margin = margin(l = 30)))
heat2

ggsave("Supplemental_heatmap.svg", heat, device = "svg", dpi = 1200, width = 12, height = 18)

#Random Forest
#Randomly dividing data using 2/3 split.
#Healthy
training_healthy = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "healthy") %>% sample_frac(2/3)
training_healthy$MostMalignantPolypType = as.factor(training_healthy$MostMalignantPolypType)
training_healthy$MostMalignantPolypType = droplevels(training_healthy$MostMalignantPolypType)

test_healthy = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "healthy")
test_healthy = test_healthy[!rownames(test_healthy) %in% row.names(training_healthy),]
test_healthy$MostMalignantPolypType = as.factor(test_healthy$MostMalignantPolypType)
test_healthy$MostMalignantPolypType = droplevels(test_healthy$MostMalignantPolypType)

#Tubular adenomas
training_TA = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "TA") %>% sample_frac(2/3)
training_TA$MostMalignantPolypType = as.factor(training_TA$MostMalignantPolypType)
training_TA$MostMalignantPolypType = droplevels(training_TA$MostMalignantPolypType)

test_TA = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "TA")
test_TA = test_TA[!rownames(test_TA) %in% row.names(training_TA),]
test_TA$MostMalignantPolypType = as.factor(test_TA$MostMalignantPolypType)
test_TA$MostMalignantPolypType = droplevels(test_TA$MostMalignantPolypType)

#Serrated polyps
training_serr = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "serrated") %>% sample_frac(2/3)
training_serr$MostMalignantPolypType = as.factor(training_serr$MostMalignantPolypType)
training_serr$MostMalignantPolypType = droplevels(training_serr$MostMalignantPolypType)

test_serr = pathway_merged %>% filter(SampleType == "aspirate") %>% dplyr::select(MostMalignantPolypType, !heatmap_filter$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "serrated")
test_serr = test_serr[!rownames(test_serr) %in% row.names(training_serr),]
test_serr$MostMalignantPolypType = as.factor(test_serr$MostMalignantPolypType)
test_serr$MostMalignantPolypType = droplevels(test_serr$MostMalignantPolypType)

#Random forest training and testing
#Healthy to TA
tmp1 = bind_rows(training_healthy, training_TA) %>% janitor::clean_names()
tmp2 = rbind(test_healthy, test_TA) %>% janitor::clean_names()
model_healthy_TA = rfPermute(formula = most_malignant_polyp_type ~ ., data = tmp1, 
                             proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_healthy_TA = as.data.frame(predict(model_healthy_TA, newdata = tmp2, type = "prob"))
roc_healthy_TA = roc(tmp2[,1], predict_healthy_TA[,2], ci=TRUE, ci.alpha=0.9, 
                     stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

#Healthy to serrated polyps
tmp3 = bind_rows(training_healthy, training_serr) %>% janitor::clean_names()
tmp4 = rbind(test_healthy, test_serr) %>% janitor::clean_names()
model_healthy_serr = rfPermute(formula = most_malignant_polyp_type ~ ., data = tmp3, 
                               proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_healthy_serr = as.data.frame(predict(model_healthy_serr, newdata = tmp4, type = "prob"))
roc_healthy_serr = roc(tmp4[,1], predict_healthy_serr[,2], ci=TRUE, ci.alpha=0.9, 
                       stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

#Serrated to TA
tmp5 = bind_rows(training_TA, training_serr) %>% janitor::clean_names()
tmp6 = rbind(test_TA, test_serr) %>% janitor::clean_names()
model_serr_TA = rfPermute(formula = most_malignant_polyp_type ~ ., data = tmp5, 
                          proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_serr_TA = as.data.frame(predict(model_serr_TA, newdata = tmp6, type = "prob"))
roc_serr_TA = roc(tmp6[,1], predict_serr_TA[,2], ci=TRUE, ci.alpha=0.9, 
                  stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

svglite("pathway_ROC_curve.svg", width = 4, height = 4)
roc_curve = NULL
roc_curve = plot(roc_healthy_TA, print.auc = T, col = "orange3", print.auc.y = .3, print.auc.x = .8, main = "ROC curve")
roc_curve = plot(roc_healthy_serr, print.auc = T, col = "cyan4", print.auc.y = .2, print.auc.x = .8, add = TRUE)
roc_curve = plot(roc_serr_TA, print.auc = T, col = "purple2", print.auc.y = .1, print.auc.x = .8, add = TRUE)
dev.off()

VIP1 = as.data.frame(importance(model_healthy_TA)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP1$rowname = stringr::str_replace_all(VIP1$rowname, "_", " ") %>% stringr::str_wrap(width = 40)
VIP1$rowname = factor(VIP1$rowname, unique(VIP1$rowname))

VIP1_plot = ggplot(data = VIP1) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "orange3", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP1$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Tubular adenoma") +
  theme_bw()

VIP2 = as.data.frame(importance(model_healthy_serr)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP2$rowname = stringr::str_replace_all(VIP2$rowname, "_", " ") %>% stringr::str_wrap(width = 40)
VIP2$rowname = factor(VIP2$rowname, unique(VIP2$rowname))

VIP2_plot = ggplot(data = VIP2) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "cyan4", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP2$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Serrated polyp") +
  theme_bw()

VIP3 = as.data.frame(importance(model_serr_TA)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP3$rowname = stringr::str_replace_all(VIP3$rowname, "_", " ") %>% stringr::str_wrap(width = 40)
VIP3$rowname = factor(VIP3$rowname, unique(VIP3$rowname))

VIP3_plot = ggplot(data = VIP3) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "purple2", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy+.5, y = rowname), hjust = 0, label = round(VIP3$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Tubular adenoma \n vs. Serrated polyp") +
  theme_bw()

pathway_VIP_plot = plot_grid(VIP1_plot, VIP2_plot, VIP3_plot, nrow = 1, labels = c("a.", "b.", "c."))

#Supplemental figure for random forest
VIPs_to_plot = pathway_merged %>% filter(SampleType == "aspirate") %>% janitor::clean_names() %>%
  dplyr::select(VIP1$rowname, VIP2$rowname, VIP3$rowname) %>% rownames_to_column()
VIPs_to_plot2 = reshape2::melt(VIPs_to_plot) %>% merge(., metadata, by.x = "rowname", by.y = "row.names")

VIPs_to_plot2$variable = stringr::str_replace_all(VIPs_to_plot2$variable, "_", " ") %>% stringr::str_wrap(width = 40)

VIP_abundances = ggplot(data = VIPs_to_plot2) +
  aes(x = variable, y = value, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1), size = 0.075, alpha = 1/3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Reads per million", fill = "Subject type") +
  scale_fill_manual(labels = c("Polyp free", "Serrated", "Tubular adenoma"), 
                    values = c("forestgreen", "steelblue3", "firebrick3"))

pathway_VIP_plot2 = plot_grid(pathway_VIP_plot, VIP_abundances, ncol = 1, labels = c("", "d."))

ggsave("pathway_VIP_abundances.png", device = "png", dpi = 600, width = 14, height = 10, units = "in")

#E lenta abundance plots
elenta = eggnog_table %>% filter(str_detect(eggnog_table$V2, "Elen"), str_detect(eggnog_table$V5, "Corio")) %>% 
  column_to_rownames(var = "Row.names") %>% select(!V2:V22) %>% merge(CAZyme_annotations, ., by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% select(!`EC#`:`#ofTools`) %>% t(.) %>% 
  merge(metadata, ., by = "row.names") %>% filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), SampleType == "aspirate") %>%
  select(!Patient:InFinalAnalysis)

elenta_melt = reshape2::melt(elenta) %>% merge(., CAZyme_annotations, by.x = "variable", by.y = "row.names") %>% 
  group_by(Row.names, DIAMOND) %>% summarise(RPKG = sum(value)) %>% merge(., metadata, by.x = "Row.names", by.y = "row.names")

fig5d = ggplot(data = elenta_melt) +
  aes(x = DIAMOND, y = RPKG + 0.001, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1), size = .1) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("forestgreen", "steelblue", "firebrick3"), 
                    labels = c("Polyp free (65)", "Serrated polyp (45)", "Tubular adenoma (47)")) +
  labs(x = NULL, y = "RPKG", fill = "Subject type", title = expression(italic("E. lenta") ~ "CAZymes"), subtitle = "Mucosal aspirates only")

fig5b_d = plot_grid(fig5b, fig5d, ncol = 1, rel_heights = c(1,1), labels = c("b.", "d."))
fig5bcd = plot_grid(fig5b_d, plot5c, nrow = 1, rel_widths = c(1,1.25), labels = c("", "c."))

fig5 = plot_grid(heat2, fig5bcd, ncol = 1, rel_heights = c(1.5, 1), labels = c("a.", ""))
fig5

ggsave(filename = "Figure_5.svg", plot = fig5, device = "svg", dpi = 600, width = 14, height = 12, units = "in")
