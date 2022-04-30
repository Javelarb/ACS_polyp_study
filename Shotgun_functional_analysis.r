#library(DESeq2)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(rstatix)
library(vegan)
library(IHW)
library(ggpubr)
library(lemon)

setwd("/media/julio/Storage/CRC/github/")
register(MulticoreParam(4)) #Use four cores

metadata <- read.delim("Shotgun_metadata.tsv", row.names=1, comment.char="#", check.names = F)
contig_table <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/contig_table.txt", row.names=1, check.names = F)
anots <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/eggnog.emapper.annotations", comment.char = "#", header=FALSE, row.names = 1)

metadata$Patient <- as.factor(metadata$Patient)
contig_table <- contig_table[rowSums(contig_table) >= 10,] #Get rid of genes with low read counts

#gene lengths
gene_lengths = read.delim("gene_lengths.txt", header=FALSE, comment.char="#") %>% .[,1:2]
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

#PCoA
#http://r-sig-ecology.471788.n2.nabble.com/Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html

eggnog_bray = vegdist(t(eggnog_table[,23:ncol(eggnog_table)]), method = "bray")
eggnog_pcoa = cmdscale(eggnog_bray, eig = T, k = nrow(t(eggnog_table[,23:ncol(eggnog_table)]))-1, add = T)
eggnog_eig = eigenvals(eggnog_pcoa)
pcoa_var = eggnog_eig/sum(eggnog_eig)
pcoa_var

eggnog_points = as.data.frame(eggnog_pcoa$points[,1:2]) %>% merge(., metadata, by = "row.names") %>% 
  dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE")

plot5b = ggplot(data = eggnog_points) +
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
plot5b

#Permanova
eggnog_merged = merge(metadata, as.matrix(eggnog_bray), by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))

adonis(formula = eggnog_merged[,14:ncol(eggnog_merged)] ~ as.numeric(eggnog_merged$BMI) + 
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
#Filter OTUs that are not present across 20% of individuals.
#RPKG2 = RPKG2[1:50000,] #Delete this when you want to run the full dataset (takes forever)

RPKG2_merged = as.matrix(t(RPKG2)) %>% merge(metadata, ., by = "row.names") 
Asps_only = RPKG2_merged %>% filter(SampleType == "aspirate") %>% column_to_rownames(var = "Row.names") %>% select(!Patient:InFinalAnalysis)

variance_filter = as.matrix(Asps_only) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(Asps_only)*0.2))

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

IHW_healthy_v_TA = ihw(pvalues = pvals_KW_healthy_v_TA2$pval, covariates = pvals_KW_healthy_v_TA2$rpkg_sum, alpha = 0.05)
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

IHW_healthy_v_ser = ihw(pvalues = pvals_KW_healthy_v_ser2$pval, covariates = pvals_KW_healthy_v_ser2$rpkg_sum, alpha = 0.05)
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

IHW_TA_v_ser = ihw(pvalues = pvals_KW_TA_v_ser2$pval, covariates = pvals_KW_TA_v_ser2$rpkg_sum, alpha = 0.05)
IHW_TA_v_ser@df$taxa = rownames(pvals_KW_TA_v_ser2)
IHW_TA_v_ser@df$group1 = "TA"
IHW_TA_v_ser@df$group2 = "serrated"

#Visualizing
Polyp_type_DAb = bind_rows(IHW_healthy_v_TA@df, IHW_healthy_v_ser@df, IHW_TA_v_ser@df) %>% 
  mutate(final_pval = weighted_pvalue*3) #%>% filter(final_pval <= 0.05)

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
volc3

plot5c = ggarrange(volc1, volc2, volc3, nrow = 1, labels = c("C.", "D.", "E."))

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

# Healthy v serrateed polyp
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

plot5d = ggarrange(fig5d1, fig5d2, fig5d3, nrow = 1, labels = c("F.", "G.", "H."))

#### Humann3 ####

pathway_table = read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/RPK_pathabundance.tsv", check.names = F, row.names = 1)
pathway_table = as.data.frame(t(pathway_table[!grepl("\\|", rownames(pathway_table)),] %>% 
                         rename_with(~ gsub('.merged_Abundance', '', .x))))

pathway_merged = merge(metadata, pathway_table, by = "row.names") %>% 
  filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE", SampleType == "aspirate")

pathway_merged = pathway_merged %>% group_by(as.factor(Patient)) %>% 
  summarise(across((ncol(metadata)+2):ncol(pathway_merged), mean)) %>% 
  remove_rownames() %>% merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% 
  .[!duplicated(.),]

#Healthy to tubular adenoma comparison
pathway_healthy_TA = pathway_merged %>% filter(!MostMalignantPolypType == "serrated")

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

pathway_IHW_healthy_v_TA = ihw(pvalues = pathway_pval_healthy_v_TA2$pval, covariates = pathway_pval_healthy_v_TA2$rpm_sum, alpha = 0.05)
pathway_IHW_healthy_v_TA@df$pathway = rownames(pathway_pval_healthy_v_TA2)
pathway_IHW_healthy_v_TA@df$group1 = "healthy"
pathway_IHW_healthy_v_TA@df$group2 = "TA"

#Healthy to serrated comparison
pathway_healthy_serrated = pathway_merged %>% filter(!MostMalignantPolypType == "TA")

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

pathway_IHW_healthy_v_serrated = ihw(pvalues = pathway_pval_healthy_v_serrated2$pval, covariates = pathway_pval_healthy_v_serrated2$rpm_sum, alpha = 0.05)
pathway_IHW_healthy_v_serrated@df$pathway = rownames(pathway_pval_healthy_v_serrated2)
pathway_IHW_healthy_v_serrated@df$group1 = "healthy"
pathway_IHW_healthy_v_serrated@df$group2 = "serrated"

#Serrated to tubular adenoma comparison
pathway_serrated_TA = pathway_merged %>% filter(!MostMalignantPolypType == "healthy")

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

pathway_IHW_serrated_v_TA = ihw(pvalues = pathway_pval_serrated_v_TA2$pval, covariates = pathway_pval_serrated_v_TA2$rpm_sum, alpha = 0.05)
pathway_IHW_serrated_v_TA@df$pathway = rownames(pathway_pval_serrated_v_TA2)
pathway_IHW_serrated_v_TA@df$group1 = "serrated"
pathway_IHW_serrated_v_TA@df$group2 = "TA"

#Notes: Only three pathways were differentially abundant among the TA vs healthy comparison. The rest were not significantly different
#PWY-6953: dTDP-3-acetamido-&alpha;-D-fucose biosynthesis	- present in 6 individuals
#PWY-7688: dTDP-&alpha;-D-ravidosamine and dTDP-4-acetyl-&alpha;-D-ravidosamine biosynthesis	- 9 individuals
#PWY-6992: 1,5-anhydrofructose degradation - present in most individuals

fig5a = ggplot(data = pathway_merged) +
  aes(x = MostMalignantPolypType, y = pathway_merged$`PWY-6992: 1,5-anhydrofructose degradation` + 1, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  geom_jitter(width = 1/4) +
  labs(x = NULL, y = "Reads per million", 
       title = "1,5-anhydrofructose \n degradation") +
  scale_x_discrete(label = c("Polyp free (37)", "SP (29)", "TA (31)")) +
  scale_fill_manual(values = c("forestgreen", "steelblue3", "firebrick3")) +
  scale_y_log10(limits = c(1,300)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

plot5ab = ggarrange(fig5a, plot5b, widths = c(1,3), labels = c("A.", "B."))
plot5ab

fig5 = ggarrange(plot5ab, plot5c, plot5d, ncol = 1)
fig5

ggsave("Figure_5.svg", fig5, device = "svg", dpi = 300, height = 12, width = 10)

#Supplement
suppl_table = bind_rows(FC_healthy_TA, FC_healthy_serrated, FC_TA_serrated) %>% filter(final_pval < 0.05) %>% 
  select(taxa, group1, group2, final_pval, Log2FoldChange) %>% 
  merge(., anots[,c(4,5,21)], by.x = "taxa", by.y = "row.names") %>% 
  rename(ORF = taxa, p_adj = final_pval, Taxa = V5, Gene = V6, Function = V22)
write.table(suppl_table, file = "Differentially_abundant_genes.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#Heatmap
heatmap_df = merge(metadata, pathway_table, by = "row.names") %>% 
  filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"), !Row.names == "OLD-SAMPLE") %>%
  select(!Patient:UNINTEGRATED) %>% column_to_rownames(var = "Row.names")

heatmap_df_melt = reshape2::melt(as.matrix(heatmap_df))

clust_order = hclust(dist(heatmap_df, method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt$Var1 = factor(heatmap_df_melt$Var1, levels = rownames(heatmap_df)[clust_order$order])

clust_order2 = hclust(dist(t(heatmap_df), method = "euclidean")) #Create hierarchical clustering
heatmap_df_melt$Var2 = factor(heatmap_df_melt$Var2, levels = colnames(heatmap_df)[clust_order2$order])

heat = ggplot(data = heatmap_df_melt) +
  geom_tile(aes(x = Var1, y = Var2, fill = log10(value+0.001))) +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "Log(Reads per million)") +
  scale_fill_gradient2(low="gold2", high="navy") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, size = 2, color = "black"), axis.text.y = element_text(size = 2, color = "black"))

ggsave("Supplemental_heatmap.png", test, device = "png", dpi = 600, width = 12, height = 14)