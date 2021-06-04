  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(rstatix)
  library(vegan)
    
  setwd("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data")
  register(MulticoreParam(4)) #Use four cores
    
  metadata <- read.delim("/media/julio/Storage/CRC/CC_docs/Shotgun_metadata_final.tsv", row.names=1, comment.char="#", check.names = F)
  contig_table <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/contig_table.txt", row.names=1, check.names = F)
  anots <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/functional/eggnog.emapper.annotations", comment.char = "#", header=FALSE, row.names = 1)
  
  metadata$Patient <- as.factor(metadata$Patient)
  contig_table <- contig_table[rowSums(contig_table) >= 10,] #Get rid of genes with low read counts
  
  #Just grab any random sample for the gene lengths to divide by
  gene_lengths = read.delim("functional/pile_up/0054_norm.txt", header=FALSE, comment.char="#") %>% .[,1:2]
  colnames(gene_lengths) = c("ORF", "gene_length")
  gene_lengths$gene_length = gene_lengths$gene_length/1000 #in Kb
  
  #From microbe census
  genome_equivs = read.table("mic_cense_genome_equivs_merged.txt")
  
  #Divide reads hitting to ORF by gene length & genome equivalents.
  RPKG = merge(gene_lengths, contig_table, by.x = "ORF", by.y = "row.names") %>% column_to_rownames(var = "ORF")
  RPKG = RPKG[,-1]/RPKG$gene_length
  RPKG2 = merge(genome_equivs, t(RPKG), by.x = "V1", by = "row.names") %>% column_to_rownames(var = "V1")
  RPKG2 = as.data.frame(t(RPKG2[,-1]/RPKG2$V2))
  
  eggnog_table = merge(anots, RPKG2, by = "row.names")
  
  #### Bar plot ####
  brite <- read.delim("/media/julio/Storage/DBs/brite.txt", header=FALSE)
  
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  brite_list <- tidyr::separate(anots, col = V14, into = c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10", 
                                                           "L11","L12","L13", "L14","L15","L16","L17","L18","L19",
                                                           "L20","L21","L22"), sep = "\\,", remove = T, extra = "drop")  %>% .[,13:(13+21)]
  brite_list <- as.data.frame(apply(brite_list, 1, lastValue))
  brite_list <- subset(brite_list, !(brite_list == "")) 
  
  colnames(brite_list) <- "V1"
  brite_list$Row.names <- rownames(brite_list)
  brite_list <- merge(brite_list, brite, by = "V1")
  
  contig_relab = t(RPKG2)/rowSums(t(RPKG2))
  contig_relab <- contig_relab %>% reshape2::melt()
  contig_relab <- contig_relab[!(contig_relab$value == 0),]
  
  plot_df <- subset(contig_relab, contig_relab$Var2 %in% brite_list$Row.names)
  plot_df2 <- merge(plot_df,brite_list[,c("Row.names","V2")], by.x = "Var2", by.y = "Row.names") %>% merge(., metadata, by.x = "Var1", by.y = "row.names")
  plot_df2$V2 = as.character(plot_df2$V2)
  
  #Take top 10 genes.
  top_genes <- group_by(plot_df2, V2) %>% summarise(., top_genes_tmp = sum(value)) %>% arrange(., desc(top_genes_tmp)) %>% slice(., 1:10)
  high_abundance <- split(top_genes$V2, 1:NROW(top_genes))
  
  #Change non top hits to other.
  plot_df2$V2[plot_df2$V2 %in% high_abundance != "TRUE"] <- "Other"
  plot_df2 <- plot_df2[order(plot_df2$V2),] #Re order
  plot_df2 <- rbind(plot_df2[!(plot_df2$V2 == "Other"),],plot_df2[(plot_df2$V2 == "Other"),]) #Move other to bottom
  plot_df2$V2 <- factor(plot_df2$V2, levels = unique(plot_df2$V2)) #Fix the order
  
  #IF you remove enzyme then you need this next line.
  relab <- aggregate(plot_df2$value, by=list(Var1=plot_df2$Var1), FUN=sum)
  plot_df2 <- merge(plot_df2, relab, by = "Var1")
  plot_df2$relab <- plot_df2$value/plot_df2$x
  plot_df2 = plot_df2 %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))
  
  sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947", "gray")
  
  brite_plot = ggplot(data = plot_df2, aes(x = Var1, weight = relab, fill = V2)) +
    geom_bar(width = 1, color = "black", size = .2) +
    theme_classic(base_size = 16) +
    facet_wrap(SampleType~MostMalignantPolypType, scales = "free", labeller = 
                 labeller(SampleType = c(`aspirate` = "Aspirate", `fecal` = "Fecal", `lavage` = "Lavage"), 
                          MostMalignantPolypType = c(`healthy`="Healthy", `TA`="TA", `serrated`="Serrated"))) +
    scale_fill_manual(values = sarah_color) +
    theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 10), 
          strip.background = element_rect(fill="lightblue")) +
    labs(x = NULL,
         y = "Relative abundance", fill = "KEGG Brite pathway", title = "Shotgun - Individuals: 104")
  brite_plot
  
  ggsave("Brite_plot.png", plot = brite_plot, device = "png", units = "in", dpi = 300, height = 10, width = 14)
  
  #PCoA
  #http://r-sig-ecology.471788.n2.nabble.com/Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html
  
  eggnog_bray = vegdist(t(eggnog_table[,23:ncol(eggnog_table)]), method = "bray")
  eggnog_pcoa = cmdscale(eggnog_bray, eig = T, k = nrow(t(eggnog_table[,23:ncol(eggnog_table)]))-1, add = T)
  eggnog_eig = eigenvals(eggnog_pcoa)
  pcoa_var = eggnog_eig/sum(eggnog_eig)
  pcoa_var
  
  eggnog_points = as.data.frame(eggnog_pcoa$points[,1:2]) %>% merge(., metadata, by = "row.names") %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))
  
  plot5a = ggplot(data = eggnog_points) +
    aes(x = V1, y = V2, color = SampleType) +
    theme_bw() +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2, aes(group = SampleType), show.legend = F, color = "black") +
    #annotate("text", x = .55, y = -.5, size = 4, label = bquote("Stress ="~.(round(eggnog_nmds$stress, digits = 2)))) +
    scale_color_manual(values=c("gold", "darkorange4", "plum2"), labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
    labs(color = "Sample type", x = bquote("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) "), y = bquote("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) "))
  plot5a
  
  #Permanova
  eggnog_merged = merge(metadata, as.matrix(eggnog_bray), by = "row.names") %>% 
    column_to_rownames(var = "Row.names") %>% dplyr::filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))
  
  adonis(formula = eggnog_merged[,14:ncol(eggnog_merged)] ~ as.numeric(eggnog_merged$BMI) + 
           as.numeric(eggnog_merged$Age) + eggnog_merged$Ethnicity + eggnog_merged$Gender + 
           (eggnog_merged$MostMalignantPolypType /as.character(eggnog_merged$Patient) / eggnog_merged$SampleType), 
         data = eggnog_merged, method = "bray", permutations = 999, parallel = 32, strata = as.factor(eggnog_merged$Plate))
  
  #### Differential abundance ####
  merged_cont_table <- merge(metadata, t(contig_table), by = "row.names")
  merged_cont_table <- merged_cont_table[!(merged_cont_table$MostMalignantPolypType %in% c("adenocarcinoma", "unknown")),]
  merged_cont_table <- merged_cont_table[!(row.names(merged_cont_table) %in% c("0837")),]
  
  rownames(merged_cont_table) <- merged_cont_table$Row.names
  
  l_asps <- merged_cont_table[(merged_cont_table$SampleType == "aspirate"),] %>% .[(.[,6] == "left"),]
  r_asps <- merged_cont_table[(merged_cont_table$SampleType == "aspirate"),] %>% .[(.[,6] == "right"),]
  fecal <- merged_cont_table[(merged_cont_table$SampleType == "fecal"),]
  #sampletype <- rbind(r_asps[!(r_asps$Patient %in% fecal$Patient),], fecal) 
  
  l_asps_deseq <- DESeqDataSetFromMatrix(countData = t(l_asps[,(ncol(metadata)+2):ncol(l_asps)]), 
                                colData = l_asps[,2:(ncol(metadata)+1)],
                                design = ~ MostMalignantPolypType)
  
  r_asps_deseq <- DESeqDataSetFromMatrix(countData = t(r_asps[,(ncol(metadata)+2):ncol(r_asps)]), 
                                         colData = r_asps[,2:(ncol(metadata)+1)],
                                         design = ~ MostMalignantPolypType)
  
  fecal_deseq <- DESeqDataSetFromMatrix(countData = t(fecal[,(ncol(metadata)+2):ncol(fecal)]), 
                                         colData = fecal[,2:(ncol(metadata)+1)],
                                         design = ~ MostMalignantPolypType)
  
  #sample_deseq <- DESeqDataSetFromMatrix(countData = t(sampletype[,(ncol(metadata)+2):ncol(sampletype)]), 
  #                                       colData = sampletype[,2:(ncol(metadata)+1)],
  #                                       design = ~ SampleType)
    
  l_asps_deseq_analysis <- DESeq(l_asps_deseq, parallel = T)
  r_asps_deseq_analysis <- DESeq(r_asps_deseq, parallel = T)
  fecal_deseq_analysis <- DESeq(fecal_deseq, parallel = T)
  #sample_deseq_analysis <- DESeq(fecal_deseq, parallel = T)
  
  l_asps_results1 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
  l_asps_results2 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
  l_asps_results3 <- results(l_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)
  
  r_asps_results1 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
  r_asps_results2 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
  r_asps_results3 <- results(r_asps_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)
  
  fecal_results1 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "healthy")) %>% subset(., padj < 0.05)
  fecal_results2 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "serrated", "healthy")) %>% subset(., padj < 0.05)
  fecal_results3 <- results(fecal_deseq_analysis, alpha = 0.05, contrast = c("MostMalignantPolypType", "TA", "serrated")) %>% subset(., padj < 0.05)
  
  #sample_deseq_results <- results(sample_deseq_analysis, alpha = 0.05) %>% subset(., padj < 0.05)
  
  summary(l_asps_results1)
  summary(l_asps_results2)
  summary(l_asps_results3)
  
  summary(r_asps_results1)
  summary(r_asps_results2)
  summary(r_asps_results3)
  
  summary(fecal_results1)
  summary(fecal_results2)
  summary(fecal_results3)
  
  summary(sample_deseq_results)
  
  l_asps_summary1 <- bind_rows(l_asps_results1@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results1@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  l_asps_summary2 <- bind_rows(l_asps_results2@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results2@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  l_asps_summary3 <- bind_rows(l_asps_results3@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(l_asps_results3@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  
  r_asps_summary1 <- bind_rows(r_asps_results1@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results1@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  r_asps_summary2 <- bind_rows(r_asps_results2@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results2@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  r_asps_summary3 <- bind_rows(r_asps_results3@listData) %>% dplyr::select("log2FoldChange", "padj") %>% magrittr::set_rownames(r_asps_results3@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("aspirate") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  
  fecal_summary1 <- bind_rows(fecal_results1@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results1@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to healthy") %>% tibble::add_column("fecal") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  fecal_summary2 <- bind_rows(fecal_results2@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results2@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("serrated to healthy") %>% tibble::add_column("fecal") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  fecal_summary3 <- bind_rows(fecal_results3@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(fecal_results3@rownames) %>% 
    merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("TA to serrated") %>% tibble::add_column("fecal") %>%
    dplyr::rename(comparison = 7, SampleType = 8)
  
  #sample_summary <- bind_rows(sample_deseq_results@listData) %>% select("log2FoldChange", "padj") %>% magrittr::set_rownames(sample_deseq_results@rownames) %>% 
  #  merge(., anots %>% select(4, 5, 21), by = "row.names") %>% tibble::add_column("aspirate to fecal") %>% rename(comparison = 6)
  
  master_summary <- bind_rows(l_asps_summary1, l_asps_summary2, l_asps_summary3, r_asps_summary1, r_asps_summary2, r_asps_summary3)
  master_summary <- master_summary[!master_summary$V5 == "Planctomycetes",]
  
  #### diff. ab. plots ####
  master_summary2 <- subset(master_summary, master_summary$comparison == "serrated to healthy")
  master_summary2$color <- master_summary2$log2FoldChange > 0
  labels1 <- master_summary2 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]
  
  ggplot(data = master_summary2) +
    theme_classic(base_size = 14, base_line_size = 1) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), size = 3, pch = 21, alpha = .5) +
    geom_text_repel(data = master_summary2, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
    geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(fill ="Subject type") +
    scale_fill_manual(values = c("forestgreen", "firebrick3"), labels = c("Healthy", "Tubular adenoma"))
  plot5c1
  
  master_summary3 <- subset(master_summary, master_summary$comparison == "serrated to healthy")
  master_summary3$color <- master_summary3$log2FoldChange > 0
  labels2 <- master_summary3 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]
  
  ggplot(data = master_summary3) +
    theme_classic(base_size = 14, base_line_size = 1) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), pch = 21, size = 3, alpha = .5) +
    geom_text_repel(data = master_summary3, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
    geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(fill ="Subject type") +
    scale_fill_manual(values = c("forestgreen", "steelblue3"), labels =  c("Healthy", "Serrated"))
  #  theme(legend.position = "none")
  
  plot_grid(plot5c1, plot5c2, rel_widths = c(.65,1))
  #ggsave("figure_5c.svg", plot = plot5c3, device = "svg", units = "in", dpi = 1000, height = 3.5, width = 5)
  
  master_summary4 <- subset(master_summary, master_summary$comparison == "TA to serrated")
  master_summary4$color <- master_summary4$log2FoldChange > 0
  labels3 <- master_summary4 %>% arrange(desc(abs(-log10(padj)))) %>% .[1:6,]
  
  ggplot(data = master_summary4) +
    theme_classic(base_size = 14, base_line_size = 1) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = color), pch = 21, size = 3, alpha = .5) +
    geom_text_repel(data = master_summary4, aes(x = log2FoldChange, y = -log10(padj), label = V6), size = 4, vjust = -1.5, box.padding = 1) +
    geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(fill ="Subject type") +
    scale_fill_manual(values = c("steelblue3", "firebrick3"), labels = c("Serrated", "Tubular adenoma"))
  
  #Checked gene abundances to eliminate false positives.
  #Should be present in at least 10% of subjects.
  #k111_458716_2 (tonB), k111_2791185_4 (GH31), k111_1237370_2 (sdaAA)
  
  for (i in unique(master_summary$Row.names)) {
    assign(paste0("L_", i), plotCounts(l_asps_deseq_analysis, gene = i, intgroup = "MostMalignantPolypType", pc = 0, returnData = T))
    assign(paste0("R_", i), plotCounts(r_asps_deseq_analysis, gene = i, intgroup = "MostMalignantPolypType", pc = 0, returnData = T))
    assign(paste0("Both_", i), rbind(get(paste0("L_", i)), get(paste0("R_", i))) %>% merge(., metadata, by = "row.names"))
    
    assign(paste0("plot_", i), ggplot(data = get(paste0("Both_", i))) +
      aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha = .5) +
      labs(title = i, x = NULL, y = "DESeq2 normalized counts") +
      scale_x_discrete(labels = c("Healthy (63)", "Serrated (45)", "TA (49)")) +
      scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
      theme_classic() +
      theme(legend.position = "none"))
  }
  
  tonB <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_458716_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
                plotCounts(r_asps_deseq_analysis, gene = "k111_458716_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
                merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
  tonB <- tonB[!duplicated(tonB),]
  dunn_test(data = tonB, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")
  
  plot5d <- ggplot(data = Both_k111_458716_2) +
    aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
    geom_boxplot(outlier.shape = NA, lwd = 1.5) +
    geom_jitter(width = .25) +
    labs(title = "", subtitle = expression(italic("tonB")), x = NULL, y = NULL) +
    scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
    scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_log10() +
    theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))
  
  GH31 <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_2791185_4", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
                plotCounts(r_asps_deseq_analysis, gene = "k111_2791185_4", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
                merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
  GH31 <- GH31[!duplicated(GH31),]
  dunn_test(data = GH31, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")
  
  plot5c <- ggplot(data = Both_k111_2791185_4) +
    aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
    geom_boxplot(outlier.shape = NA, lwd = 1.5) +
    geom_jitter(width = .25) +
    labs(title = "", subtitle = "GH31", x = NULL, y = NULL) +
    scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
    scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_log10() +
    theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))
  
  sdaAA <- rbind(plotCounts(l_asps_deseq_analysis, gene = "k111_1237370_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T),
                plotCounts(r_asps_deseq_analysis, gene = "k111_1237370_2", intgroup = "MostMalignantPolypType", pc = 0, returnData = T)) %>% 
    merge(., metadata, by = "row.names") %>% group_by(Patient) %>% summarise(avg = mean(count)) %>% merge(., metadata[,c(1,10)], by = "Patient") 
  sdaAA <- sdaAA[!duplicated(sdaAA),]
  dunn_test(data = sdaAA, formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")
  
  plot5b <- ggplot(data = Both_k111_1237370_2) +
      aes(y = count, x = MostMalignantPolypType.x, fill = MostMalignantPolypType.x) +
      geom_boxplot(outlier.shape = NA, lwd = 1.5) +
      geom_jitter(width = .25) +
      labs(title = "Mucosal aspirates only", subtitle = expression(italic("sdaA")), x = NULL, y = "DESeq2 normalized reads") +
      scale_x_discrete(labels = c("Healthy (65)", "Serrated (45)", "TA (47)")) +
      scale_fill_manual(values=c("forestgreen", "steelblue3", "firebrick3")) +
      theme_classic() +
      theme(legend.position = "none") +
      scale_y_log10() +
      theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(colour = "black"))
  
  plot5b_d <- plot_grid(plot5b,plot5c,plot5d, nrow = 1, rel_widths = c(1,.95,.95))
  plot5b_d
  
  #ggsave(filename = "fig5.svg", plot = plot_grid(plot5a, plot5b_d, ncol = 1, rel_heights = c(.6,.4)), device = "svg", dpi = 300, height = 8, width = 8.5)
