library(tidyverse)
library(vegan)
library(nlme)
library(rfPermute)
library(rgl)
library(ggrepel)
library(cluster)
library(car)
library(MASS)
library(cowplot)
library(ape)
library(rstatix)
library(IHW)
library(pROC)
library(svglite)

source("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/ANCOM_V2.1.R")
source("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/CustomRFPermute.R")

options(scipen=10000)
set.seed(seed = 999)

setwd("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data")
metadata <- read.delim("/media/julio/Storage/CRC/CC_docs/Shotgun_metadata_final.tsv", row.names=1, comment.char="#", check.names = F)
metadata$Patient <- as.factor(metadata$Patient)

OTU_table <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/iggsearch_taxonomy/iggsearch_2_merged/total_mapped_reads.tsv", row.names=1, comment.char="#", check.names = F)
species <- read.delim("/media/julio/Storage/Software/IGGsearch/iggdb_v1.0.0/iggdb_v1.0.0.species", row.names=1, comment.char="#", check.names = F)

#### Read counts####
read_counts_1 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/raw_read_counts.txt", row.names=1, quote="\"", comment.char="", check.names = F)
read_counts_2 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/NH_read_counts.txt", row.names=1, quote="\"", comment.char="", check.names = F)
read_counts_3 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/QF_read_counts.txt", row.names = 1, check.names = F)

read_counts_all <- merge(read_counts_1, read_counts_3, by = "row.names")
read_counts_all <- merge(read_counts_all, read_counts_2, by.x = "Row.names", by.y = "row.names")
names(read_counts_all) <- c("SampleID", "Raw", "QF", "QF+Decon")
read_counts_all <- reshape2::melt(read_counts_all) %>% merge(., metadata, by.x = "SampleID", by.y = "row.names")

read_counts_all <- read_counts_all[!(read_counts_all$MostMalignantPolypType %in% c("adenocarcinoma", "unknown")),]
#read_counts_all <- read_counts_all[!(read_counts_all$SampleType %in% "lavage"),]

ggplot(data = read_counts_all) +
  aes(x = SampleType, y = value, fill = variable) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = 'Sample Type',
       y = 'PE150 read count', title = "Sample set 2 - Shotgun reads per sample") +
  theme_bw(base_size = 14) +
  #facet_wrap(.~SampleType) +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = 0.25, size = 1)

#### Making OTU tables ####
OTU_table_plus_taxonomy <- merge(t(OTU_table), species, by = "row.names")

#Put the species name in place of OTU ID but keep unique
row.names(OTU_table_plus_taxonomy) <- make.names(OTU_table_plus_taxonomy$species_name, unique = T) 

#Make entries integers.
OTU_table_plus_taxonomy <- as.data.frame(apply(t(OTU_table_plus_taxonomy[,2:(1+(nrow(OTU_table)))]), MARGIN = c(1,2), as.integer))

#Check to see if there are any contaminants. Only 10 taxa should be there. See:
# https://files.zymoresearch.com/protocols/_d6305_d6306_zymobiomics_microbial_community_dna_standard.pdf
pos_ctrl <-OTU_table_plus_taxonomy[rownames(OTU_table_plus_taxonomy) == "COMM-STD",] %>% .[,!(colSums(.) == 0)] #No contaminants to filter in my case, but also no fungi.
OTU_table_plus_taxonomy <- OTU_table_plus_taxonomy[!(rownames(OTU_table_plus_taxonomy) %in% "COMM-STD") & rowSums(OTU_table_plus_taxonomy) > 1000,] 

merged_OTU_table <- merge(metadata, OTU_table_plus_taxonomy, by = "row.names")
merged_OTU_table <- merged_OTU_table[!(merged_OTU_table$MostMalignantPolypType %in% c("unknown","adenocarcinoma")),]
row.names(merged_OTU_table) <- merged_OTU_table$Row.names

fecal_OTU <- merged_OTU_table[grepl("fecal", merged_OTU_table$SampleType),]
asp_OTU <- subset(merged_OTU_table, merged_OTU_table$SampleType == "aspirate")
lavage_OTU <- merged_OTU_table[grepl("lavage", merged_OTU_table$SampleType),]
asp_L_OTU <- asp_OTU[grepl("left", asp_OTU$ColonLocation),]
asp_R_OTU <- asp_OTU[grepl("right", asp_OTU$ColonLocation),]

#### alpha diversity ####
alpha_div <- as.data.frame(diversity(merged_OTU_table[,15:ncol(merged_OTU_table)], index = "shannon"))
alpha_div <- merge(alpha_div, metadata, by = "row.names")
#alpha_div <- alpha_div[!(alpha_div$SampleType == "lavage"),]

fig2a1 <- ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = `diversity(merged_OTU_table[, 15:ncol(merged_OTU_table)], index = "shannon")`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("Shotgun - Individuals:"~.(length(unique(alpha_div$Patient))))) +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
fig2a1

qqp(alpha_div$`diversity(merged_OTU_table[, 15:ncol(merged_OTU_table)], index = "shannon")`, "norm") #Pretty good fit for lm
alpha_lm <- NULL
alpha_lm$x <- as.numeric(alpha_div$`diversity(merged_OTU_table[, 15:ncol(merged_OTU_table)], index = "shannon")`)
alpha_lm$y <- as.factor(alpha_div$SampleType)
alpha_lm$z <- as.factor(alpha_div$MostMalignantPolypType)
alpha_lm$p <- as.factor(alpha_div$Plate)
alpha_lm$i <- as.factor(alpha_div$Patient)
alpha_lm <- as.data.frame(alpha_lm)
alpha_lm <- within(alpha_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))

#species counts instead of shannon
specno <- as.data.frame(specnumber(merged_OTU_table[,15:ncol(merged_OTU_table)]))
specno <- merge(specno, metadata, by = "row.names")
#specno <- specno[!(specno$SampleType == "lavage"),]

fig2a2 <- ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(merged_OTU_table[, 15:ncol(merged_OTU_table)])`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.text = element_text(size = 10))
fig2a2

qqp(specno$`specnumber(merged_OTU_table[, 15:ncol(merged_OTU_table)])`, "norm")

specno_lm <- NULL
specno_lm$x <- as.numeric(specno$`specnumber(merged_OTU_table[, 15:ncol(merged_OTU_table)])`)
specno_lm$y <- as.factor(specno$SampleType)
specno_lm$z <- as.factor(specno$MostMalignantPolypType)
specno_lm$i <- as.factor(specno$Patient)
specno_lm$p <- as.factor(specno$Plate)
specno_lm <- as.data.frame(specno_lm)
specno_lm <- within(specno_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = specno_lm, random = list(p=~1, i=~1)))

fig2a <- plot_grid(fig2a1, fig2a2, rel_heights = c(1,1), ncol = 1)
fig2a
#ggsave("Figure_3a.svg", plot = fig3a, device = "svg", units = "in", dpi = 1000, height = 4, width = 7.5)

alpha_supp1 = ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = alpha_div$`diversity(merged_OTU_table[, 15:ncol(merged_OTU_table)], index = "shannon")`, 
      fill = as.character(MostMalignantPolypType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = "Sample set 2 - shotgun data") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c(" (78)", "Tubular adenoma (75)", "Serrated polyp (58)")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp2 = ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(merged_OTU_table[, 15:ncol(merged_OTU_table)])`, 
      fill = as.character(MostMalignantPolypType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type', title = "") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c(" (78)", "Tubular adenoma (75)", "Serrated polyp (58)")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp = plot_grid(alpha_supp1, alpha_supp2, rel_widths = c(1,1.65))
#ggsave("Shotgun_alpha_polyp.png", plot = alpha_supp, device = "png", units = "in", dpi = 300, height = 4, width = 10)

#### Beta div####
NMDS <- metaMDS(merged_OTU_table[,-(1:(1+ncol(metadata)))], trymax = 999, parallel = 32)
NMDS_points <- NMDS$points %>% merge(., metadata, by = "row.names") %>% filter(!MostMalignantPolypType %in% c("unknown", "adenocarcinoma"))

asp_dm <- as.data.frame(as.matrix(vegdist(asp_OTU[,15:ncol(asp_OTU)], method = "bray")))
asp_pca <- pcoa(asp_dm)
merged_asp_pca <- merge(asp_pca$vectors[,1:2], metadata, by = "row.names")

fecal_dm <- as.data.frame(as.matrix(vegdist(fecal_OTU[,15:ncol(fecal_OTU)], method = "bray")))
fecal_pca <- pcoa(fecal_dm)
merged_fecal_pca <- merge(fecal_pca$vectors[,1:2], metadata, by = "row.names")

lavage_dm <- as.data.frame(as.matrix(vegdist(lavage_OTU[,15:ncol(lavage_OTU)], method = "bray")))
lavage_pca <- pcoa(lavage_dm)
merged_lavage_pca <- merge(lavage_pca$vectors[,1:2], metadata, by = "row.names")

#samp_type_perma  
perma1 <- adonis(formula = merged_OTU_table[,15:ncol(merged_OTU_table)] ~ as.numeric(merged_OTU_table$BMI) + 
          as.numeric(merged_OTU_table$Age) + merged_OTU_table$Ethnicity + merged_OTU_table$Gender + 
          (merged_OTU_table$MostMalignantPolypType /as.character(merged_OTU_table$Patient) / merged_OTU_table$SampleType), 
          data = merged_OTU_table, method = "bray", permutations = 999, parallel = 32, strata = as.factor(merged_OTU_table$Plate))
perma1
#capture.output(perma1, file = "perma_out.txt")

coef1 <- coefficients(perma1)["merged_OTU_table$MostMalignantPolypType1",]
top.coef1 <- coef1[rev(order(abs(coef1)))[1:12]]
par(mar=c(2,12,1,1)+.1)
barplot(sort(top.coef1), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=1)

#asp perma
perma2 <- adonis(formula = asp_OTU[,15:ncol(asp_OTU)] ~ as.numeric(asp_OTU$BMI) + as.numeric(asp_OTU$Age) + asp_OTU$Ethnicity 
          + asp_OTU$Gender + asp_OTU$ColonLocation + as.character(asp_OTU$PrepType) + asp_OTU$MostMalignantPolypType / as.character(asp_OTU$Patient), 
          data = asp_OTU, method = "bray", permutations = 999, parallel = 32, strata = as.factor(asp_OTU$Plate))

coef2 <- coefficients(perma2)["asp_OTU$MostMalignantPolypType2",]
top.coef2 <- as.data.frame(coef2[rev(order(abs(coef2)))[1:10]])
top.coef2$Taxa <- rownames(top.coef2)
top.coef2 <- top.coef2 %>% mutate(Taxa = stringr::str_replace_all(top.coef2$Taxa, "\\.", " ")) %>% 
  rename(value = `coef2[rev(order(abs(coef2)))[1:10]]`) %>% arrange(., value)
top.coef2$fill[(top.coef2$value > 0)] <- "Tubular adenoma"
top.coef2$fill[(top.coef2$value < 0)] <- "Healthy"
top.coef2$Taxa <- factor(top.coef2$Taxa, levels = top.coef2$Taxa)

ggplot(data = top.coef2) +
  aes(x = value, y = Taxa, fill = fill) +
  geom_col(color = "black", size = 1) +
  labs(x = "PERMANOVA coefficient value", y = NULL, fill = NULL, title = "Aspirate PERMANOVA coefficients",
       subtitle = bquote("Individuals:"~.(length(unique(asp_OTU$Patient))))) +
  scale_fill_manual(values=c("forestgreen", "firebrick3")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), 
        title = element_text(size = 12), legend.text = element_text(size = 12, color = "black"))

#fecal perma
perma3 <- adonis(formula = fecal_OTU[,15:ncol(fecal_OTU)] ~ fecal_OTU$MostMalignantPolypType 
                 + as.numeric(fecal_OTU$Age) + as.numeric(fecal_OTU$BMI) + fecal_OTU$Gender 
                 + fecal_OTU$Ethnicity, data = fecal_OTU, method = "bray", permutations = 999, parallel = 32, 
                 strata = as.factor(fecal_OTU$Plate))

coef3 <- coefficients(perma3)["fecal_OTU$MostMalignantPolypType2",]
top.coef3 <- coef3[rev(order(abs(coef3)))[1:12]]
barplot(sort(top.coef3), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=1)

#Lavage
perma4 <- adonis(formula = lavage_OTU[,15:ncol(lavage_OTU)] ~ lavage_OTU$MostMalignantPolypType + lavage_OTU$PrepType + 
                  as.numeric(lavage_OTU$Age) + as.numeric(lavage_OTU$BMI) + lavage_OTU$Gender + 
                  lavage_OTU$Ethnicity, data = lavage_OTU, method = "bray", permutations = 999, parallel = 32,
                  strata = as.factor(lavage_OTU$Plate))

coef4 <- coefficients(perma4)["lavage_OTU$MostMalignantPolypType2",]
top.coef4 <- coef4[rev(order(abs(coef4)))[1:12]]
barplot(sort(top.coef4), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=1)

#General plot
ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = SampleType, fill = MostMalignantPolypType), size = 6, alpha = 0.6) + 
  labs(x= "MDS1", y = "MDS2", subtitle = bquote("Individuals:"~.(length(unique(NMDS_points$Patient))))) +
  scale_shape_manual(values = c(21,22,24), name = "Sample type", labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  stat_ellipse(linetype = 2, aes(group = SampleType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), 
                    labels = c(" (78)", "Serrated polyp (58)", "Tubular adenoma (75)")) +
  geom_text(label = NMDS_points$Patient, size = 3, color = "black")

fig2b <- ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_bw() +
  geom_point(aes(pch = MostMalignantPolypType, fill = SampleType), size = 5, alpha = 0.5) + 
  labs(x= "MDS1", y = "MDS2") +
  scale_shape_manual(values = c(21,22,24), name = "Subject type", labels = c("Polyp free (78)", "Serrated polyp (58)", "Tubular adenoma (75)")) +
  #stat_ellipse(linetype = 2, aes(group = SampleType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Sample type", values=c("gold", "darkorange4", "plum2"), 
                    labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  geom_text(label = NMDS_points$Patient, size = 2, color = "black") +
  annotate("text", x = 1.25, y = -1.2, label = bquote("Stress ="~.(round(NMDS$stress, digits = 2))))
fig2b

ggsave("Figure_2b.svg", plot = plot_grid(fig2a, fig2b, nrow = 1, rel_widths = c(.35,.65), labels = c("E.", "F.")), 
       device = "svg", units = "in", dpi = 300, height = 5, width = 11)

#aspirate only
ggplot(data = merged_asp_pca) +
  aes(x = Axis.1, y = Axis.2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = ColonLocation, fill = MostMalignantPolypType), size = 8) + 
  labs(x= "PCoA 1", y = "PCoA 2") +
  scale_shape_manual(values = c(23,22), name = "Colon side") +
  #ggtitle(bquote("Patient:"~R^2~"= 0.94, p = 0.001"), subtitle = ~ atop("Colon location:"~R^2~"= 0.001, p = 0.006",
   #                                                          "Subject type:" ~R^2 ~ "= 0.02, p = 0.001")) +
  stat_ellipse(linetype = 2, aes(color = MostMalignantPolypType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_asp_pca$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#fecal
ggplot(data = merged_fecal_pca) +
  aes(x = Axis.1, y = Axis.2, fill = MostMalignantPolypType) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 8) + 
  #ggtitle(bquote("Patient:"~R^2~"= 0.94, p = 0.99"), subtitle = bquote("Subject type:"~R^2~"= 0.07, p = 0.99")) +
  stat_ellipse(linetype = 2, aes(color = MostMalignantPolypType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_fecal_pca$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#Lavage
ggplot(data = merged_lavage_pca) +
  aes(x = Axis.1, y = Axis.2, fill = MostMalignantPolypType) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 8) + 
  labs(x= "PC1 (14%)", y = "PC2 (12%)") +
  #ggtitle(bquote("Patient:"~R^2~"= 0.90, p = 0.99"), subtitle = bquote("Subject type:"~R^2~"= 0.10, p = 0.99")) +
  stat_ellipse(linetype = 2, aes(color = MostMalignantPolypType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_lavage_pca$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#### ANCOM v2.1 ####
#Notes: ANCOM is necessary for sample type comparisons due to repeated measurements, despite its reduced sensitivity.

ancom_otu1 <- merged_OTU_table[!(merged_OTU_table$SampleType == "lavage"),] #Compare fecal with asps
ancom_otu1 <- as.data.frame(t(ancom_otu1[,-(1:(1+ncol(metadata)))])) # Transpose otu table.

ancom_meta1 <- metadata
ancom_meta1$Patient <- as.factor(ancom_meta1$Patient)
ancom_meta1$Row.names <- rownames(metadata)
rownames(ancom_meta1) <- NULL

ancom_pre1 <- feature_table_pre_process(feature_table = ancom_otu1, meta_data = ancom_meta1, sample_var = "Row.names",
                                       group_var = "SampleType", lib_cut = 0, neg_lb = T)

ancom_res1 <- ANCOM(feature_table = ancom_pre1$feature_table, meta_data = ancom_pre1$meta_data, struc_zero = ancom_pre1$structure_zeros, 
                   main_var = "SampleType", rand_formula = "~ 1 | Patient", p_adj_method = "fdr")

ancom_out <- ancom_res1$out %>% filter(!(W == "Inf"), detected_0.7 == TRUE) %>% arrange(desc(W))
#write.table(ancom_out, "Fecal_v_lav_DAb.txt", sep = "\t", quote = F, row.names = F)

#Plot the abundance of the top ancom results.
diff_ab_df <- merged_OTU_table[,colnames(merged_OTU_table) %in% ancom_out$taxa_id] %>% merge(metadata, ., by = "row.names")

plot_data_column = function (column) {
  ggplot(diff_ab_df) +
    aes(x = SampleType, y = diff_ab_df[,column]) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    labs(y = colnames(diff_ab_df)[column])
}

DAb_taxa <- lapply((ncol(metadata)+2):ncol(diff_ab_df), plot_data_column)

#### Taxa barplots ####
#Turn to relative abundance
relab_OTU = as.data.frame(t(OTU_table)) %>% select(!c("COMM-STD", "0725", "0728")) %>% 
  merge(., species, by = "row.names") %>% column_to_rownames(var = "Row.names")

relab_OTU[,1:(nrow(OTU_table)-3)] = decostand(x = relab_OTU[,1:(nrow(OTU_table)-3)], method = "total", MARGIN = 2)

#Melt and separate taxonomy, remove whitespaces
barplot_df <- reshape2::melt(as.matrix(relab_OTU[,1:(nrow(OTU_table)-3)])) %>% 
  merge(., species[,c(2,9)], by.x = "Var1", by.y = "row.names") %>% mutate(gtdb_taxonomy = gsub(".__", "", .$gtdb_taxonomy)) %>% 
  tidyr::separate(., col = gtdb_taxonomy, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "\\;", remove = T, extra = "drop") %>% mutate_all(na_if,"")
barplot_df$value <- as.numeric(barplot_df$value)
barplot_df <- rename(barplot_df, taxonomy = L5)

#Take top 10
top_taxa <- group_by(barplot_df, taxonomy) %>% summarise(., top_taxa_tmp = sum(value)) %>% arrange(., desc(top_taxa_tmp)) %>% slice(., 1:7)
high_abundance <- split(top_taxa$taxonomy, 1:NROW(top_taxa))
high_abundance <- high_abundance[!is.na(high_abundance)]

#Replace not top 10 with other.
barplot_df$taxonomy[barplot_df$taxonomy %in% high_abundance != "TRUE"] <- "Other"
barplot_df2 <- aggregate(barplot_df$value, by=list(taxonomy=barplot_df$taxonomy, Var2 = barplot_df$Var2), FUN=sum) %>% merge(., metadata, by.x = "Var2", by.y = "row.names")

barplot_df2 <- barplot_df2[order(barplot_df2$taxonomy),] #Re order
barplot_df2 <- rbind(barplot_df2[!(barplot_df2$taxonomy == "Other"),],barplot_df2[(barplot_df2$taxonomy == "Other"),]) #Move other to bottom
barplot_df2$taxonomy <- factor(barplot_df2$taxonomy, levels = unique(barplot_df2$taxonomy)) #Fix the order

#Remove metadata that I dont want to plot
barplot_df2 <- barplot_df2[!(barplot_df2$MostMalignantPolypType %in% c("unknown", "adenocarcinoma")),]
#barplot_df2 <- barplot_df2[!(barplot_df2$SampleType %in% "lavage"),]

#Custom color pallette.
sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947" ,"gray")

ggplot(data = barplot_df2, aes(x = Var2, weight = x, fill = taxonomy)) +
  geom_bar(width = 1, color = "black", size = .2) +
  theme_classic(base_size = 16) +
  facet_grid(.~SampleType, space = "free", scales = "free", labeller = 
               labeller(SampleType = c(`aspirate` = "Aspirate", `fecal` = "Fecal", `lavage` = "Lavage"))) + 
  #                      MostMalignantPolypType = c(`healthy`="Polyp free", `non-serrated`="Tubular Adenoma", `serrated`="Serrated"))) +
  scale_fill_manual(values = sarah_color) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 8), 
        strip.background = element_rect(fill="lightblue")) +
  labs(x = NULL,
       y = "Relative abundance", fill = "Family")

fig4a <- ggplot(data = barplot_df2) +
  aes(y = x, x = taxonomy, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), labels = c("Polyp free (78)", "Serrated polyp (58)", "Tubular adenoma (75)")) +
  geom_jitter(position = position_jitterdodge(), size = .5, alpha = .2) +
  facet_wrap(~SampleType, labeller = labeller(SampleType = c(`aspirate` = "Aspirate (156)", `fecal` = "Fecal (35)", `lavage` = "Lavage (20)"))) +
  labs(y = "Relative Abundance", x = " ", subtitle = bquote("Shotgun - Individuals:"~.(length(unique(barplot_df2$Patient))))) +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white")) +
  scale_x_discrete(limits = rev(levels(barplot_df2$taxonomy)))
fig4a

ggsave("Fig4a.svg", plot = fig4a, device = "svg", dpi = 300, height = 4, width = 12)

sfig4 <- ggplot(data = barplot_df2) +
  aes(y = x, x = taxonomy, fill = SampleType) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(name = "Sample type", values=c("gold", "darkorange4", "plum2"), labels = c("Aspirates (156)", "Fecal (35)", "Lavage (20)")) +
  geom_jitter(position = position_jitterdodge(), size = .5, alpha = .2) +
  labs(y = "Relative Abundance", x = "Phylum", subtitle = bquote("Shotgun - Individuals:"~.(length(unique(barplot_df2$Patient))))) +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white")) +
  scale_x_discrete(limits = rev(levels(barplot_df2$taxonomy)))
sfig4

#ggsave("supp_fig_4.png", sfig4, device = "png", dpi = 300, units = "in", height = 4, width = 6)

#Calculation for mean % of each taxon across sample type
mean_calc <- barplot_df %>% group_by(Var2,taxonomy) %>% summarise(x=sum(value)) %>% 
  merge(metadata, ., by.x = "row.names", by.y = "Var2") %>% group_by(SampleType,taxonomy) %>% summarise(y=mean(x))

#### Correlations (show network analysis) ####
#sparcc_corr <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/sparcc_corr.txt", check.names = F, row.names = 1)
#sparcc_pval <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/sparcc_pvals.txt", check.names = F, row.names = 1)

#heatmap_df = reshape2::melt(as.matrix(sparcc_corr))
#heatmap_df2 = reshape2::melt(as.matrix(sparcc_pval))
#heatmap_df$p_val = heatmap_df2$value #Combine p values and correlation values
#heatmap_df = heatmap_df %>% mutate(sig_p = ifelse(p_val <= .05, T, F))
#heatmap_df4 = heatmap_df %>% group_by(Var1) %>% summarise(n_p = sum(sig_p)) %>% filter(!n_p == 0) #Count the number of significant p values.
#heatmap_df = heatmap_df %>% filter(Var1 %in% heatmap_df4$Var1) #Remove the taxa from the original data frame that don't have any significant correlations.

#Unmelt filtered dataframe for hierarchical clustering
#heatmap_df5 = reshape2::dcast(data = heatmap_df[,1:3], formula = Var1 ~ Var2) %>% column_to_rownames(var = "Var1")
#heatmap_df6 = hclust(dist(heatmap_df5, method = "euclidean")) #Create hierarchical clustering
#heatmap_df$Var1 = factor(heatmap_df$Var1, levels = rownames(heatmap_df5)[heatmap_df6$order])

#ggplot(data = heatmap_df) +
#  geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  #geom_point(data = heatmap_df[heatmap_df$sig_p == T,], aes(x = Var2, y = Var1), shape = 8, size = 1) +
#  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
#  theme_classic() +
#  labs(x = NULL, y = NULL, fill = "SparCC correlation") +
#  theme(axis.text.x = element_text(angle = 90))

#elenta_cor <- as.data.frame(sparcc_corr$Eggerthella.lenta)
#elenta_pval <- as.data.frame(sparcc_pval$Eggerthella.lenta)

#elenta <- cbind(elenta_cor, elenta_pval)
#elenta$Taxon <- rownames(sparcc_corr)

#elenta2 <- elenta[elenta$`sparcc_pval$Eggerthella.lenta` <= 0.05,]
#elenta2 <- as.data.frame(elenta2[rev(order(abs(elenta2$`sparcc_corr$Eggerthella.lenta`)))[1:10],])

#elenta2 <- arrange(elenta2, elenta2$`sparcc_corr$Eggerthella.lenta`)
#elenta2$Taxon <- factor(elenta2$Taxon, levels = elenta2$Taxon)

#elenta_plot <- ggplot(data = elenta2) +
#  aes(x = `sparcc_corr$Eggerthella.lenta`, y = Taxon) +
#  geom_col(color = "black", size = .5) +
#  labs(title = "E. lenta", y = NULL, x = "Significant SparCC correlations") +
#  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
#        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), 
#        title = element_text(size = 12), legend.text = element_text(size = 12, color = "black"),
#        plot.title = element_text(face = "italic"))

#cscinds_cor <- as.data.frame(sparcc_corr$Clostridium.scindens)
#cscinds_pval <- as.data.frame(sparcc_pval$Clostridium.scindens)

#cscinds <- cbind(cscinds_cor, cscinds_pval)
#cscinds$Taxon <- rownames(sparcc_corr)

#cscinds2 <- cscinds[cscinds$`sparcc_pval$Clostridium.scindens` <= 0.05,]
#cscinds2 <- as.data.frame(cscinds2[rev(order(abs(cscinds2$`sparcc_corr$Clostridium.scindens`)))[1:10],])

#cscinds2 <- arrange(cscinds2, cscinds2$`sparcc_corr$Clostridium.scindens`)
#cscinds2$Taxon <- factor(cscinds2$Taxon, levels = cscinds2$Taxon)

#cscinds_plot <- ggplot(data = cscinds2) +
#  aes(x = `sparcc_corr$Clostridium.scindens`, y = Taxon) +
#  geom_col(color = "black", size = .5) +
#  labs(title = "C. scindens", y = NULL, x = "Significant SparCC correlations") +
#  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
#        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), 
#        title = element_text(size = 12), legend.text = element_text(size = 12, color = "black"), 
#        plot.title = element_text(face = "italic"))

blank <- ggplot(data = OTU_table) +
  geom_blank() +
  theme_classic() +
  theme(axis.line = element_blank())

#plot4e <- plot_grid(blank, cscinds_plot, elenta_plot, blank, nrow = 1, rel_widths = c(.2,1,1,.2))

#fig4 <- plot_grid(fig4a, plot4b_d, rf_taxa_plot, plot4e, ncol = 1, rel_heights = c(1.5,1,1,.8))
#ggsave("Figure_4.svg", plot = fig4, device = "svg", units = "in", dpi = 300, height = 18, width = 14)

#### Reviewer rebuttal analysis ####

#Make a list of low relative abundance taxa to filter.
#filter_list = decostand(merged_OTU_table[,-(1:(ncol(metadata)+1))], method = "total") %>% colMeans(.) %>% 
#  as.data.frame(.) %>% rownames_to_column() 
#names(filter_list)[2] = "relab"
#filter_list = filter_list %>% filter(relab < 0.001) #Get rid of microbes with an average relative abundance less than 0.1%

#Filter by sample prevelance (20%)
filter_list_tmp = asp_OTU %>% remove_rownames() %>% column_to_rownames(var = "Row.names") %>% 
  select(!Patient:PrepType)

filter_list = as.matrix(filter_list_tmp) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(filter_list_tmp)*0.2))

#Differential abundance analysis will be performed using kruskal wallis test with independent hypothesis weighting

#Healthy to tubular adenoma comparison
healthy_v_TA = asp_OTU %>% filter(!MostMalignantPolypType == "serrated") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(15:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
healthy_v_TA[,3:ncol(healthy_v_TA)] = decostand(healthy_v_TA[,3:ncol(healthy_v_TA)], method = "total")
healthy_v_TA = healthy_v_TA %>% select(Patient, MostMalignantPolypType, !filter_list$Var2)

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
names(pvals_KW_healthy_v_TA2)[1] = "relab_sum"
pvals_KW_healthy_v_TA2$pval = as.numeric(pvals_KW_healthy_v_TA2$pval)

IHW_healthy_v_TA = ihw(pvalues = pvals_KW_healthy_v_TA2$pval, covariates = pvals_KW_healthy_v_TA2$relab_sum, alpha = 0.05)
IHW_healthy_v_TA@df$taxa = rownames(pvals_KW_healthy_v_TA2)
IHW_healthy_v_TA@df$group1 = "healthy"
IHW_healthy_v_TA@df$group2 = "TA"

#Healthy to serrated comparison
healthy_v_ser = asp_OTU %>% filter(!MostMalignantPolypType == "TA") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(15:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
healthy_v_ser[,3:ncol(healthy_v_ser)] = decostand(healthy_v_ser[,3:ncol(healthy_v_ser)], method = "total")
healthy_v_ser = healthy_v_ser %>% select(Patient, MostMalignantPolypType, !filter_list$Var2)

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
names(pvals_KW_healthy_v_ser2)[1] = "relab_sum"
pvals_KW_healthy_v_ser2$pval = as.numeric(pvals_KW_healthy_v_ser2$pval)

IHW_healthy_v_ser = ihw(pvalues = pvals_KW_healthy_v_ser2$pval, covariates = pvals_KW_healthy_v_ser2$relab_sum, alpha = 0.05)
IHW_healthy_v_ser@df$taxa = rownames(pvals_KW_healthy_v_ser2)
IHW_healthy_v_ser@df$group1 = "healthy"
IHW_healthy_v_ser@df$group2 = "serrated"

#TA to serrated comparison
TA_v_ser = asp_OTU %>% filter(!MostMalignantPolypType == "healthy") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(15:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
TA_v_ser[,3:ncol(TA_v_ser)] = decostand(TA_v_ser[,3:ncol(TA_v_ser)], method = "total")
TA_v_ser = TA_v_ser %>% select(Patient, MostMalignantPolypType, !filter_list$Var2)

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
names(pvals_KW_TA_v_ser2)[1] = "relab_sum"
pvals_KW_TA_v_ser2$pval = as.numeric(pvals_KW_TA_v_ser2$pval)

IHW_TA_v_ser = ihw(pvalues = pvals_KW_TA_v_ser2$pval, covariates = pvals_KW_TA_v_ser2$relab_sum, alpha = 0.05)
IHW_TA_v_ser@df$taxa = rownames(pvals_KW_TA_v_ser2)
IHW_TA_v_ser@df$group1 = "TA"
IHW_TA_v_ser@df$group2 = "serrated"

#Visualizing
Polyp_type_DAb = bind_rows(IHW_healthy_v_TA@df, IHW_healthy_v_ser@df, IHW_TA_v_ser@df) %>% 
  mutate(final_pval = weighted_pvalue*3) %>% filter(final_pval <= 0.05)

relab_OTU_table = asp_OTU %>% group_by(as.factor(Patient)) %>% 
  summarise(across(`Methanosphaera.cuniculi`:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
relab_OTU_table[,3:ncol(relab_OTU_table)] = decostand(relab_OTU_table[,3:ncol(relab_OTU_table)], method = "total")
relab_OTU_table = relab_OTU_table %>% select(Patient, MostMalignantPolypType, unique(Polyp_type_DAb$taxa))
DAb_plot_df = reshape2::melt(relab_OTU_table)

fig4b = ggplot(data = DAb_plot_df) +
  aes(x = variable, y = value + 0.0001, fill = MostMalignantPolypType) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, size = 0.5) +
  scale_y_log10(limits = c(0.0001, 0.1)) +
  labs(x = NULL, y = "Relative abundance", fill = "Subject type", title = "Mucosal aspirates only") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  #ggpubr::stat_pvalue_manual(data = custom_pvals, y.position = 0.1, label = "final_pval") +
  scale_fill_manual(values=c("forestgreen", "steelblue", "firebrick3"), labels = c("Polyp free (35)", "Serrated polyp (30)", "Tubular adenoma (28)"))

#### Random forest ####
#Randomly dividing data using 2/3 split. Also removes taxa under 0.001 relative abundance.

#Healthy
training_healthy = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "healthy") %>% sample_frac(2/3)
training_healthy$MostMalignantPolypType = droplevels(training_healthy$MostMalignantPolypType)

test_healthy = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "healthy")
test_healthy = test_healthy[!rownames(test_healthy) %in% rownames(training_healthy),]
test_healthy$MostMalignantPolypType = droplevels(test_healthy$MostMalignantPolypType)

#Tubular adenomas
training_TA = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "TA") %>% sample_frac(2/3)
training_TA$MostMalignantPolypType = droplevels(training_TA$MostMalignantPolypType)

test_TA = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "TA")
test_TA = test_TA[!rownames(test_TA) %in% rownames(training_TA),]
test_TA$MostMalignantPolypType = droplevels(test_TA$MostMalignantPolypType)

#Serrated polyps
training_serr = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "serrated") %>% sample_frac(2/3)
training_serr$MostMalignantPolypType = droplevels(training_serr$MostMalignantPolypType)

test_serr = asp_OTU %>% select(MostMalignantPolypType, !filter_list$Var2) %>% 
  select(!Row.names:PrepType) %>% filter(MostMalignantPolypType == "serrated")
test_serr = test_serr[!rownames(test_serr) %in% rownames(training_serr),]
test_serr$MostMalignantPolypType = droplevels(test_serr$MostMalignantPolypType)

#Random forest training and testing
#Healthy to TA
model_healthy_TA = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_healthy, training_TA), 
                        proximity = T, importance = T, ntree = 501, num.cores = 32)
predict_healthy_TA = as.data.frame(predict(model_healthy_TA, newdata = rbind(test_healthy, test_TA), type = "prob"))
roc_healthy_TA = roc(rbind(test_healthy, test_TA)[,1], predict_healthy_TA[,2], ci=TRUE, ci.alpha=0.9, 
                     stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

#Healthy to serrated polyps
model_healthy_serr = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_healthy, training_serr), 
                             proximity = T, importance = F, ntree = 501, num.cores = 32)
predict_healthy_serr = as.data.frame(predict(model_healthy_serr, newdata = rbind(test_healthy, test_serr), type = "prob"))
roc_healthy_serr = roc(rbind(test_healthy, test_serr)[,1], predict_healthy_serr[,2], ci=TRUE, ci.alpha=0.9, 
                     stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

#Serrated to TA
model_serr_TA = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_serr, training_TA), 
                             proximity = T, importance = F, ntree = 501, num.cores = 32)
predict_serr_TA = as.data.frame(predict(model_serr_TA, newdata = rbind(test_serr, test_TA), type = "prob"))
roc_serr_TA = roc(rbind(test_serr, test_TA)[,1], predict_serr_TA[,2], ci=TRUE, ci.alpha=0.9, 
                     stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

svglite("ROC_curve.svg", width = 4, height = 4)
roc_curve = NULL
roc_curve = plot(roc_healthy_TA, print.auc = T, col = "orange3", print.auc.y = .3, print.auc.x = .8, main = "ROC curve")
roc_curve = plot(roc_healthy_serr, print.auc = T, col = "cyan4", print.auc.y = .2, print.auc.x = .8, add = TRUE)
roc_curve = plot(roc_serr_TA, print.auc = T, col = "purple2", print.auc.y = .1, print.auc.x = .8, add = TRUE)
dev.off()

VIP1 = as.data.frame(varImpPlot(model_healthy_TA, n.var = 10, type = 1)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10)
VIP1$rowname = factor(VIP1$rowname, unique(VIP1$rowname))
VIP1_plot = ggplot(data = VIP1) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "orange3", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP1$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Tubular adenoma") +
  theme_bw()

VIP2 = as.data.frame(varImpPlot(model_healthy_serr, n.var = 10, type = 1)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10)
VIP2$rowname = factor(VIP2$rowname, unique(VIP2$rowname))
VIP2_plot = ggplot(data = VIP2) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "cyan4", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP2$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Serrated polyp") +
  theme_bw()

VIP3 = as.data.frame(varImpPlot(model_serr_TA, n.var = 10, type = 1)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10)
VIP3$rowname = factor(VIP3$rowname, unique(VIP3$rowname))
VIP3_plot = ggplot(data = VIP3) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "purple2", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy+.5, y = rowname), hjust = 0, label = round(VIP3$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Tubular adenoma \n vs. Serrated polyp") +
  theme_bw()

fig4d = plot_grid(VIP1_plot, VIP2_plot, VIP3_plot, nrow = 1, labels = c("D.", "E.", "F."))
fig4bc = plot_grid(fig4b, blank, nrow = 1, labels = c("B.", "C."), rel_widths = c(1.25,1))
fig4 = plot_grid(fig4a, fig4bc, fig4d, ncol = 1, rel_heights = c(1.5,1,.8), labels = c("A.", "", ""))
ggsave("Figure_4.svg", fig4, device = "svg", dpi = 300, width = 12, height = 12)

#Supplemental figure for random forest
asp_OTU[,15:ncol(asp_OTU)] = decostand(x = asp_OTU[,15:ncol(asp_OTU)], method = "total", MARGIN = 1)
test = asp_OTU %>% select(VIP1$rowname, VIP2$rowname, VIP3$rowname) %>% rownames_to_column()
test2 = reshape2::melt(test) %>% merge(., metadata, by.x = "rowname", by.y = "row.names")

ggplot(data = test2) +
  aes(x = variable, y = value + 0.0001, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 0.25, alpha = 1/2) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Relative abundance", fill = "Subject type") +
  scale_fill_manual(labels = c("Polyp free", "Serrated", "Tubular adenoma"), 
                    values = c("forestgreen", "steelblue3", "firebrick3"))
ggsave("VIP_abundances.png", device = "png", dpi = 300, height = 4, width = 8)

### old random forest analysis ###
#Dashes in the OTU number messing things up. Writing a table and reading it fixes it in a lazy way.
#write.table(asp_OTU[!asp_OTU$MostMalignantPolypType == "unknown",], "/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", quote = F, sep = "\t")
#Randomforest <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", row.names=1, comment.char="#", check.names = F)

#rf_OTU_1 <- Randomforest[!Randomforest$MostMalignantPolypType == "serrated", 15:ncol(Randomforest)]
#rf_OTU_1 <- rf_OTU_1/rowSums(rf_OTU_1)
#rf_OTU_1 <- rf_OTU_1[,colMeans(rf_OTU_1) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
#rf_OTU_1 <- merge(rf_OTU_1, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
#rf_OTU_1$MostMalignantPolypType <- droplevels(rf_OTU_1$MostMalignantPolypType)

#rf_OTU$metadata <- sample(1:3, 154, replace = T)
#tuneRF(rf_OTU[,!(colnames(rf_OTU) == "metadata")], rf_OTU$metadata)
#rf_out_1 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_1, proximity = T, importance = F, ntree = 601, num.cores = 32)

#conf4b <- plotConfMat2(rf_out_1)
#plot4b <- proximityPlot2(rf_out_1, class.cols = c("forestgreen", "firebrick3"))
#plot4b <- plot4b$g + annotate("text", x = .4, y = -.45, label = conf4b[["labels"]][["title"]])
#varimps1 <- as.data.frame(varImpPlot(rf_out_1, type = 1, n.var = 5))
#top10_VIP.1 <- varimps1 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

#rf_OTU_2 <- Randomforest[!Randomforest$MostMalignantPolypType == "TA", 15:ncol(Randomforest)]
#rf_OTU_2 <- rf_OTU_2/rowSums(rf_OTU_2)
#rf_OTU_2 <- rf_OTU_2[,colMeans(rf_OTU_2) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
#rf_OTU_2 <- merge(rf_OTU_2, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
#rf_OTU_2$MostMalignantPolypType <- droplevels(rf_OTU_2$MostMalignantPolypType)

#rf_out_2 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_2, proximity = T, importance = F, ntree = 601, num.cores = 32)

#conf4c <- plotConfMat2(rf_out_2)
#impHeatmap(rf_out_2, alpha = 0.05, n = 10, ranks = F)
#plot4c <- proximityPlot2(rf_out_2, class.cols = c("forestgreen", "steelblue3"))
#plot4c <- plot4c$g + annotate("text", x = -.5, y = -.5, label = conf4c[["labels"]][["title"]]) +
#  coord_cartesian(clip = 'off')
#varimps2 <- as.data.frame(varImpPlot(rf_out_2, type = 1, n.var = 5))
#top10_VIP.2 <- varimps2 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

#rf_OTU_3 <- Randomforest[!Randomforest$MostMalignantPolypType == "healthy", 15:ncol(Randomforest)]
#rf_OTU_3 <- rf_OTU_3/rowSums(rf_OTU_3)
#rf_OTU_3 <- rf_OTU_3[,colMeans(rf_OTU_3) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
#rf_OTU_3 <- merge(rf_OTU_3, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
#rf_OTU_3$MostMalignantPolypType <- droplevels(rf_OTU_3$MostMalignantPolypType)

#rf_out_3 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_3, proximity = T, importance = F, ntree = 601, num.cores = 32)

#conf4d <- plotConfMat2(rf_out_3)
#impHeatmap(rf_out_3, alpha = 0.05, n = 10, ranks = F)
#plot4d <- proximityPlot2(rf_out_3, class.cols = c("steelblue3", "firebrick3"))
#plot4d <- plot4d$g + annotate("text", x = .5, y = -.45, label = conf4d[["labels"]][["title"]]) +
#  coord_cartesian(clip = 'off')
#varimps3 <- as.data.frame(varImpPlot(rf_out_3, type = 1, n.var = 5))
#top10_VIP.3 <- varimps3 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

#plot4b_d <- plot_grid(plot4b, plot4c, plot4d, nrow = 1)

#write.table(asp_OTU[!asp_OTU$MostMalignantPolypType %in% c("unknown"),], "/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", quote = F, sep = "\t")
#Randomforest <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", row.names=1, comment.char="#", check.names = F)
#rf_OTU <- Randomforest[,15:ncol(Randomforest)]
#rf_OTU <- rf_OTU/rowSums(rf_OTU)
#rf_OTU <- rf_OTU[,colMeans(rf_OTU) >= 0.002]
#top10_VIP <- rbind(top10_VIP.1, top10_VIP.2)

#for (i in 1:10) {
#  assign(paste0("tax_of_int", i), as.data.frame(rf_OTU[,colnames(rf_OTU) == print(top10_VIP$rowname[[i]])]))
#}

#tax_of_int_list <- lapply(seq(10), function(i) get(paste0("tax_of_int", i)))
#tax_of_int_list <- lapply(tax_of_int_list, function(x) {colnames(x)[1] <- "Taxon";   x})
#tax_of_int_list <- lapply(tax_of_int_list, function(x) {rownames(x) <- rownames(rf_OTU);   x})
#tax_of_int_list <- lapply(tax_of_int_list, function(x,y) {merge(x, y, by = "row.names")}, metadata)

#for (i in 1:10) {
#  assign(paste0("taxa_plot", i), 
#         ggplot(data = tax_of_int_list[[i]]) +
#           aes(x= MostMalignantPolypType, y = as.numeric(Taxon), fill = MostMalignantPolypType) +
#           geom_boxplot(outlier.shape = NA) +
#           theme_classic() +
#           geom_jitter(alpha = 0.5, width = .15) +
#           scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3")) +
#           scale_x_discrete(labels = c("Healthy", "Serrated", "TA")) +
#           theme(plot.title = element_text(face = "italic"), legend.position = "none") +
#           labs(title = print(top10_VIP$rowname[[i]]), y = "Relative abundance", x = NULL))
#}

#plot_list <- lapply(seq(10), function(i) get(paste0("taxa_plot", i)))
#rf_taxa_plot <- plot_grid(plotlist = plot_list, ncol = 5, nrow = 2)
#rf_taxa_plot

#Dunnets test, multiple samples per individual are averaged.
#dunns <- lapply(seq(10), function(i) {group_by(tax_of_int_list[[i]], Patient) %>% 
#    summarise(avg = mean(Taxon)) %>% merge(., metadata[,c(1,10)], by = "Patient") %>% .[!duplicated(.),] %>% 
#    dunn_test(data = ., formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")})