library(reshape2)
library(ggplot2)
library(vegan)
library(nlme)
library(dplyr)
library(rfPermute)
library(rgl)
library(EcolUtils)
library(ggrepel)
library(cluster)
library(car)
library(MASS)
library(cowplot)
library(ape)
library(rstatix)

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
read_counts_all <- melt(read_counts_all) %>% merge(., metadata, by.x = "SampleID", by.y = "row.names")

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
OTU_table_plus_taxonomy <- OTU_table_plus_taxonomy[!(rownames(OTU_table_plus_taxonomy) %in% "COMM-STD"),] #Take it out

#For picking rarefaction depth.
sort(rowSums(OTU_table_plus_taxonomy))

#Permute
rared_OTU <- as.data.frame((rrarefy.perm(OTU_table_plus_taxonomy, sample = 1500, n = 10, round.out = T))) #Change 'sample =' parameter to rarefaction depth.
rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= 1500-(1500*.1), colSums(rared_OTU) >= 2])
#write.table(rared_OTU, "rared_OTU_table.txt", sep = "\t", quote = F)

merged_rared_table <- merge(metadata, rared_OTU, by = "row.names")
merged_rared_table <- merged_rared_table[!(merged_rared_table$MostMalignantPolypType %in% c("unknown","adenocarcinoma")),]
row.names(merged_rared_table) <- merged_rared_table$Row.names

fecal_OTU <- merged_rared_table[grepl("fecal", merged_rared_table$SampleType),]
asp_OTU <- subset(merged_rared_table, merged_rared_table$SampleType == "aspirate")
lavage_OTU <- merged_rared_table[grepl("lavage", merged_rared_table$SampleType),]
asp_L_OTU <- asp_OTU[grepl("left", asp_OTU$ColonLocation),]
asp_R_OTU <- asp_OTU[grepl("right", asp_OTU$ColonLocation),]
#write.table(t(asp_OTU[,-(1:14)]), file = "rared_OTU_table.txt", sep = '\t', quote = FALSE)

#### alpha diversity ####
alpha_div <- as.data.frame(diversity(merged_rared_table[,15:ncol(merged_rared_table)], index = "shannon"))
alpha_div <- merge(alpha_div, metadata, by = "row.names")
#alpha_div <- alpha_div[!(alpha_div$SampleType == "lavage"),]

fig2a1 <- ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = `diversity(merged_rared_table[, 15:ncol(merged_rared_table)], index = "shannon")`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("Shotgun - Individuals:"~.(length(unique(alpha_div$Patient))))) +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("Healthy", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))
fig2a1

qqp(alpha_div$`diversity(merged_rared_table[, 15:ncol(merged_rared_table)], index = "shannon")`, "norm") #Pretty good fit for lm
alpha_lm <- NULL
alpha_lm$x <- as.numeric(alpha_div$`diversity(merged_rared_table[, 15:ncol(merged_rared_table)], index = "shannon")`)
alpha_lm$y <- as.factor(alpha_div$SampleType)
alpha_lm$z <- as.factor(alpha_div$MostMalignantPolypType)
alpha_lm$p <- as.factor(alpha_div$Plate)
alpha_lm$i <- as.factor(alpha_div$Patient)
alpha_lm <- as.data.frame(alpha_lm)
alpha_lm <- within(alpha_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))

#species counts instead of shannon
specno <- as.data.frame(specnumber(merged_rared_table[,15:ncol(merged_rared_table)]))
specno <- merge(specno, metadata, by = "row.names")
#specno <- specno[!(specno$SampleType == "lavage"),]

fig2a2 <- ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(merged_rared_table[, 15:ncol(merged_rared_table)])`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type', title = "") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("Healthy", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))
fig2a2

qqp(specno$`specnumber(merged_rared_table[, 15:ncol(merged_rared_table)])`, "norm")

specno_lm <- NULL
specno_lm$x <- as.numeric(specno$`specnumber(merged_rared_table[, 15:ncol(merged_rared_table)])`)
specno_lm$y <- as.factor(specno$SampleType)
specno_lm$z <- as.factor(specno$MostMalignantPolypType)
specno_lm$i <- as.factor(specno$Patient)
specno_lm$p <- as.factor(specno$Plate)
specno_lm <- as.data.frame(specno_lm)
specno_lm <- within(specno_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = specno_lm, random = list(p=~1, i=~1)))

fig2a <- plot_grid(fig2a1, fig2a2, rel_widths = c(1,1))
fig2a
#ggsave("Figure_3a.svg", plot = fig3a, device = "svg", units = "in", dpi = 1000, height = 4, width = 7.5)

alpha_supp1 = ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = alpha_div$`diversity(merged_rared_table[, 15:ncol(merged_rared_table)], index = "shannon")`, 
      fill = as.character(MostMalignantPolypType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = "Sample set 2 - shotgun data") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Healthy (75)", "Tubular adenoma (78)", "Serrated polyp (58)")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp2 = ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(merged_rared_table[, 15:ncol(merged_rared_table)])`, 
      fill = as.character(MostMalignantPolypType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type', title = "") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Healthy (75)", "Tubular adenoma (78)", "Serrated polyp (58)")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp = plot_grid(alpha_supp1, alpha_supp2, rel_widths = c(1,1.65))
#ggsave("Shotgun_alpha_polyp.png", plot = alpha_supp, device = "png", units = "in", dpi = 300, height = 4, width = 10)

#### Beta div####
NMDS <- metaMDS(merged_rared_table[,-(1:(1+ncol(metadata)))], trymax = 999)
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
perma1 <- adonis(formula = merged_rared_table[,15:ncol(merged_rared_table)] ~ as.numeric(merged_rared_table$BMI) + 
          as.numeric(merged_rared_table$Age) + merged_rared_table$Ethnicity + merged_rared_table$Gender + 
          (merged_rared_table$MostMalignantPolypType /as.character(merged_rared_table$Patient) / merged_rared_table$SampleType), 
          data = merged_rared_table, method = "bray", permutations = 999, parallel = 32, strata = as.factor(merged_rared_table$Plate))
perma1
#capture.output(perma1, file = "perma_out.txt")

coef1 <- coefficients(perma1)["merged_rared_table$MostMalignantPolypType1",]
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
                    labels = c("Healthy (75)", "Serrated polyp (58)", "Tubular adenoma (78)")) +
  geom_text(label = NMDS_points$Patient, size = 3, color = "black")

fig2b <- ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_bw() +
  geom_point(aes(pch = MostMalignantPolypType, fill = SampleType), size = 5, alpha = 0.5) + 
  labs(x= "MDS1", y = "MDS2") +
  scale_shape_manual(values = c(21,22,24), name = "Subject type", labels = c("Healthy (75)", "Serrated polyp (58)", "Tubular adenoma (78)")) +
  #stat_ellipse(linetype = 2, aes(group = SampleType), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Sample type", values=c("gold", "darkorange4", "plum2"), 
                    labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  geom_text(label = NMDS_points$Patient, size = 2, color = "black") +
  annotate("text", x = 1.5, y = -1.2, label = bquote("Stress ="~.(round(NMDS$stress, digits = 2))))
fig2b

ggsave("Figure_2b.svg", plot = plot_grid(fig2a, fig2b, ncol = 1, rel_heights = c(.35,.65)), 
       device = "svg", units = "in", dpi = 300, height = 8, width = 7.5)

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
ancom_otu1 <- merged_rared_table[!(merged_rared_table$SampleType == "lavage"),] #Compare fecal with asps
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
diff_ab_df <- merged_rared_table[,colnames(merged_rared_table) %in% ancom_out$taxa_id] %>% merge(metadata, ., by = "row.names")

plot_data_column = function (column) {
  ggplot(diff_ab_df) +
    aes(x = SampleType, y = diff_ab_df[,column]) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    labs(y = colnames(diff_ab_df)[column])
}
DAb_taxa <- lapply((ncol(metadata)+2):ncol(diff_ab_df), plot_data_column)

# Performing on aspirates only. Compare between healthy and serr.
ancom_otu <- asp_OTU[!asp_OTU$MostMalignantPolypType == "TA",]
ancom_otu <- as.data.frame(t(ancom_otu[,-(1:(1+ncol(metadata)))])) # Transpose otu table.

ancom_meta <- metadata
ancom_meta$Patient <- as.factor(ancom_meta$Patient)
ancom_meta$Row.names <- rownames(metadata)
rownames(ancom_meta) <- NULL

ancom_pre <- feature_table_pre_process(feature_table = ancom_otu, meta_data = ancom_meta, sample_var = "Row.names",
                                  group_var = "MostMalignantPolypType", lib_cut = 0, neg_lb = T)

ancom_res <- ANCOM(feature_table = ancom_pre$feature_table, meta_data = ancom_pre$meta_data, struc_zero = ancom_pre$structure_zeros, 
              main_var = "MostMalignantPolypType", rand_formula = "~ 1 | Patient", p_adj_method = "fdr")

ancom_out1 <- ancom_res$out %>% filter(!(W == "Inf"), detected_0.7 == TRUE) %>% arrange(desc(W))

#### Random forest ####
#Dashes in the OTU number messing things up. Writing a table and reading it fixes it in a lazy way.
write.table(asp_OTU[!asp_OTU$MostMalignantPolypType == "unknown",], "/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", quote = F, sep = "\t")
Randomforest <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", row.names=1, comment.char="#", check.names = F)

rf_OTU_1 <- Randomforest[!Randomforest$MostMalignantPolypType == "serrated", 15:ncol(Randomforest)]
rf_OTU_1 <- rf_OTU_1/rowSums(rf_OTU_1)
rf_OTU_1 <- rf_OTU_1[,colMeans(rf_OTU_1) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
rf_OTU_1 <- merge(rf_OTU_1, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
rf_OTU_1$MostMalignantPolypType <- droplevels(rf_OTU_1$MostMalignantPolypType)

#rf_OTU$metadata <- sample(1:3, 154, replace = T)
#tuneRF(rf_OTU[,!(colnames(rf_OTU) == "metadata")], rf_OTU$metadata)
rf_out_1 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_1, proximity = T, importance = F, ntree = 601, num.cores = 32)

conf4b <- plotConfMat2(rf_out_1)
plot4b <- proximityPlot2(rf_out_1, class.cols = c("forestgreen", "firebrick3"))
plot4b <- plot4b$g + annotate("text", x = .4, y = -.45, label = conf4b[["labels"]][["title"]])
varimps1 <- as.data.frame(varImpPlot(rf_out_1, type = 1, n.var = 5))
top10_VIP.1 <- varimps1 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

rf_OTU_2 <- Randomforest[!Randomforest$MostMalignantPolypType == "TA", 15:ncol(Randomforest)]
rf_OTU_2 <- rf_OTU_2/rowSums(rf_OTU_2)
rf_OTU_2 <- rf_OTU_2[,colMeans(rf_OTU_2) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
rf_OTU_2 <- merge(rf_OTU_2, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
rf_OTU_2$MostMalignantPolypType <- droplevels(rf_OTU_2$MostMalignantPolypType)

rf_out_2 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_2, proximity = T, importance = F, ntree = 601, num.cores = 32)

conf4c <- plotConfMat2(rf_out_2)
impHeatmap(rf_out_2, alpha = 0.05, n = 10, ranks = F)
plot4c <- proximityPlot2(rf_out_2, class.cols = c("forestgreen", "steelblue3"))
plot4c <- plot4c$g + annotate("text", x = -.5, y = -.5, label = conf4c[["labels"]][["title"]]) +
  coord_cartesian(clip = 'off')
varimps2 <- as.data.frame(varImpPlot(rf_out_2, type = 1, n.var = 5))
top10_VIP.2 <- varimps2 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

rf_OTU_3 <- Randomforest[!Randomforest$MostMalignantPolypType == "healthy", 15:ncol(Randomforest)]
rf_OTU_3 <- rf_OTU_3/rowSums(rf_OTU_3)
rf_OTU_3 <- rf_OTU_3[,colMeans(rf_OTU_3) >= 0.002] #RF doesnt like having more features than samples, so remove low abundance taxa.
rf_OTU_3 <- merge(rf_OTU_3, metadata[,c(1,10)], by = "row.names") %>% tibble::column_to_rownames(., var = "Row.names") %>% dplyr::select(-Patient)
rf_OTU_3$MostMalignantPolypType <- droplevels(rf_OTU_3$MostMalignantPolypType)

rf_out_3 <- rfPermute(formula = MostMalignantPolypType ~ ., data = rf_OTU_3, proximity = T, importance = F, ntree = 601, num.cores = 32)

conf4d <- plotConfMat2(rf_out_3)
impHeatmap(rf_out_3, alpha = 0.05, n = 10, ranks = F)
plot4d <- proximityPlot2(rf_out_3, class.cols = c("steelblue3", "firebrick3"))
plot4d <- plot4d$g + annotate("text", x = .5, y = -.45, label = conf4d[["labels"]][["title"]]) +
  coord_cartesian(clip = 'off')
varimps3 <- as.data.frame(varImpPlot(rf_out_3, type = 1, n.var = 5))
top10_VIP.3 <- varimps3 %>% tibble::rownames_to_column() %>% arrange(desc(MeanDecreaseAccuracy)) %>% .[1:5,]

plot4b_d <- plot_grid(plot4b, plot4c, plot4d, nrow = 1)

#### Adding taxonomy ####
write.table(asp_OTU[!asp_OTU$MostMalignantPolypType %in% c("unknown"),], "/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", quote = F, sep = "\t")
Randomforest <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Randomforest.txt", row.names=1, comment.char="#", check.names = F)
rf_OTU <- Randomforest[,15:ncol(Randomforest)]
rf_OTU <- rf_OTU/rowSums(rf_OTU)
rf_OTU <- rf_OTU[,colMeans(rf_OTU) >= 0.002]
top10_VIP <- rbind(top10_VIP.1, top10_VIP.2)

for (i in 1:10) {
  assign(paste0("tax_of_int", i), as.data.frame(rf_OTU[,colnames(rf_OTU) == print(top10_VIP$rowname[[i]])]))
}

tax_of_int_list <- lapply(seq(10), function(i) get(paste0("tax_of_int", i)))
tax_of_int_list <- lapply(tax_of_int_list, function(x) {colnames(x)[1] <- "Taxon";   x})
tax_of_int_list <- lapply(tax_of_int_list, function(x) {rownames(x) <- rownames(rf_OTU);   x})
tax_of_int_list <- lapply(tax_of_int_list, function(x,y) {merge(x, y, by = "row.names")}, metadata)

for (i in 1:10) {
  assign(paste0("taxa_plot", i), 
         ggplot(data = tax_of_int_list[[i]]) +
           aes(x= MostMalignantPolypType, y = as.numeric(Taxon), fill = MostMalignantPolypType) +
           geom_boxplot(outlier.shape = NA) +
           theme_classic() +
           geom_jitter(alpha = 0.5, width = .15) +
           scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3")) +
           scale_x_discrete(labels = c("Healthy", "Serrated", "TA")) +
           theme(plot.title = element_text(face = "italic"), legend.position = "none") +
           labs(title = print(top10_VIP$rowname[[i]]), y = "Relative abundance", x = NULL))
}

plot_list <- lapply(seq(10), function(i) get(paste0("taxa_plot", i)))
rf_taxa_plot <- plot_grid(plotlist = plot_list, ncol = 5, nrow = 2)
rf_taxa_plot

#Dunnets test, multiple samples per individual are averaged.
dunns <- lapply(seq(10), function(i) {group_by(tax_of_int_list[[i]], Patient) %>% 
    summarise(avg = mean(Taxon)) %>% merge(., metadata[,c(1,10)], by = "Patient") %>% .[!duplicated(.),] %>% 
    dunn_test(data = ., formula = avg ~ MostMalignantPolypType, p.adjust.method = "fdr")})
  
#### Taxa barplots ####
#Turn to relative abundance
OTU_table3 <- OTU_table[!(rownames(OTU_table) %in% c("COMM-STD", "0725", "0728")),]
rared_OTU2 <- as.data.frame((rrarefy.perm(OTU_table3, sample = 1500, n = 10, round.out = T)))
rared_OTU2 <- as.data.frame(rared_OTU2[rowSums(rared_OTU2) >= 1500-(1500*.1), colSums(rared_OTU2) >= 2])
OTU_table_taxa <- merge(t(rared_OTU2), species, by = "row.names")
row.names(OTU_table_taxa) <- OTU_table_taxa$Row.names
OTU_table_taxa[,2:(nrow(rared_OTU2)+1)] <- OTU_table_taxa[,2:(nrow(rared_OTU2)+1)]/colSums(OTU_table_taxa[,2:(nrow(rared_OTU2)+1)])

#Melt and separate taxonomy, remove whitespaces
barplot_df <- reshape2::melt(as.matrix(OTU_table_taxa[,2:(nrow(rared_OTU2)+1)])) %>% 
  merge(., species[,c(2,9)], by.x = "Var1", by.y = "row.names") %>% mutate(gtdb_taxonomy = gsub(".__", "", .$gtdb_taxonomy)) %>% 
  tidyr::separate(., col = gtdb_taxonomy, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "\\;", remove = T, extra = "drop") %>% mutate_all(na_if,"")
barplot_df$value <- as.numeric(barplot_df$value)
barplot_df <- rename(barplot_df, taxonomy = L2)

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
  #                      MostMalignantPolypType = c(`healthy`="Healthy", `non-serrated`="Tubular Adenoma", `serrated`="Serrated"))) +
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
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "steelblue3", "firebrick3"), labels = c("Healthy (75)", "Serrated polyp (58)", "Tubular adenoma (78)")) +
  geom_jitter(position = position_jitterdodge(), size = .5, alpha = .2) +
  facet_wrap(~SampleType, labeller = labeller(SampleType = c(`aspirate` = "Aspirate (156)", `fecal` = "Fecal (35)", `lavage` = "Lavage (20)"))) +
  labs(y = "Relative Abundance", x = "Family", subtitle = bquote("Shotgun - Individuals:"~.(length(unique(barplot_df2$Patient))))) +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white")) +
  scale_x_discrete(limits = rev(levels(barplot_df2$taxonomy)))
fig4a

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

ggsave("supp_fig_4.png", sfig4, device = "png", dpi = 300, units = "in", height = 4, width = 6)

#Calculation for mean % of each taxon across sample type
mean_calc <- barplot_df %>% group_by(Var2,taxonomy) %>% summarise(x=sum(value)) %>% 
  merge(metadata, ., by.x = "row.names", by.y = "Var2") %>% group_by(SampleType,taxonomy) %>% summarise(y=mean(x))

#### Correlations ####
sparcc_corr <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/sparcc_corr.txt", check.names = F, row.names = 1)
sparcc_pval <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/sparcc_pvals.txt", check.names = F, row.names = 1)

elenta_cor <- as.data.frame(sparcc_corr$Eggerthella.lenta)
elenta_pval <- as.data.frame(sparcc_pval$Eggerthella.lenta)

elenta <- cbind(elenta_cor, elenta_pval)
elenta$Taxon <- rownames(sparcc_corr)

elenta2 <- elenta[elenta$`sparcc_pval$Eggerthella.lenta` <= 0.05,]
elenta2 <- as.data.frame(elenta2[rev(order(abs(elenta2$`sparcc_corr$Eggerthella.lenta`)))[1:10],])

elenta2 <- arrange(elenta2, elenta2$`sparcc_corr$Eggerthella.lenta`)
elenta2$Taxon <- factor(elenta2$Taxon, levels = elenta2$Taxon)

elenta_plot <- ggplot(data = elenta2) +
  aes(x = `sparcc_corr$Eggerthella.lenta`, y = Taxon) +
  geom_col(color = "black", size = .5) +
  labs(title = "E. lenta", y = NULL, x = "Significant SparCC correlations") +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), 
        title = element_text(size = 12), legend.text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "italic"))

cscinds_cor <- as.data.frame(sparcc_corr$Clostridium.scindens)
cscinds_pval <- as.data.frame(sparcc_pval$Clostridium.scindens)

cscinds <- cbind(cscinds_cor, cscinds_pval)
cscinds$Taxon <- rownames(sparcc_corr)

cscinds2 <- cscinds[cscinds$`sparcc_pval$Clostridium.scindens` <= 0.05,]
cscinds2 <- as.data.frame(cscinds2[rev(order(abs(cscinds2$`sparcc_corr$Clostridium.scindens`)))[1:10],])

cscinds2 <- arrange(cscinds2, cscinds2$`sparcc_corr$Clostridium.scindens`)
cscinds2$Taxon <- factor(cscinds2$Taxon, levels = cscinds2$Taxon)

cscinds_plot <- ggplot(data = cscinds2) +
  aes(x = `sparcc_corr$Clostridium.scindens`, y = Taxon) +
  geom_col(color = "black", size = .5) +
  labs(title = "C. scindens", y = NULL, x = "Significant SparCC correlations") +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), 
        title = element_text(size = 12), legend.text = element_text(size = 12, color = "black"), 
        plot.title = element_text(face = "italic"))

blank <- ggplot(data = OTU_table) +
  geom_blank() +
  theme_classic() +
  theme(axis.line = element_blank())

plot4e <- plot_grid(blank, cscinds_plot, elenta_plot, blank, nrow = 1, rel_widths = c(.2,1,1,.2))

fig4 <- plot_grid(fig4a, plot4b_d, rf_taxa_plot, plot4e, ncol = 1, rel_heights = c(1.5,1,1,.8))
#ggsave("Figure_4.svg", plot = fig4, device = "svg", units = "in", dpi = 300, height = 18, width = 14)

#### Taylor's power law ####
#merged_rared_table <- merge(metadata, rared_OTU, by = "row.names")
#df1 <- merged_rared_table[grepl("left", merged_rared_table$ColonLocation),] #Stuff I care about
#df2 <- merged_rared_table[grepl("right", merged_rared_table$ColonLocation),] #Stuff I care about pt 2
#df4 <- merged_rared_table[grepl("fecal", merged_rared_table$ColonLocation),]
#df2 <- merged_rared_table[merged_rared_table$Patient %in% df1$Patient,] #Stuff I care about pt 2
#No significant differences between polyp and healthy samples in any combination.
#df3 <- rbind(df1, df2)
#df3 <- df3 %>% group_by(Patient) %>% filter(n()>1) #Filter entries with only 1 per patient
#taylor_out <- NULL #create empty variable for loop

#for (i in df3$Patient) {
#  assign(paste0("OTU_table_", i), df3[df3$Patient == i,]) #will create lots of variables
#  log.m <- log(apply(get(paste0("OTU_table_", i))[,-(1:(1+ncol(metadata)))], 2, mean))
#  log.v <- log(apply(get(paste0("OTU_table_", i))[,-(1:(1+ncol(metadata)))], 2, var))
#  log.m <- log.m[is.finite(log.m)]
#  log.v <- log.v[is.finite(log.v)]
#  log.m_v <- merge(log.m, log.v, by = "row.names")
#  fit <- lm(log.m_v$y ~ log.m_v$x)
#  c.value <- as.numeric(coef(fit)[1]) #constant
#  z.value <- as.numeric(coef(fit)[2]) #var of interest, z > 1 = aggregated, z = 1 random, z < 1 uniform
#  individual <- i
#  individual <- cbind(individual, z.value)
#  taylor_out <- rbind(taylor_out, individual)
#}

#taylor_out <- as.data.frame(unique(taylor_out))
#taylor_final <- merge(taylor_out, metadata, by.x = "individual", by.y = "Patient")
#taylor_final <- taylor_final[,c(1:2,11)]
#taylor_final <- unique(taylor_final)
#taylor_final <- taylor_final[!grepl("unknown", taylor_final$MostMalignantPolypType),]

#ggplot(data = taylor_final) +
#  aes(x = taylor_final$MostMalignantPolypType, y = taylor_final$z.value, fill = taylor_final$MostMalignantPolypType) +
#  geom_boxplot(outlier.shape = NA, lwd = 1) +
#  theme_classic(base_size = 14, base_line_size = 1) +
#  labs(x = 'Subject status', y = 'Taylor`s exponent' , fill = 'Subject status', 
#       title = 'Proximal + distal') +
#  geom_point(position = position_jitterdodge(jitter.width = .5)) +
#  theme(legend.position = "none")
#TukeyHSD(aov(formula = taylor_final$z.value ~ taylor_final$MostMalignantPolypType))
#count <- taylor_final[grepl("serrated", taylor_final$MostMalignantPolypType),]