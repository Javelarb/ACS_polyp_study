library(tidyverse)
library(vegan)
library(nlme)
library(rfPermute)
library(rgl)
library(cluster)
library(car)
library(IHW)
library(MASS)
library(cowplot)
library(ape)
library(pROC)
library(svglite)

source("/media/julio/Storage/Software/ANCOM-master/programs/ancom.R")
source("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/CustomRFPermute.R")

options(scipen=10000, expressions = 5e5)
set.seed(seed = 999)

setwd("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data")

metadata <- read.delim("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/Shotgun_metadata.tsv", row.names=1, comment.char="#", check.names = F)
metadata$Patient <- as.factor(metadata$Patient)
metadata[metadata == "unknown"] <- NA

OTU_table <- read.delim("/media/julio/Storage/CRC/github/Shotgun_OTU_table.tsv", row.names=1, comment.char="#", check.names = F)
species <- read.delim("/media/julio/Storage/Software/IGGsearch/iggdb_v1.0.0/iggdb_v1.0.0.species", row.names=1, comment.char="#", check.names = F)

#### Read counts####
#read_counts_1 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/raw_read_counts.txt", row.names=1, quote="\"", comment.char="", check.names = F)
#read_counts_2 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/NH_read_counts.txt", row.names=1, quote="\"", comment.char="", check.names = F)
#read_counts_3 <- read.table("/media/julio/Storage/CRC/Summer_2019_colon_cancer_data/QF_read_counts.txt", row.names = 1, check.names = F)

#read_counts_all <- merge(read_counts_1, read_counts_3, by = "row.names")
#read_counts_all <- merge(read_counts_all, read_counts_2, by.x = "Row.names", by.y = "row.names")
#names(read_counts_all) <- c("SampleID", "Raw", "QF", "QF+Decon")
#read_counts_all <- reshape2::melt(read_counts_all) %>% merge(., metadata, by.x = "SampleID", by.y = "row.names")

#read_counts_all <- read_counts_all[!(read_counts_all$MostMalignantPolypType %in% c("adenocarcinoma", NA)),]
#read_counts_all <- read_counts_all[!(read_counts_all$SampleType %in% "lavage"),]

#ggplot(data = read_counts_all) +
#  aes(x = SampleType, y = value, fill = variable) +
#  geom_boxplot(outlier.shape = NA) +
#  labs(x = 'Sample Type',
#       y = 'PE150 read count', title = "Sample set 2 - Shotgun reads per sample") +
#  theme_bw(base_size = 14) +
  #facet_wrap(.~SampleType) +
#  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = 0.25, size = 1)

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
merged_OTU_table <- merged_OTU_table[!(merged_OTU_table$MostMalignantPolypType %in% c(NA,"adenocarcinoma")),]
row.names(merged_OTU_table) <- merged_OTU_table$Row.names

fecal_OTU <- merged_OTU_table[grepl("fecal", merged_OTU_table$SampleType),]
asp_OTU <- subset(merged_OTU_table, merged_OTU_table$SampleType == "aspirate")
lavage_OTU <- merged_OTU_table[grepl("lavage", merged_OTU_table$SampleType),]
asp_L_OTU <- asp_OTU[grepl("left", asp_OTU$ColonLocation),]
asp_R_OTU <- asp_OTU[grepl("right", asp_OTU$ColonLocation),]

#### alpha diversity ####
alpha_div <- as.data.frame(diversity(merged_OTU_table[,16:ncol(merged_OTU_table)], index = "shannon"))
alpha_div <- merge(alpha_div, metadata, by = "row.names")
#alpha_div <- alpha_div[!(alpha_div$SampleType == "lavage"),]

fig2a1 <- ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = `diversity(merged_OTU_table[, 16:ncol(merged_OTU_table)], index = "shannon")`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("Shotgun - Individuals:"~.(length(unique(alpha_div$Patient))))) +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
fig2a1

qqp(alpha_div$`diversity(merged_OTU_table[, 16:ncol(merged_OTU_table)], index = "shannon")`, "norm") #Pretty good fit for lm
alpha_lm <- NULL
alpha_lm$x <- as.numeric(alpha_div$`diversity(merged_OTU_table[, 16:ncol(merged_OTU_table)], index = "shannon")`)
alpha_lm$y <- as.factor(alpha_div$SampleType)
alpha_lm$z <- as.factor(alpha_div$MostMalignantPolypType)
alpha_lm$p <- as.factor(alpha_div$Plate)
alpha_lm$i <- as.factor(alpha_div$Patient)
alpha_lm <- as.data.frame(alpha_lm)
alpha_lm <- within(alpha_lm, y <- relevel(y, "fecal"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))

#species counts instead of shannon
specno <- as.data.frame(specnumber(merged_OTU_table[,16:ncol(merged_OTU_table)]))
specno <- merge(specno, metadata, by = "row.names")
#specno <- specno[!(specno$SampleType == "lavage"),]

fig2a2 <- ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(merged_OTU_table[, 16:ncol(merged_OTU_table)])`, 
      fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "darkorange4", "plum2"), labels = c("", "Tubular adenoma", "Serrated polyp")) +
  scale_x_discrete(labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  theme(legend.position = "none", axis.text = element_text(size = 10))
fig2a2

qqp(specno$`specnumber(merged_OTU_table[, 16:ncol(merged_OTU_table)])`, "norm")

specno_lm <- NULL
specno_lm$x <- as.numeric(specno$`specnumber(merged_OTU_table[, 16:ncol(merged_OTU_table)])`)
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
  aes(x = as.character(SampleType), y = alpha_div$`diversity(merged_OTU_table[, 16:ncol(merged_OTU_table)], index = "shannon")`, 
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
  aes(x = as.character(SampleType), y = `specnumber(merged_OTU_table[, 16:ncol(merged_OTU_table)])`, 
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
NMDS_points <- NMDS$points %>% merge(., metadata, by = "row.names") %>% filter(!MostMalignantPolypType %in% c(NA, "adenocarcinoma"))

asp_dm <- as.data.frame(as.matrix(vegdist(asp_OTU[,16:ncol(asp_OTU)], method = "bray")))
asp_pca <- pcoa(asp_dm)
merged_asp_pca <- merge(asp_pca$vectors[,1:2], metadata, by = "row.names")

fecal_dm <- as.data.frame(as.matrix(vegdist(fecal_OTU[,16:ncol(fecal_OTU)], method = "bray")))
fecal_pca <- pcoa(fecal_dm)
merged_fecal_pca <- merge(fecal_pca$vectors[,1:2], metadata, by = "row.names")

lavage_dm <- as.data.frame(as.matrix(vegdist(lavage_OTU[,16:ncol(lavage_OTU)], method = "bray")))
lavage_pca <- pcoa(lavage_dm)
merged_lavage_pca <- merge(lavage_pca$vectors[,1:2], metadata, by = "row.names")

#samp_type_perma
perma1_df = merged_OTU_table[complete.cases(merged_OTU_table),]
perma1 <- adonis2(formula = perma1_df[,16:ncol(perma1_df)] ~ as.numeric(perma1_df$BMI) + 
                   as.numeric(perma1_df$Age) + perma1_df$Ethnicity + perma1_df$Gender + 
                   (perma1_df$MostMalignantPolypType /as.character(perma1_df$Patient) / perma1_df$SampleType), 
                 data = perma1_df, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma1_df$Plate))
perma1

#asp perma
perma2_df = asp_OTU[complete.cases(asp_OTU),]
perma2 <- adonis2(formula = perma2_df[,16:ncol(perma2_df)] ~ as.numeric(perma2_df$BMI) + as.numeric(perma2_df$Age) + perma2_df$Ethnicity 
                 + perma2_df$Gender + perma2_df$ColonLocation + as.character(perma2_df$PrepType) + perma2_df$MostMalignantPolypType / as.character(perma2_df$Patient), 
                 data = perma2_df, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df$Plate))
perma2

#asp simper test
asp_simper = summary(simper(perma2_df[,16:ncol(perma2_df)], group = perma2_df$MostMalignantPolypType, permutations = how()))

#asp mantel test
dist2 = read.delim("/media/julio/Storage/CRC/github/pathway_BC_mat.txt", row.names=1, check.names = F) %>% 
  filter(row.names(.) %in% row.names(asp_dm)) %>% select_if(colnames(.) %in% colnames(asp_dm))
dist3 = read.delim("/media/julio/Storage/CRC/github/gene_BC_mat.txt", row.names=1, check.names = F)  %>% 
  filter(row.names(.) %in% row.names(asp_dm)) %>% select_if(colnames(.) %in% colnames(asp_dm))

mantel(asp_dm, dist2) #p = 0.001 and r = 0.33
mantel(asp_dm, dist3) #p = 0.001 and r = 0.70

#Pairwise permanovas with aspirates only
#Healthy vs. TA
perma2_df_hta = perma2_df %>% filter(!MostMalignantPolypType == "serrated")
pairwise_perma1 = adonis2(formula = perma2_df_hta[,16:ncol(perma2_df_hta)] ~ as.numeric(perma2_df_hta$BMI) + as.numeric(perma2_df_hta$Age) + perma2_df_hta$Ethnicity 
        + perma2_df_hta$Gender + perma2_df_hta$ColonLocation + as.character(perma2_df_hta$PrepType) + perma2_df_hta$MostMalignantPolypType / as.character(perma2_df_hta$Patient), 
        data = perma2_df_hta, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_hta$Plate))
#Polyp type R2 = 1.9%

#Healthy vs. Serrated
perma2_df_hs = perma2_df %>% filter(!MostMalignantPolypType == "TA")
pairwise_perma2 =adonis2(formula = perma2_df_hs[,16:ncol(perma2_df_hs)] ~ as.numeric(perma2_df_hs$BMI) + as.numeric(perma2_df_hs$Age) + perma2_df_hs$Ethnicity 
        + perma2_df_hs$Gender + perma2_df_hs$ColonLocation + as.character(perma2_df_hs$PrepType) + perma2_df_hs$MostMalignantPolypType / as.character(perma2_df_hs$Patient), 
        data = perma2_df_hs, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_hs$Plate))
#Polyp type R2 = 1.5%

#TA vs. Serrated
perma2_df_tas = perma2_df %>% filter(!MostMalignantPolypType == "healthy")
pairwise_perma3 = adonis2(formula = perma2_df_tas[,16:ncol(perma2_df_tas)] ~ as.numeric(perma2_df_tas$BMI) + as.numeric(perma2_df_tas$Age) + perma2_df_tas$Ethnicity 
        + perma2_df_tas$Gender + perma2_df_tas$ColonLocation + as.character(perma2_df_tas$PrepType) + perma2_df_tas$MostMalignantPolypType / as.character(perma2_df_tas$Patient), 
        data = perma2_df_tas, method = "bray", permutations = 999, parallel = 32, strata = as.factor(perma2_df_tas$Plate))
#Polyp type R2 = 2.7%

#Serrated vs healthy has the smallest effect size. TA vs serrated has strongest.
#write.table(bind_rows(pairwise_perma1, pairwise_perma2, pairwise_perma3), file = "Pairwise_permanovas.txt", sep = "\t", quote = F, row.names = T, col.names = T)

#fecal perma
perma3_df = fecal_OTU[complete.cases(fecal_OTU),]
perma3 <- adonis2(formula = perma3_df[,16:ncol(perma3_df)] ~ perma3_df$MostMalignantPolypType 
                 + as.numeric(perma3_df$Age) + as.numeric(perma3_df$BMI) + perma3_df$Gender 
                 + perma3_df$Ethnicity, data = perma3_df, method = "bray", permutations = 999, parallel = 32, 
                 strata = as.factor(perma3_df$Plate))
perma3

#Lavage
perma4_df = lavage_OTU[complete.cases(lavage_OTU),]
perma4 <- adonis2(formula = perma4_df[,16:ncol(perma4_df)] ~ perma4_df$MostMalignantPolypType + perma4_df$PrepType + 
                   as.numeric(perma4_df$Age) + as.numeric(perma4_df$BMI) + perma4_df$Gender + 
                   perma4_df$Ethnicity, data = perma4_df, method = "bray", permutations = 999, parallel = 32,
                 strata = as.factor(perma4_df$Plate))
perma4

#General plot
ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = SampleType, fill = MostMalignantPolypType), size = 6, alpha = 0.6) + 
  labs(x= "MDS1", y = "MDS2", subtitle = bquote("Individuals:"~.(length(unique(NMDS_points$Patient))))) +
  scale_shape_manual(values = c(21,22,24), name = "Sample type", labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  stat_ellipse(linetype = 2, aes(group = SampleType), size = 1) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
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
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Sample type", values=c("gold", "darkorange4", "plum2"), 
                    labels = c("Aspirate (156)", "Fecal (35)", "Lavage (20)")) +
  geom_text(label = NMDS_points$Patient, size = 2, color = "black") +
  annotate("text", x = 1.25, y = -1.2, label = bquote("Stress ="~.(round(NMDS$stress, digits = 2)))) 
fig2b

#ggsave("Figure_2b.svg", plot = plot_grid(fig2a, fig2b, nrow = 1, rel_widths = c(.35,.65), labels = c("E.", "F.")), 
#       device = "svg", units = "in", dpi = 300, height = 5, width = 11)

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
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
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
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
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
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_lavage_pca$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#### ANCOM v2.1 ####
#Notes: ANCOM is necessary for sample type comparisons due to repeated measurements.

ancom_otu1 <- merged_OTU_table[!(merged_OTU_table$SampleType == "lavage"),] #Compare fecal with asps
ancom_otu1 <- as.data.frame(t(ancom_otu1[,-(1:(1+ncol(metadata)))])) # Transpose otu table.

ancom_meta1 <- metadata
ancom_meta1$Patient <- as.factor(ancom_meta1$Patient)
ancom_meta1$Row.names <- row.names(metadata)
rownames(ancom_meta1) <- NULL

ancom_pre1 <- feature_table_pre_process(feature_table = ancom_otu1, meta_data = ancom_meta1, sample_var = "Row.names",
                                        group_var = "SampleType", lib_cut = 0, neg_lb = T)

ancom_res1 <- ANCOM(feature_table = ancom_pre1$feature_table, meta_data = ancom_pre1$meta_data, struc_zero = ancom_pre1$structure_zeros, 
                    main_var = "SampleType", rand_formula = "~ 1 | Patient", p_adj_method = "BH")

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

#### Taxa plots ####
#Turn to relative abundance
relab_OTU = as.data.frame(t(OTU_table[rowSums(OTU_table) > 1000,])) %>% dplyr::select(!c("COMM-STD", "0725", "0728")) %>% 
  merge(., species, by = "row.names") %>% column_to_rownames(var = "Row.names")

relab_OTU[,1:(ncol(relab_OTU)-10)] = decostand(x = relab_OTU[,1:(ncol(relab_OTU)-10)], method = "total", MARGIN = 2)

#Melt and separate taxonomy, remove whitespaces
barplot_df <- reshape2::melt(as.matrix(relab_OTU[,1:(ncol(relab_OTU)-10)])) %>% 
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
barplot_df2 <- barplot_df2[!(barplot_df2$MostMalignantPolypType %in% c("adenocarcinoma", NA)),]
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
  geom_jitter(position = position_jitterdodge(jitter.width = .12), size = .5, alpha = .2) +
  facet_wrap(~SampleType, labeller = labeller(SampleType = c(`aspirate` = "Aspirate (156)", `fecal` = "Fecal (35)", `lavage` = "Lavage (20)"))) +
  labs(y = "Relative Abundance", x = " ", subtitle = bquote("Shotgun - Individuals:"~.(length(unique(barplot_df2$Patient))))) +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white")) +
  scale_x_discrete(limits = rev(levels(barplot_df2$taxonomy)))
fig4a

#ggsave("Fig4a.svg", plot = fig4a, device = "svg", dpi = 300, height = 4, width = 12)

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

blank <- ggplot(data = OTU_table) +
  geom_blank() +
  theme_classic() +
  theme(axis.line = element_blank())

## Plot of individual species
merged_relab = decostand(merged_OTU_table[,(ncol(metadata)+2):ncol(merged_OTU_table)], method = "total")
merged_relab = merged_relab[,!colSums(merged_relab) == 0]

indiv_relabs_df = reshape2::melt(as.matrix(merged_relab)) %>% merge(., metadata, by.x = "Var1", by.y = "row.names")

indiv_relabs_plot = ggplot(data = indiv_relabs_df) +
  aes(x = SampleType, y = value + 0.0001, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1), size = .1, alpha = .2) +
  facet_wrap(.~Var2, ncol = 35)+
  theme_bw() +
  scale_y_log10() +
  labs(x = NULL, y = "Relative abundance") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 4.5)) +
  scale_fill_manual(values = c("forestgreen", "steelblue3", "firebrick3"))
#ggsave(filename = "Individual_rel_abs.png", plot = indiv_relabs_plot, device = "png", width = 40, height = 49, limitsize = T)

#### Reviewer rebuttal analysis ####

#Filter by sample prevalence (1/3)
filter_list_tmp = asp_OTU %>% remove_rownames() %>% column_to_rownames(var = "Row.names") %>% 
  dplyr::select(!Patient:InFinalAnalysis)

filter_list = as.matrix(filter_list_tmp) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(filter_list_tmp)*(1/3)))

#Differential abundance analysis will be performed using kruskal wallis test with independent hypothesis weighting

#Healthy to tubular adenoma comparison
healthy_v_TA = asp_OTU %>% filter(!MostMalignantPolypType == "serrated") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(16:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
healthy_v_TA[,3:ncol(healthy_v_TA)] = decostand(healthy_v_TA[,3:ncol(healthy_v_TA)], method = "total")
healthy_v_TA = healthy_v_TA %>% dplyr::select(Patient, MostMalignantPolypType, !filter_list$Var2)

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

IHW_healthy_v_TA = ihw(pvalues = pvals_KW_healthy_v_TA2$pval, covariates = pvals_KW_healthy_v_TA2$relab_sum, alpha = 0.05, nbins = 6, nfolds = 5, null_proportion = T)

IHW_healthy_v_TA@df$taxa = row.names(pvals_KW_healthy_v_TA2)
IHW_healthy_v_TA@df$group1 = "healthy"
IHW_healthy_v_TA@df$group2 = "TA"

#Healthy to serrated comparison
healthy_v_ser = asp_OTU %>% filter(!MostMalignantPolypType == "TA") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(16:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
healthy_v_ser[,3:ncol(healthy_v_ser)] = decostand(healthy_v_ser[,3:ncol(healthy_v_ser)], method = "total")
healthy_v_ser = healthy_v_ser %>% dplyr::select(Patient, MostMalignantPolypType, !filter_list$Var2)

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

IHW_healthy_v_ser = ihw(pvalues = pvals_KW_healthy_v_ser2$pval, covariates = pvals_KW_healthy_v_ser2$relab_sum, alpha = 0.05, nbins = 6, nfolds = 5, null_proportion = T)

IHW_healthy_v_ser@df$taxa = row.names(pvals_KW_healthy_v_ser2)
IHW_healthy_v_ser@df$group1 = "healthy"
IHW_healthy_v_ser@df$group2 = "serrated"

#TA to serrated comparison
TA_v_ser = asp_OTU %>% filter(!MostMalignantPolypType == "healthy") %>% group_by(as.factor(Patient)) %>% 
  summarise(across(16:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
TA_v_ser[,3:ncol(TA_v_ser)] = decostand(TA_v_ser[,3:ncol(TA_v_ser)], method = "total")
TA_v_ser = TA_v_ser %>% dplyr::select(Patient, MostMalignantPolypType, !filter_list$Var2)

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

IHW_TA_v_ser = ihw(pvalues = pvals_KW_TA_v_ser2$pval, covariates = pvals_KW_TA_v_ser2$relab_sum, alpha = 0.05, nbins = 6, nfolds = 5, null_proportion = T)

IHW_TA_v_ser@df$taxa = row.names(pvals_KW_TA_v_ser2)
IHW_TA_v_ser@df$group1 = "TA"
IHW_TA_v_ser@df$group2 = "serrated"

#Visualizing
Polyp_type_DAb = bind_rows(IHW_healthy_v_TA@df, IHW_healthy_v_ser@df, IHW_TA_v_ser@df) %>% filter(adj_pvalue <= 0.25)

relab_OTU_table = asp_OTU %>% group_by(as.factor(Patient)) %>% 
  summarise(across(`Methanosphaera.cuniculi`:ncol(asp_OTU), mean)) %>% remove_rownames() %>% 
  merge(metadata[,c(1,10)], ., by.y = "as.factor(Patient)", by.x = "Patient") %>% .[!duplicated(.),]
relab_OTU_table[,3:ncol(relab_OTU_table)] = decostand(relab_OTU_table[,3:ncol(relab_OTU_table)], method = "total") #Round about way of doing Dunn's post-hoc test
relab_OTU_table = relab_OTU_table %>% dplyr::select(Patient, MostMalignantPolypType, unique(Polyp_type_DAb$taxa))
DAb_plot_df = reshape2::melt(relab_OTU_table)

fig4b = ggplot(data = DAb_plot_df) +
  aes(x = variable, y = value + 0.0001, fill = MostMalignantPolypType) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .15), alpha = 0.5, size = 0.5) +
  scale_y_log10(limits = c(0.0001, 0.1)) +
  labs(x = NULL, y = "Relative abundance", fill = "Subject type", title = "Mucosal aspirates only") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  #ggpubr::stat_pvalue_manual(data = custom_pvals, y.position = 0.1, label = "final_pval") +
  scale_fill_manual(values=c("forestgreen", "steelblue", "firebrick3"), labels = c("Polyp free (35)", "Serrated polyp (30)", "Tubular adenoma (28)"))
fig4b

#### Random forest ####
#Filter by sample prevelance (20%)
filter_list_tmp = asp_OTU %>% remove_rownames() %>% column_to_rownames(var = "Row.names") %>% 
  dplyr::select(!Patient:InFinalAnalysis)

filter_list = as.matrix(filter_list_tmp) %>% reshape2::melt(.) %>% mutate(not_zero = ifelse(value != 0, T, F)) %>% 
  group_by(Var2) %>% summarise(n_not_zero = sum(not_zero)) %>% filter(n_not_zero <= round(nrow(filter_list_tmp)*(1/3)))

#Randomly dividing data using 2/3 split.
#Healthy
training_healthy = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "healthy") %>% sample_frac(2/3)
training_healthy$MostMalignantPolypType = as.factor(training_healthy$MostMalignantPolypType)
training_healthy$MostMalignantPolypType = droplevels(training_healthy$MostMalignantPolypType)

test_healthy = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "healthy")
test_healthy = test_healthy[!rownames(test_healthy) %in% row.names(training_healthy),]
test_healthy$MostMalignantPolypType = as.factor(test_healthy$MostMalignantPolypType)
test_healthy$MostMalignantPolypType = droplevels(test_healthy$MostMalignantPolypType)

#Tubular adenomas
training_TA = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "TA") %>% sample_frac(2/3)
training_TA$MostMalignantPolypType = as.factor(training_TA$MostMalignantPolypType)
training_TA$MostMalignantPolypType = droplevels(training_TA$MostMalignantPolypType)

test_TA = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "TA")
test_TA = test_TA[!rownames(test_TA) %in% row.names(training_TA),]
test_TA$MostMalignantPolypType = as.factor(test_TA$MostMalignantPolypType)
test_TA$MostMalignantPolypType = droplevels(test_TA$MostMalignantPolypType)

#Serrated polyps
training_serr = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "serrated") %>% sample_frac(2/3)
training_serr$MostMalignantPolypType = as.factor(training_serr$MostMalignantPolypType)
training_serr$MostMalignantPolypType = droplevels(training_serr$MostMalignantPolypType)

test_serr = asp_OTU %>% dplyr::select(MostMalignantPolypType, !filter_list$Var2) %>% 
  dplyr::select(!Row.names:InFinalAnalysis) %>% filter(MostMalignantPolypType == "serrated")
test_serr = test_serr[!rownames(test_serr) %in% row.names(training_serr),]
test_serr$MostMalignantPolypType = as.factor(test_serr$MostMalignantPolypType)
test_serr$MostMalignantPolypType = droplevels(test_serr$MostMalignantPolypType)

#Random forest training and testing
#Healthy to TA
model_healthy_TA = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_healthy, training_TA), 
                             proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_healthy_TA = as.data.frame(predict(model_healthy_TA, newdata = rbind(test_healthy, test_TA), type = "prob"))
roc_healthy_TA = roc(rbind(test_healthy, test_TA)[,1], predict_healthy_TA[,2], ci=TRUE, ci.alpha=0.9, 
                     stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

#Healthy to serrated polyps
model_healthy_serr = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_healthy, training_serr), 
                               proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_healthy_serr = as.data.frame(predict(model_healthy_serr, newdata = rbind(test_healthy, test_serr), type = "prob"))
roc_healthy_serr = roc(rbind(test_healthy, test_serr)[,1], predict_healthy_serr[,2], ci=TRUE, ci.alpha=0.9, 
                       stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

#Serrated to TA
model_serr_TA = rfPermute(formula = MostMalignantPolypType ~ ., data = bind_rows(training_serr, training_TA), 
                          proximity = T, importance = T, ntree = 501, num.cores = 32, num.rep = 999)
predict_serr_TA = as.data.frame(predict(model_serr_TA, newdata = rbind(test_serr, test_TA), type = "prob"))
roc_serr_TA = roc(rbind(test_serr, test_TA)[,1], predict_serr_TA[,2], ci=TRUE, ci.alpha=0.9, 
                  stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)

ROC_df1 = tibble(Specificity = roc_healthy_TA[["specificities"]], Sensitivity = roc_healthy_TA[["sensitivities"]], Comparison = "Polyp-free vs. TA")
ROC_df2 = tibble(Specificity = roc_healthy_serr[["specificities"]], Sensitivity = roc_healthy_serr[["sensitivities"]], Comparison = "Polyp-free vs. SP")
ROC_df3 = tibble(Specificity = roc_serr_TA[["specificities"]], Sensitivity = roc_serr_TA[["sensitivities"]], Comparison = "TA vs. SP")
ROC_df = rbind(ROC_df1, ROC_df2, ROC_df3)

ROC_df$Comparison = factor(x = ROC_df$Comparison, levels = unique(ROC_df$Comparison))

roc_curve = ggplot(data = ROC_df) +
  aes(x = Specificity, y = Sensitivity, color = Comparison) +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "gray") +
  geom_path() +
  scale_x_reverse() + 
  theme_bw() +
  scale_color_manual(values = c("orange3", "cyan4", "purple2")) +
  labs(title = "ROC Curve", subtitle = expression(italic("Mucosal aspirates only"))) +
  annotate("text", x = 0.7, y = .3, hjust = 0, color = "orange3", label = bquote("AUC:"~.(round(roc_healthy_TA[["auc"]][1], digits = 3))~.
                                                                                ("(")~.
                                                                                (round(roc_healthy_TA[["ci"]][1], digits = 2))~.
                                                                                ("-")~.
                                                                                (round(roc_healthy_TA[["ci"]][3], digits = 2))~.
                                                                               (")"))) +
  annotate("text", x = 0.7, y = .225, hjust = 0, color = "cyan4", label = bquote("AUC:"~.(round(roc_healthy_serr[["auc"]][1], digits = 3))~.
                                                                               ("(")~.
                                                                               (round(roc_healthy_serr[["ci"]][1], digits = 2))~.
                                                                               ("-")~.
                                                                               (round(roc_healthy_serr[["ci"]][3], digits = 2))~.
                                                                               (")"))) +
  annotate("text", x = 0.7, y = .15, hjust = 0, color = "purple2", label = bquote("AUC:"~.(round(roc_serr_TA[["auc"]][1], digits = 3))~.
                                                                               ("(")~.
                                                                               (round(roc_serr_TA[["ci"]][1], digits = 2))~.
                                                                               ("-")~.
                                                                               (round(roc_serr_TA[["ci"]][3], digits = 2))~.
                                                                               (")")))

VIP1 = as.data.frame(importance(model_healthy_TA)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP1$rowname = factor(VIP1$rowname, unique(VIP1$rowname))
VIP1_plot = ggplot(data = VIP1) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "orange3", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP1$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Tubular adenoma") +
  theme_bw()

VIP2 = as.data.frame(importance(model_healthy_serr)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP2$rowname = factor(VIP2$rowname, unique(VIP2$rowname))
VIP2_plot = ggplot(data = VIP2) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "cyan4", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy + .5, y = rowname), hjust = 0, label = round(VIP2$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Polyp free \n vs. Serrated polyp") +
  theme_bw()

VIP3 = as.data.frame(importance(model_serr_TA)) %>% rownames_to_column() %>% 
  arrange(MeanDecreaseAccuracy) %>% top_n(10, wt = MeanDecreaseAccuracy)
VIP3$rowname = factor(VIP3$rowname, unique(VIP3$rowname))
VIP3_plot = ggplot(data = VIP3) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = rowname, yend = rowname), lty = 2) +
  geom_point(aes(x = MeanDecreaseAccuracy, y = rowname), pch = 21, fill = "purple2", size = 3) +
  #geom_text(aes(x = MeanDecreaseAccuracy+.5, y = rowname), hjust = 0, label = round(VIP3$MeanDecreaseAccuracy, digits = 1)) +
  labs(x = "Mean Decrease Accuracy", y = NULL, title = "Tubular adenoma \n vs. Serrated polyp") +
  theme_bw()

#Supplemental figure for random forest
asp_OTU[,16:ncol(asp_OTU)] = decostand(x = asp_OTU[,16:ncol(asp_OTU)], method = "total", MARGIN = 1)
VIP_df = asp_OTU %>% dplyr::select(VIP1$rowname, VIP2$rowname, VIP3$rowname) %>% rownames_to_column()
VIP_df2 = reshape2::melt(VIP_df) %>% merge(., metadata, by.x = "rowname", by.y = "row.names")

VIP_abundances = ggplot(data = VIP_df2) +
  aes(x = variable, y = value + 0.0001, fill = MostMalignantPolypType) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .1), size = 0.075, alpha = 1/3) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Relative abundance", fill = "Subject type", title = "Random Forest VIPs", subtitle = expression(italic("Mucosal aspirates only"))) +
  scale_fill_manual(labels = c("Polyp free (64)", "Serrated (45)", "Tubular adenoma (47)"), 
                    values = c("forestgreen", "steelblue3", "firebrick3"))

fig4bd = plot_grid(roc_curve, VIP1_plot, VIP2_plot, nrow = 1, labels = c("b.", "c.", "d."), rel_widths = c(1.4,1.15,1))
fig4ef = plot_grid(VIP3_plot, VIP_abundances, nrow = 1, labels = c("e.", "f."), rel_widths = c(1,2.5))
fig4 = plot_grid(fig4a, fig4bd, fig4ef, ncol = 1, rel_heights = c(1.4,1,1), labels = c("a.", "", ""))

ggsave("Figure_4.svg", device = "svg", dpi = 600, height = 12, width = 13)
