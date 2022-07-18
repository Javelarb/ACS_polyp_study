library(nlme)
library(tidyverse)
library(compositions)
library(vegan)
library(rfPermute)
library(ggrepel)
library(car)
library(MASS)
library(reshape2)
library(cowplot)
library(ape)

options(scipen=10000)
set.seed(seed = 999)

setwd("/media/julio/Storage/CRC/github/")

metadata <- read.delim("16S_metadata.tsv", stringsAsFactors=TRUE, row.names = 1)
OTU_table <- read.delim("16S_OTU_table.tsv", row.names=1, check.names = FALSE)
OTU_taxonomy <- read.delim("16S_taxonomy.tsv", row.names=1)
  
###### Making a boxplot of read counts for all samples ################
#Read_counts <- read.delim("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/Read_counts.tsv", stringsAsFactors = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
#Read_counts <- Read_counts[,-(2:4)]
#Read_counts <- reshape2::melt(as.matrix(Read_counts))

#read_counts_plus_metadata <- merge(Read_counts, metadata, by.x = "Var1", by.y = "row.names")

#ggplot(data = read_counts_plus_metadata) +
#  aes(x = SampleType, y = value, fill = Var2) +
#  geom_boxplot() +
#  labs(x = 'Sample Type',
#       y = 'PE300 read count', fill = "QF step", title = "Sample set 1 - 16S") +
#  theme_classic(base_size = 14) +
#  geom_point(position = position_jitterdodge(jitter.width = .2))

#TukeyHSD(aov(formula = read_counts_plus_metadata$value ~ as.character(read_counts_plus_metadata$SampleType)))
  
###### filtering  #####

#First, filter out mock community standards by name.
pos_ctrl <- OTU_table[,(names(OTU_table) %in% c("529","673"))] %>% .[!(rowSums(.) == 0),] %>% merge(., OTU_taxonomy, by = "row.names")

#Contaminants as determind by pos controls
OTU_table <- OTU_table[(!rownames(OTU_table) %in% c("f579c4bf7c2f8d98427bbd26e4b1d1cd", 
                                                  "7aca7ad71ebba7564f8e1fdf194cde9f", 
                                                  "7718b068b6862fe80b5fc06cc8f5a849", 
                                                  "4bf2b38d96348196aa1543cbae3a2b05", 
                                                  "227a951ca93c62c13a7870185fd2f116")),]

#remove mocks
OTU_no_mock <- OTU_table[,!(names(OTU_table) %in% c("529","673"))]

#We must merge them to filter out unassigned sequence IDs
OTU_table_plus_taxonomy <- as.data.frame(merge(OTU_taxonomy, OTU_no_mock, by.x = "row.names", by.y = "row.names"))

#Put the taxonomy name in place of ASV ID but keep unique
row.names(OTU_table_plus_taxonomy) <- make.names(OTU_table_plus_taxonomy$Taxon, unique = T)

#Filter all rows that contain mitochondria in their taxonomy.
OTU_table_no_mitochondia <- OTU_table_plus_taxonomy[!grepl("mitochondria", OTU_table_plus_taxonomy$`Taxon`),]
OTU_table_no_mitochondia <- OTU_table_no_mitochondia[,c(2,4:ncol(OTU_table_no_mitochondia))]
write.table(t(OTU_table_no_mitochondia), "OTU_clean_16S.txt", quote = F, sep = "\t", col.names = F)

OTU_clean <- read.delim("OTU_clean_16S.txt", row.names = 1, check.names = F)

#Generate rarefaction curve
rare_curve = rarecurve(OTU_clean, step = 1000, label = F, sample = 2500, xlab = "Read Depth", col = "orange")

Nmax <- sapply(rare_curve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rare_curve, max)

plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Read depth",
     ylab = "ASVs", type = "n", main = "16S")
for (i in seq_along(rare_curve)) {
  N <- attr(rare_curve[[i]], "Subsample")
  lines(N, rare_curve[[i]], col = "gray")
}
abline(v = 2500, lty = 2, col = "red")

#Didn't actually rarefy, just dropped samples under 2,500 sequences.
#Variable name is misleading post reviewer comments.
rared_OTU = OTU_clean[!rowSums(OTU_clean) <= 2500,]

#Performing the actual rarefaction
#rared_OTU <- as.data.frame((rrarefy.perm(OTU_clean, sample = rd, n = 10, round.out = T))) #Change 'sample =' parameter to rarefaction depth.

#This only keeps the samples that meet the rarefaction cutoff.
#rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= rd-(rd*.1), colSums(rared_OTU) >= 2]) #Change'x' in X-(X*.1) where X is rarefaction depth.

merged_OTU_table <- merge(metadata, rared_OTU, by = "row.names")
row.names(merged_OTU_table) <- merged_OTU_table$Row.names

lavage_OTU <- merged_OTU_table[grepl("lavage", merged_OTU_table$SampleType),]
asp_OTU <- merged_OTU_table[grepl("aspirate", merged_OTU_table$SampleType),]
brush_OTU <- merged_OTU_table[grepl("brush", merged_OTU_table$SampleType),]
asp_L_OTU <- asp_OTU[grepl("left", asp_OTU$ColonLocation),]
asp_R_OTU <- asp_OTU[grepl("right", asp_OTU$ColonLocation),]
  
###### alpha diversity ####
alpha_div <- as.data.frame(diversity(rared_OTU, index = "shannon"))
alpha_div <- merge(alpha_div, metadata, by = "row.names")
alpha_div <- alpha_div[!grepl("unknown", alpha_div$LabDiagnosis),]
  
fig2a1 <- ggplot(data = alpha_div) +
    aes(x = as.character(SampleType), y = `diversity(rared_OTU, index = "shannon")`, fill = as.character(SampleType)) +
    geom_boxplot(outlier.shape = NA, lwd = 1) +
    theme_classic(base_size = 14, base_line_size = 1) +
    labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("16S - Individuals:"~.(length(unique(alpha_div$Patient))))) +
    geom_point(position = position_jitterdodge(jitter.width = .2), alpha = 1/2) +
    scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Polyp free (39)", "Tubular adenoma (70)", "Serrated polyp (36)")) +
    scale_x_discrete(labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
    scale_y_continuous(breaks = seq(0,5, by = .5)) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
fig2a1
  
ggplot(data = subset(alpha_div, alpha_div$SampleType == "brush")) +
    aes(x = as.character(SampleType), y = `diversity(rared_OTU, index = "shannon")`, fill = SampleSite) +
    geom_boxplot(outlier.shape = NA, lwd = 1) +
    theme_classic(base_size = 14, base_line_size = 1) +
    labs(x = 'Sample type', y = 'Shannon index', fill = 'Sample site') +
    geom_point(position = position_jitterdodge(jitter.width = .1))

alpha_lm <- NULL
alpha_lm$x <- as.numeric(alpha_div$`diversity(rared_OTU, index = "shannon")`)
alpha_lm$y <- as.factor(alpha_div$SampleType)
alpha_lm$z <- as.factor(alpha_div$LabDiagnosis)
alpha_lm$p <- as.factor(alpha_div$Plate)
alpha_lm$i <- as.factor(alpha_div$Patient)
alpha_lm <- as.data.frame(alpha_lm)
alpha_lm <- within(alpha_lm, z <- relevel(z, "non-serrated"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(p=~1, i=~1)))

#species counts instead of shannon
specno <- as.data.frame(specnumber(rared_OTU))
specno <- merge(specno, metadata, by = "row.names")
specno <- specno[!grepl("unknown", specno$LabDiagnosis),]

fig2a2 <- ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(rared_OTU)`, fill = as.character(SampleType)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = 1/2) +
  scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Polyp free (39)", "Tubular adenoma (70)", "Serrated polyp (36)")) +
  scale_x_discrete(labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
  scale_y_continuous(breaks = seq(0,300, by = 50)) +
  theme(legend.position = "none", axis.text = element_text(size = 10)) 
fig2a2

ggplot(data = subset(specno, specno$SampleType == "brush")) +
  aes(x = as.character(SampleType), y = `specnumber(rared_OTU)`, fill = SampleSite) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = 'Sample type', y = 'ASV Richness', fill = 'Sample site') +
  geom_point(position = position_jitterdodge(jitter.width = .1))

specno_lm <- NULL
specno_lm$x <- as.numeric(specno$`specnumber(rared_OTU)`)
specno_lm$y <- as.factor(specno$SampleType)
specno_lm$z <- as.factor(specno$LabDiagnosis)
specno_lm$p <- as.factor(specno$Patient)
specno_lm <- as.data.frame(specno_lm)
summary(specno_lm <- lme(x ~ y * z, data = specno_lm, random = ~1|p))

fig2a <- plot_grid(fig2a1, fig2a2, rel_widths = c(1, 1), ncol = 1)
fig2a

alpha_supp1 = ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = `diversity(rared_OTU, index = "shannon")`, 
      fill = as.character(LabDiagnosis)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = "Sample set 1 - 16S data") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Polyp free (39)", "Tubular adenoma (70)", "Serrated polyp (36)")) +
  scale_x_discrete(labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp2 = ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(rared_OTU)`, 
      fill = as.character(LabDiagnosis)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type', title = "") +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Polyp free (39)", "Tubular adenoma (71)", "Serrated polyp (37)")) +
  scale_x_discrete(labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))

alpha_supp = plot_grid(alpha_supp1, alpha_supp2, rel_widths = c(1,1.65))

#ggsave("16S_alpha_polyp.png", plot = alpha_supp, device = "png", units = "in", dpi = 300, height = 4, width = 10)

###### Beta div ######
NMDS <- metaMDS(rared_OTU, trymax = 99)
NMDS_points <- as.data.frame(NMDS$points) %>% merge(., metadata, by = "row.names") %>% filter(!LabDiagnosis %in% c("unknown", "adenocarcinoma"))

asp_mds <- metaMDS(asp_OTU[,(ncol(metadata)+2):ncol(asp_OTU)])
merged_asp_mds <- merge(asp_mds$points[,1:2], metadata, by = "row.names")

brush_mds <- metaMDS(brush_OTU[,(ncol(metadata)+2):ncol(brush_OTU)])
merged_brush_mds <- merge(brush_mds$points[,1:2], metadata, by = "row.names")

lavage_mds <- metaMDS(lavage_OTU[,(ncol(metadata)+2):ncol(lavage_OTU)])
merged_lavage_mds <- merge(lavage_mds$points[,1:2], metadata, by = "row.names")

#merged_OTU_table <- merged_OTU_table[!merged_OTU_table$SampleType == "lavage",]
#samp_type_perma. sig diff between aspirates and brushes. but not lavage and brushes, and aspirate and lavage
perma1 <- adonis(formula = merged_OTU_table[,(ncol(metadata)+2):ncol(merged_OTU_table)] ~ as.numeric(merged_OTU_table$BMI) + 
                     as.numeric(merged_OTU_table$Age) + merged_OTU_table$Ethnicity + merged_OTU_table$Gender + 
                     (merged_OTU_table$LabDiagnosis /as.character(merged_OTU_table$Patient) / merged_OTU_table$SampleType), 
                   data = merged_OTU_table, method = "bray", permutations = 999, parallel = 32, strata = as.factor(merged_OTU_table$Plate))

#asp perma
perma2 <- adonis(formula = asp_OTU[,(ncol(metadata)+2):ncol(asp_OTU)] ~ as.numeric(asp_OTU$BMI) + as.numeric(asp_OTU$Age) + asp_OTU$Ethnicity 
               + asp_OTU$Gender + asp_OTU$ColonLocation + asp_OTU$LabDiagnosis / as.character(asp_OTU$Patient), 
                 data = asp_OTU, method = "bray", permutations = 999, parallel = 32, strata = as.factor(asp_OTU$Plate))

#Brush perma
perma3 <- adonis(formula = brush_OTU[,(ncol(metadata)+2):ncol(brush_OTU)] ~ as.numeric(brush_OTU$BMI) + as.numeric(brush_OTU$Age) + brush_OTU$Ethnicity 
                 + brush_OTU$Gender + brush_OTU$ColonLocation + brush_OTU$LabDiagnosis / as.character(brush_OTU$Patient) / brush_OTU$SampleSite, 
                 data = brush_OTU, method = "bray", permutations = 999, parallel = 32, strata = as.factor(brush_OTU$Plate))

#Lavage perma
perma4 <- adonis(formula = lavage_OTU[,(ncol(metadata)+2):ncol(lavage_OTU)] ~ lavage_OTU$LabDiagnosis +
          as.numeric(lavage_OTU$Age) + as.numeric(lavage_OTU$BMI) + lavage_OTU$Gender + 
          lavage_OTU$Ethnicity + as.character(lavage_OTU$Patient), data = lavage_OTU, 
          method = "bray", permutations = 999, parallel = 32, strata = as.factor(lavage_OTU$Plate))

#General plot
ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_classic() +
  geom_point(aes(pch = SampleType, fill = LabDiagnosis), size = 6, alpha = .6) + 
  labs(x= "MDS1", y = "MDS2", title = bquote("16S - Individuals:"~.(length(unique(NMDS_points$Patient))))) +
  scale_shape_manual(values = c(21,22,24), name = "Sample type", labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
  #stat_ellipse(linetype = 2, aes(group = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue3"), labels = c("Polyp free (39)", "Tubular adenoma (71)", "Serrated polyp (37)")) +
  geom_text(label = NMDS_points$Patient, size = 3, color = "black")

fig2b <- ggplot(data = NMDS_points) +
  aes(x = MDS1, y = MDS2) +
  theme_bw() +
  geom_point(aes(pch = LabDiagnosis, fill = SampleType), size = 5, alpha = .5) + 
  labs(x= "MDS1", y = "MDS2") +
  scale_shape_manual(values = c(21,22,24), name = "Subject type", 
                     labels = c("Polyp free (39)", "Serrated polyp (37)", "Tubular adenoma (71)"),
                     limits = c("healthy", "serrated", "non-serrated")) +
  #stat_ellipse(linetype = 2, aes(group = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(order = 1, override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Sample type", values=c("gold", "lightsalmon", "plum2"), labels = c("Aspirate (52)", "Brush (64)", "Lavage (31)")) +
  geom_text(label = NMDS_points$Patient, size = 2, color = "black") +
  annotate("text", x = 1.5, y = 3, label = bquote("Stress ="~.(round(NMDS$stress, digits = 2))))
fig2b

#ggsave("Figure_2a.svg", plot = plot_grid(fig2a, fig2b, nrow = 1, rel_widths = c(.35,.65), labels = c("A.", "B.")), 
#       device = "svg", units = "in", dpi = 300, height = 5, width = 11)
  
#aspirate only
ggplot(data = merged_asp_mds) +
  aes(x = MDS1, y = MDS2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = ColonLocation, fill = LabDiagnosis), size = 8) + 
  scale_shape_manual(values = c(23,22), name = "Colon side") +
  #ggtitle(bquote("Patient:"~R^2~"= 0.86, p = 0.001"), subtitle = ~ atop("Colon location:"~R^2~"= 0.001, p = 0.07",
  #                                                                     "Subject type:" ~R^2 ~ "= 0.06, p = 0.001")) +
  stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_asp_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#Brush only
ggplot(data = merged_brush_mds) +
  aes(x = MDS1, y = MDS2, fill = SampleSite) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 6, alpha = .5) +
  #ggtitle(bquote("Patient:"~R^2~"= 0.79, p = 0.001"), subtitle = ~ atop("Hyperlocal sampling site:"~R^2~"= 0.06, p = 0.79",
   #                                                                     "Subject type:" ~R^2 ~ "= 0.04, p = 0.001")) +
  #stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Hyperlocal Site") +
  geom_text(label = merged_brush_mds$Patient, size = 3) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))
  
#Lavages only
ggplot(data = merged_lavage_mds) +
  aes(x = MDS1, y = MDS2, fill = LabDiagnosis) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 8) +
  #ggtitle(bquote("Patient:"~R^2~"= 0.88, p = 0.68"), subtitle = bquote("Subject type:" ~R^2 ~ "= 0.08, p = 0.44")) +
  stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_lavage_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))
  
###### 16S Brushes versus opposite wall ##############
OTU_brush <- OTU_clean[(row.names(OTU_clean) %in% c("114","115","505","506","591","591-5","601","602", "614", "615", "686","686-5", "685","685-5")),]

#Alpha
brush_alpha_div <- as.data.frame(diversity(OTU_brush, index = "shannon"))
brush_alpha_div <- merge(brush_alpha_div, metadata, by = "row.names")
brush_alpha_div$SampleSite <- as.character(brush_alpha_div$SampleSite)
brush_alpha_div$SampleSite[!(brush_alpha_div$SampleSite == "healthy")] <- "polyp"

brush_lme <- NULL
brush_lme$x <- as.numeric(brush_alpha_div$`diversity(OTU_brush, index = "shannon")`)
brush_lme$y <- as.factor(brush_alpha_div$ColonLocation)
brush_lme$z <- as.factor(brush_alpha_div$SampleSite)
brush_lme$p <- as.factor(brush_alpha_div$Plate)
brush_lme$i <- as.factor(brush_alpha_div$Patient)
brush_lme <- as.data.frame(brush_lme)
summary(lme(x ~ z * y, data = brush_lme, random = list(p=~1, i=~1)))

fig3b1 <- ggplot(data = brush_alpha_div) +
  aes(x = SampleSite, y = `diversity(OTU_brush, index = "shannon")`, fill = SampleSite) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("Individuals:"~.(length(unique(brush_alpha_div$Patient))))) +
  geom_path(aes(group = Patient), linetype = 2) +
  geom_point(aes(pch = ColonLocation), fill = "white", size = 6, alpha = .8) + 
  scale_shape_manual(values = c(21,22), name = "Colon location", labels = c("Left", "Right")) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Tissue type", values=c("forestgreen", "firebrick3"), labels = c("Polyp free", "Polyp")) +
  geom_text(label = brush_alpha_div$Patient, size = 3, color = "black") +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "none")
fig3b1

#species counts instead of shannon
brush_specno <- as.data.frame(specnumber(OTU_brush))
brush_specno <- merge(brush_specno, metadata, by = "row.names")
brush_specno$SampleSite <- as.character(brush_specno$SampleSite)
brush_specno$SampleSite[!(brush_specno$SampleSite == "healthy")] <- "polyp"

brush_lme2 <- NULL
brush_lme2$x <- as.numeric(brush_specno$`specnumber(OTU_brush)`)
brush_lme2$y <- as.factor(brush_specno$ColonLocation)
brush_lme2$z <- as.factor(brush_specno$SampleSite)
brush_lme2$p <- as.factor(brush_specno$Plate)
brush_lme2$i <- as.factor(brush_specno$Patient)
brush_lme2 <- as.data.frame(brush_lme2)
summary(lme(x ~ z * y, data = brush_lme2, random = list(p=~1, i=~1)))

fig3b2 <- ggplot(data = brush_specno) +
  aes(x = SampleSite, y = `specnumber(OTU_brush)`, fill = SampleSite) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type', title = "") +
  geom_path(aes(group = Patient), linetype = 2) +
  geom_point(aes(pch = ColonLocation), fill = "white", size = 6, alpha = .8) + 
  scale_shape_manual(values = c(21,22), name = "Colon location", labels = c("Left", "Right")) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Tissue type", values=c("forestgreen", "firebrick3"), labels = c("Polyp free", "Polyp")) +
  geom_text(label = brush_specno$Patient, size = 3, color = "black") +
  theme(axis.line = element_line(color = "black"), panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12), legend.position = "none")
fig3b2

fig3b <- plot_grid(fig3b1, fig3b2)
#ggsave("Figure_2b.svg", plot = fig3b, device = "svg", units = "in", dpi = 300, height = 4, width = 6)

#Beta
brush_mds <- metaMDS(OTU_brush, trymax = 999)
brush_mds_merged <- merge(brush_mds$points[,1:2], metadata, by.x = "row.names", by.y = "row.names")

brush_mds_merged$SampleSite <- as.character(brush_mds_merged$SampleSite)
brush_mds_merged$SampleSite[!(brush_mds_merged$SampleSite == "healthy")] <- "polyp"
  
fig3c <- ggplot(data = brush_mds_merged) +
  aes(x = MDS1, y = MDS2) +
  theme_bw() +
  geom_path(aes(group = Patient), linetype = 2) +
  geom_point(aes(pch = ColonLocation, fill = SampleSite), size = 7, alpha = .8) + 
  labs(title = "") +
  scale_shape_manual(values = c(21,22), name = "Colon location", labels = c("Left", "Right")) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Tissue type", values=c("forestgreen", "firebrick3"), labels = c("Polyp free", "Polyp")) +
  geom_text(label = brush_mds_merged$Patient, size = 4, color = "black") +
  annotate("text", x = -2.1, y = 1.7, label = bquote("Stress ="~.(round(brush_mds$stress, digits = 2))), hjust = 0)
fig3c

#ggsave("Figure_2c.svg", plot = fig3c, device = "svg", units = "in", dpi = 300, height = 4, width = 5.5)

#permanova
OTU_brush2 <- merge(metadata, OTU_brush, by = "row.names")
OTU_brush2$SampleSite <- as.character(OTU_brush2$SampleSite)
OTU_brush2$SampleSite[!(OTU_brush2$SampleSite == "healthy")] <- "polyp"

adonis(formula = OTU_brush2[,15:ncol(OTU_brush2)] ~ OTU_brush2$ColonLocation + OTU_brush2$LabDiagnosis / 
      as.character(OTU_brush2$Patient) / as.character(OTU_brush2$SampleSite), 
       data = OTU_brush2[,15:ncol(OTU_brush2)], method = "bray", permutations = 999)

#Barplot
brush_relab_table <- decostand(OTU_brush, method = "total")

#Melt and separate taxonomy, remove whitespaces
brush_barplot_df <- reshape2::melt(as.matrix(brush_relab_table)) %>% mutate(Var2 = gsub(".__", "", .$Var2)) %>% 
  tidyr::separate(., col = Var2, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop") %>% mutate_all(na_if,"")
brush_barplot_df$value <- as.numeric(brush_barplot_df$value)

#Taxonomy level of interest
brush_barplot_df <- rename(brush_barplot_df, taxonomy = L6)

#Take top 10
top_taxa <- group_by(brush_barplot_df, taxonomy) %>% summarise(., top_taxa_tmp = sum(value)) %>% arrange(., desc(top_taxa_tmp)) %>% slice(., 1:11)
high_abundance <- split(top_taxa$taxonomy, 1:NROW(top_taxa)) 
high_abundance <- high_abundance[!is.na(high_abundance)]

#Replace not top 10 with other. 
brush_barplot_df$taxonomy[brush_barplot_df$taxonomy %in% high_abundance != "TRUE"] <- "Other"
brush_barplot_df2 <- aggregate(brush_barplot_df$value, by=list(taxonomy=brush_barplot_df$taxonomy, Var1 = brush_barplot_df$Var1), FUN=sum) %>% merge(., metadata, by.x = "Var1", by.y = "row.names")

brush_barplot_df2 <- brush_barplot_df2[order(brush_barplot_df2$taxonomy),] #Re order
brush_barplot_df2 <- rbind(brush_barplot_df2[!(brush_barplot_df2$taxonomy == "Other"),],brush_barplot_df2[(brush_barplot_df2$taxonomy == "Other"),]) #Move other to bottom
brush_barplot_df2$taxonomy <- factor(brush_barplot_df2$taxonomy, levels = unique(brush_barplot_df2$taxonomy)) #Fix the order

#Custom color pallette.
sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947" ,"gray")
polyplabels1 <- c("TA","healthy","HPP","healthy","TA","healthy","SSP","healthy","TA","healthy","TA","healthy","TA","healthy")

brush_barplot_df3 <- brush_barplot_df2[!brush_barplot_df2$Patient == 49,]
brush_barplot1 <- ggplot(data = brush_barplot_df3) +
  aes(x = as.character(Var1), weight = x, fill = taxonomy) +
  geom_bar(color = "black", size = .5) +
  theme_classic() +
  facet_grid(.~as.factor(Patient), scales = "free") +
  scale_fill_manual(values = sarah_color) +
  scale_x_discrete(labels = c("           ",""))+
  geom_text(label = brush_barplot_df3$SampleSite, y = -.13, size = 4) +
  geom_text(label = brush_barplot_df3$ColonLocation, y = 1.045, size = 4) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none", axis.title.y = element_text(size = 14)) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  coord_cartesian(clip = 'off') +
  ylim(0,1.04)

brush_barplot_df4 <- brush_barplot_df2[brush_barplot_df2$Patient == 49,]
brush_barplot2 <- ggplot(data = brush_barplot_df4, aes(x = as.character(Var1), weight = x, fill = taxonomy)) +
  geom_bar(color = "black", size = .5) +
  theme_classic() +
  facet_grid(.~as.factor(Patient), scales = "free") +
  scale_fill_manual(values = sarah_color) +
  scale_x_discrete(labels = c("           ","", "", "")) +
  geom_text(label = brush_barplot_df4$SampleSite, y = -.13, size = 4) +
  geom_text(label = brush_barplot_df4$ColonLocation, y = 1.045, size = 4) +
  theme(axis.text.x=element_text(angle = 90), axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(x = NULL, y = NULL, fill = "Genus") +
  coord_cartesian(clip = 'off') +
  ylim(0,1.04)

genus_barplot <- plot_grid(brush_barplot1, brush_barplot2, nrow = 1, rel_widths = c(1.6,1))
#ggsave("genus_barplot.svg", genus_barplot, device = "svg", dpi = 300, height = 5, width = 7)

#Blank plot
blank <- ggplot(data = brush_mds_merged) +
  geom_blank() +
  theme_classic() +
  theme(axis.line = element_blank())

fig3a_c <- plot_grid(blank, fig3b, fig3c, rel_widths = c(2,2.5,2.75), nrow = 1, labels = c("A.", "B.", "C."))
fig3 <- plot_grid(fig3a_c, genus_barplot, ncol = 1, rel_heights = c(.7,1), labels = c("", "D."))

#ggsave("Figure_3.svg", device = "svg", dpi = 600, height = 8, width = 12)
  
###### F. nucleatum investigation #########
#fuso_OTU <- as.data.frame(OTU_no_mock[(rownames(OTU_no_mock) %in% c("784613cec54042f06ac530324696f6bd", "3b9d5137ed8c9088b160e93286595664", "b44db83ce82cbbd0253f5b8131bf35a1",
#                                                             "57ca1171911f37a20631c6f0c04bd109", "ee74a96f7669f1a201ceb1705ec2fd5f")),])
#fuso_OTU <- as.data.frame(colSums(fuso_OTU))
#fuso_OTU <- (fuso_OTU/rd)
#fuso_OTU <- merge(fuso_OTU, metadata, by = "row.names")

#ggplot(data = fuso_OTU) +
#  aes(x = as.character(fuso_OTU$SampleType), y = fuso_OTU$`colSums(fuso_OTU)`, fill = as.character(fuso_OTU$LabDiagnosis)) +
#  geom_boxplot(outlier.shape = NA, lwd = 1.2) +
#  theme_classic(base_size = 14, base_line_size = 1) +
#  geom_point(position = position_jitterdodge(jitter.width = .1)) +
#  labs(title =  expression(paste(italic("F. nucleatum"), " abundance"))
#       , x = "Sample Type", y = expression(paste("Relative abundance of ", italic("F. nucleatum"))), fill = 'Polyp pathology') 

#No significant differences in fuso relative ab.
#kruskal.test(fuso_OTU$`colSums(fuso_OTU)` ~ fuso_OTU$LabDiagnosis)

#ggplot(data = fuso_OTU) +
#  aes(x = as.character(fuso_OTU$SampleSite), y = fuso_OTU$`colSums(fuso_OTU)`, fill = as.character(fuso_OTU$SampleSite)) +
#  geom_boxplot(outlier.shape = NA) +
#  geom_jitter(width = .1) +
#  labs(title = 'F. nucleatum abundance', subtitle = "No significant differences",
#       x = 'Sample Site', y = 'Number of F. nucleatum hits', fill = 'Tissue Site')

###### E lenta ############
taxofint <- as.data.frame(rared_OTU[, grep('lenta', colnames(rared_OTU))])
row.names(taxofint) <- row.names(rared_OTU)
taxofint2 <- as.data.frame(rowSums(taxofint))
taxofint3 <- taxofint2/rowSums(rared_OTU)
taxofint4 <- merge(taxofint3, metadata, by = "row.names")
taxofint4 <- taxofint4[!grepl("unknown", taxofint4$LabDiagnosis),]

ggplot(data = taxofint4[taxofint4$SampleType == "aspirate",]) +
  aes(x = as.character(LabDiagnosis), y = `rowSums(taxofint)`+0.0001, fill = as.character(LabDiagnosis)) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_jitter(size = 2, alpha = .9) +
  labs(x = NULL, y = 'Relative abundance', fill = 'Subject type', title = "Eggerthella lenta") +
  #geom_point(position = position_jitterdodge(jitter.width = .1), alpha = 0.3) +
  scale_fill_manual(values=c("forestgreen", "firebrick3", "steelblue3")) +
  theme(plot.title = element_text(face = "italic"), legend.position = "none") +
  scale_x_discrete(labels = c("Polyp free (17)", "TA (24)", "Serrated (11)")) +
  scale_y_log10()

#ggsave("16S_elenta.png", device = "png", dpi = 300, width = 4, height = 3)

###### Taxa barplots ####
#Turn to relative abundance
relab_table <- rared_OTU/rowSums(rared_OTU)
  
#Melt and separate taxonomy, remove whitespaces
barplot_df <- reshape2::melt(as.matrix(relab_table)) %>% mutate(Var2 = gsub(".__", "", .$Var2)) %>% 
  tidyr::separate(., col = Var2, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop") %>% mutate_all(na_if,"")
barplot_df$value <- as.numeric(barplot_df$value)
  
#Taxonomy level of interest
barplot_df <- rename(barplot_df, taxonomy = L5)
  
#Take top 10
top_taxa <- group_by(barplot_df, taxonomy) %>% summarise(., top_taxa_tmp = sum(value)) %>% arrange(., desc(top_taxa_tmp)) %>% slice(., 1:11)
high_abundance <- split(top_taxa$taxonomy, 1:NROW(top_taxa)) 
high_abundance <- high_abundance[!is.na(high_abundance)]

#Replace not top 10 with other. 
barplot_df$taxonomy[barplot_df$taxonomy %in% high_abundance != "TRUE"] <- "Other"
barplot_df2 <- aggregate(barplot_df$value, by=list(taxonomy=barplot_df$taxonomy, Var1 = barplot_df$Var1), FUN=sum) %>% merge(., metadata, by.x = "Var1", by.y = "row.names")

barplot_df2 <- barplot_df2[order(barplot_df2$taxonomy),] #Re order
barplot_df2 <- rbind(barplot_df2[!(barplot_df2$taxonomy == "Other"),],barplot_df2[(barplot_df2$taxonomy == "Other"),]) #Move other to bottom
barplot_df2$taxonomy <- factor(barplot_df2$taxonomy, levels = unique(barplot_df2$taxonomy)) #Fix the order

#Remove metadata that I dont want to plot
barplot_df2 <- barplot_df2[!(barplot_df2$LabDiagnosis %in% c("unknown", "adenocarcinoma")),]

#Custom color pallette.
sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947" ,"gray")

ggplot(data = barplot_df2, aes(x = Var1, weight = x, fill = taxonomy)) +
  geom_bar(width = 1, color = "black", size = .2) +
  theme_classic(base_size = 16) +
  facet_grid(.~SampleType+LabDiagnosis, space = "free", scales = "free", labeller = 
               labeller(SampleType = c(`aspirate` = "Aspirate", `brush` = "Brush", `lavage` = "Lavage"), 
                        LabDiagnosis = c(`healthy`="Healthy", `non-serrated`="Tubular Adenoma", `serrated`="Serrated"))) +
  scale_fill_manual(values = sarah_color) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 8), 
        strip.background = element_rect(fill="lightblue")) +
  labs(x = NULL, y = "Relative abundance", fill = "Family", title = bquote("Individuals:"~.(length(unique(barplot_df2$Patient))))) 

#group_by(barplot_df2, Var1) %>% filter(., MostMalignantPolypType == "serrated") %>% count(.) #Used to count sample numbers
ggplot(data = barplot_df2) +
  aes(y = x, x = taxonomy, fill = as.factor(LabDiagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  theme_bw(base_size = 14, base_line_size = 1) +
  scale_fill_manual(name = "Subject type", values=c("forestgreen", "firebrick3", "steelblue3"), labels = c("Polyp free (44)", "Tubular adenoma (55)", "Serrated polyp (20)")) +
  geom_jitter(position = position_jitterdodge(), size = .5) +
  facet_wrap(~SampleType, labeller = labeller(SampleType = c(`aspirate` = "Aspirate (51)", `brush` = "Brush (37)", `lavage` = "Lavage (31)"))) +
  labs(y = "Relative Abundance", x = "Family", subtitle = bquote("Individuals:"~.(length(unique(barplot_df2$Patient))))) +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="lightblue")) +
  scale_x_discrete(limits = rev(levels(barplot_df2$taxonomy)))

###### ANCOM ####

source("/media/julio/Storage/Software/ANCOM-master/programs/ancom.R")

ancom_otu <- merged_OTU_table[!merged_OTU_table$SampleType == "lavage",]
ancom_otu <- as.matrix(t(ancom_otu[,-(1:(1+ncol(metadata)))]))
  
ancom_meta <- metadata
ancom_meta$Patient <- as.factor(ancom_meta$Patient)
ancom_meta$Row.names <- rownames(metadata)
rownames(ancom_meta) <- NULL
  
ancom_pre <- feature_table_pre_process(feature_table = ancom_otu, meta_data = ancom_meta, sample_var = "Row.names",
                                         group_var = "SampleType", lib_cut = 0, neg_lb = T)
  
ancom_res <- ANCOM(feature_table = ancom_pre$feature_table, meta_data = ancom_pre$meta_data, struc_zero = ancom_pre$structure_zeros, 
                     main_var = "SampleType", rand_formula = "~ 1 | Patient", p_adj_method = "fdr")
  
#Riketsialles is mitochondria according to BLAST
ancom_out <- ancom_res$out %>% filter(!(W == "Inf"), detected_0.7 == TRUE) %>% arrange(desc(W))
View(ancom_out)

#Removing rarefaction did not change the number of differentially abundant taxa.

difab_df = merged_OTU_table %>% filter(!SampleType == "lavage")
difab_df[,(ncol(metadata)+2):ncol(difab_df)] = decostand(difab_df[,(ncol(metadata)+2):ncol(difab_df)], method = "total", MARGIN = 1)

difab1 <- ggplot(data = difab_df) +
  aes(x = SampleType, y = `k__Bacteria; p__Firmicutes; c__Bacilli; o__Gemellales; f__Gemellaceae` + 0.0001) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = .5, size = 1) +
  scale_x_discrete(label = c("Aspirate (52)", "Brush (64)")) +
  labs(y = "Relative abundance", title = "f__Gemellaceae", x = NULL) +
  scale_y_log10() +
  theme_bw()

#Determined to be mitochondria by BLAST.
difab2 <- ggplot(data = difab_df) +
  aes(x = SampleType, y = `k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales`+0.0001) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = .5, size = 1) +
  labs(y = "Relative abundance", title = "o__Rickettsiales", x = NULL) +
  scale_y_log10() +
  scale_x_discrete(label = c("Aspirate (52)", "Brush (64)")) +
  theme_bw()

difab3 <- ggplot(data = difab_df) +
  aes(x = SampleType, y = `k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__.17`+0.0001) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(label = c("Aspirate (52)", "Brush (64)")) +
  geom_jitter(alpha = .5, size = 1) +
  labs(y = "Relative abundance", title = "g__Streptococcus; s__.17", x = NULL) +
  scale_y_log10() +
  theme_bw()

difab4 <- ggplot(data = difab_df) +
  aes(x = SampleType, y = `k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__` + 0.0001) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = .5, size = 1) +
  labs(y = "Relative abundance", title = "g__Streptococcus; s__", x = NULL) +
  scale_x_discrete(label = c("Aspirate (52)", "Brush (64)")) +
  scale_y_log10() +
  theme_bw()

difab <- plot_grid(difab1, difab3, difab4, ncol = 3, nrow = 1)

#ggsave("brush_v_aspirates.png", difab, device = "png", dpi = 300, width = 9, height = 3)