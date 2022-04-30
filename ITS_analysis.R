library(reshape2)
library(tidyverse)
library(EcolUtils)
library(vegan)
library(nlme)
library(labdsv)
library(RColorBrewer)
library(rfPermute)
library(cowplot)

setwd("/media/julio/Storage/CRC/github/")

#Issue with taxonomy, where the finest resolution is fungi. Unite database.
ITS_OTU_table <- read.delim("ITS_OTU_table.tsv", row.names=1, check.names = FALSE)
ITS_OTU_taxonomy <- read.delim("ITS_taxonomy.tsv", row.names=1)
ITS_metadata <- read.delim("ITS_metadata.tsv", row.names=1, check.names = FALSE)

#### ITS read counts and filtering ####
#ITS_read_count <- read.delim("/qiime2/ITS/dada2/7d41928f-e4eb-422d-a057-f423b907314d/data/metadata.tsv", row.names=1, comment.char="#")
#ITS_read_count <- ITS_read_count[,-(2:4)] #Remove filtered, denoised, and merged columns
#melted_read_counts <- reshape2::melt(as.matrix(ITS_read_count))
#read_counts_plus_metadata <- merge(melted_read_counts, ITS_metadata, by.x = "Var1", by.y = "row.names")
#row.names(read_counts_plus_metadata) <- read_counts_plus_metadata$Row.names

#Read count by sample type:
#ggplot(data = read_counts_plus_metadata) +
#  aes(x = SampleType, y = value, fill = Var2) +
#  geom_boxplot(outlier.shape = NA, lwd = 1) +
#  theme_classic(base_size = 14, base_family = "", base_line_size = 1) +
#  scale_y_continuous(breaks = seq(0,200000,25000)) +
#  labs(x = NULL, y = 'Read count', fill = 'Denoising step') +
#  geom_point(position = position_jitterdodge(jitter.width = .1)) +
#  scale_y_log10()

#Test if there are sig dif between sample types, no repeated measurements taken into account but it is not significant regardless.
#TukeyHSD(aov(formula = read_counts_plus_metadata$value ~ as.character(read_counts_plus_metadata$SampleType)))

#Filtering
ITS_OTU_no_mock <- ITS_OTU_table[,!(names(ITS_OTU_table) %in% c("227-ITS","363-ITS","406-ITS", "591-5-ITS"))] #Remove samples with Mock spike in
ITS_OTU_table_plus_taxonomy <- as.data.frame(merge(ITS_OTU_taxonomy, ITS_OTU_no_mock, by = "row.names"))
ITS_OTU_table_no_unassigned <- ITS_OTU_table_plus_taxonomy[!grepl("Unassigned", ITS_OTU_table_plus_taxonomy$`Taxon`),]
row.names(ITS_OTU_table_no_unassigned) <- ITS_OTU_table_no_unassigned$Row.names
ITS_OTU_clean <- as.data.frame(t(ITS_OTU_table_no_unassigned[,4:ncol(ITS_OTU_table_no_unassigned)]))

rare_curve = rarecurve(ITS_OTU_clean, step = 250, label = F, sample = 1000, xlab = "Read Depth", col = "orange")

Nmax <- sapply(rare_curve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rare_curve, max)

plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Read depth",
     ylab = "ASVs", type = "n", main = "ITS")
for (i in seq_along(rare_curve)) {
  N <- attr(rare_curve[[i]], "Subsample")
  lines(N, rare_curve[[i]], col = "gray")
}
abline(v = 1000, lty = 2, col = "red")

#rd_ITS <- 1000
ITS_rared_OTU = ITS_OTU_clean[!rowSums(ITS_OTU_clean) <= 1000,]

#ITS_rared_OTU <- as.data.frame((rrarefy.perm(ITS_OTU_clean, sample = rd_ITS, n = 10, round.out = T)))
#ITS_rared_OTU <- as.data.frame(ITS_rared_OTU[rowSums(ITS_rared_OTU) >= rd_ITS-(rd_ITS*.1), colSums(ITS_rared_OTU) >= 1])
merged_OTU_table <- merge(ITS_metadata, ITS_rared_OTU, by = "row.names")
row.names(merged_OTU_table) <- merged_OTU_table$Row.names

lavage_OTU <- merged_OTU_table[grepl("lavage", merged_OTU_table$SampleType),]
asp_OTU <- merged_OTU_table[grepl("aspirate", merged_OTU_table$SampleType),]
brush_OTU <- merged_OTU_table[grepl("brush", merged_OTU_table$SampleType),]
asp_L_OTU <- asp_OTU[grepl("left", asp_OTU$ColonLocation),]
asp_R_OTU <- asp_OTU[grepl("right", asp_OTU$ColonLocation),]

#### alpha diversity ####
alpha_div <- as.data.frame(diversity(ITS_rared_OTU, index = "shannon"))
alpha_div <- merge(alpha_div, ITS_metadata, by = "row.names")
alpha_div <- alpha_div[!grepl("unknown", alpha_div$LabDiagnosis),]

ITS_shannon <- ggplot(data = alpha_div) +
  aes(x = as.character(SampleType), y = `diversity(ITS_rared_OTU, index = "shannon")`, fill = SampleType) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Shannon index', fill = 'Subject type', title = bquote("ITS - Individuals:"~.(length(unique(alpha_div$Patient))))) +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = 1/2) +
  scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Polyp free (32)", "Tubular adenoma (41)", "Serrated polyp (25)")) +
  scale_x_discrete(labels = c("Aspirate (41)", "Brush (41)", "Lavage (16)")) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
ITS_shannon

ITS_shannon_polyp <- ggplot(data = alpha_div) +
  aes(x = SampleType, y = `diversity(ITS_rared_OTU, index = "shannon")`, fill = LabDiagnosis) +
  geom_boxplot() +
  theme_classic() +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = .5) +
  scale_fill_manual(values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Polyp free (32)", "Tubular adenoma (41)", "Serrated polyp (25)")) +
  labs(x = NULL, y = "Shannon index", title = "Sample set 1 - ITS data", fill = "Subject type") +
  scale_x_discrete(labels = c("Aspirate (41)", "Brush (41)", "Lavage (16)")) + 
  theme(legend.position = "none")
ITS_shannon_polyp

alpha_lm <- NULL
alpha_lm$x <- as.numeric(alpha_div$`diversity(ITS_rared_OTU, index = "shannon")`)
alpha_lm$y <- as.factor(alpha_div$SampleType)
alpha_lm$z <- as.factor(alpha_div$LabDiagnosis)
alpha_lm$i <- as.factor(alpha_div$Patient)
alpha_lm <- as.data.frame(alpha_lm)
#alpha_lm <- within(alpha_lm, z <- relevel(z, "healthy"))
summary(lme(x ~ z * y, data = alpha_lm, random = list(i=~1)))
#Looks like tubular adenomas have significantly less shannon fungus diversity compared to healthy.

#species counts instead of shannon
specno <- as.data.frame(specnumber(ITS_rared_OTU))
specno <- merge(specno, ITS_metadata, by = "row.names")
specno <- specno[!grepl("unknown", specno$LabDiagnosis),]

ITS_rich <- ggplot(data = specno) +
  aes(x = as.character(SampleType), y = `specnumber(ITS_rared_OTU)`, fill = SampleType) +
    geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'Richness', fill = 'Subject type') +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = 1/2) +
  scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Polyp free (32)", "Tubular adenoma (41)", "Serrated polyp (25)")) +
  scale_x_discrete(labels = c("Aspirate (41)", "Brush (41)", "Lavage (16)")) +
  theme(legend.position = "none", axis.text = element_text(size = 10))
ITS_rich

ITS_rich_polyp = ggplot(data = specno) +
  aes(x = SampleType, y = `specnumber(ITS_rared_OTU)`, fill = LabDiagnosis) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = .5) +
  scale_x_discrete(labels = c("Aspirate (41)", "Brush (41)", "Lavage (16)")) +
  theme_classic() +
  scale_fill_manual(values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Polyp free (32)", "Tubular adenoma (41)", "Serrated polyp (25)")) +
  labs(x = NULL, y = "Richness", title = " ", fill = "Subject type")
ITS_rich_polyp

ITS_alpha <- plot_grid(ITS_shannon, ITS_rich, ncol = 1, rel_widths = c(1,1))
ITS_alpha_polyp <- plot_grid(ITS_shannon_polyp, ITS_rich_polyp, rel_widths = c(1, 1.55))

ggsave("ITS_alpha.png", plot = ITS_alpha, device = "png", units = "in", dpi = 300, height = 5, width = 8)
ggsave("ITS_alpha_polyp.png", plot = ITS_alpha_polyp, device = "png", units = "in", dpi = 300, height = 4, width = 10)

specno_lm <- NULL
specno_lm$x <- as.numeric(specno$`specnumber(ITS_rared_OTU)`)
specno_lm$y <- as.factor(specno$SampleType)
specno_lm$z <- as.factor(specno$LabDiagnosis)
specno_lm$p <- as.factor(specno$Patient)
specno_lm <- as.data.frame(specno_lm)
#specno_lm <- within(specno_lm, z <- relevel(z, "healthy"))
summary(lme(x ~ y * z, data = specno_lm, random = ~1|p))

#Test sig dif in eveness without plotting.
specno$evenness <- alpha_div$`diversity(ITS_rared_OTU, index = "shannon")`/log(specno$`specnumber(ITS_rared_OTU)`)
specno_lm$e <- as.numeric(specno$evenness)
summary(lme(e ~ y * z, data = specno_lm, random = ~1|p))

#### Beta div ####
rared_mds <- metaMDS(ITS_rared_OTU, distance = "bray", k = 3, trymax = 999, parallel = 32)
merged_rared_mds <- merge(rared_mds$points[,1:2], ITS_metadata, by = "row.names")

asp_mds <- metaMDS(asp_OTU[,(ncol(ITS_metadata)+2):ncol(asp_OTU)])
merged_asp_mds <- merge(asp_mds$points[,1:2], ITS_metadata, by = "row.names")

brush_mds <- metaMDS(brush_OTU[,(ncol(ITS_metadata)+2):ncol(brush_OTU)])
merged_brush_mds <- merge(brush_mds$points[,1:2], ITS_metadata, by = "row.names")

lavage_mds <- metaMDS(lavage_OTU[,(ncol(ITS_metadata)+2):ncol(lavage_OTU)])
merged_lavage_mds <- merge(lavage_mds$points[,1:2], ITS_metadata, by = "row.names")

#samp_type_perma
merged_OTU_table2 <- merged_OTU_table[!merged_OTU_table$Ethnicity == "unknown",]
perma1 <- adonis(formula = merged_OTU_table2[,(ncol(ITS_metadata)+2):ncol(merged_OTU_table2)] ~ as.numeric(merged_OTU_table2$BMI) + 
                   as.numeric(merged_OTU_table2$Age) + merged_OTU_table2$Gender + merged_OTU_table2$Ethnicity +
                   (merged_OTU_table2$LabDiagnosis /as.character(merged_OTU_table2$Patient) / merged_OTU_table2$SampleType), 
                 data = merged_OTU_table2, method = "bray", permutations = 999, parallel = 32)
perma1

#asp perma
perma2 <- adonis(formula = asp_OTU[,(ncol(ITS_metadata)+2):ncol(asp_OTU)] ~ asp_OTU$LabDiagnosis + as.character(asp_OTU$Patient) + asp_OTU$ColonLocation, 
       data = asp_OTU, method = "bray", permutations = 999, parallel = 16)

#Brush perma
adonis(formula = brush_OTU[,(ncol(ITS_metadata)+2):ncol(brush_OTU)] ~ brush_OTU$LabDiagnosis / as.character(brush_OTU$Patient)/ as.character(brush_OTU$SampleSite), 
       data = brush_OTU, method = "bray", permutations = 999, parallel = 16)

#Lavage perma
adonis(formula = lavage_OTU[,(ncol(ITS_metadata)+2):ncol(lavage_OTU)] ~ lavage_OTU$LabDiagnosis / as.character(lavage_OTU$Patient), 
       data = lavage_OTU, method = "bray", permutations = 999, parallel = 16)

merged_rared_mds = merged_rared_mds %>% mutate(LabDiagnosis = gsub("non-serrated", "TA", .$LabDiagnosis))

#General plot
ITS_beta = ggplot(data = merged_rared_mds) +
  aes(x = MDS1, y = MDS2) +
  theme_bw() +
  geom_point(aes(pch = LabDiagnosis, fill = SampleType), size = 5, alpha = .5) + 
  labs(x= "MDS1", y = "MDS2") +
  scale_shape_manual(values = c(21,22,24), name = "Subject type", 
                     labels = c("Polyp free (32)", "Serrated polyp (25)", "Tubular adenoma (41)"),
                     limits = c("healthy", "serrated", "TA")) +
  #stat_ellipse(linetype = 2, aes(group = LabDiagnosis), size = 1) +
  guides(color = "none", fill = guide_legend(order = 1, override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Sample type", values=c("gold", "lightsalmon", "plum2"), 
                    labels = c("Aspirate (41)", "Brush (41)", "Lavage (16)")) +
  geom_text(label = merged_rared_mds$Patient, size = 2, color = "black") +
  annotate("text", x = .25, y = -.35, hjust = 1, label = bquote("Stress ="~.(round(rared_mds$stress, digits = 2))~.(", k = 3")))
ITS_beta

#ggsave("Figure_2ITS.svg", plot = plot_grid(ITS_alpha, ITS_beta, nrow = 1, rel_widths = c(.35,.65), labels = c("C.", "D.")), 
#       device = "svg", units = "in", dpi = 300, height = 5, width = 11)

#aspirate only
ggplot(data = merged_asp_mds) +
  aes(x = MDS1, y = MDS2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = ColonLocation, fill = LabDiagnosis), size = 8) + 
  scale_shape_manual(values = c(23,22), name = "Colon side") +
  stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_asp_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#Brush only
ggplot(data = merged_brush_mds) +
  aes(x = MDS1, y = MDS2, fill = SampleSite) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 8) +
  #stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Hyperlocal Site") +
  geom_text(label = merged_brush_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))

#Lavages only
ggplot(data = merged_lavage_mds) +
  aes(x = MDS1, y = MDS2, fill = LabDiagnosis) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(pch = 21, size = 8) +
  ggtitle(bquote("Patient:"~R^2~"= 0.78, p = 0.90"), subtitle = bquote("Subject type:" ~R^2 ~ "= 0.13, p = 0.83")) +
  stat_ellipse(linetype = 2, aes(color = LabDiagnosis), size = 1) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_discrete(name = "Subject Type") +
  geom_text(label = merged_lavage_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15))