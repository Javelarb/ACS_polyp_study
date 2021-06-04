library(reshape2)
library(tidyverse)
library(EcolUtils)
library(vegan)
library(nlme)
library(labdsv)
library(RColorBrewer)
library(rfPermute)
library(cowplot)

setwd("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/R_analysis")

ITS_OTU_table <- read.delim("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/qiime2/ITS/dada2/OTU_table", row.names=1, check.names = FALSE)

#Issue with taxonomy, where the finest resolution is fungi. Unite database.
ITS_OTU_taxonomy <- read.delim("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/qiime2/ITS/3ae24d7d-e98d-4893-a28d-1663d80dc47d/data/taxonomy.tsv", row.names=1)
ITS_metadata <- read.delim("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/qiime2/ITS_metadata2.tsv", row.names=1, check.names = FALSE)

#### ITS read counts and filtering ####
ITS_read_count <- read.delim("/media/julio/Storage/CRC/Winter_2019_amplicon_CC_data/qiime2/ITS/dada2/7d41928f-e4eb-422d-a057-f423b907314d/data/metadata.tsv", row.names=1, comment.char="#")
ITS_read_count <- ITS_read_count[,-(2:4)] #Remove filtered, denoised, and merged columns
melted_read_counts <- reshape2::melt(as.matrix(ITS_read_count))
read_counts_plus_metadata <- merge(melted_read_counts, ITS_metadata, by.x = "Var1", by.y = "row.names")
row.names(read_counts_plus_metadata) <- read_counts_plus_metadata$Row.names

#Read count by sample type:
ggplot(data = read_counts_plus_metadata) +
  aes(x = SampleType, y = value, fill = Var2) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_family = "", base_line_size = 1) +
  scale_y_continuous(breaks = seq(0,200000,25000)) +
  labs(x = NULL, y = 'Read count', fill = 'Denoising step') +
  geom_point(position = position_jitterdodge(jitter.width = .1)) +
  scale_y_log10()

#Test if there are sig dif between sample types, no repeated measurements taken into account but it is not significant regardless.
TukeyHSD(aov(formula = read_counts_plus_metadata$value ~ as.character(read_counts_plus_metadata$SampleType)))

#Filtering
ITS_OTU_no_mock <- ITS_OTU_table[,!(names(ITS_OTU_table) %in% c("227-ITS","363-ITS","406-ITS", "591-5-ITS"))] #Remove samples with Mock spike in
ITS_OTU_table_plus_taxonomy <- as.data.frame(merge(ITS_OTU_taxonomy, ITS_OTU_no_mock, by = "row.names"))
ITS_OTU_table_no_unassigned <- ITS_OTU_table_plus_taxonomy[!grepl("Unassigned", ITS_OTU_table_plus_taxonomy$`Taxon`),]
row.names(ITS_OTU_table_no_unassigned) <- ITS_OTU_table_no_unassigned$Row.names
ITS_OTU_clean <- as.data.frame(t(ITS_OTU_table_no_unassigned[,4:ncol(ITS_OTU_table_no_unassigned)]))

sort(rowSums(ITS_OTU_clean))
rd_ITS <- 1000

ITS_rared_OTU <- as.data.frame((rrarefy.perm(ITS_OTU_clean, sample = rd_ITS, n = 10, round.out = T)))
ITS_rared_OTU <- as.data.frame(ITS_rared_OTU[rowSums(ITS_rared_OTU) >= rd_ITS-(rd_ITS*.1), colSums(ITS_rared_OTU) >= 1])
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
  geom_point(position = position_jitterdodge(jitter.width = .8), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Healthy (33)", "Tubular adenoma (43)", "Serrated polyp (28)")) +
  scale_x_discrete(labels = c("Aspirate (42)", "Brush (46)", "Lavage (16)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))
ITS_shannon

ITS_shannon_polyp <- ggplot(data = alpha_div) +
  aes(x = SampleType, y = alpha_div$`diversity(ITS_rared_OTU, index = "shannon")`, fill = LabDiagnosis) +
  geom_boxplot() +
  theme_classic() +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = .5) +
  scale_fill_manual(values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Healthy (33)", "Tubular adenoma (43)", "Serrated polyp (28)")) +
  labs(x = NULL, y = "Shannon index", title = "Sample set 1 - ITS data", fill = "Subject type") +
  scale_x_discrete(labels = c("Aspirate (42)", "Brush (46)", "Lavage (16)")) + 
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
  aes(x = as.character(SampleType), y = specno$`specnumber(ITS_rared_OTU)`, fill = SampleType) +
    geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = NULL, y = 'ASV richness', fill = 'Subject type', title = " ") +
  geom_point(position = position_jitterdodge(jitter.width = .8), alpha = .5) +
  scale_fill_manual(name = "Subject type", values=c("gold", "lightsalmon", "plum2"), labels = c("Healthy (33)", "Tubular adenoma (43)", "Serrated polyp (28)")) +
  scale_x_discrete(labels = c("Aspirate (42)", "Brush (46)", "Lavage (16)")) +
  theme(legend.position = "none", axis.line = element_line(color = "black"), panel.background = element_blank(), 
        axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12), title = element_text(size = 12))
ITS_rich

ITS_rich_polyp = ggplot(data = specno) +
  aes(x = SampleType, y = specno$`specnumber(ITS_rared_OTU)`, fill = LabDiagnosis) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = .4), alpha = .5) +
  scale_x_discrete(labels = c("Aspirate (42)", "Brush (46)", "Lavage (16)")) +
  theme_classic() +
  scale_fill_manual(values=c("forestgreen", "firebrick3", "steelblue"), labels = c("Healthy (33)", "Tubular adenoma (43)", "Serrated polyp (28)")) +
  labs(x = NULL, y = "Richness", title = " ", fill = "Subject type")
ITS_rich_polyp

ITS_alpha <- plot_grid(ITS_shannon, ITS_rich, rel_widths = c(1, 1))
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
rared_mds <- metaMDS(ITS_rared_OTU, distance = "bray", k = 2, trymax = 999, parallel = 32)
merged_rared_mds <- merge(rared_mds$points[,1:2], ITS_metadata, by = "row.names")

asp_mds <- metaMDS(asp_OTU[,15:ncol(asp_OTU)])
merged_asp_mds <- merge(asp_mds$points[,1:2], ITS_metadata, by = "row.names")

brush_mds <- metaMDS(brush_OTU[,15:ncol(brush_OTU)])
merged_brush_mds <- merge(brush_mds$points[,1:2], ITS_metadata, by = "row.names")
#plot3d(merged_pcoa[,2:4], col = as.integer(merged_pcoa$SampleType), size = 10)

lavage_mds <- metaMDS(lavage_OTU[,15:ncol(lavage_OTU)])
merged_lavage_mds <- merge(lavage_mds$points[,1:2], ITS_metadata, by = "row.names")

#samp_type_perma
merged_OTU_table <- merged_OTU_table[!merged_OTU_table$Ethnicity == "unknown",]
perma1 <- adonis(formula = merged_OTU_table[,14:ncol(merged_OTU_table)] ~ as.numeric(merged_OTU_table$BMI) + 
                   as.numeric(merged_OTU_table$Age) + merged_OTU_table$Gender + merged_OTU_table$Ethnicity +
                   (merged_OTU_table$LabDiagnosis /as.character(merged_OTU_table$Patient) / merged_OTU_table$SampleType), 
                 data = merged_OTU_table[,14:ncol(merged_OTU_table)], method = "bray", permutations = 999, parallel = 32)

#Some left over code to pull out PERMANOVA coefficients.
#coef1 <- coefficients(perma1)["merged_OTU_table$LabDiagnosis2",]
#top.coef1 <- coef1[rev(order(abs(coef1)))[1:10]]
#par(mar=c(2,16,1,1)+.1)
#barplot(sort(top.coef1), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=1)

#asp perma
perma2 <- adonis(formula = asp_OTU[,15:ncol(asp_OTU)] ~ asp_OTU$LabDiagnosis + as.character(asp_OTU$Patient) + asp_OTU$ColonLocation, 
       data = asp_OTU[,15:ncol(asp_OTU)], method = "bray", permutations = 999, parallel = 16)
coef2 <- coefficients(perma2)["asp_OTU$LabDiagnosis1",]
top.coef2 <- coef2[rev(order(abs(coef2)))[1:10]]
barplot(sort(top.coef2), horiz = T, las = 1, main = "Top taxa", col = "deepskyblue3", cex.axis=1, cex.names=1)

#Brush perma
adonis(formula = brush_OTU[,15:ncol(brush_OTU)] ~ brush_OTU$LabDiagnosis / as.character(brush_OTU$Patient)/ as.character(brush_OTU$SampleSite), 
       data = brush_OTU[,15:ncol(brush_OTU)], method = "bray", permutations = 999, parallel = 16)

#Lavage perma
adonis(formula = lavage_OTU[,15:ncol(lavage_OTU)] ~ lavage_OTU$LabDiagnosis / as.character(lavage_OTU$Patient), 
       data = lavage_OTU[,15:ncol(lavage_OTU)], method = "bray", permutations = 999, parallel = 16)

#General plot
ggplot(data = merged_rared_mds) +
  aes(x = MDS1, y = MDS2) +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_point(aes(pch = LabDiagnosis, fill = SampleType), size = 8, alpha = .4) + 
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text(label = merged_rared_mds$Patient, size = 4) +
  theme(plot.title = element_text(size = 15, vjust = -1), plot.subtitle = element_text(size = 15)) +
  scale_shape_manual(values = c(21,22,24), name = "Subject type", labels = c("Healthy (33)", "Tubular adenoma (43)", "Serrated polyp (28)")) +
  scale_fill_manual(name = "Sample type", values=c("gold", "lightsalmon", "plum2"), 
                    labels = c("Aspirate (42)", "Brush (46)", "Lavage (16)"))

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

###### TPL #######
df3 <- rbind(asp_R_OTU, asp_L_OTU)
df3 <- df3 %>% group_by(Patient) %>% filter(n()>1)
taylor_out = NULL
for (i in df3$Patient) {
  assign(paste0("OTU_table_", i), df3[df3$Patient == i,]) #will create lots of variables
  log.m <- log(apply(get(paste0("OTU_table_", i))[,-(1:(1+ncol(ITS_metadata)))], 2, mean))
  log.v <- log(apply(get(paste0("OTU_table_", i))[,-(1:(1+ncol(ITS_metadata)))], 2, var))
  log.m <- log.m[is.finite(log.m)]
  log.v <- log.v[is.finite(log.v)]
  log.m_v <- merge(log.m, log.v, by = "row.names")
  fit <- lm(log.m_v$y ~ log.m_v$x)
  c.value <- as.numeric(coef(fit)[1]) #constant
  z.value <- as.numeric(coef(fit)[2]) #var of interest, z > 1 = aggregated, z = 1 random, z < 1 uniform
  individual <- i
  individual <- cbind(individual, z.value)
  taylor_out <- rbind(taylor_out, individual)
}

taylor_out <- as.data.frame(unique(taylor_out))
taylor_final <- merge(taylor_out, ITS_metadata, by.x = "individual", by.y = "Patient")
taylor_final <- taylor_final[,c(1:2,11)]
taylor_final <- unique(taylor_final)

ggplot(data = taylor_final) +
  aes(x = taylor_final$LabDiagnosis, y = taylor_final$z.value, fill = taylor_final$LabDiagnosis) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  theme_classic(base_size = 14, base_line_size = 1) +
  labs(x = 'Subject status', y = 'Taylor`s exponent' , fill = 'Subject status', 
       title = 'Proximal + distal') +
  geom_point(position = position_jitterdodge(jitter.width = .5)) +
  theme(legend.position = "none")
TukeyHSD(aov(formula = taylor_final$z.value ~ taylor_final$LabDiagnosis))