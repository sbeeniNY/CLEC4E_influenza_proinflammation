library(dplyr)
library(tidyr)
library(ggplot2)
library(DESeq2) 
library(ggrepel)

output_dir <- 'figures/'

# Fig5b : weight loss
surv <- read.table('data/survival.txt', header = TRUE, sep='\t')

df_long <- surv %>%
  group_by(Sample) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  pivot_longer(
    cols = c("D.1", "D0", "D1", "D2", "D3", "D4", "D5"),
    names_to = "day",
    values_to = "value"
  ) %>%
  select(Sample, day, replicate, value)

df_long <- df_long[df_long$Sample %in% c('PBS+PBS','PBS+IVR180','Piceatannol+IVR180','Piceatannol+PBS'),]
df_long$Sample <- factor(df_long$Sample, levels = c('PBS+PBS','Piceatannol+PBS','PBS+IVR180','Piceatannol+IVR180'))
df_long <- df_long[df_long$day != "D.1",]
df_long$day <- factor(df_long$day, levels = c("D0","D1","D2","D3","D4","D5"))

df_median <- df_long %>%
  group_by(Sample, day) %>%
  summarize(median_value = median(value), 
            se_value = sd(value) / sqrt(n()),  
            .groups = "drop")

ggplot(df_median, aes(x = day, y = median_value, color = Sample, group = Sample)) +
  geom_errorbar(aes(ymin = median_value - se_value, ymax = median_value + se_value), 
                width = 0.2, position = position_dodge(width = 0)) +  
  geom_point(size = 2, position = position_dodge(width = 0)) +  
  geom_line(data = df_median, aes(x = day, y = median_value, color = Sample, group = Sample),
            position = position_dodge(width = 0), linewidth = 1) +
  theme_classic() +
  labs(x = "Day", y = "Value") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("PBS+PBS" = "black", 
                                 "PBS+IVR180" = "#f37f88", 
                                 "Piceatannol+IVR180" = "#8bd1f4", 
                                 "Piceatannol+PBS" = "#aaaaab")) +  
  theme(legend.position = "top")  
ggsave(paste0(output_dir, "Fig5b.pdf"), width = 6, height = 4)

# Fig5c : viral load

df_long <- surv %>%
  group_by(Sample) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  pivot_longer(
    cols = c("plaqueAssayTiter"),
    names_to = "viral_load",
    values_to = "value"
  ) %>%
  select(Sample, viral_load, replicate, value)

df_long <- df_long[df_long$Sample %in% c('PBS+PBS','PBS+IVR180','Piceatannol+IVR180','Piceatannol+PBS'),]
df_long$Sample <- factor(df_long$Sample, levels = c('PBS+PBS','Piceatannol+PBS','PBS+IVR180','Piceatannol+IVR180'))

ggplot(df_long, aes(x = Sample, y = value, color = Sample, group = Sample)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Viral Load", y = "Value") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("PBS+PBS" = "black", 
                                 "PBS+IVR180" = "#f37f88", 
                                 "Piceatannol+IVR180" = "#8bd1f4", 
                                 "Piceatannol+PBS" = "#aaaaab")) + 
  theme(legend.position = "top")  
ggsave(paste0(output_dir, "Fig5c.pdf"), width = 4, height = 4)

# Fig5d : Fold change - Fold change

bulk_counts <- read.table("data/count_w_geneName.tsv", sep = "\t", header = TRUE, row.names = 1)
sample_metadata <- data.frame(SampleID = colnames(bulk_counts), Group = gsub("S([1-6]).*", "S\\1", colnames(bulk_counts))  )

sample_metadata[sample_metadata$Group=="S1","Treat"] <- "PBS"
sample_metadata[sample_metadata$Group=="S2","Treat"] <- "Piceatannol"
sample_metadata[sample_metadata$Group=="S3","Treat"] <- "TriacsinC"
sample_metadata[sample_metadata$Group=="S4","Treat"] <- "PBS"
sample_metadata[sample_metadata$Group=="S5","Treat"] <- "Piceatannol"
sample_metadata[sample_metadata$Group=="S6","Treat"] <- "TriacsinC"
sample_metadata$Treat <- factor(sample_metadata$Treat)

sample_metadata[sample_metadata$Group=="S1","Infection"] <- "IVR180"
sample_metadata[sample_metadata$Group=="S2","Infection"] <- "IVR180"
sample_metadata[sample_metadata$Group=="S3","Infection"] <- "IVR180"
sample_metadata[sample_metadata$Group=="S4","Infection"] <- "PBS"
sample_metadata[sample_metadata$Group=="S5","Infection"] <- "PBS"
sample_metadata[sample_metadata$Group=="S6","Infection"] <- "PBS"
sample_metadata$Infection <- factor(sample_metadata$Infection, levels = c("PBS", "IVR180"))

bulk_counts <- round(bulk_counts)
# low expression filter
keep_genes <- rowSums(bulk_counts >= 10) >= (ncol(bulk_counts) / 2)
bulk_counts_filtered <- bulk_counts[keep_genes, ] 

# high expression filter
gene_mean <- rowMeans(bulk_counts_filtered)
log_gene_mean <- log10(gene_mean + 1)
z_scores <- (log_gene_mean - mean(log_gene_mean)) / sd(log_gene_mean)
outlier_genes <- names(z_scores[z_scores > 5])
bulk_counts_filtered <- bulk_counts_filtered[!rownames(bulk_counts_filtered) %in% outlier_genes, ] 

#1. Infected untreated vs uninfected untreated → effects of infection
sample_metadata_filtered <- sample_metadata %>% filter(Treat == "PBS")  
bulk_counts_pbs <- bulk_counts_filtered[, sample_metadata_filtered$SampleID]

dds <- DESeqDataSetFromMatrix(
  countData = bulk_counts_pbs,
  colData = sample_metadata_filtered,
  design = ~ Infection 
)
dds <- DESeq(dds)
res_comp1 <- results(dds, contrast = c("Infection", "IVR180", "PBS")) # IVR180 - PBS

#2. Infected treated vs infected untreated→ effects of treatment
sample_metadata_filtered <- sample_metadata %>% filter(Infection == "IVR180")
bulk_counts_ivr180 <- bulk_counts_filtered[, sample_metadata_filtered$SampleID]

dds <- DESeqDataSetFromMatrix(
  countData = bulk_counts_ivr180,
  colData = sample_metadata_filtered,
  design = ~ Treat 
)
dds <- DESeq(dds)
res_comp2 <- results(dds, contrast = c("Treat", "Piceatannol", "PBS")) # Piceatannol - PBS

#3. Check enrichment of signature II in signature I.
comp1_deg_filtered <- as.data.frame(res_comp1) %>% filter(pvalue < 0.05, abs(log2FoldChange) > 0.5)
comp2_deg_filtered <- as.data.frame(res_comp2) %>% filter(pvalue < 0.05, abs(log2FoldChange) > 0.5)

comp1_up <- comp1_deg_filtered %>% filter(log2FoldChange > 0) %>% row.names()
comp1_down <- comp1_deg_filtered %>% filter(log2FoldChange < 0) %>% row.names()
comp2_up <- comp2_deg_filtered %>% filter(log2FoldChange > 0) %>% row.names()
comp2_down <- comp2_deg_filtered %>% filter(log2FoldChange < 0) %>% row.names()
comp1_up_comp2_dn <- intersect(comp1_up, comp2_down)
comp1_down_comp2_up <- intersect(comp1_down, comp2_up)
comp1_up_only <- setdiff(comp1_up,comp1_up_comp2_dn)
comp2_down_only <- setdiff(comp2_down, comp1_up_comp2_dn)
comp1_down_only <- setdiff(comp1_down, comp1_down_comp2_up)
comp2_up_only <- setdiff(comp2_up, comp1_down_comp2_up)

res <- merge(data.frame(res_comp2), data.frame(res_comp1), by = "row.names", suffixes = c("_res2", "_res1"))
rownames(res) <- res$Row.names
res$Row.names <- NULL 

res$label <- ifelse(rownames(res) %in% comp1_up_comp2_dn, "comp1_up_comp2_dn",
              ifelse(rownames(res) %in% comp1_down_comp2_up, "comp1_down_comp2_up",
                ifelse(rownames(res) %in% comp1_up_only, "comp1_up_only",
                  ifelse(rownames(res) %in% comp2_down_only, "comp2_down_only",
                    ifelse(rownames(res) %in% comp1_down_only, "comp1_down_only",
                      ifelse(rownames(res) %in% comp2_up_only, "comp2_up_only", "not significant"))))))

red_points <- res[c('Ifi30','Itga5','Mafb','Marco'),]

ggplot(res, aes(x = log2FoldChange_res1, y = log2FoldChange_res2)) +
  geom_point(aes(color = label, alpha = 0.7)) + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +  
  theme_classic() +
  labs(y = "log2FoldChange_infected treated vs infected untreated", x = "log2FoldChange_infected untreated vs uninfected untreated") +
  scale_color_manual(values = c("comp1_up_comp2_dn" = "#f98089", "comp1_down_comp2_up" = "#89d2f7", 
                                "comp1_up_only" = "#ffcd87", "comp2_down_only" = "#b7e4a5", 
                                "comp1_down_only" = "#c193b9", "comp2_up_only" = "#f9b7c5", 
                                "not significant" = "grey")) +
  theme(legend.position = "bottom") +
  geom_text_repel(data = red_points, aes(label = rownames(red_points)), 
                   color = "black", size = 3, 
                   segment.color = 'black')  
ggsave(paste0(output_dir, "Fig5d.pdf"), width = 10, height = 8.5)

