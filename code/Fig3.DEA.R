library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

output_dir <- 'figures/'

# Fig3a: Volcano plot of Monocytes from DEA (Differential Expression Analysis)

Mono1_DEGs <- read.table('data/Clusterwise_DEG_DESeq2/clusterwise_DESeq2.Mono.1.txt',sep="\t",header=T)
Mono2_DEGs <- read.table('data/Clusterwise_DEG_DESeq2/clusterwise_DESeq2.Mono.2.txt',sep="\t",header=T)
Mono3_DEGs <- read.table('data/Clusterwise_DEG_DESeq2/clusterwise_DESeq2.Mono.3.txt',sep="\t",header=T)

Mono1_Influenza <- Mono1_DEGs[Mono1_DEGs$contrast == "Influenza vs control",] 
Mono1_COVID <- Mono1_DEGs[Mono1_DEGs$contrast == "COVID-19 vs control",] 

Mono2_Influenza <- Mono2_DEGs[Mono2_DEGs$contrast == "Influenza vs control",] 
Mono2_COVID <- Mono2_DEGs[Mono2_DEGs$contrast == "COVID-19 vs control",] 

Mono3_Influenza <- Mono3_DEGs[Mono3_DEGs$contrast == "Influenza vs control",] 
Mono3_COVID <- Mono3_DEGs[Mono3_DEGs$contrast == "COVID-19 vs control",] 

cell_type_map <- c("Mono1" = "Classical_Monocyte1", "Mono2" = "Classical_Monocyte2", "Mono3" = "Non-classical_Monocyte")

padj_threshold <- 0.05
logfc_threshold <- log2(1.5)

for (cell_prefix in names(cell_type_map)) {
    data_flu <- get(paste0(cell_prefix, "_Influenza"))
    data_covid <- get(paste0(cell_prefix, "_COVID"))
    data_flu$comparison <- "Influenza\nvs control"
    data_covid$comparison <- "COVID-19\nvs control"
    combined_data <- rbind(data_flu, data_covid)
    combined_data$comparison <- factor(combined_data$comparison, levels = c("Influenza\nvs control", "COVID-19\nvs control"))

    processed_data <- combined_data %>%
        filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
        mutate(
            neg_log10_padj = -log10(padj),
            significance = case_when(log2FoldChange > logfc_threshold & padj < padj_threshold  ~ "Up-regulated",
                                    log2FoldChange < -logfc_threshold & padj < padj_threshold ~ "Down-regulated",
                                    TRUE ~ "Not significant"))

    top_genes <- processed_data %>% filter(significance != "Not significant") %>%
                group_by(comparison) %>%
                arrange(padj, desc(abs(log2FoldChange))) %>%
                slice_head(n = 10) %>% ungroup()

    combined_plot <- ggplot(processed_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
        geom_point(aes(color = significance), alpha = 0.7, size = 3) +
        scale_color_manual(values = c("Up-regulated" = "#E64B35",
                                    "Down-regulated" = "#3C5488",
                                    "Not significant" = "grey")) +
        geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed") +
        geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
        geom_text_repel(data = top_genes,
                        aes(label = gene.name),
                        size = 4,
                        box.padding = 0.5,
                        point.padding = 0.5,
                        segment.color = 'grey50',
                        max.overlaps = Inf) +
        facet_wrap(~ comparison, ncol = 1, strip.position = "right", scales = "free_y") +
        labs(title = cell_type_map[cell_prefix],
            x = "log2(Fold Change)",
            y = "-log10(Adjusted P-value)") +
        theme_minimal(base_size = 14) +
        theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey90", linetype = "blank")
        )
  ggsave(paste0(output_dir, "Fig3a.", cell_prefix, ".pdf"), plot = combined_plot, width = 6, height = 8)
}


# Fig3d: Volcano plot of mDC DEGs : influenza only

DEGs <- read.table('data/Clusterwise_DEG_DESeq2/clusterwise_DESeq2.mDC.txt',sep="\t",header=T)
mDC_Influenza <- DEGs[DEGs$contrast == "Influenza vs control",] 

processed_data <- mDC_Influenza %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      log2FoldChange > logfc_threshold & padj < padj_threshold  ~ "Up-regulated",
      log2FoldChange < -logfc_threshold & padj < padj_threshold ~ "Down-regulated",
      TRUE                                                     ~ "Not significant"
    ))

if (any(is.infinite(processed_data$neg_log10_padj))) {
  max_val <- max(processed_data$neg_log10_padj[is.finite(processed_data$neg_log10_padj)], na.rm = TRUE)
  processed_data$neg_log10_padj <- ifelse(is.infinite(processed_data$neg_log10_padj), max_val * 1.1, processed_data$neg_log10_padj)
}

top_genes <- processed_data[processed_data$gene.name %in% c("CXCL8","CXCL2","CLEC4E","STSL","SIRPA","TLR2","HMGN1","SNHG32"),]

mDC_plot <- ggplot(processed_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("Up-regulated" = "#E64B35",
                                "Down-regulated" = "#3C5488",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
  geom_text_repel(data = top_genes,
                  aes(label = gene.name),
                  size = 4,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = 'grey50',
                  max.overlaps = Inf) +
  labs(title = "Myelo.Dend.",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "none"
  )
  
ggsave(paste0(output_dir, "Fig3d.pdf"), plot = mDC_plot, width = 6, height = 4)


