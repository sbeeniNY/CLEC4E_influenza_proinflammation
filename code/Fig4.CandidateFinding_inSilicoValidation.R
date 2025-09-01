
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

output_dir <- 'figures/'

# Fig4a: Heatmap of candidate genes
load("data/targets_heatmap.RData")

targets <- c("ACSL1", "ALDH2", "SDCBP", "SLC7A7", "NAMPT", "AQP9", "CLEC10A", "UPP1", "VCAN", "CLEC4E")

col_KD <- circlize::colorRamp2(c(0,5), c("#E0E0E0","#DA4453")) 
col_fun <- circlize::colorRamp2(c(1,0,-1), c("#DA4453","white","#4A89DC")) 
col_corr <- circlize::colorRamp2(c(.5,0,-.5), c("#d73027","#ffffbf","#4575b4")) 

h1 <- Heatmap(as.matrix(tD[,-1]), name="log2(FC.DEGs)", col = col_fun, cluster_rows = F)
h2 <- Heatmap(as.matrix(-log10(tmp[,iKD])), name="-lg(P(KDs))", col = col_KD, cluster_rows = F)
h3 <- Heatmap(as.matrix(t3g[,-1]), name="corr", col = col_corr, cluster_rows = F)
h4 <- Heatmap(as.matrix(log10(tmp[,2:5]+.1)), name="log10(#Lit)", col = c("#fdfdd7","#225ea8"), cluster_rows = F)
h5 <- Heatmap(tmp[,6:9], name="Tally", col = c("#ffffe4","#993404"), 
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%i", tmp[i, j+5]), x, y, gp = gpar(fontsize = 10))},
              cluster_rows = F)

ht.list <- h1 + h2 + h3 + h4 + h5 

pdf(paste0(output_dir, "Fig4a.pdf"), width = 20, height = 8)
draw(ht.list, main_heatmap=4)
dev.off()

# Fig4b: in silico validation of candidate genes with bulk RNA-seq data
load('data/Mouse_Validation.RData')

dot_data <- data.frame(gene = rep(rownames(mouse.FC),n=ncol(mouse.FC)),
                      sample = rep(colnames(mouse.FC), each = nrow(mouse.FC)))
dot_data$log2FC <- 0
dot_data$qval <- 0
for(i in rownames(mouse.FC)){
    for(j in colnames(mouse.FC)){
        dot_data[dot_data$gene == i & dot_data$sample == j,]$log2FC <- mouse.FC[i,j]
        dot_data[dot_data$gene == i & dot_data$sample == j,]$qval <- mouse.qval[i,j]
}}
dot_data$gene <- factor(dot_data$gene, levels = rev(c("Clec4e","Vcan","Upp1","Clec10a","Aqp9","Nampt","Slc7a7","Sdcbp","Aldh2","Acsl1")))

pdf(paste0(output_dir, "Fig4b.pdf"), width = 7, height = 4)
ggplot(dot_data, aes(x = sample, y = gene, color = log2FC, size = qval)) + 
  geom_point() + scale_size_continuous(limits = c(0, 2.5), range = c(1, 10)) +
  scale_color_gradient2(low = "#4A89DC", mid = "white", high = "#DA4453", 
                        midpoint = 0, limits = c(-2, 4.2)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Fig4c: mDC subtype
Lee_mDC <- readRDS('data/mDC_sub.RDS')

Lee_mDC_colors <- c(
  "mDC.1" = "#F06292",
  "mDC.2" = "#A5D6A7",
  "mDC.3" = "#66BB6A",
  "mDC.4" = "#FFB74D",
  "mDC.5" = "#FFF176",
  "mDC.6" = "#4FC3F7",
  "mDC.7" = "#BA68C8"
)

DimPlot(Lee_mDC, group.by = 'newCluster', raster = FALSE) + scale_color_manual(values = Lee_mDC_colors) + NoLegend()
ggsave(paste0(output_dir, "Fig4c.pdf"), height = 4, width = 5)

dot_data <- data.frame(
    clusters = rep(unique(Lee_mDC$newCluster), each = length(unique(Lee_mDC$subject.group.name))),
    condition = rep(unique(Lee_mDC$subject.group.name), times = length(unique(Lee_mDC$newCluster))),
    percent_expressed = 0,
    average_expression = 0
)

for (i in 1:length(unique(Lee_mDC$newCluster))) {
    for (j in 1:length(unique(Lee_mDC$subject.group.name))) {
        tryCatch({
            sub_mDC <- subset(Lee_mDC, subset = newCluster == unique(Lee_mDC$newCluster)[i] & subject.group.name == unique(Lee_mDC$subject.group.name)[j])
        }, error = function(e) {
            message(paste0("i: ", i, " j: ", j))
            return(NULL)
        })
        if (ncol(sub_mDC)>1) {
        a <- DotPlot(object = sub_mDC, features = "CLEC4E")
        dot_data$percent_expressed[dot_data$clusters == unique(mDC$newCluster)[i] & dot_data$condition == unique(mDC$subject.group.name)[j]] <- a$data$pct.exp
        dot_data$average_expression[dot_data$clusters == unique(mDC$newCluster)[i] & dot_data$condition == unique(mDC$subject.group.name)[j]] <- a$data$avg.exp
        } else {
            next
        }
    }}

dot_data$condition <- factor(dot_data$condition, levels = c('healthy control','Influenza patient','Asymptomatic case of COVID-19 patient','mild COVID-19 patient','severe COVID-19 patient'))
dot_data <- dot_data[dot_data$condition %in% c('healthy control','Influenza patient'),]

ggplot(dot_data, aes(x = condition, y = clusters)) +
    geom_point(aes(size = percent_expressed, color = average_expression)) +  
    scale_size(range = c(1, 10), name = "Percent Expressed") +  # 크기 조정
    scale_color_gradient(low = "grey", high = "red") +  # 색상 그라디언트 설정
    labs(x = "Condition", y = "Clusters") +     
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # x축 레이블 45도 기울이기
ggsave(paste0(output_dir, "Fig4c_dot.pdf"), width = 4, height = 4)

# Fig4d: in silico validation of candidate genes with public single cell data

public_mDC <- readRDS("data/Flu_mDC_sub.RDS")

public_mDC_colors <- c(
  "0" = "#FFAB91",
  "1" = "#FFD54F",
  "2" = "#E6EE9C",
  "3" = "#8BC34A",
  "4" = "#4DB6AC",
  "5" = "#B39DDB"
)

DimPlot(public_mDC, reduction = "umap.cca",group.by = 'newCluster') + scale_color_manual(values = public_mDC_colors) + NoLegend()
ggsave(paste0(output_dir, "Fig4d.pdf"), height = 5, width = 5) 

dot_data <- data.frame(
    clusters = rep(unique(public_mDC$newCluster), each = length(unique(public_mDC$Infection))),
    condition = rep(unique(public_mDC$Infection), times = length(unique(public_mDC$newCluster))),
    percent_expressed = 0,
    average_expression = 0
)

for (i in 1:length(unique(public_mDC$newCluster))) {
    for (j in 1:length(unique(public_mDC$Infection))) {
        tryCatch({
            sub_public_mDC <- subset(public_mDC, subset = newCluster == unique(public_mDC$newCluster)[i] & Infection == unique(public_mDC$Infection)[j])
        }, error = function(e) {
            message(paste0("i: ", i, " j: ", j))
            return(NULL)    
        })
        a <- DotPlot(object = sub_public_mDC, features = "CLEC4E")
        dot_data$percent_expressed[dot_data$clusters == unique(public_mDC$newCluster)[i] & dot_data$condition == unique(public_mDC$Infection)[j]] <- a$data$pct.exp
        dot_data$average_expression[dot_data$clusters == unique(public_mDC$newCluster)[i] & dot_data$condition == unique(public_mDC$Infection)[j]] <- a$data$avg.exp
    }}

dot_data$condition <- factor(dot_data$condition, levels = c('Control','Infected'))

ggplot(dot_data, aes(x = condition, y = clusters)) +
    geom_point(aes(size = percent_expressed, color = average_expression)) +  
    scale_size(range = c(1, 10), name = "Percent Expressed") +  # 크기 조정
    scale_color_gradient(low = "grey", high = "red") +  # 색상 그라디언트 설정
    labs(x = "Condition", y = "Clusters") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # x축 레이블 45도 기울이기
ggsave(paste0(output_dir, "Fig4d_dot.pdf"), width = 4, height = 4)

# Fig4e: Module score, finding influenza-responsive mDC subtypes in public single cell data 
Idents(Lee_mDC) <- 'subject.group'
Lee_mDC <- JoinLayers(Lee_mDC)
Flu.Marker <- FindMarkers(subset(Lee_mDC, subset = newCluster %in% c('mDC.2','mDC.3')), ident.1='Influenza patient', ident.2='healthy control', only.pos = TRUE)
Flu.marker.up <- Flu.Marker %>% filter(p_val_adj < 0.05) %>% dplyr::arrange(avg_log2FC) 
Flu.marker.up <- Flu.marker.up[Flu.marker.up$pct.1>0.7 & Flu.marker.up$pct.2<0.3,] %>% dplyr::arrange(desc(avg_log2FC)) 

public_mDC <- AddModuleScore(public_mDC, features = rownames(Flu.marker.up), name="Lee_mDC_Flu_marker_Score")
FeaturePlot(public_mDC, features = 'Lee_mDC_Flu_marker_Score', raster = FALSE) + scale_color_gradientn(colors = c("#DEDEDE","#A5D6A7","#66BB6A","#448461","#448461"))
ggsave('Lee_mDC_Flu_marker_Score.pdf', width = 5, height = 5)
