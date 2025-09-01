
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggbeeswarm)

output_dir <- 'figures/'

Lee_dataset <- readRDS('data/integrated.annotated.Lee_dataset.RDS')

# Fig2a: UMAP plot of all cells colored by celltype cluster
custom_colors <- c(
    "Mono.1" = "#D5BB21FF", "Mono.2" = "#F8B620FF", "Mono.3" = "#F89217FF", "Mono.4" = "#F06719FF", "Mono.5" = "#E03426FF",  
    "mDC.1" = "#f8b8c4", "mDC.2" = "#f596a3", "mDC.3" = "#f17e96", "mDC.4" = "#e35177", "mDC.5" = "#ec4269", "mDC.6" = "#ad3f5c", "mDC.7" = "#a7213a",
    "Myelo.Dend." = "#F64971FF", "Dend." = "#FC719EFF",
    "T1" = "#7873C0FF", "T2" = "#4F7CBAFF", "T3" = "#1BA3C6FF", "T4" = "#2CB5C0FF", "T5" = "#30BCADFF", "T6" = "#21B087FF", "Prolif." = "#33A65CFF",
    "B1" = "#57A337FF", "B2" = "#A2B627FF",
    "NK" = "#A26DC2FF",
    "RBC" = "#A9A9A9", "PL" = "#696969"
)

pdf(paste0(output_dir, "Fig2a.pdf"), width = 10, height = 9.5)
DimPlot(Lee_dataset, group.by = 'celltype.cluster', shuffle = TRUE) + scale_color_manual(values = custom_colors)
dev.off()

# Fig2b: UMAP plot of all cells colored by subject type
infect_colors <- c(
  "Asymptomatic case of COVID-19 patient" = "#C5A1C5",
  "healthy control" = "#D3D3D3",
  "Influenza patient" = "#A3D9A5",
  "mild COVID-19 patient" = "#C5A1C5",
  "severe COVID-19 patient" = "#8B5E8B"
)

pdf(paste0(output_dir, "Fig2b.pdf"), width = 7, height = 8)
DimPlot(Lee_dataset, group.by = 'subject.type', raster = FALSE, shuffle = TRUE) + scale_color_manual(values = infect_colors) + NoLegend()
dev.off()

# Fig2c: Dot plot of celltype-specific markers
markers <- c("CD14", "FCGR3A","CD68","CD163","CSF1R","CD1C","NRP1","CLEC4C",
             "CD3D", "CD3E", "CD4","CD8A", "IFNG","NCAM1","ZNF683", "MKI67","CD79A", "MS4A1")

pdf(paste0(output_dir, "Fig2c.pdf"), width = 16, height = 5)
DotPlot(Lee_dataset, features = rev(markers), group.by = "celltype.cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient(low = "lightgrey", high = "#C80815") + coord_flip()
dev.off()

# Fig2d: Stacked bar plot of composition of celltypes of each cell cluster by diease group
x <- rep(levels(Lee_dataset$celltype.cluster), length(levels(Lee_dataset$subject.type)))
y <- rep(levels(Lee_dataset$subject.type), each = length(levels(Lee_dataset$celltype.cluster)))
data <- data.frame(x,y)
data$value <- 0
data$x <- factor(data$x, levels=levels(Lee_dataset$celltype.cluster)) 
data$y <- factor(data$y, levels=levels(Lee_dataset$subject.type))

for(i in levels(Lee_dataset$celltype.cluster)){
	for( j in levels(Lee_dataset$subject.type)){
		data[data$x==i&data$y==j,"value"] <- table(subset(Lee_dataset,subset=celltype.cluster==i)$subject.type)[j]
	}}

ggplot(data, aes(fill = y, y = value, x = x)) + 
    geom_bar(position = "fill", stat = "identity", colour = "white") +
    scale_fill_manual(values = infect_colors) +  # 색상 지정
    theme_classic() + NoLegend() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(output_dir, "Fig2d.pdf"), width = 20) 

# Fig2e: Odds ratio of celltype-specific markers between disease groups


make.Pairwise.Tables <- function(geneSets1,geneSets2,background){
    mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
    mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

    d <- t(mem1) %*% mem2;
    b <- abs(t(mem1) %*% (mem2-1))
    c <- abs(t(mem1-1) %*% (mem2))
    a <- t(mem1-1) %*% (mem2-1);

    ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

    pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
    names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
    return(pairwise.tables)
}

do.FisherExactTest <- function(table.count,or = 1,alternative = "greater",N_bg = NULL){
    if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
    out <- fisher.test(x = table.count,or = or,alternative = alternative)
    odds.ratio <- out$estimate
    p.value <- out$p.value;
    geneSet1.count <- rowSums(table.count)[2]
    geneSet2.count <- colSums(table.count)[2]
    expected.count <- geneSet1.count/N_bg * geneSet2.count
    overlap.count <- table.count[2,2];
    fold.change <- overlap.count/expected.count
 
    out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
    names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
    return(out)
}

perform.AllPairs.FET <- function(geneSets1,geneSets2,background, or = 1,alternative = "greater", adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL){
    pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
    if (do.multicore){
        cl <- parallel::makeCluster(n.cores)
        registerDoParallel(cl)
        cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
        fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
        fact <- factor(fact[1:length(pairwise.tables)])
        split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
        output <- foreach (tbl = split.tables,.combine = 'c',.export = c("do.FisherExactTest")) %dopar% {
            out <- lapply(tbl,function(x,y,z) do.FisherExactTest(table.count = x,or = y,alternative = z),y = or,z = alternative) 
		    return(out)
        }
    }else{
        output <- lapply(pairwise.tables,function(x,y,z) do.FisherExactTest(table.count = x,or = y,alternative = z),y = or,z = alternative)
    }
    output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
    int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
    int.element <- do.call(c,int.element)
    if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
    return(output)
}

full.name.map = c("Mono." = "Monocytes","Myelo.Dend." = "Myeloid Dendritic","Dend." = "Dendritic",
                    "RBC" = "Erythrocyte","PL" = "Platelet","B" = "B-cell","T" = "T-cell","NK" = "NK","Prolif." = "Ki67+ Lymphocytes")
vec1 = Lee_dataset@meta.data$celltype.cluster;
vec2 = Lee_dataset@meta.data$sample
names(vec1) = names(vec2) = rownames(Lee_dataset@meta.data)
  
mod1 = split(names(vec1),factor(vec1))
mod2 = split(names(vec2),factor(vec2))
res = perform.AllPairs.FET(mod1,mod2,background = names(vec1))
res$set1_proportion = res$actual.overlap/res$set1_size
res$set2_proportion = res$actual.overlap/res$set2_size
res$disease = Lee_dataset@meta.data$subject.group.name[match(res$set2_Name,Lee_dataset@meta.data$sample)]
res$cell.type = full.name.map[gsub("[0-9]$","",res$set1_Name)]
  
cls.id = unique(res$set1_Name)
wout.matrix = matrix(0,nrow = 0,ncol = 3);
colnames(wout.matrix) = c("mu.x","mu.y","p.value")
dis.vec = cls.vec = c()
dis.names = setdiff(unique(res$disease),c("healthy control","Asymptomatic case of COVID-19 patient"))
for (i in 1:length(cls.id)){
    yval = subset(res,set1_Name == cls.id[i] & disease == "healthy control")$odds.ratio
    for (j in 1:length(dis.names)){
      xval = subset(res,set1_Name == cls.id[i] & disease == dis.names[j])$odds.ratio
      wout = t.test(x = xval,y = yval)
      out = c("mu.x" = mean(xval),mu.y = mean(yval),p.value = wout$p.value)
      wout.matrix = rbind(wout.matrix,out)
      dis.vec = c(dis.vec,dis.names[j])
      cls.vec = c(cls.vec,cls.id[i]) 
    }
  }
  
wout.data = data.frame(cluster.id = cls.vec,case.id = dis.vec,ctrl.id = rep("healthy",nrow(wout.matrix)),as.data.frame(wout.matrix))
wout.data$FC = wout.data$mu.x/wout.data$mu.y
wout.data$cell.type = full.name.map[gsub("[0-9]$","",wout.data$cluster.id)]
  
dis.colmap = infect_colors
res <- res[res$disease != "Asymptomatic case of COVID-19 patient",]
res$set1_Name <- factor(res$set1_Name, levels = c( "Mono.1", "Mono.2", "Mono.3", "Mono.4", "Mono.5", "Myelo.Dend.", "Dend.","T1", "T2", "T3", "T4", "T5", "T6", "NK", "Prolif.", "B1", "B2",  "RBC", "PL"))
ggplot(data = res) + 
    geom_quasirandom(aes(x = set1_Name,y = odds.ratio,colour = disease),alpha = 0.7,dodge.width = 0.5) + 
    geom_boxplot(aes(x = set1_Name,y = odds.ratio,colour = disease,fill = disease),alpha = 0.7) + 
    geom_text(data = subset(wout.data,FC > 1 & p.value < 0.2 & case.id == "Influenza patient"),aes(x = cluster.id),y = 5,label = "*",size = 8,colour = infect_colors["Influenza patient"]) + 
    geom_text(data = subset(wout.data,FC > 1 & p.value < 0.2 & case.id == "mild COVID-19 patient"),aes(x = cluster.id),y = 5.5,label = "*",size = 8,colour = infect_colors["mild COVID-19 patient"]) + 
    geom_text(data = subset(wout.data,FC > 1 & p.value < 0.2 & case.id == "severe COVID-19 patient"),aes(x = cluster.id),y = 6,label = "*",size = 8,colour = infect_colors["severe COVID-19 patient"]) + 
    labs(x = "Cluster",y = "OR") + 
    guides(colour = "none",fill = guide_legend(title = "Disease",ncol = 3)) + 
    scale_colour_manual(values = dis.colmap) + 
    scale_fill_manual(values = dis.colmap) + 
    scale_y_continuous(limits = c(0,7.5)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 16,angle = 45,vjust= 1,hjust =1 ),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19),
          strip.text = element_text(angle = 90,hjust = 0.5,size = 17),
          legend.position = "bottom",legend.direction = "horizontal",legend.title = element_text(size = 18),
          legend.text = element_text(size = 15))
ggsave(paste0(output_dir, "Fig2e.pdf"), width = 20) 
