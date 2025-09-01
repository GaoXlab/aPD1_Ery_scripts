library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(SingleR)
library(scater)
library(ExperimentHub)
library(AnnotationHub)
library(AnnotationDbi)
library(celldex)
library(reshape2)
library(ggsci)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

source('./stat.r') # aPD1_Ery_scripts/scRNAseq_process/Mus
workpath <- "./3.combination/Spleen/" # aPD1_Ery_scripts/scRNAseq_process/Mus/
setwd(workpath)

dir.create('Tsubset')
rdsfile <- str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.rds")

immune.combined <- readRDS(rdsfile)
DefaultAssay(immune.combined) <- "RNA"
sample_list <- c('IgG_mEry','aPD1','aPD1_mEry')
immune.combined$sample_label <- factor(immune.combined$sample_label, levels=sample_list)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 8) + NoLegend()
ggsave(str_c('./Tsubset/supFig1.Spleen_SeuratClusters.',dtVar,'.checkinput.pdf'), p1, width=5.5, height=5.5)

# 1. T细胞亚群细分
Idents(immune.combined) <- immune.combined$celltype
pbmc_T <- subset(immune.combined, idents=c('CD4 T cell', 'CD8 T cell', 'NK T cell'))#))
DefaultAssay(pbmc_T) <- "integrated"
# Run the standard workflow for visualization and clustering
pbmc_T <- ScaleData(pbmc_T, verbose = FALSE)
pbmc_T <- RunPCA(pbmc_T, npcs = 20, verbose = FALSE)
pbmc_T <- RunTSNE(pbmc_T, reduction = "pca", dims = 1:20)
pbmc_T <- FindNeighbors(pbmc_T, reduction = "pca", dims = 1:20)
pbmc_T <- FindClusters(pbmc_T, resolution = 0.5)
p2 <- DimPlot(pbmc_T, reduction = "tsne", label = TRUE, repel = TRUE) + theme_bw()
ggsave(str_c('./Tsubset/supFig1.Spleen_T_SeuratClusters.',dtVar,'.pdf'), p2, width=5.5, height=4.5)

# 2. annotation
DefaultAssay(pbmc_T) <- "RNA"
Idents(pbmc_T) <- factor(pbmc_T$seurat_clusters, levels=c(0,11,4,10,5,1,2,3,6,7,8,12,9,13))
markers <- c('Cd3e','Cd4','Cd8a','Cd44','Gzmk','Tbx21','Ccl5','Cxcr3','Klrb1','Foxp3','Il2ra','Ctla4','Ikzf2','Ly6c2','Cd69','Sell','Tcf7','Il7r','Bcl2','Il6ra','Klrb1c','Rorc','Ccr6','Cd24a')
p1 <- DotPlot(pbmc_T, features = markers, dot.scale = 3, col.min=-1, cols=c('white','red','red4')) + RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('./Tsubset/supFig1.Spleen_T_SeuratClusters.dotplot.',dtVar,'.pdf'), p1, width=8, height=4)

table(pbmc_T$seurat_clusters)
pbmc_T$celltype = dplyr::case_when(
  pbmc_T$seurat_clusters %in% c(1,2) ~ "CD4 Naive",
  # pbmc_T$seurat_clusters %in% c(2) ~ "CD4 Tem",
  pbmc_T$seurat_clusters %in% c(6) ~ "Central Treg",
  pbmc_T$seurat_clusters %in% c(3) ~ "Effector Treg",
  pbmc_T$seurat_clusters %in% c(0,11) ~ "CD8 Naive",
  pbmc_T$seurat_clusters %in% c(4) ~ "CD8 Tem",
  pbmc_T$seurat_clusters %in% c(5) ~ "CD8 Te",
  pbmc_T$seurat_clusters %in% c(10) ~ "CD8 Tcm",
  pbmc_T$seurat_clusters %in% c(9,13) ~ "NKT",
  pbmc_T$seurat_clusters %in% c(7) ~ "DN T",
  pbmc_T$seurat_clusters %in% c(8) ~ "DP T",
  pbmc_T$seurat_clusters %in% c(12) ~ "Tgd")
  # pbmc_T$seurat_clusters %in% c(7) ~ "Unknown")
table(pbmc_T$celltype)

marjor_cluster <- c("CD4 Naive","CD4 Tem","Central Treg","Effector Treg","CD8 Naive","CD8 Tcm","CD8 Tem","CD8 Te","NKT","Tgd","DN T","DP T")
pbmc_T$celltype <- factor(pbmc_T$celltype, levels = marjor_cluster)
Idents(pbmc_T) = 'celltype'
p1 <- DimPlot(pbmc_T, reduction = "tsne", label = TRUE, label.size = 4, repel=TRUE, cols=c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5','grey'))+theme_bw() + theme(legend.position="bottom")#+ NoLegend()
ggsave(str_c('./Tsubset/supFig2.Spleen_T_Anno1_SeuratClusters.',dtVar,'.pdf'), p1, width=4.5, height=5)
p1 <- DotPlot(pbmc_T, features = markers, dot.scale = 3, col.min=-2, cols=c('white','red','red4')) + RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('./Tsubset/supFig2.Spleen_T_Anno1_SeuratClusters.dotplot.',dtVar,'.pdf'), p1, width=8, height=4)


# 3. Calculating Cell numbers and ratio
rds = pbmc_T
rds2 = immune.combined
colors = c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5','grey')
group_colors = c('#A9B8C6', '#0E5F98FF', 'red4')
sample_labels = c('IgG_mEry','aPD1','aPD1_mEry')
marjor_cluster <- c("CD4 Naive","Central Treg","Effector Treg","CD8 Naive","CD8 Tcm","CD8 Tem","CD8 Te","NKT","Tgd","DN T","DP T")
cellnumber_summary(rds, rds2, sample_labels,sample_labels, marjor_cluster, colors,group_colors, str_c('Tsubset/supFig2.Spleen_inTotal_',dtVar))
rds2 = pbmc_T
cellnumber_summary(rds, rds2, sample_labels,sample_labels, marjor_cluster, colors, group_colors,str_c('Tsubset/supFig2.Spleen_inTsubset_',dtVar))

### 保存文件
saveRDS(pbmc_T, file = str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.Tsubset.rds"))



# 4. cellcycle figure
pbmc_T = readRDS(str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.Tsubset.rds"))
DefaultAssay(pbmc_T) <- "RNA"
Idents(pbmc_T) = 'celltype'

phase_df <- as.data.frame(cbind(as.character(pbmc_T$celltype), as.character(pbmc_T$sample_label), as.character(pbmc_T$seurat_clusters), 
                                as.character(pbmc_T$Phase), as.numeric(pbmc_T$S.Score), as.numeric(pbmc_T$G2M.Score)))
colnames(phase_df) <- c('celltype', 'samplelabel', 'seurat_clusters', 'Phase', 'S.Score', 'G2M.Score')
phase_df$samplelabel = factor(phase_df$samplelabel, levels=c('IgG_mEry','aPD1','aPD1_mEry'))
phase_df$seurat_clusters = factor(phase_df$seurat_clusters, levels=as.character(sort(as.numeric(names(table(phase_df$seurat_clusters))))))
phase_df$Phase = factor(phase_df$Phase, levels=c('G1', 'G2M', 'S'))
phase_df$S.Score = as.numeric(as.character(phase_df$S.Score))
phase_df$G2M.Score = as.numeric(as.character(phase_df$G2M.Score))

cellcount <- aggregate(phase_df$Phase, by=list(type=phase_df$samplelabel, celltype=phase_df$celltype), table)
total_cells <- table(phase_df$samplelabel)
cellcount_perc <- c()
for(sample in c('IgG_mEry','aPD1','aPD1_mEry')){
    print(sum(cellcount[cellcount$type==sample,-c(1,2)]))
    cellcount_perc <- rbind(cellcount_perc, cellcount[cellcount$type==sample, -c(1,2)]/total_cells[sample])
}
cellcount_perc <- as.data.frame(cellcount_perc)
cellcount_perc$celltype <- rep(cellcount[cellcount$type=='IgG_mEry','celltype'], 3)
cellcount_perc$samplelabel <- c(rep('IgG_mEry', length(cellcount[cellcount$type=='IgG_mEry','celltype'])),rep('aPD1', length(cellcount[cellcount$type=='aPD1','celltype'])),rep('aPD1_mEry', length(cellcount[cellcount$type=='aPD1_mEry','celltype'])))
cellcount_perc_melt <- melt(cellcount_perc)
marjor_cluster <- c("CD4 Naive","CD4 Tem","Effector Treg","Central Treg","CD8 Naive","CD8 Tcm","CD8 Tem","NKT","DN T", "Tgd", "Unknown")
cellcount_perc_melt$celltype <- factor(cellcount_perc_melt$celltype, levels=marjor_cluster)

cellcount = as.data.frame(cellcount)
cellcount$sum = apply(cellcount[,-c(1,2)],1,sum)
nor <- c()
for(i in c(1:3)){
    nor <- cbind(nor, cellcount$x[,i]/cellcount$sum)
}
cellcount$G1_perc = nor[,1]
cellcount$G2M_perc = nor[,2]
cellcount$S_perc = nor[,3]

cellcount_perc <- as.data.frame(cellcount[, setdiff(colnames(cellcount),'x')])
cellcount_perc_melt <- melt(cellcount_perc)

cellcount_perc_melt$variable <- gsub('_perc','',cellcount_perc_melt$variable)
cellcount_perc_melt$variable <- gsub('_perc','',cellcount_perc_melt$variable)
cellcount_perc_melt$type <- factor(cellcount_perc_melt$type, levels=c('IgG_mEry','aPD1','aPD1_mEry'))
cellcount_perc_melt$value <- cellcount_perc_melt$value*100
p3_0 <- ggplot(data=cellcount_perc_melt[(cellcount_perc_melt$variable!='sum')&(cellcount_perc_melt$celltype=='CD8 Te'), ], aes(x=type, y=value, fill=variable)) +
  geom_bar(stat="identity")+scale_fill_brewer(palette="Paired")+theme_bw()+facet_grid(.~celltype) +ylab('Proportion (%)')+xlab('')+
  theme(strip.background=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
ggsave('Tsubset/CD8TE_cellcycle.pdf', p3_0, width=2.4, height=3)


my_comparisons=list(c('IgG_mEry','aPD1'),c('aPD1_mEry','aPD1'),c('IgG_mEry','aPD1_mEry'))
p3_1 <- ggplot(data = phase_df[phase_df$celltype=='CD8 Te', ], aes(x = samplelabel, y = S.Score, fill=samplelabel)) + #
  geom_boxplot(outlier.colour = NA)+theme_bw()+scale_fill_manual(values=c('#A9B8C6', '#0E5F98FF', 'red4'))+facet_grid(. ~ celltype)+ 
  theme(strip.background=element_blank(), legend.position='right',legend.title=element_blank(),axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", vjust=0.4, size=3)

p3_2 <- ggplot(data = phase_df[phase_df$celltype=='CD8 Te', ], aes(x = samplelabel, y = G2M.Score, fill=samplelabel)) + #
  geom_boxplot(outlier.colour = NA)+theme_bw()+scale_fill_manual(values=c('#A9B8C6', '#0E5F98FF', 'red4'))+ facet_grid(. ~ celltype)+ 
  theme(strip.background=element_blank(), legend.position='right',legend.title=element_blank(),axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", vjust=0.4, size=3)
pdf('Tsubset/CD8TE_cellcycle2.pdf', width=5, height=3)
plot_grid(p3_1, p3_2, ncol=2)
dev.off()


sessionInfo()