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
library(org.Mm.eg.db)
library(rlist)
library(clusterProfiler)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")
source('./stat.r') # aPD1_Ery_scripts/scRNAseq_process/Mus
workpath <- "./3.combination/Spleen/" # aPD1_Ery_scripts/scRNAseq_process/Mus/
setwd(workpath)

dir.create('Myeloidsubset')
rdsfile <- str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.rds")

immune.combined <- readRDS(rdsfile)
DefaultAssay(immune.combined) <- "RNA"
sample_list <- c('IgG_mEry','aPD1','aPD1_mEry')
immune.combined$sample_label <- factor(immune.combined$sample_label, levels=sample_list)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4) + NoLegend()
ggsave(str_c('./Myeloidsubset/supFig1.Spleen_SeuratClusters.',dtVar,'.checkinput.pdf'), p1, width=5.5, height=5.5)

# 1. Myeloid细胞亚群细分
Idents(immune.combined) = immune.combined$celltype
pbmc_Myeloid <- subset(immune.combined, idents=c('Basophil','DC','Macrophage','Monocyte','Neutrophil'))#))
DefaultAssay(pbmc_Myeloid) <- "integrated"
# Run the standard workflow for visualization and clustering
pbmc_Myeloid <- ScaleData(pbmc_Myeloid, verbose = FALSE)
pbmc_Myeloid <- RunPCA(pbmc_Myeloid, npcs = 20, verbose = FALSE)
pbmc_Myeloid <- RunTSNE(pbmc_Myeloid, reduction = "pca", dims = 1:20)
pbmc_Myeloid <- FindNeighbors(pbmc_Myeloid, reduction = "pca", dims = 1:20)
pbmc_Myeloid <- FindClusters(pbmc_Myeloid, resolution = 0.6)
p2 <- DimPlot(pbmc_Myeloid, reduction = "tsne", label = TRUE, repel = TRUE) + theme_bw()
ggsave(str_c('./Myeloidsubset/supFig1.Spleen_Myeloid_SeuratClusters.',dtVar,'.pdf'), p2, width=5.5, height=4.5)

# 2. annotation
DefaultAssay(pbmc_Myeloid) <- "RNA"
Idents(pbmc_Myeloid) <- pbmc_Myeloid$seurat_clusters
Idents(pbmc_Myeloid) <- factor(Idents(pbmc_Myeloid), levels=c(11,0,3,7,5,9,8,10,13,2,14,1,6,4,12))
# macrophages / dc / eosi / baso
myeliod_markers <- unique(c('Adgre1','Spic','Mrc1','Cd86','Cd163','Cd68','Cd80','Itgax','Xcr1','Clec9a','Sirpa','Siglech','Bst2','Cd83','Cd74','Itgam','Ly6g','Cxcr2','Cd14','Ccr2','Csf1r','Arg2','Il1b','S100a8','S100a9','Fcgr3','Cd84','Ctsd','Stat1','Stat3','Tgfb1','Wfdc17','Clec4e','Il1f9','Fcer1a','Cpa3','Fcer1g','Cst3','Ccr3','Siglecf'))
myeliod_markers <- unique(c('Adgre1','Spic','Cd68','Mrc1','Itgax','Xcr1','Clec9a','Sirpa','Siglech','Bst2','Ly6d','Cd83','Cd274','Fscn1', 'Cacnb3', 'Ccr7', 'Cd40', 'Tmem123','Stat3',
                'Itgam','Ly6g','Cxcr2','Cd14','Arg2','Il1f9','Il1b','Fcgr3','Cd84','Ctsd','Kit','Cd34','Sox4','Ccr2','Csf1r','Lyz2','S100a4','S100a8','Ccl9',
                'Cd36','Cd300e','Stat1','Stat3','Fcer1a','Cpa3','Fcer1g','Cst3','Cd74'))#,'Siglecf','Itga2'))
                # 'S100a8','S100a9','Tgfb1','Wfdc17','Clec4e',,'Ccr3','Siglecf'))
p2 <- DotPlot(pbmc_Myeloid, features =myeliod_markers, cols = c("blue", "red", "grey"), dot.scale = 3) + RotatedAxis()
p1 <- DotPlot(pbmc_Myeloid, features = myeliod_markers, dot.scale = 3, col.min=-2, cols=c('white','red','red4')) + RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('Myeloidsubset/supFig1.Non-Lymphocytes.Spleen_SeuratClusters.dotplot.',dtVar,'.pdf'), p1, width=10, height=4)

sum(table(pbmc_Myeloid$seurat_clusters))
pbmc_Myeloid$celltype = dplyr::case_when(
  pbmc_Myeloid$seurat_clusters %in% c(0,3,11) ~ "Macrophage",
  pbmc_Myeloid$seurat_clusters %in% c(5,7,9) ~ "cDC",
  pbmc_Myeloid$seurat_clusters %in% c(8,10,13) ~ "pDC",
  pbmc_Myeloid$seurat_clusters %in% c(2) ~ "mDC",
  pbmc_Myeloid$seurat_clusters %in% c(14) ~ "Neutrophil",
  pbmc_Myeloid$seurat_clusters %in% c(1) ~ "PMN-MDSC",
  pbmc_Myeloid$seurat_clusters %in% c(4) ~ "M-MDSC",
  pbmc_Myeloid$seurat_clusters %in% c(6) ~ "Monocyte",
  pbmc_Myeloid$seurat_clusters %in% c(12) ~ "Basophil")
table(pbmc_Myeloid$celltype)
marjor_cluster = c("Macrophage","cDC","pDC","mDC","PMN-MDSC","M-MDSC","Neutrophil","Monocyte","Basophil")
pbmc_Myeloid$celltype <- factor(pbmc_Myeloid$celltype, levels=marjor_cluster)
Idents(pbmc_Myeloid) = 'celltype'
p1 <- DimPlot(pbmc_Myeloid, reduction = "tsne", label = TRUE, label.size = 4, repel=TRUE, cols=c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5','grey'))+theme_bw() + theme(legend.position="bottom")
ggsave(str_c('./Myeloidsubset/Fig1.Spleen_Myeloid_Anno1_SeuratClusters.',dtVar,'.pdf'), p1, width=4.5, height=5)
p1 <- DotPlot(pbmc_Myeloid, features = myeliod_markers, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis()+theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('./Myeloidsubset/supFig2.Spleen_Myeloid_Anno1_SeuratClusters.dotplot.',dtVar,'.pdf'), p1, width=9, height=4)


# 3. Calculating Cell numbers and ratio
rds = pbmc_Myeloid
rds2 = immune.combined
sample_labels = c('IgG_mEry','aPD1','aPD1_mEry')
group_colors = c('#A9B8C6', '#0E5F98FF', 'red4')
colors = c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5','grey')
cellnumber_summary(rds, rds2, sample_labels,sample_labels, marjor_cluster, colors,group_colors,  str_c('./Myeloidsubset/supFig2.Spleen_inTotal_',dtVar))
# rds2 = pbmc_Myeloid
# cellnumber_summary(rds, rds2, sample_labels,sample_labels, marjor_cluster, colors,group_colors,  str_c('./Myeloidsubset/supFig2.Spleen_inMyeloidsubset_',dtVar))

### 保存文件
saveRDS(pbmc_Myeloid, file = str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.Myeloidsubset.rds"))



# 4. cellcycle figure
pbmc_Myeloid = readRDS(str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.Myeloidsubset.rds"))
DefaultAssay(pbmc_Myeloid) <- "RNA"
Idents(pbmc_Myeloid) = 'celltype'

phase_df <- as.data.frame(cbind(as.character(pbmc_Myeloid$celltype), as.character(pbmc_Myeloid$sample_label), as.character(pbmc_Myeloid$seurat_clusters), 
                                as.character(pbmc_Myeloid$Phase), as.numeric(pbmc_Myeloid$S.Score), as.numeric(pbmc_Myeloid$G2M.Score)))
colnames(phase_df) <- c('celltype', 'samplelabel', 'seurat_clusters', 'Phase', 'S.Score', 'G2M.Score')
phase_df$samplelabel = factor(phase_df$samplelabel, levels=c('IgG_mEry','aPD1','aPD1_mEry'))
phase_df$seurat_clusters = factor(phase_df$seurat_clusters, levels=as.character(sort(as.numeric(names(table(phase_df$seurat_clusters))))))
phase_df$Phase = factor(phase_df$Phase, levels=c('G1', 'G2M', 'S'))
phase_df$S.Score = as.numeric(as.character(phase_df$S.Score))
phase_df$G2M.Score = as.numeric(as.character(phase_df$G2M.Score))

phase_df[phase_df$celltype=='M-MDSC', 'celltype'] = 'MDSC'
phase_df[phase_df$celltype=='PMN-MDSC', 'celltype'] = 'MDSC'

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
marjor_cluster <- c("Basophil","cDC","mDC","pDC","Macrophage","Monocyte","Neutrophil","MDSC")#,"M-MDSC","PMN-MDSC")
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
p3_0 <- ggplot(data=cellcount_perc_melt[(cellcount_perc_melt$variable!='sum')&(cellcount_perc_melt$celltype=='MDSC'), ], aes(x=type, y=value, fill=variable)) +
  geom_bar(stat="identity")+scale_fill_brewer(palette="Paired")+theme_bw()+facet_grid(.~celltype) +ylab('Proportion (%)')+xlab('')+
  theme(strip.background=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
ggsave('Myeloidsubset/MDSC_cellcycle.pdf', p3_0, width=2.4, height=3)

my_comparisons=list(c('IgG_mEry','aPD1'),c('aPD1_mEry','aPD1'),c('IgG_mEry','aPD1_mEry'))
p3_1 <- ggplot(data = phase_df[phase_df$celltype=='MDSC', ], aes(x = samplelabel, y = S.Score, fill=samplelabel)) + #
  geom_boxplot(outlier.colour = NA)+theme_bw()+scale_fill_manual(values=c('#A9B8C6', '#0E5F98FF', 'red4'))+facet_grid(. ~ celltype)+
  theme(strip.background=element_blank(), legend.position='right',legend.title=element_blank(),axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", vjust=0.4, size=3)

p3_2 <- ggplot(data = phase_df[phase_df$celltype=='MDSC', ], aes(x = samplelabel, y = G2M.Score, fill=samplelabel)) + #
  geom_boxplot(outlier.colour = NA)+theme_bw()+scale_fill_manual(values=c('#A9B8C6', '#0E5F98FF', 'red4'))+ facet_grid(. ~ celltype)+ 
  theme(strip.background=element_blank(), legend.position='right',legend.title=element_blank(),axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", vjust=0.4, size=3)
pdf('Myeloidsubset/MDSC_cellcycle2.pdf', width=5, height=3)
plot_grid(p3_1, p3_2, ncol=2)
dev.off()


sessionInfo()