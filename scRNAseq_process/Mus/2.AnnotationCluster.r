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

rdsfile <- str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.rds")
immune.combined <- readRDS(rdsfile)
DefaultAssay(immune.combined) <- "RNA"
sample_list <- c('IgG_mEry','aPD1','aPD1_mEry')
immune.combined$sample_label <- factor(immune.combined$sample_label, levels=sample_list)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4) + NoLegend()
ggsave(str_c('supFig1.Spleen_SeuratClusters.',dtVar,'.pdf'), p1, width=5.5, height=5.5)


# 1. annotation
Idents(immune.combined) <- factor(Idents(immune.combined), levels=c(1,5,3,6,8,0,2,4,17,7,12,9,13,21,14,10,19,16,11,20,18,15))
epc_marker <- c('Gata1','Klf1','Sox6','Epor','Hemgn','Tfrc',
            'Gypa','Bpgm','Hbb-bs','Hbb-bt','Hba-a1','Hba-a2',
            'Gata2','Klf2','Cebpa','Irf8','Pu.1','Ly6c2','Fcgr1','Csf1r','Csf2ra','Csf3r','C5ar1','Cebpe','Wfdc21')
marker <- c('Cd3e','Tcf7','Cd4','Cd8a','Cd19','Cd79a','Klrb1c','Lyz2','Ccr2','Csf1r','Ly6g','Fcgr3','Cxcr2','Arg2','Il1b','Adgre1','Spic','Mrc1','Cd86','Cd80','Cd163','Cd68','Sirpa','Itgax','Siglech','Bst2','Clec9a','Xcr1','Cd83','Fscn1','Fcer1a','Cpa3','Fcer1g','Cst3','Kit','Gata1','Klf1','Sox6','Epor','Hemgn','Tfrc','Hbb-bt','Sdc1','Tnfrsf17','Bst1','Mcam')
p1 <- DotPlot(immune.combined, features = unique(c(marker)), dot.scale = 3, col.min=-0.5, cols=c('white','red','red4')) + RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('supFig1.Spleen_SeuratClusters.dotplot.',dtVar,'.original.pdf'), p1, width=12, height=6)
p1 <- DotPlot(immune.combined, features = unique(c(marker, c('Kit','Gata1','Klf1','Sox6','Epor','Hemgn','Tfrc'))), cols = c("blue", "red", "grey"), dot.scale = 3) + RotatedAxis() + theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
bm_markers <- c('Kit','Cd34','Hlf','Mpo','Elane','Ms4a3','Fcgr3','Fcgr2b','Ly6a','Spi1','Camp','Ltf','Itgam','Ly6g','S100a8','S100a9','S100a6','Ly6c2','S100a4','Csf1r','Ly86','Adgre1','C1qb','Mrc1',#'Ms4a2','Lmo4',#,'Prtn3','Ctsg','Ms4a3''Fcer1a','Prss34',
'Siglech','Ms4a2','Prss34','Pax5','Vpreb1', 'Ebf1','Cd79a','Il7r','Cd3g','Cd3e','Klrd1','Ptprc',
'Apoe','Car2','Car1','Zfpm1','Gata1','Klf1','Tfrc','Hbb-bt','Hbb-bs','Fos','Fosb','Jun','Junb','Jund',
'Cdh5','Col1a1','Kitl','Spp1')
p1 <- DotPlot(immune.combined, features = bm_markers, cols = c("blue", "red", "grey"), dot.scale = 3) + RotatedAxis() + theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave(str_c('supFig1.Spleen_SeuratClusters.dotplot.',dtVar,'.original.bm.pdf'), p1, width=20, height=6)

table(immune.combined$seurat_clusters)
table(immune.combined$seurat_clusters)
immune.combined$celltype = dplyr::case_when(
  immune.combined$seurat_clusters %in% c(0,2,4,17,7) ~ "B cell",
  immune.combined$seurat_clusters %in% c(15) ~ "Plasma cell",
  immune.combined$seurat_clusters %in% c(1,5) ~ "CD4 T cell",
  immune.combined$seurat_clusters %in% c(3,6,8) ~ "CD8 T cell",
  immune.combined$seurat_clusters %in% c(12) ~ "NK T cell",  
  immune.combined$seurat_clusters %in% c(9) ~ "NK cell",
  immune.combined$seurat_clusters %in% c(13) ~ "Monocyte",
  immune.combined$seurat_clusters %in% c(14,21) ~ "Neutrophil",
  immune.combined$seurat_clusters %in% c(10) ~ "Macrophage",
  immune.combined$seurat_clusters %in% c(11,16,19) ~ "DC",
  immune.combined$seurat_clusters %in% c(20) ~ "Basophil",
  immune.combined$seurat_clusters %in% c(18) ~ "CD34+Kit+ Progenitor cell")
table(immune.combined$celltype)

marjor_cluster <- c('CD4 T cell','CD8 T cell', 'B cell', 'NK cell','NK T cell', 'Monocyte','Neutrophil','Macrophage','DC', 'Basophil', "CD34+Kit+ Progenitor cell","Plasma cell")
immune.combined$celltype <- factor(immune.combined$celltype, levels = marjor_cluster)

Idents(immune.combined) = 'celltype'
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4, repel=TRUE, cols=c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5', pal_npg(alpha=0.2)(2), 'grey'))+theme_bw() +
      theme(legend.position="bottom")#+ NoLegend()
ggsave(str_c('supFig2.Spleen_Anno1_SeuratClusters.',dtVar,'.pdf'), p1, width=5.5, height=5.5)
p1 <- DotPlot(immune.combined, features = c('Cd3e','Cd4','Cd8a','Cd19','Cd79a','Nkg7','Klrb1c','Ccr2','Lyz2','Csf1r','Ly6g','Fcgr3','Cxcr2','Il1b','Arg2','Adgre1','Spic','Mrc1','Itgax','Siglech','Bst2','Clec9a','Sirpa','Fcer1a','Cpa3','Fcer1g','Cst3','Kit','Cd34','Mpo','Hlf','Gata1','Klf1','Sox6','Epor','Hemgn','Tfrc','Ptprc','Sdc1','Tnfrsf17'), dot.scale = 3, col.min=-1, cols=c('white','red','red4')) + 
      RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5), panel.grid=element_blank())
ggsave(str_c('supFig2.Spleen_Anno1_SeuratClusters.dotplot_removeRedu.',dtVar,'.pdf'), p1, width=10, height=4)

# 2. Calculating Cell numbers and ratio
# source('./stat.r')
rds = immune.combined
rds2 = immune.combined
sample_labels = c('IgG_mEry','aPD1','aPD1_mEry')
colors = c(pal_nejm(alpha = 0.2)(8), '#7E6148E5', '#B09C85E5', pal_npg(alpha=0.2)(2), 'grey')
group_colors = c('#A9B8C6', '#0E5F98FF', 'red4')
cellnumber_summary(rds, rds2, sample_labels, sample_labels, marjor_cluster, colors, group_colors, str_c('supFig2.Spleen_',dtVar))

# 3. cellcycle
library(Hmisc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- c('Cenpu', capitalize(tolower(cc.genes$s.genes)))
g2m.genes <- c('Pimreg','Jpt1', capitalize(tolower(cc.genes$g2m.genes)))
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
p1 <- DimPlot(immune.combined,pt.size=0.2)#, split.by="Phase")
ggsave(str_c('supFig2.Spleen_Anno1_SeuratClusters.cellcycle.',dtVar,'.pdf'), p1, width=10, height=4)


### 保存结果
Idents(immune.combined) = 'celltype'
saveRDS(immune.combined, file = str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.rds"))


sessionInfo()