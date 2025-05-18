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
library(ggsci)
library(stringr)
library(RColorBrewer)
library(scRepertoire)
library(scales)
library(org.Hs.eg.db)
library(rlist)
library(clusterProfiler)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

source('aPD1_Ery_scripts/scRNAseq_process/Human/stat.r')

workpath <- "aPD1_Ery_scripts/scRNAseq_process/Human/"
setwd(workpath)

# 1. subset myeloid 
immune.combined <- readRDS(str_c("./P1_3_t-SNE_30PCA_0.6Resolution/P1_3_t-SNE_30PCA_0.6Resolution.AnnoManual.rds"))
Idents(immune.combined) <- immune.combined$celltype 
DefaultAssay(immune.combined) <- "integrated"
APC <- subset(immune.combined, idents=c('Monocyte','cDC','pDC'))
DefaultAssay(APC) <- "integrated"
# Run the standard workflow for visualization and clustering
APC <- ScaleData(APC, verbose = FALSE)
APC <- RunPCA(APC, npcs = 30, verbose = FALSE)
APC <- RunTSNE(APC, reduction = "pca", dims = 1:30)
APC <- FindNeighbors(APC, reduction = "pca", dims = 1:30)
APC <- FindClusters(APC, resolution = 0.5)
p2 <- DimPlot(APC, reduction = "tsne", label = TRUE, repel = TRUE) + theme_bw()
ggsave(str_c('./supFig.P1_3_Myeloid_SeuratClusters.',dtVar,'.pdf'), p2, width=5.5, height=4.5)

# 2. annotation 
DefaultAssay(APC) <- "RNA"
Idents(APC) <- factor(APC$seurat_clusters, levels=c(12,10,8,4,14,11,7,5,13,3,6,2,1,0,9))
APCmajor_cluster <- c('IGKC','GZMB','JCHAIN','LILRA4','TCF4','IRF4','MZB1','XBP1','CD1C','FCER1A','CLEC9A','CST3','CD86','CD163','PPBP','PF4','MYL9','MARCO','CD3G','CD8A','TRDC','NKG7','KLRB1','ITGAM','ITGAX','FCGR3A','CD68', 'IFITM2','SIGLEC10','LILRB1','CX3CR1','CD14','FUT4','FCN1','VCAN','S100A8','S100A9','STAT1','NFKB1','CD63','MSR1', 'CFS1', 'CCL2', 'CCR2','TGFB1','TGFBR2','TGFBR1','TGFBR3','RSAD2','ISG15','IRF7','AIF1','APOE','HLA-DRA','CD74','CTSB','CTSD','FCGR2B','CEACAM8','CD274','NOS2','IL1B','TCL1A','CD79A','CD79B')
p1 <- DotPlot(APC, features = APCmajor_cluster, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave(str_c('./supFig.P1_3_Myeloid_SeuratClusters.dotplot.',dtVar,'.pdf'), p1, width=12.5, height=4.5)

APC$celltype = dplyr::case_when(
  APC$seurat_clusters %in% c(9) ~ "B",
  APC$seurat_clusters %in% c(3) ~ "CD16+ Mono",
  APC$seurat_clusters %in% c(6,13) ~ "CD14+ NFKB1+ Mono",
  APC$seurat_clusters %in% c(0,1,2) ~ "CD14+ APC",
  APC$seurat_clusters %in% c(10,12,8) ~ "pDC",
  APC$seurat_clusters %in% c(4,14) ~ "cDC",  
  APC$seurat_clusters %in% c(11) ~ "Platelet",
  APC$seurat_clusters %in% c(5,7) ~ "NK")
APC_celltype_levels <- c('pDC','cDC','NK','Platelet',"B","CD16+ Mono","CD14+ APC","CD14+ NFKB1+ Mono")
APC_celltype_colors <- c((pal_gsea("default", n = 11, alpha = 0.2, reverse = TRUE)(11)[1:4]), "#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#08519C")
APC$celltype <- factor(APC$celltype, levels=APC_celltype_levels)
Idents(APC) <- APC$celltype 
DefaultAssay(APC) <- "RNA"
p1 <- DotPlot(APC, features = APCmajor_cluster, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) 
ggsave('APC_markers.dotplot.pdf', p1, width=11,height=5)
APC$celltype <- factor(APC$celltype, levels=APC_celltype_levels)
Idents(APC) <- factor(APC$celltype)
p1 <- DimPlot(APC, reduction = "tsne", label = TRUE, repel=TRUE, label.size = 3, cols=APC_celltype_colors) + NoLegend() + theme_bw()
ggsave('APC.pdf', p1, width=5.5, height=3.5)
p1 <- DimPlot(APC, reduction = "tsne", label = FALSE, label.size = 2.5, repel=FALSE, cols=APC_celltype_colors)+theme_bw() +
      theme(legend.position="bottom", axis.text=element_blank())+ NoLegend()+ylab('')+xlab('')
ggsave(str_c('APC.',dtVar,'.no.pdf'), p1, width=3, height=3)

p1 <- DimPlot(APC, reduction = "tsne", label = TRUE, repel=TRUE, label.size = 3, cols=c(pal_gsea("default", n = 11, alpha = 0.4, reverse = TRUE)(11)[1:4], "#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#08519C")) + NoLegend() + theme_bw()
ggsave('APC_legend.pdf', p1, width=5.5, height=3.5)


### 3. MDSC, markers from (PMID: 33293583)
DefaultAssay(APC) <- "RNA"
mdsc_gene <- c('CCR7','CD1A','CD1B','CD1C','CD207','CD209','CD4','CD40','CD80','CD83','CD86','CMKLR1','HLA-DOA','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-DRB6','ITGA4','ITGAX','ITGAM','LY75','NPR1','PDCDILG2')
gene <- list()
gene$gene <- mdsc_gene#il18_mediated_signaling_pathway_mm10
APC <- AddModuleScore(object=APC, features=gene, ctrl=100, name='CD_features')
colnames(APC@meta.data)[length(colnames(APC@meta.data))] <- 'MDSC'
Idents(APC) = APC$celltype
# Idents(APC) = factor(APC$seurat_clusters, levels=c(12,6,13,3,9,4,11,5,8,14,2,10,7,1,0))
a <- VlnPlot(APC, features='MDSC', pt.size=0.03)
Idents(APC) = "sample_label"
APC_patient = subset(APC, idents=c("P1_3_baseline","P1_3_C3D1"))
Idents(APC_patient) = "celltype"
APC_mono = subset(APC_patient, idents=c("CD16+ Mono","CD14+ APC","CD14+ NFKB1+ Mono"))#,8,4,13,0,1,3,6,14
Idents(APC_mono) = factor(APC_mono$celltype, levels=c("CD14+ APC","CD16+ Mono","CD14+ NFKB1+ Mono"))
a_cd14 <- VlnPlot(APC_mono, features='MDSC', pt.size=0)+ scale_fill_manual(values=c("#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#08519C")) #scale_fill_brewer(palette="Blues")
ggsave('MDSC_scoring.pdf',a_cd14, width=4.5, height=4)


### 4. cellcycle
if(all(colnames(APC)==colnames(immune.combined$Phase[colnames(APC)]))){
  print('same')
  APC$Phase = immune.combined$Phase[colnames(APC)]
  APC$S.Score = immune.combined$S.Score[colnames(APC)]
  APC$G2M.Score = immune.combined$G2M.Score[colnames(APC)]
}
### 保存结果
saveRDS(APC, file = str_c("./P1_3_t-SNE_30PCA_0.6Resolution/P1_3_t-SNE_30PCA_0.6Resolution.Myeloid.AnnoManual.rds"))


### 5. cell population summary
APC <- readRDS(str_c("./P1_3_t-SNE_30PCA_0.6Resolution/P1_3_t-SNE_30PCA_0.6Resolution.Myeloid.AnnoManual.rds"))
cellnumber_summary(APC, names(table(APC$sample_label)), APC_celltype_levels, APC_celltype_colors, 'APC')
