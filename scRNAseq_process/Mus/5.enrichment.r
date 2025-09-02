library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(scater)
library(celldex)
library(reshape2)
library(ggsci)
library(clusterProfiler)
library(GOplot)
library (Hmisc) 
library(org.Mm.eg.db)
library(rlist)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")
source('../Rscript/GOHelper.R') # source('aPD1_Ery_scripts/scRNAseq_process/Rscript/GOHelper.R')
source('../Rscript/GOCluster.R') # source('aPD1_Ery_scripts/scRNAseq_process/Rscript/GOCluster.R')

workpath <- "./3.combination/Spleen/" # aPD1_Ery_scripts/scRNAseq_process/Mus/
setwd(workpath)

pbmc_T = readRDS(str_c("./Spleen_t-SNE_20PCA_0.6Resolution/Spleen_t-SNE_20PCA_0.6Resolution.AnnoManual.Tsubset.rds"))
DefaultAssay(pbmc_T) <- "RNA"
Idents(pbmc_T) = 'celltype'
pbmc_T_CD8TE <- subset(pbmc_T, ident='CD8 Te')
pbmc_T_CD8TE$celltype.stim <- paste(Idents(pbmc_T_CD8TE), pbmc_T_CD8TE$sample_label, sep = "_")
Idents(pbmc_T_CD8TE) <- "celltype.stim"

df <- FindMarkers(pbmc_T_CD8TE, ident.1 = str_c('CD8 Te_aPD1_mEry'), ident.2 = NULL, verbose = FALSE, only.pos = TRUE)#
write.table(df, str_c('Tsubset/aPD1_mEryvsAllOtherCell_in_CD8Te_posDEGs.log'), sep='\t')
transID <- bitr(rownames(df[df$p_val<0.05, ]), fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Mm.eg.db')
ego <- enrichGO(gene= transID$ENTREZID, OrgDb = 'org.Mm.eg.db',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 1,readable= TRUE)
# p1 <- dotplot(ego, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 30, label_format = 70)

sample="aPD1_mEry"
markers_all <- FindMarkers(pbmc_T_CD8TE, ident.1 = str_c('CD8 Te_aPD1_mEry'), ident.2 = NULL, verbose = FALSE, only.pos = TRUE)#
pos_markers <- markers_all[(markers_all$avg_log2FC>0.25)&(markers_all$p_val<0.05),]
trans.markers <- bitr(rownames(pos_markers), fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Mm.eg.db')
ego <- enrichGO(gene= trans.markers$ENTREZID,OrgDb = org.Mm.eg.db,ont= "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 1,readable= TRUE)
tmp = ego@result
colnames(tmp) <- c('category','ID','term','GeneRatio','BgRatio','pvalue','adj_pval','qvalue','genes','count')
tmp$genes <- gsub('/', ',', tmp$genes)
pos_markers$ID <- rownames(pos_markers)
colnames(pos_markers)[2] <- 'logFC'
pos_circ = circle_dat(tmp, pos_markers)
pos_circ$genes = capitalize (tolower (pos_circ$genes)) 
pos_p1 <- dotplot(ego, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 10, label_format = 70)

markers <- rbind(pos_circ)#, neg_circ)
d3 <- head(ego$Description,15)#c(pos_p1$data$Description)#, neg_p1$data$Description)
markers$genes <- rownames(markers)
chord = chord_dat(pos_circ,pos_markers[,c('ID','logFC')],d3)
pdf('Tsubset/chord_cd8te_RvsRest_onlypos.pdf',width=15,height=7)
GOChord(chord, title="GOChord plot",#标题设置
        space = 0.02, #GO term处间隔大小设置
        gene.order = 'logFC', gene.space = 0.25, gene.size=4,#基因排序，间隔，名字大小设置
        lfc.col=c('#F16913','#FDD0A2','#FFF5EB','white'), border.size=0.01,#,##上调下调颜色设置
        ribbon.col=c(rev(colorRampPalette(brewer.pal(9, "Reds"))(15))),limit = c(0,1))#c(pal_nejm(alpha = 0.6)(5))#, pal_lancet ("lanonc",alpha = 0.6) (4))#GO term 颜色设置
dev.off()

sessionInfo()