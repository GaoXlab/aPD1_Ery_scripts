library(ggplot2)
library(flowCore)
library(diffcyt)
library(HDCytoData)
library(CATALYST)
library(cytofWorkflow)
library(readxl)
library(stringr)
library(ggthemes)
library(ggsci)
library(ggrepel)
library(data.table)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(scales)
library(dplyr)
packinfo <- installed.packages(fields = c("Package", "Version"))
for(pk in c('flowCore','diffcyt','HDCytoData','CATALYST','cytofWorkflow')){
    print(packinfo[pk,c("Package", "Version")])
}
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

set.seed(2024)
setwd('./cytof_process/')
panel_plt_filename <- "CyTOF_panel.xlsx"

p1 = './Data/' ### rawdta
fs1=list.files(p1,'*fcs' )
sample <- read.flowSet(files = fs1,path = p1)

md_plt <- "Spleen_metadata.xlsx"
md_plt <- read_excel(str_c('./', md_plt))

panel_plt <- read_excel(str_c('./', panel_plt_filename))
panel_plt$Antigen <- gsub("/", "_", panel_plt$Antigen)
panel_plt$Antigen <- gsub("-", "_", panel_plt$Antigen)
panel_plt$fcs_colname <- str_c(panel_plt$Metal,panel_plt$Isotope,'Di')
all(panel_plt$fcs_colname %in% colnames(sample))
panel_plt <- panel_plt[panel_plt$fcs_colname %in% intersect(panel_plt$fcs_colname, colnames(sample)),]
panel_plt$marker_class <- "type"
panel_plt <- panel_plt[, c("fcs_colname", "Antigen", "marker_class" )]
colnames(panel_plt) <- c("fcs_colname", "antigen", "marker_class")

sce <- prepData(sample, panel_plt, md_plt, features = panel_plt$fcs_colname)
sce$condition <- factor(sce$condition, levels=c('IgG_mEry', 'aPD1', 'aPD1_mEry'))
sce$sample_id <- factor(sce$sample_id, levels=c("IgG_mEry", "aPD1", "aPD1_mEry"))

### 
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
ggsave(str_c('supFig1.Cytof_marksDensity.',dtVar,'.pdf'), p, width=10,height=7)

### We call ConsensusClusterPlus() with maximum number of clusters maxK = 40.
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 40, seed = 1234)
### reference markers: 
T <- c('CD45','CD3','TCR_b chain','CD4','CD8','CD44','CD183','CD62L','CD127','CD107a(LAMP_1)','CD69','CD25','FoxP3','CTLA4','ICOS','Granzyme B','CTLA4','TIM3','LAG3','TIGIT(VSTM3)')
B <- c('CD19.2','Gd158Di','CD40','CD38','CD185','ST2','CD196')
NK <- c('NK1.1')
MONO <- c('CD11b')
DC <- c('CD11c','MHC II','CD103')
MACROPHAGES <- c('F4_80','CD206','CD163')
M_MDSC <- c('Ly6C')
G_MDSC <- c('Ly6G', 'Gr_1')
APC <- c('CD80','CD86','PD_L1')
nouse <- c('PD1')
markers <- c(T, B, NK, MONO, DC, MACROPHAGES, M_MDSC, G_MDSC, APC, nouse)
order_id <- c()
for(i in markers){
    order_id <- c(order_id, which(rownames(sce)==i))
}
cluste_levels = c(13,5,11,31,12,9,8,17,25,32,33,14,18,37,26,21,10,6,7,2,3,15,19,22,23,27,28,34,38,30,35,39,20,40,36,4,29,1,16,24)
sce@ metadata$ cluster_codes$ meta40 = factor(sce@ metadata$ cluster_codes$ meta40, levels=cluste_levels)

a=plotExprHeatmap(sce, 
                by = "cluster_id", k = "meta40",
                row_clust = FALSE,
                col_clust = FALSE,
                # bars = TRUE, perc = TRUE)
                features = rownames(sce)[order_id],
                bars = TRUE, perc = TRUE)
pdf(str_c('SupFig2.cluster_plotExprHeatmap_row_col_F.',dtVar,'.pdf'), width=14,height=7)
a
dev.off()

#####  改n-cells显示log10(n_cells)
source('Rscript/utils-differential.R') ### 改的是.anno_counts中的rowAnnotation
source('Rscript/plotExprHeatmap.R')
source('Rscript/validityChecks.R')
source('Rscript/utils-differential.R')
source('Rscript/plotCounts.R')
a=plotExprHeatmap(sce, 
                by = "cluster_id", k = "meta40",
                row_clust = FALSE,
                col_clust = FALSE,fun='median',scale="first",
                features = rownames(sce)[order_id],
                bars = TRUE, perc = TRUE)

pdf(str_c('SupFig2.cluster_plotExprHeatmap_row_col_F.',dtVar,'.paperUsed.log10.pdf'), width=14,height=7)
a
dev.off()


p <- plotAbundances(sce, k = "meta40", by = "cluster_id", shape_by = "patient_id") #scale_fill_manual(values=c('grey', 'steelblue1', 'red'))
ggsave(str_c('SupFig3.plotAbundances_boxplot.',dtVar,'.pdf'), p, width=15,height=13)

df_all <- p$data
df_all$patient_id <- factor(df_all$patient_id, levels=c('M1', 'M2', 'M3'))
write.table(df_all[order(df_all$cluster_id, df_all$patient_id), ], str_c('proportion.',dtVar,'log'), sep='\t')



sce <- runDR(sce, "TSNE",features = "type") 
p1 <- plotDR(sce, "TSNE", color_by = "meta40") + theme_few() + xlab('')+ylab('')+ggplot2::theme(plot.title=element_blank())
ggsave( str_c('SupFig4.ClustersTSNE.',dtVar,'.pdf'), p1, width=5.6, height=3.8)

p1_data <- p1$data
p1_data$celltype = dplyr::case_when(
  p1_data$meta40 %in% c(21,10,6,7,2,3,15,19,22,23) ~ "B",
  p1_data$meta40 %in% c(13) ~ "CD4 Naïve",
  p1_data$meta40 %in% c(5) ~ "CD4 Tcm",
  p1_data$meta40 %in% c(11,31) ~ "CD4 Tem",
  p1_data$meta40 %in% c(12) ~ "Central Treg",
  p1_data$meta40 %in% c(9) ~ "Effector Treg",
  p1_data$meta40 %in% c(8) ~ "CD4 Te",
  p1_data$meta40 %in% c(17) ~ "DP T",
  p1_data$meta40 %in% c(25,32) ~ "CD8 Naïve",
  p1_data$meta40 %in% c(33) ~ "CD8 Tem",
  p1_data$meta40 %in% c(14,18,37) ~ "CD8 Tcm",
  p1_data$meta40 %in% c(26) ~ "CD8 Te",
  p1_data$meta40 %in% c(27,28,34) ~ "NK",
  p1_data$meta40 %in% c(30,38) ~ "DC",
  p1_data$meta40 %in% c(35,39,20) ~ "Monocyte",
  p1_data$meta40 %in% c(40) ~ "PMN-MDSC",
  p1_data$meta40 %in% c(36,4) ~ "M-MDSC",
  p1_data$meta40 %in% c(29) ~ "RP Macrophage",
  p1_data$meta40 %in% c(1,16,24) ~ "Unknown")
table(p1_data$celltype)
p1_data$celltype=factor(p1_data$celltype, levels=c("B","CD4 Naïve","CD4 Te","CD4 Tem","CD4 Tcm","Central Treg","Effector Treg","CD8 Naïve","CD8 Te","CD8 Tem","CD8 Tcm","RP Macrophage","Monocyte","DC","M-MDSC","PMN-MDSC","NK","Unknown","DP T"))
mini_mpg <- c(); for(i in names(table(p1_data$celltype))){ mini_mpg <- rbind(mini_mpg, head(p1_data[p1_data$celltype==i,],1))}
p <- ggplot(p1_data, aes(x, y, color = celltype))+geom_point(size=0.03)+
    scale_color_manual(values=c(rev(pal_nejm(alpha = 0.04)(7)), pal_npg("nrc", alpha = 0.04)(9), '#7E6148E5', '#B09C85E5', 'grey'))+theme_bw()+
    ggrepel::geom_text_repel(data = mini_mpg, aes(label = celltype), colour='black', size=3.5)
ggsave(str_c('Fig1.Anno_ClustersTSNE.',dtVar,'.pdf'), p, width=5.6, height=3.8)
write.table(table(p1_data$celltype,p1_data$sample_id),str_c('cytof_cellsinsample',dtVar,'.csv'),sep=',')



sessionInfo()
