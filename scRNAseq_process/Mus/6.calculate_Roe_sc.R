library(ggplot2)
library(reshape2)
source('aPD1_Ery_scripts/scRNAseq_process/Mus/stat.r')

Tsubset = read.table("aPD1_Ery_scripts/scRNAseq_process/Mus/Tsubset/supFig2.Spleen_inTotal_2024-08-30cells_statistics.log", sep='\t', head=T)
Tsubset = Tsubset[,c('CellNumIgG_mEry','CellNumaPD1','CellNumaPD1_mEry')]
Myesubset = read.table("aPD1_Ery_scripts/scRNAseq_process/Mus/Myeloidsubset/supFig2.Spleen_inTotal_2024-08-30cells_statistics.log", sep='\t', head=T)
Myesubset = Myesubset[,c('CellNumIgG_mEry','CellNumaPD1','CellNumaPD1_mEry')]
total_cell = read.table("aPD1_Ery_scripts/scRNAseq_process/Mus/supFig2.Spleen_2024-08-30cells_statistics.log", sep='\t', head=T)
total_cell = total_cell[,c('CellNumIgG_mEry','CellNumaPD1','CellNumaPD1_mEry')]

total_cells = apply(total_cell, 2, sum); total = sum(total_cells)
total_clusters = apply(total_cell, 1, sum)
total_T_clusters = apply(Tsubset, 1, sum)
total_Mye_clusters = apply(Myesubset, 1, sum)

chisq_T_result2 = calculate_roe(total_T_clusters, Tsubset, total)
chisq_Mye_result2 = calculate_roe(total_Mye_clusters, Myesubset, total)
chisq_Total_result2 = calculate_roe(total_clusters, total_cell, total)

symp <- symnum(chisq_T_result2$CvK, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_T_result2$Signif_CvK = symp
symp <- symnum(chisq_T_result2$KvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_T_result2$Signif_KvR = symp
symp <- symnum(chisq_T_result2$CvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_T_result2$Signif_CvR = symp

symp <- symnum(chisq_Mye_result2$CvK, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Mye_result2$Signif_CvK = symp
symp <- symnum(chisq_Mye_result2$KvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Mye_result2$Signif_KvR = symp
symp <- symnum(chisq_Mye_result2$CvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Mye_result2$Signif_CvR = symp

symp <- symnum(chisq_Total_result2$CvK, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Total_result2$Signif_CvK = symp
symp <- symnum(chisq_Total_result2$KvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Total_result2$Signif_KvR = symp
symp <- symnum(chisq_Total_result2$CvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_Total_result2$Signif_CvR = symp
library(openxlsx)
write.xlsx(list('chisq_T'=chisq_T_result2,'chisq_Mye'=chisq_Mye_result2,'chisq_Total'=chisq_Total_result2), 'Roe.xlsx' ,sep='\t')
