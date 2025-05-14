library(stringr)
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

cytof = read.table('./cytof_cellsinsample.csv', sep=',', head=TRUE, fill=TRUE)

total_cells = apply(cytof, 2, sum)
total = sum(total_cells)
total_clusters = apply(cytof, 1, sum)
expected_cells = c()
for(cl in names(total_clusters)){
   expected_cells =  rbind(expected_cells, total_cells*as.numeric(total_clusters[cl]) / as.numeric(total))
}
expected_cells = as.data.frame(round(expected_cells),0)
rownames(expected_cells) = names(total_clusters)

Roe = as.data.frame(cbind(cytof[,'IgG_mEry']/expected_cells[,'IgG_mEry'], 
            cytof[,'aPD1']/expected_cells[,'aPD1'],
            cytof[,'aPD1_mEry']/expected_cells[,'aPD1_mEry']))
rownames(Roe) = names(total_clusters)
colnames(Roe) = names(total_cells)


CD4T = apply(cytof[grep('CD4|Treg', rownames(cytof)), ], 2, sum); CD4T_sum = sum(CD4T)
CD8T = apply(cytof[grep('CD8', rownames(cytof)), ], 2, sum); CD8T_sum = sum(CD8T)
MDSC = apply(cytof[grep('MDSC', rownames(cytof)), ], 2, sum); MDSC_sum = sum(MDSC)
Mye = apply(cytof[grep('RP Macrophage|Monocyte|DC|MDSC', rownames(cytof)), ], 2, sum); Mye_sum = sum(Mye)
cytof_inLarge = as.data.frame(rbind(CD4T, CD8T, MDSC, Mye)); rownames(cytof_inLarge) <- c('CD4T_sum', 'CD8T_sum', 'MDSC_sum', 'Mye_sum')
total_clusters_inLarge = data.frame(CD4T_sum, CD8T_sum, MDSC_sum, Mye_sum)
expected_cells_inLarge = c()
for(cl in names(total_clusters_inLarge)){
   expected_cells_inLarge =  rbind(expected_cells_inLarge, total_cells*as.numeric(total_clusters_inLarge[cl]) / as.numeric(total))
}
expected_cells_inLarge = as.data.frame(round(expected_cells_inLarge),0)
rownames(expected_cells_inLarge) = names(total_clusters_inLarge)

Roe_inLarge = as.data.frame(cbind(cytof_inLarge[,'IgG_mEry']/expected_cells_inLarge[,'IgG_mEry'], 
            cytof_inLarge[,'aPD1']/expected_cells_inLarge[,'aPD1'],
            cytof_inLarge[,'aPD1_mEry']/expected_cells_inLarge[,'aPD1_mEry']))
rownames(Roe_inLarge) = names(total_clusters_inLarge)
colnames(Roe_inLarge) = names(total_cells)


expected_cells_1 = as.data.frame(rbind(expected_cells, expected_cells_inLarge)); colnames(expected_cells_1) <- str_c('expected_',colnames(expected_cells))
cytof_1 = as.data.frame(rbind(cytof, cytof_inLarge)); colnames(cytof_1) <- str_c('real_',colnames(cytof))
Roe_1 = as.data.frame(rbind(Roe, Roe_inLarge)); colnames(Roe_1) <- str_c('Roe_',colnames(Roe_1))
all(rownames(expected_cells_1)==rownames(Roe_1))
df_2_print <- cbind(cytof_1, expected_cells_1, Roe_1)

options(digits = 2)
#### 测试chisq.test
chisq_result_KvR = c()
chisq_result_KvC = c()
chisq_result_RvC = c()
total_clusters_names = rownames(df_2_print)
colnames(cytof_1) = gsub('real_','',colnames(cytof_1))
colnames(expected_cells_1) = gsub('expected_','',colnames(expected_cells_1))
for(cl in total_clusters_names){
    mydata = as.data.frame(rbind(cytof_1[cl, c('aPD1_mEry','aPD1')],  expected_cells_1[cl, c('aPD1_mEry','aPD1')])); 
    chisq_result_KvR = c(chisq_result_KvR, chisq.test(mydata)$p.value)
    mydata = as.data.frame(rbind(cytof_1[cl, c('IgG_mEry','aPD1')],  expected_cells_1[cl, c('IgG_mEry','aPD1')])); 
    chisq_result_KvC = c(chisq_result_KvC, chisq.test(mydata)$p.value)
    mydata = as.data.frame(rbind(cytof_1[cl, c('IgG_mEry','aPD1_mEry')],  expected_cells_1[cl, c('IgG_mEry','aPD1_mEry')])); 
    chisq_result_RvC = c(chisq_result_RvC, chisq.test(mydata)$p.value)
}
chisq_result_KvR = as.data.frame(chisq_result_KvR)
chisq_result_KvR$cl = total_clusters_names
symp <- symnum(chisq_result_KvR$chisq_result_KvR, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_result_KvR$Signif_KvR = symp

chisq_result_KvC = as.data.frame(chisq_result_KvC)
chisq_result_KvC$cl = total_clusters_names
symp <- symnum(chisq_result_KvC$chisq_result_KvC, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_result_KvC$Signif_KvC = symp

chisq_result_RvC = as.data.frame(chisq_result_RvC)
chisq_result_RvC$cl = total_clusters_names
symp <- symnum(chisq_result_RvC$chisq_result_RvC, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","."," "))
chisq_result_RvC$Signif_RvC = symp

chisq_result_m1 <- merge(chisq_result_KvR, chisq_result_KvC, by='cl')
chisq_result_merge <- merge(chisq_result_m1, chisq_result_RvC, by='cl')

df_2_print$cl <- rownames(df_2_print)
df_2_print_merge <- merge(df_2_print, chisq_result_merge, by='cl')

write.table(df_2_print_merge,str_c('./chisq_for_cellpopulation.Roe_for_plot.',dtVar,'.log'), sep='\t', row.names=FALSE)
