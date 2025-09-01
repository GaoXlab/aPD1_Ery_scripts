library(reshape2)
library(dplyr)
library(plyr)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
cellnumber_summary <- function(rds, rds2, sample_labels, sample_groups, celltypes, colors, group_colors, name){
    # Idents(rds) <- factor(rds$celltype, levels=celltypes)
    cell_num_inclusters <- c()
    cellnumber <- c()
    cellratio_in_sample <- c()
    sample_label_list <- names(table(rds$sample_label))
    for(clusters_i in celltypes){
    	cellnum_in_cluster <- WhichCells(rds, idents = clusters_i)
    	sample_label_cell_num <- c()
    	cellratio_in_sample_label <- c()
    	for(sample_label_i in sample_label_list){
            sample_label_cell_num <- c(sample_label_cell_num, length(grep(str_c('^', sample_label_i, '_[A-Z0-9]'), cellnum_in_cluster)))
    		cellratio_in_sample_label <- c(cellratio_in_sample_label, length(grep(str_c('^', sample_label_i, '_[A-Z0-9]'), cellnum_in_cluster))/table(rds2$sample_label)[sample_label_i])
      }
    	cellnumber <- rbind(cellnumber, sample_label_cell_num)
    	cellratio_in_sample <- rbind(cellratio_in_sample, cellratio_in_sample_label)
    	cell_num_inclusters <- c(cell_num_inclusters, length(cellnum_in_cluster))
    }
    cells_statistics <- as.data.frame(cbind(cellnumber, cellnumber/cell_num_inclusters, cellratio_in_sample))
    colnames(cells_statistics) <- c(str_c('CellNum',sample_label_list), str_c('CellRatio_in_Clusters',sample_label_list), str_c('CellRatio_in_Samples',sample_label_list))
    rownames(cells_statistics) <- celltypes
    write.table(cells_statistics, str_c(name, 'cells_statistics.log'), sep='\t')
    # cells_statistics <- read.table('BoneMarrow_cellsratio.csv', sep=',')
    df3 <- melt(cbind(rownames(cells_statistics), cells_statistics[grep('CellRatio_in_Samples', colnames(cells_statistics))]))
    df3$label <- gsub("CellRatio_in_Samples", "", as.character(df3$variable))
    colnames(df3) <- c('celltype', 'variable', 'value', 'label')
    df3$label <- factor(df3$label, levels=sample_labels)
    df3$celltype <- factor(df3$celltype, levels=celltypes)
    df3$value100 <- df3$value*100
    df3$group <- gsub('_1$|_2$|_3$|_4$','', df3$label)
    df3$group <- factor(df3$group, sample_groups)
    p1 <- ggplot(data=df3, aes(x=label, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+facet_grid(.~celltype)+
          theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5))+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='bottom')
    p2 <- ggplot(data=df3, aes(x=label, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='right')
    p3 <- ggplot(data=df3, aes(x=group, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+ylab('Percentage (%)')+
          scale_fill_manual(values=colors) + theme(legend.position='right')
    print(head(df3))
    df3_trans <- data_summary(df3, varname="value100", groupnames=c("celltype", "group"))
    # Convert dose to a factor variable
    df3_trans$celltype=as.factor(df3_trans$celltype)
    print(head(df3_trans))
    g <- plot_grid(p1, p2, p3, rel_widths=c(1,1,0.4), ncol=3)
    ggsave(str_c("./", name, "cellsratio_InAnnotationClusters.barplot2.pdf"), g, width=16, height=5)
    # colors=c('grey','#4c85bb','red','#FB6A4A')
    p4 <- ggplot(df3_trans, aes(x=celltype, y=value100, fill=group)) + 
      geom_bar(stat="identity", position=position_dodge()) +
      # geom_errorbar(aes(ymin=ifelse(value100-sd<0, 0 ,value100-sd ), ymax=value100+sd), width=.2, position=position_dodge(.9))+
      scale_fill_manual(values=group_colors)+ theme(legend.position='bottom')+
      theme_classic()+theme(legend.position='right', axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + ylab('Percentage (%)')
    # p4 <- ggplot(data=df3, aes(x=group, y=value100, fill=celltype)) + geom_bar(stat="identity")+theme_classic()+ylab('Percentage (%)')+
    #       scale_fill_manual(values=colors) + theme(legend.position='right')
    ggsave(str_c("./", name, "cellsratio_InAnnotationClusters.barplot3.pdf"), p4, width=10, height=4)
}


calculate_roe <- function(calculate_df, subset, total){
    expected_tmp_cells = c()
    for(cl in names(calculate_df)){
        expected_tmp_cells =  rbind(expected_tmp_cells, total_cells*as.numeric(calculate_df[cl]) / as.numeric(total))
    }
    expected_tmp_cells = as.data.frame(expected_tmp_cells)
    rownames(expected_tmp_cells) = names(calculate_df)
    Roe_tmp = as.data.frame(cbind(subset[,'CellNumIgG_mEry']/expected_tmp_cells[,'CellNumIgG_mEry'], 
                subset[,'CellNumaPD1']/expected_tmp_cells[,'CellNumaPD1'],
                subset[,'CellNumaPD1_mEry']/expected_tmp_cells[,'CellNumaPD1_mEry']))
    rownames(Roe_tmp) = names(calculate_df)
    colnames(Roe_tmp) = names(total_cells)
    Roe_tmp$cl = rownames(Roe_tmp)
    #### 测试chisq.test
    chisq_tmp_result = c()
    for(cl in names(calculate_df)){
        ck = as.data.frame(rbind(subset[cl, c('CellNumIgG_mEry','CellNumaPD1')],  round(expected_tmp_cells[cl, c('CellNumIgG_mEry','CellNumaPD1')],0))); 
        kr = as.data.frame(rbind(subset[cl, c('CellNumaPD1','CellNumaPD1_mEry')],  round(expected_tmp_cells[cl, c('CellNumaPD1','CellNumaPD1_mEry')],0))); 
        cr = as.data.frame(rbind(subset[cl, c('CellNumIgG_mEry','CellNumaPD1_mEry')],  round(expected_tmp_cells[cl, c('CellNumIgG_mEry','CellNumaPD1_mEry')],0))); 
        chisq_tmp_result = rbind(chisq_tmp_result, c(chisq.test(ck)$p.value, chisq.test(kr)$p.value, chisq.test(cr)$p.value))
    }
    chisq_tmp_result = as.data.frame(round(chisq_tmp_result,6)); colnames(chisq_tmp_result) = c('CvK','KvR','CvR')
    chisq_tmp_result$cl = names(calculate_df)
    Roe_tmp = merge(Roe_tmp, chisq_tmp_result, by='cl')
    return(Roe_tmp)
}