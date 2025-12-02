# Paper Related Code and Data

## scRNAseq scipts  
### 1. preprocess (cellranger to QC)
sc_datapreprocess.human.sh  
sc_datapreprocess.mus.sh
### 2. cell type annotation and differential analysis (R version 4.3.1)
Enter either directory to run the pipeline; the results will be stored in the corresponding output directory.
- ./Human
- ./Mus

## cytof scripts
/usr/local/bin/R CMD BATCH cytof_process/1.cytof_process.r
/usr/local/bin/R CMD BATCH cytof_process/2.calculate_Roe_cytof.r
