#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(Seurat)
library(stringr)
library(DoubletFinder)
library(cowplot)
library(ggplot2)

option_list = list(
    make_option(c("-w", "--work_path"), type = "character", default = "/", help = "File path (path of output dataset)"),
    make_option(c("-l", "--sample_label_list"), type = "character", default = "",  help = "The list of sample_label using for every sample"),
    make_option(c("-f", "--input_file_list"),  type = "character", default = "xx_raw.rds",  help = "Path and full name of inputfile for disease"),
    make_option(c("-u", "--nFeature_RNA_upper_limit_list"), type = "character", default = "",  help = "Upper limit of nFeature_RNA"),
    make_option(c("-m", "--mt_upper_limit_list"), type = "character", default = "",  help = "Upper limit of percent.mt"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "Combined", help = "Prefix of output filename")
    )
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

sample_label_list <- as.character(opt$sample_label_list)
inputfile_list <- as.character(opt$input_file_list)
workpath <- as.character(opt$work_path)
output_prename <- as.character(opt$output_prefix)
nFeature_RNA_upper_limit_list <- as.character(opt$nFeature_RNA_upper_limit_list)
percent.mt_upper_limit_list <- as.character(opt$mt_upper_limit_list)

# 反馈读参结果
print(str_c("-w: ", workpath))
print(str_c("-o: ", output_prename))
print(str_c("-f: ", inputfile_list))
print(str_c("-u: ", nFeature_RNA_upper_limit_list))
print(str_c("-m: ", percent.mt_upper_limit_list))


############################################################
# 函数： filter_object 
single_sample <- function(inputfile, sample_label, nFeature_RNA_upper_limit, percent.mt_upper_limit){
	################# object ##################
	# Load the raw rds dataset
	pbmc <- readRDS(str_c(inputfile))
	pbmc$sample_label <- sample_label
	Idents(pbmc) <- sample_label
	pbmc@project.name <- sample_label
	pbmc@ meta.data$ orig.ident <- sample_label
	# Screening cells 
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_RNA_upper_limit & percent.mt < percent.mt_upper_limit)
	# print(str_c("### The original barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- RenameCells(object = pbmc, add.cell.id = sample_label)
	# print(str_c("### The new barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- NormalizeData(pbmc, verbose = FALSE)
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	# print("The pbmc object info: ")
	return(pbmc)
}

single_sample_dedoublet <- function(inputfile, sample_label, nFeature_RNA_upper_limit, percent.mt_upper_limit){
	################# object ##################
	# Load the raw rds dataset
	pbmc <- readRDS(str_c(inputfile))
	pbmc$sample_label <- sample_label
	Idents(pbmc) <- sample_label
	pbmc@project.name <- sample_label
	pbmc@ meta.data$ orig.ident <- sample_label
	# Screening cells (这之后细胞总个数不会再发生改变)
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_RNA_upper_limit & percent.mt < percent.mt_upper_limit)
	# print(str_c("### The original barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- RenameCells(object = pbmc, add.cell.id = sample_label)
	# print(str_c("### The new barcodes of pbmc object:"))
	# head(x = colnames(x = pbmc))
	pbmc <- NormalizeData(pbmc, verbose = FALSE)
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	# print("The pbmc object info: ")
    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc)
    pbmc <- RunTSNE(pbmc, dims = 1:20)

    multiplet_rate <- as.data.frame(rbind(c(0.004,800),c(0.008,1600),c(0.016,3200),c(0.023,4800),c(0.031,6400),c(0.039,8000),c(0.046,9600),c(0.054,11200),c(0.061,12800),c(0.069,14400),c(0.076,16000)))
    colnames(multiplet_rate) <- c('rate','cells')

    rate = multiplet_rate [(which(multiplet_rate$cells>ncol(pbmc))[1]-1), 'rate']

    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- pbmc@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- pbmc@meta.data$ClusteringResults
    nExp_poi <- round(rate*nrow(pbmc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	# filer_DFname = colnames(pbmc@meta.data)[grep('^pANN_', colnames(pbmc@meta.data))]
    # pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = filer_DFname, sct = FALSE)
    name = colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]
    p0 <- DimPlot(pbmc, reduction = "tsne", group.by = name, pt.size = 0.9)

    #通过小提琴图来展示结果
    p1 <- VlnPlot(pbmc, features = c("nCount_RNA"), group.by = name, log = T)
    p2 <- VlnPlot(pbmc, features = c("nFeature_RNA"), group.by = name, log = F)

    g <- plot_grid(p0,p1,p2,ncol=3,rel_widths=c(1.3,1,1))
    ggsave(str_c(sample_label,'.doublet.pdf'), g, width=12, height=4)

    Idents(pbmc) = name
    pbmc = subset(pbmc, idents='Singlet') 

    return(pbmc)
}

outputpath <- paste0(workpath, "/", output_prename, collapse = "")


################################################################
############################ Begin #############################
################################################################
setwd(workpath)
################# object input ##################
pbmc_list <- c()
for(i in 1:(length(strsplit(inputfile_list, split=",")[[1]]))){
	pbmc <- single_sample_dedoublet(strsplit(inputfile_list, split=",")[[1]][i], as.character(strsplit(sample_label_list, split=",")[[1]][i]), as.numeric(strsplit(nFeature_RNA_upper_limit_list, split=",")[[1]][i]), as.numeric(strsplit(percent.mt_upper_limit_list, split=",")[[1]][i]))
	print(str_c("===== pbmc ", i, " object is input! ====="))
	print(as.character(strsplit(sample_label_list, split=",")[[1]][i]))
	pbmc_list <- c(pbmc_list, pbmc)
}


################## Standard Workflow ###################
# Run the Standard workflow 
pbmc_remain <- pbmc_list
pbmc_remain[[1]] <- NULL
cca_merge <- merge(x = pbmc_list[[1]], y = pbmc_remain)
merge.list <- SplitObject(cca_merge, split.by = "sample_label")
merge.list <- lapply(X = merge.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merge.list)
# perform intergration
immune.anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features, dims = 1:10)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:10)
DefaultAssay(immune.combined) <- "integrated"

print("The immune.combined object info:")
print(immune.combined)
saveRDS(immune.combined, file = str_c(outputpath,"_Origin_Integrated.rds"))
print(str_c(output_prename, " origin rds has saved!"))

# visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

pdf(str_c(outputpath,"_3.1.DimPlot.pdf")) 
DimPlot(immune.combined, reduction = "pca")
dev.off()
pdf(str_c(outputpath,"_3.2.DimHeatmap.pdf"))
DimHeatmap(immune.combined, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()
pdf(str_c(outputpath,"_3.3.ElbowPlot.pdf"))
ElbowPlot(immune.combined)
dev.off()


sessionInfo()