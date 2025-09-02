#!/bin/bash
#SBATCH -J BM4
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH -c 40
#SBATCH --mem=400G

module load cellranger/6.1.1
module load R/4.2.1
module load gcc/11.2.0
module load fftw/3.3.10

### 1.preprocess
### for only single-cell RNA-seq:
rawdata_path="$(pwd)"
sample="IgG_mEry aPD1 aPD1_mEry"
for s in $sample
do
cellranger count --id=$s --fastqs=$rawdata_path/$s --sample=$s --transcriptome=refdata-gex-mm10-2020-A/ --localcores=40
done

### 2.quality control (single sample)
R_script="$(pwd)/Rscript/"
sample="IgG_mEry aPD1 aPD1_mEry"
for s in $sample
do
subpath="${s}/outs/filtered_feature_bc_matrix/"
mkdir $s
Rscript ${R_script}/1.Seurat3_QualityControl.r -i 1.rawdata/${subpath} -w 2.cellanno/$s/ -o ${s} -s Mus > 2.cellanno/${s}/1.QualityControl.log
done


### 3.combined
rawdata_path="$(pwd)/Mus/"
combined="Spleen"
mkdir $combined
cd $combined
ln -s ../../2.cellanno/IgG_mEry/*_raw.rds IgG_mEry_raw.rds
ln -s ../../2.cellanno/aPD1/*_raw.rds aPD1_raw.rds
ln -s ../../2.cellanno/aPD1_mEry/*_raw.rds aPD1_mEry_raw.rds
cd ..

sample_names="IgG_mEry,aPD1,aPD1_mEry"
rds_files="IgG_mEry_raw.rds,aPD1_raw.rds,aPD1_mEry_raw.rds"
filter_fea="6000,6000,6000"
filter_mt="8,8,8"
# Because the integrated dataset contained clusters annotated as “unknown,” we adopted a more relaxed strategy by removing doublets based on the theoretical expectation (nExp_poi) without applying homotypic adjustment. Validation analyses indicated that this approach did not affect the overall proportion of major cell populations. 
Rscript ${R_script}/2.IntegratedSamples.Seurat3_PCAselection.dedoublet.r -w ${rawdata_path}/3.combination/$combined -l $sample_names -f $rds_files -u $filter_fea -m $filter_mt -o ${combined} > ${rawdata_path}/3.combination/$combined/2.PCAselection.log
Rscript ${R_script}/CCA_3.Seurat3_tSNEorUMAP.r -w ${rawdata_path}/3.combination/$combined -p 20 -f ${combined}_Origin_Integrated.rds -o ${combined} > ${rawdata_path}/3.combination/$combined/3.CCA_p20.log

# annotation
cd $rawdata_path
Rscript 2.AnnotationCluster.r
Rscript 3.Myeloid.r
Rscript 4.Tcell.r  