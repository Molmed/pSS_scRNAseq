#!/bin/bash -l
#SBATCH -A sens2020501
#SBATCH -t 18:00:00
#SBATCH -p core
#SBATCH -n 16
#SBATCH -J process_gex

#210818 G.Arvidsson

#ml R/4.0.0
#ml R_packages/4.0.0 #Seurat_3.2.0 harmony_1.0

ml R/4.0.4
ml R_packages/4.0.4 #SeuratObject_4.0.0 Seurat_4.0.1 harmony_1.0

cd /castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/code

echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SNIC_TMP: ${SNIC_TMP}"

echo "start R script:"
date
Rscript /castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/code/Analysis_master_bianca.R
