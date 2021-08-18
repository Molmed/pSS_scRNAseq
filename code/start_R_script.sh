#!/bin/bash -l
#SBATCH -A sens2020501
#SBATCH -t 18:00:00
#SBATCH -p node
#SBATCH -J process_gex

ml R/4.0.4
ml R_packages/4.0.4

cd /castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/code

echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SNIC_TMP: ${SNIC_TMP}"

echo "start R script:"
date
Rscript /castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/code/Analysis_master_bianca.R
