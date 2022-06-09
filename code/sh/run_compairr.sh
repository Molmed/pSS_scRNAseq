#!/bin/bash -l

# 211110
# Gustav Arvidsson
# run compairr to cluster VDJ clonotypes
# https://github.com/uio-bmi/compairr

wd=/Users/gusarv/Documents/projekt/SjS/data/pss_bcells_scRNAseq

cd $wd

compairr --version > "${wd}/results/VDJ_compairr_version.txt"

compairr \
--cluster "${wd}/results/VDJ_immcantation_out.tsv" \
--indels \
--differences 1 \
--nucleotides \
--threads 4 \
--log "${wd}/results/VDJ_compairr_out_nt.log.txt" \
--output "${wd}/results/VDJ_compairr_out_nt.tsv"

compairr \
--cluster "${wd}/results/VDJ_immcantation_out.tsv" \
--differences 2 \
--threads 4 \
--log "${wd}/results/VDJ_compairr_out_aa.log.txt" \
--output "${wd}/results/VDJ_compairr_out_aa.tsv"
