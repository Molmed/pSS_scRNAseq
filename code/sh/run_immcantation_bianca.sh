#!/bin/bash -l
#SBATCH -A sens2020501
#SBATCH -t 8:00:00
#SBATCH -p core
#SBATCH -n 16
#SBATCH -J immcantation

#run immcantation on bianca for B VDJ output from cellranger
#210923

projdir=/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/data/clono_B_VDJ_tmp

cd $projdir

for i in *; do  echo $i

  cd $i

  singularity exec -B ${projdir}/${i}:/data /castor/project/proj_nobackup/sjs1/R/immcantation/immcantation_suite-4.2.0.sif \
  changeo-10x \
  -s /data/*_B_filtered_contig.fasta \
  -a /data/*_B_filtered_contig_annotations.csv \
  -o /data/immcantation_out \
  -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IG* \
  -g human \
  -n $i \
  -t ig \
  -x 0.1 \
  -f airr \
  -p 16

  cd ..

done






# ll /usr/local/share/germlines/imgt/human/vdj
# ll /usr/local/share/igblast
