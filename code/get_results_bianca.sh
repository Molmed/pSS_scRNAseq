#get_results_bianca.sh

#########################################
###---- on bianca
#########################################

#empty wharf
ll /proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/
rm -r /proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/*
ll /proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

#cp result plots
cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*27.png \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*.pdf \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*.csv \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

#cp gex
cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/GEX8* \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

ll /proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/



#########################################
####---- on mac
#########################################

cd /Users/gusarv/Documents/projekt/SjS/data/pss_bcells_scRNAseq/results
sftp -r gustava-sens2020501@bianca-sftp.uppmax.uu.se:gustava-sens2020501
get *
