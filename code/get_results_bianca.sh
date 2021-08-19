#get_results_bianca.sh

#########################################
###---- on bianca
#########################################

#empty wharf
rm -r /proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/*

#cp results
cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/* \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/



#########################################
####---- on mac
#########################################

cd /Users/gusarv/Documents/projekt/SjS/data/pss_bcells_scRNAseq/results
sftp -r gustava-sens2020501@bianca-sftp.uppmax.uu.se:gustava-sens2020501
get *
