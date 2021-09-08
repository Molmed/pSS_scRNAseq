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
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*06.png \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*03.pdf \
/proj/sens2020501/nobackup/wharf/gustava/gustava-sens2020501/

cp -r \
/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/results/*03.csv \
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


#########################################
####---- rename files
#########################################

#mac
rename 's/noRiboRegression/nFeatureMitoRegression/' *
rename 's/Markers_21/Markers_nFeatureMitoRiboRegression_21/' *
rename 's/myPars210826/myPars_nFeatureMitoRiboRegression_210826/' *
rename 's/clustRes0.3_2/clustRes0.3_nFeatureMitoRiboRegression_2/' *
rename 's/__noRegression/_noRegression/' *

#bianca
rename noRiboRegression nFeatureMitoRegression *
rename Markers_21 Markers_nFeatureMitoRiboRegression_21 *
rename myPars210826 myPars_nFeatureMitoRiboRegression_210826 *
rename clustRes0.3_2 clustRes0.3_nFeatureMitoRiboRegression_2 *
rename __noRegression _noRegression *

mmv -n GEX\?_210826.rds GEX\#1_nFeatureMitoRiboRegression_210826.rds
mmv GEX\?_210826.rds GEX\#1_nFeatureMitoRiboRegression_210826.rds

mmv -n GEX\?_dnsample.500_210826.rds GEX\#1_dnsample.500_nFeatureMitoRiboRegression_210826.rds
mmv GEX\?_dnsample.500_210826.rds GEX\#1_dnsample.500_nFeatureMitoRiboRegression_210826.rds

mmv -n GEX\?_dnsample.500_nFeatureMitoRiboRegression_210826.rds GEX\#1_dnsample.500_210826.rds
mmv GEX\?_dnsample.500_nFeatureMitoRiboRegression_210826.rds GEX\#1_dnsample.500_210826.rds

mmv -n GEX\?_nFeatureMitoRiboRegression_210826.rds GEX\#1_210826.rds
mmv GEX\?_nFeatureMitoRiboRegression_210826.rds GEX\#1_210826.rds

mmv -n GEX\?_nFeatureRegression_210902.rds GEX\#1_210902.rds
mmv GEX\?_nFeatureRegression_210902.rds GEX\#1_210902.rds
