#00_dependencies.R

#install.packages('ggthemes', dependencies = TRUE)
#reticulate::conda_export("../results/environment_R_pSSscRNAseq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

####------------------------------------------------------------------------------
####----- system
####------------------------------------------------------------------------------
#install.packages('future', dependencies = TRUE)
suppressMessages(library(future)) #for planning multiprocessess

#BiocManager::install('BiocParallel', dependencies = TRUE)
suppressMessages(library(BiocParallel))

#install.packages('pryr', dependencies = TRUE)
suppressMessages(library(pryr)) #check memory usage


####------------------------------------------------------------------------------
####----- general
####------------------------------------------------------------------------------
#install.packages('stringr', dependencies = TRUE)
suppressMessages(library(stringr)) #for string modulations, find and replace, sort etc.

#BiocManager::install('biomaRt', dependencies = TRUE)
suppressMessages(library(biomaRt))

#BiocManager::install('rtracklayer', dependencies = TRUE)
suppressMessages(library(rtracklayer))

#install.packages("dplyr", dependencies = TRUE)
suppressMessages(library(dplyr))

#install.packages("tidyr", dependencies = TRUE)
suppressMessages(library(tidyr))

#install.packages("future.apply", dependencies = TRUE)
suppressMessages(library(future.apply))

#BiocManager::install('msa', dependencies = TRUE)
suppressMessages(library(msa))

#install.packages("seqinr", dependencies = TRUE)
suppressMessages(library(seqinr))

#install.packages("Matrix", dependencies = TRUE)
suppressMessages(library(Matrix))

#install.packages("readxl", dependencies = TRUE)
suppressMessages(library(readxl))

#install.packages("nipals", dependencies = TRUE)
suppressMessages(library(nipals))

#install.packages("devtools", dependencies = TRUE)
#install.packages('remotes', dependencies = TRUE)


####------------------------------------------------------------------------------
####----- plot
####------------------------------------------------------------------------------
#install.packages("plotly", dependencies = TRUE)
suppressMessages(library(plotly))

#install.packages("ggplot2", dependencies = TRUE)
suppressMessages(library(ggplot2))

#install.packages("cowplot", dependencies = TRUE)
suppressMessages(library(cowplot))

#install.packages("rafalib", dependencies = TRUE)
suppressMessages(library(rafalib))

#install.packages("RColorBrewer", dependencies = TRUE)
suppressMessages(library(RColorBrewer))

#install.packages("ggthemes", dependencies = TRUE)
suppressMessages(library(ggthemes))

#install.packages("scales", dependencies = TRUE)
suppressMessages(library(scales))

#install.packages("VennDiagram", dependencies = TRUE)
suppressMessages(library(VennDiagram))

#install.packages("ggrepel", dependencies = TRUE)
suppressMessages(library(ggrepel))

#install.packages("UpSetR", dependencies = TRUE)
suppressMessages(library(UpSetR))

#remotes::install_github("czarnewski/niceRplots")
library(niceRplots)
#source("../../niceRplots/R/add_fig_label.R")
#source("../../niceRplots/R/helper_functions.R")
#source("../../niceRplots/R/plotting_functions.R")
#source("../../niceRplots/R/trajectory_plots.R")

#install.packages('pheatmap', dependencies = TRUE)
suppressMessages(library(pheatmap))

#install.packages('ggpubr', dependencies = TRUE)
suppressMessages(library(ggpubr))

#install.packages('gridExtra', dependencies = TRUE)
suppressMessages(library(gridExtra))

#install.packages('scales', dependencies = TRUE)
suppressMessages(library(scales)) #scales for plotting

#install.packages('lessR', dependencies = TRUE)
#suppressMessages(library(lessR))

#install.packages('gridtext', dependencies = TRUE)
suppressMessages(library(gridtext))


####------------------------------------------------------------------------------
####----- scRNAseq
####------------------------------------------------------------------------------
#install.packages('Seurat', dependencies = TRUE)
suppressMessages(library(Seurat)) #read and process scRNAseq data

#BiocManager::install("scater", dependencies = TRUE)
suppressMessages(library(scater))

#BiocManager::install("edgeR", dependencies = TRUE)
suppressMessages(library(edgeR))

#BiocManager::install("scran", dependencies = TRUE)
suppressMessages(library(scran))

#BiocManager::install("scDblFinder", dependencies = TRUE)
suppressMessages(library(scDblFinder))

#install.packages("harmony", dependencies = TRUE)
suppressMessages(library(harmony))

#devtools::install_github("powellgenomicslab/scPred")
suppressMessages(library(scPred))

#install.packages("tidyseurat")
suppressMessages(library(tidyseurat))

#BiocManager::install("glmGamPoi", dependencies = TRUE)
suppressMessages(library(glmGamPoi))


####------------------------------------------------------------------------------
####----- VDJ
####------------------------------------------------------------------------------
#BiocManager::install("ggtree", dependencies = TRUE)
#suppressMessages(library(ggtree))

#BiocManager::install("scRepertoire", dependencies = TRUE)
##devtools::install_github("ncborcherding/scRepertoire@dev")
suppressMessages(library(scRepertoire))#; packageVersion("scRepertoire") #local :[1] ???1.3.4???; bianca:[1] ???1.0.0???

#install.packages("dowser")
#suppressMessages(library(dowser))

#install.packages("phangorn", dependencies = TRUE)
#suppressMessages(library(phangorn))

#install.packages('stringdist', dependencies = TRUE)
suppressMessages(library(stringdist))


####------------------------------------------------------------------------------
#### Immcantation
####------------------------------------------------------------------------------
#install.packages('alakazam', dependencies = TRUE)
suppressMessages(library(alakazam))

#install.packages('shazam', dependencies = TRUE)
suppressMessages(library(shazam))

#install.packages('airr', dependencies = TRUE)
suppressMessages(library(airr))

#install.packages('scoper', dependencies = TRUE)
suppressMessages(library(scoper))


####------------------------------------------------------------------------------
####GSEA
####------------------------------------------------------------------------------
#install.packages('enrichR', dependencies = TRUE)
suppressMessages(library(enrichR))

#BiocManager::install("fgsea", dependencies = TRUE)
suppressMessages(library(fgsea))

#BiocManager::install("EnrichmentBrowser")
suppressMessages(library(EnrichmentBrowser))


####------------------------------------------------------------------------------
#### save sessionInfo()
####------------------------------------------------------------------------------
getwd()
setwd("/Users/gusarv/Documents/projekt/SjS/data/pss_bcells_scRNAseq/code")
thesession <- sessionInfo()
saveRDS(thesession, paste0("../results/sessionInfo", Sys.Date(), ".rds"))
