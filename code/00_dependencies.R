#00_dependencies.R

#install.packages('ggthemes', dependencies = TRUE)

#reticulate::conda_export("../results/environment_R_pSSscRNAseq")

####----- system
suppressMessages(library(future)) #for planning multiprocessess
suppressMessages(library(BiocParallel))
suppressMessages(library(pryr)) #check memory usage


####----- general
suppressMessages(library(stringr)) #for string modulations, find and replace, sort etc.
suppressMessages(library(biomaRt))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(future.apply))
suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(Matrix))

####----- plot
suppressMessages(library(plotly))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(rafalib))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggthemes))
suppressMessages(library(scales))

# remotes::install_github("czarnewski/niceRplots")
library(niceRplots)
#source("../../niceRplots/R/add_fig_label.R")
#source("../../niceRplots/R/helper_functions.R")
#source("../../niceRplots/R/plotting_functions.R")
#source("../../niceRplots/R/trajectory_plots.R")
suppressMessages(library(pheatmap))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(scales)) #scales for plotting
#suppressMessages(library(lessR))


####----- scRNAseq
suppressMessages(library(Seurat)) #read and process scRNAseq data
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(scDblFinder))
suppressMessages(library(harmony))
suppressMessages(library(scPred))
#suppressMessages(library(Nebulosa))
#devtools::install_github("ncborcherding/scRepertoire@dev")
suppressMessages(library(tidyseurat))
suppressMessages(library(edgeR))
suppressMessages(library(glmGamPoi))

####----- VDJ
suppressMessages(library(scRepertoire))#; packageVersion("scRepertoire") #local :[1] ‘1.3.4’; bianca:[1] ‘1.0.0’
#suppressMessages(library(dowser))
#install.packages("dowser")
#install.packages("ggtree") #package ‘ggtree’ is not available for this version of R
#install.packages("phangorn") #package ‘phangorn’ is not available for this version of R

# Immcantation
suppressMessages(library(alakazam))
suppressMessages(library(shazam))
suppressMessages(library(airr))
suppressMessages(library(scoper))

#GSEA
suppressMessages(library(enrichR))
suppressMessages(library(fgsea))



#save sessionInfo()
thesession <- sessionInfo()
saveRDS(thesession, paste0("../results/sessionInfo", Sys.Date(), ".rds"))



