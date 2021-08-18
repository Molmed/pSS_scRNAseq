#00_dependencies.R

####----- system
suppressMessages(library(future)) #for planning multiprocessess
suppressMessages(library(BiocParallel))
suppressMessages(library(pryr)) #check memory usage

####----- general
suppressMessages(library(stringr)) #for string modulations, find and replace etc.
suppressMessages(library(biomaRt))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(future.apply))

####----- plot
suppressMessages(library(plotly))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(rafalib))
# remotes::install_github("czarnewski/niceRplots")
library(niceRplots)
#source("../../niceRplots/R/add_fig_label.R")
#source("../../niceRplots/R/helper_functions.R")
#source("../../niceRplots/R/plotting_functions.R")
#source("../../niceRplots/R/trajectory_plots.R")

####----- scRNAseq
suppressMessages(library(Seurat)) #read and process scRNAseq data
suppressMessages(library(scater))
suppressMessages(library(scDblFinder))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(scRepertoire))

