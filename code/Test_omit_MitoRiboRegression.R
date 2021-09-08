paste0("R script started at: ", Sys.time())




## ----------------------------------------------------------------------------------------------------------------

#source("./00_dependencies.R")

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
#library(niceRplots)

source("../../niceRplots/R/add_fig_label.R")
source("../../niceRplots/R/helper_functions.R")
source("../../niceRplots/R/plotting_functions.R")
source("../../niceRplots/R/trajectory_plots.R")

####----- scRNAseq
suppressMessages(library(Seurat)) #read and process scRNAseq data
suppressMessages(library(scDblFinder))
suppressMessages(library(harmony))
suppressMessages(library(scater))
#suppressMessages(library(Nebulosa))


plan("multiprocess", workers = 16) #parallellization for seurat using future package:
options(future.globals.maxSize = 100000 * 1024 ^ 2) #resolve memory related issues related to parallellization 100000 here refers to approx 100GB memory

setwd("/castor/project/proj_nobackup/sjs1/R/pss_bcells_scRNAseq/code"); getwd()

sessionInfo()

paste0("Libraries loaded at: ", Sys.time())



gex <- readRDS("../results/GEX3_210826.rds")
Idents(object = gex) <- "orig.ident"

print("GEX3 size:")
print(object.size(gex), units = "GB")

annot <- readRDS("../suppl/annot/refdata-gex-GRCh38-2020-A.genes.gtf.annot.rds")

head(annot)

gex <- NormalizeData(gex, scale.factor = 10000)
gex <- FindVariableFeatures(gex, selection.method = "vst", nfeatures = 4000)

VariableFeatures(gex) <-
  VariableFeatures(gex)[VariableFeatures(gex) %in%
                          annot$external_gene_name[annot$gene_biotype == "protein_coding"]]

gex <- ScaleData(gex, 
                 features = VariableFeatures(gex),
                 vars.to.regress = c("nFeature_RNA"))

paste0("NormalizeData, FindVariableFeatures and ScaleData done at: ", Sys.time())


saveRDS(gex, file = paste0("../results/GEX4_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".rds"))


print("GEX4 size:")
print(object.size(gex), units = "GB")

gex <- RunPCA(
  gex,
  assay = "RNA",
  features = VariableFeatures(gex),
  npcs = 100,
  reduction.name = "pca_1",
  verbose = TRUE
)


gex@assays$RNA@scale.data <- matrix(0)
gc()

gex

paste0("PCA1 done at:", Sys.time())


saveRDS(gex, file = paste0("../results/GEX5_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".rds"))


print("GEX5 size:")
print(object.size(gex), units = "GB")

gex <- RunHarmony(
  gex,
  group.by.vars = "orig.ident",
  reduction = "pca_1",
  reduction.save = "harmony_1",
  assay = "RNA",
  project.dim = FALSE, #project.dim = FALSE needed for seurat object v4.0.0??
  verbose = TRUE
)



saveRDS(gex, file = paste0("../results/GEX6_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".rds"))

print("GEX6 size:")
print(object.size(gex), units = "GB")


gex <- RunUMAP(
  gex,
  dims = 1:50, #or 100?
  reduction = "harmony_1",
  metric = "correlation",
  reduction.name = "umap_1",
  min.dist = 0.4, #local cell separation
  spread = .5, #global cell separation
  n.neighbors = 30, #30 on bianca!
  repulsion.strength = 0.4, #~2x min.dist, repulsion to cells fr annat cluster, global separation
  negative.sample.rate = 50, #ant ggr n.neighbors celler, global distance
  n.epochs = 100,
  n.components = 2
)

paste0("UMAP1 done at: ", Sys.time())


saveRDS(gex, file = paste0("../results/GEX7_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".rds"))

print("GEX7 size:")
print(object.size(gex), units = "GB")


gex <- FindNeighbors(gex,
                     reduction = "harmony_1",
                     dims = 1:50,
                     k.param = 30,
                     #~30 on bianca!
                     verbose = TRUE
)

res <- 0.3 #Bianca RNA_nn resolution

gex <- FindClusters(gex, 
                    resolution = res, 
                    verbose = TRUE, 
                    graph.name = "RNA_nn")

table(gex@meta.data$seurat_clusters)

paste0("Clustering1 done at: ", Sys.time())

Idents(object = gex) <- "orig.ident"
saveRDS(gex, file = paste0("../results/GEX8_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".rds"))

print("GEX8 size:")
print(object.size(gex), units = "GB")

p1.u.s <- DimPlot(gex, 
                  reduction = "umap_1", 
                  pt.size = .1, 
                  group.by = "orig.ident", 
                  raster=FALSE) + 
  NoLegend() + 
  ggtitle("UMAP_1, by Sample")

p1.u.c <- DimPlot(gex, 
                  reduction = "umap_1", 
                  pt.size = .1, 
                  group.by = "seurat_clusters", 
                  label = TRUE, 
                  repel = TRUE, 
                  raster=FALSE) + 
  NoLegend() + 
  ggtitle("UMAP_1, by Cluster")

ggsave2(paste0("../results/UMAP_allCells_clustRes", res, "_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".png"), 
        plot_grid(p1.u.s, p1.u.c, ncol = 2, labels = "AUTO", align = "h"), 
        width = 30, height = 15, unit = "cm")

ig.features <- c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", 
                 "IGHG4",  "IGHD", "IGHE", "IGHM" )
b.markers <- c("CD79A", "CD79B", "MS4A1", "CD19", "CD27", 
               "IGHA1", "IGHD", "IGHM", "JCHAIN", "MME")
pbmc.markers <- c("MS4A1", "CD19", "CD27", "CD79A", "GNLY", "CD3E", "CD14", "LYZ", "CD8A")
my_pars <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb")

####-----  PBMC markers
p.pbmc.f <- FeaturePlot(gex, 
                        features = pbmc.markers, 
                        reduction = "umap_1", 
                        raster=FALSE, 
                        order = TRUE)

ggsave2(paste0("../results/FeaturePlot_allCells_pbmcMarkers_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".png"), 
        p.pbmc.f, 
        width=25, height=25, unit="cm")

####-----  B markers
p.b.markers.f <- FeaturePlot(gex, 
                             features = b.markers, 
                             reduction = "umap_1", 
                             raster=FALSE, 
                             order = TRUE)
ggsave2(paste0("../results/FeaturePlot_allCells_BMarkers_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".png"), 
        p.b.markers.f, 
        width=30, height=20, unit="cm")

####-----  IG markers
p.ig.markers.f <- FeaturePlot(gex, 
                              features = ig.features, 
                              reduction = "umap_1", 
                              raster=FALSE, 
                              order = TRUE)

ggsave2(paste0("../results/FeaturePlot_allCells_IGMarkers_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".png"), 
        p.ig.markers.f, 
        width=25, height=25, unit="cm")

####-----  my_pars
p.myPars.f <- FeaturePlot(gex, 
                          features = my_pars , 
                          reduction = "umap_1", 
                          raster=FALSE, 
                          order = TRUE)
ggsave2(paste0("../results/FeaturePlot_allCells_myPars_nFeatureRegression_", format(Sys.time(), "%y%m%d"), ".png"), 
        p.myPars.f, 
        width=25, height=25, unit="cm")


