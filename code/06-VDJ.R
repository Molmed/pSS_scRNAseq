## ----include = FALSE----------------------------------------------------------------------------------------------------

#source("./00_dependencies.R")
#devtools::install_github("ncborcherding/scRepertoire")
#devtools::install_github("ncborcherding/scRepertoire@dev")
#devtools::install_github("ncborcherding/scRepertoire@refine")
library(scRepertoire)




## ----include = FALSE----------------------------------------------------------------------------------------------------
# library(niceRplots)
# saveRDS(gex, file = paste0("../suppl/GEX14_", format(Sys.time(), "%y%m%d"), ".rds"))
# gex <- readRDS("../suppl/GEX14_210628.rds")
# sort(sapply(ls(), function(x){format(object.size(get(x)), units = "Gb")}))
# colnames(gex)
# gex@assays$RNA@scale.data <- matrix(0); gc()
# # gex@assays$RNA@scale.data <- matrix(0)
# saveRDS(gex,"../results/B-cells.rds")
#
# plot_meta(gex,"umap_2",feat = "seurat_clusters")



## -----------------------------------------------------------------------------------------------------------------------

# read 10x VDJ BCR files
fl_vdj_B <- list.files("../data/clono_VDJ", pattern = "_B_")

cl_vdj_B <- lapply(fl_vdj_B, function(x) {
  return(read.csv(paste0("../data/clono_VDJ/", x), stringsAsFactors = FALSE))
})

names(cl_vdj_B) <- sub("_.*", "", fl_vdj_B)

cl_vdj_B <- lapply(seq_along(cl_vdj_B), function(i) {
  cl_vdj_B[[i]]$barcode <-
    paste0(names(cl_vdj_B)[i], "_", gsub("-1", "", cl_vdj_B[[i]]$barcode))
  return(cl_vdj_B[[i]])
})

names(cl_vdj_B) <- sub("_.*", "", fl_vdj_B)

#handle samples where libraries were made twice
cl_vdj_B[["P007"]] <- rbind(cl_vdj_B[["P007a"]], cl_vdj_B[["P007b"]])
cl_vdj_B[["P007a"]] <- NULL
cl_vdj_B[["P007b"]] <- NULL

cl_vdj_B[["P008"]] <- rbind(cl_vdj_B[["P008a"]], cl_vdj_B[["P008b"]])
cl_vdj_B[["P008a"]] <- NULL
cl_vdj_B[["P008b"]] <- NULL

names(cl_vdj_B)
cl_vdj_B <- cl_vdj_B[order(names(cl_vdj_B))]
names(cl_vdj_B)

saveRDS(cl_vdj_B, paste0("../results/cl_vdj_B_test_",format(Sys.time(), "%y%m%d"),".rds"))

# names(cl_vdj_B[1])
# head(cl_vdj_B[[1]])
# length(cl_vdj_B)

lapply(cl_vdj_B, function(x){
  print(head(x, n = 2))
  print(dim(x))
})

#using scRepertoire to generate a combined data frame
combined <- combineBCR(cl_vdj_B,
                       samples = names(cl_vdj_B),
                       #ID = sub("_.*", "", fl_vdj_B),
                       ID = rep("", times = length(cl_vdj_B)),
                       removeNA = TRUE,
                       removeMulti = TRUE) #Error: vector memory exhausted (limit reached?)

str(combined)
head(combined)


saveRDS(combined, paste0("../results/combined_BCR_scRepertoire_test_",format(Sys.time(), "%y%m%d"),".rds"))
write.csv(combined, paste0("../results/combined_BCR_scRepertoire_test_",format(Sys.time(), "%y%m%d"),".csv"))


# # read 10x VDJ TCR files
# fl_vdj_T <- list.files("../data/clono_VDJ",pattern = "_T_")
# cl_vdj_T <- lapply(fl_vdj_T, function(x){
#   return(read.csv(paste0("../data/clono_VDJ/", x)))
# })
# names(cl_vdj_T) <- sub("_.*", "", fl_vdj_T)
# head(cl_vdj_T[[1]])
#
# combined <- combineTCR(cl_vdj_T,
#                        samples = sub("_.*", "", fl_vdj_T),
#                        #ID = sub("_.*", "", fl_vdj_T),
#                        ID = rep("", times = length(cl_vdj_T)),
#                        removeNA = FALSE,
#                        removeMulti = TRUE,
#                        cells = c("T-AB"))
#
# # read TRUST4 files
# fl_TRUST4 <- list.files("../data/clono_TRUST4")
# cl_TRUST4 <- lapply(fl_TRUST4, function(x){
#   return(read.table(paste0("../data/clono_TRUST4/", x)))
# })
# names(cl_TRUST4) <- sub("_.*", "", fl_TRUST4)
# head(cl_TRUST4[[1]])
#
# combineTRUST4(cl_TRUST4)







# library(immunarch)
#
# immdata_10x <- repLoad("../data/clono_VDJ")
# immdata_trust4 <- repLoad("../data/clono_TRUST4")
# str(immdata_10x)


#plot_meta(gex,"umap_2",feat = "seurat_clusters")

# dir.create(paste0("../data/B_clono_VDJ"))
# for(i in sub("_.*","",file_list)){
#   dir.create(paste0("../data/B_clono_VDJ/",i))
#   system(paste0("cp ../data/clono_VDJ/",grep("_B_",grep(paste0(i,"_"),file_list,value = T),value = T)," ",paste0("../data/B_clono_VDJ/",i,"/") ))
# }



## -----------------------------------------------------------------------------------------------------------------------



