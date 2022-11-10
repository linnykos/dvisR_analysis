rm(list=ls())
library(Seurat)
library(SAVER)
library(SeuratData)
library(cowplot)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# SeuratData::InstallData("bmcite")
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- 'RNA'
bm <- Seurat::NormalizeData(bm) %>% Seurat::FindVariableFeatures()

mat <- bm[["RNA"]]@counts[Seurat::VariableFeatures(bm, assay = "RNA"),]
mat <- scale(mat)

print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 4)
save(saver_res, bm,
     date_of_run, session_info,
     file = "../../../../../out/dvisR_analysis/kevin_tmp/Writeup1/Writeup1_bm_saver.RData")
