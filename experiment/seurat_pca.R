library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
bm <- SeuratData::LoadData(ds = "bmcite")
Seurat::DefaultAssay(bm) <- 'RNA'
# Just an example way to preprocess data, to demonstrate the PCA part
# You could use SCTransform into RunPCA also
bm <- Seurat::NormalizeData(bm) %>% Seurat::FindVariableFeatures() %>% Seurat::ScaleData() %>% Seurat::RunPCA()
celltype_factor <- as.numeric(as.factor(bm@meta.data$celltype.l2))

scale_mat <- t(bm[["RNA"]]@scale.data[bm[["RNA"]]@var.features,])
dim(scale_mat); round(scale_mat[1:5,1:5],2)
idx1 <- which(colnames(scale_mat) == "IGHA1")
idx2 <- which(colnames(scale_mat) == "GYPE")
plot(scale_mat[,idx1], scale_mat[,idx2], asp = T, col = celltype_factor,
     xlab = colnames(scale_mat)[idx1],
     ylab = colnames(scale_mat)[idx2])
# This is what plotting the scale.data matrix looks like
# It still looks like it's very aligned with the axes 

pca_obj <- bm[["pca"]]
pca_mat <- pca_obj@cell.embeddings %*% t(pca_obj@feature.loadings)
dim(pca_mat); round(pca_mat[1:5,1:5],2)
idx1 <- which(colnames(pca_mat) == "IGHA1")
idx2 <- which(colnames(pca_mat) == "GYPE")
plot(pca_mat[,idx1], pca_mat[,idx2], asp = T, col = celltype_factor,
     xlab = colnames(pca_mat)[idx1],
     ylab = colnames(pca_mat)[idx2])
# By plotting the low-rank matrix, some of the structure is more visible 
# You can verify that pca_mat is indeed low-rank by running
# Matrix::rankMatrix(pca_mat) (it should return 50), but this line takes a few minutes to run