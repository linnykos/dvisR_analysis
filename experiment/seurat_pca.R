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

scale_mat <- bm[["RNA"]]@scale.data[bm[["RNA"]]@var.features,]
idx1 <- which(rownames(scale_mat) == "IGHA1")
idx2 <- which(rownames(scale_mat) == "GYPE")
plot(scale_mat[idx1,], scale_mat[idx2,], asp = T, col = celltype_factor,
     xlab = rownames(scale_mat)[idx1],
     ylab = rownames(scale_mat)[idx2])
# This is what plotting the scale.data matrix looks like
# It still looks like it's very aligned with the axes 

pca_obj <- bm[["pca"]]
pca_mat <- pca_obj@cell.embeddings %*% t(pca_obj@feature.loadings)
idx1 <- which(colnames(pca_mat) == "IGHA1")
idx2 <- which(colnames(pca_mat) == "GYPE")
plot(pca_mat[,idx1], pca_mat[,idx2], asp = T, col = celltype_factor,
     xlab = colnames(pca_mat)[idx1],
     ylab = colnames(pca_mat)[idx2])
