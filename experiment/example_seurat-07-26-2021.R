library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "../filtered_gene_bc_matrices/hg19/")
class(pbmc.data)
pbmc.data[1:5,1:5]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
class(pbmc)
names(pbmc)
pbmc[["RNA"]]
pbmc[["RNA"]]@counts[1:5,1:5]
names(pbmc[["RNA"]])

###############

library(SeuratData)
SeuratData::AvailableData()
SeuratData::InstallData("pbmc3k")
data("pbmc3k")
class(pbmc3k)
pbmc3k[["RNA"]]@counts[90:100,1:10]
as.matrix(pbmc3k[["RNA"]]@counts[90:100,1:10])
zz <- as.matrix(pbmc3k[["RNA"]]@counts)
dim(zz)
length(which(zz == 0))
length(which(zz == 0))/prod(dim(zz))

tmp <- apply(zz, 1, function(x){length(which(x != 0))})
idx <- order(tmp, decreasing = T)[c(2000,2001)]

plot(jitter(zz[idx[1],]), jitter(zz[idx[2],]), asp = T, pch = 16, col = rgb(0,0,0,0.1))

###########

i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- Matrix::sparseMatrix(i, j, x = x))  

A
class(A)

B <- as.matrix(A)
