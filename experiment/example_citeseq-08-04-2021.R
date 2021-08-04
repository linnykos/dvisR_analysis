library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm)
bm <- FindVariableFeatures(bm)
bm <- ScaleData(bm)

#############

# how do I now extract the 2000-gene preprocessed dataset?
bm
class(bm)
names(bm)
bm[["RNA"]]
zz <- bm[["RNA"]]
class(zz)

zz@counts[81:90,20:30] # raw data
zz@data[81:90,20:30] # log normalized
zz@scale.data[81:90,20:30] # the data but now scaled (by the variables, in this case, the rows)

dim(zz@counts)
dim(zz@data)
dim(zz@scale.data) # we can either (option 1) take this as our dataset to use for our method
dat1 <- zz@scale.data
quantile(apply(zz@scale.data, 1, mean))
?Seurat::ScaleData # (the scale.max argument makes the means slightly not 0)

# (option 2) is to take the ONLY log-normalized data (i.e., not scaled) and extract the
# appropriate 2000 genes
dat2 <- zz@data[Seurat::VariableFeatures(bm),]
dim(dat2)
# (typically though, we want to stick with option 1)

all(rownames(dat1) %in% rownames(dat2)) # this shows us the two datasets have the same genes


#####################

# where are my celltypes?
head(bm@meta.data)
dim(bm@meta.data)
celltype <- bm@meta.data[,"celltype.l1"]
table(celltype)

####################3

tmp <- apply(dat2, 1, mean)
idx <- order(tmp, decreasing = T)[1:2]

plot(dat2[idx[1],], dat2[idx[2],],
     pch = 16, col = as.numeric(as.factor(celltype)), asp = T)

plot(dat2[100,], dat2[200,],
     pch = 16, col = as.numeric(as.factor(celltype)), asp = T)
# to see why there's so many zero's, recall how log-normalization works
# in ?Seurat::NormalizeData

plot(dat1[100,], dat1[200,],
     pch = 16, col = as.numeric(as.factor(celltype)), asp = T)

