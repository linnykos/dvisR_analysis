rm(list=ls())
set.seed(10)
sessionInfo()

## uncomment the following line to install package
# devtools::install_github("linnykos/dvisR", subdir = "dvisR")

## run the following lines for an interactive example
dat <- read.csv("../../data/Zeisel_preprocessed.csv", row.names = 1)
dat <- dat[,1:30]
cell_type <- read.table("../../data/Zeisel_cell_info.txt", sep = "\t", header = 1)
cluster_labels <- as.numeric(as.factor(cell_type$level1class))

set.seed(10)
res <- dvisR::dvisR_system(dat, cluster_labels = cluster_labels,
                    system_options = system_options_default(minimum_instances_first_phase = 4,
                                                            ntrials = 10, classifier = classifier_xgboost_closure(verbose = 0)), 
                    pch = 16, asp = T)

names(res)
res
head(res$df)
dvisR_embedding(res, embedding_method = embedding_umap)

par(mar = c(0.5, 0.5, 0.5, 0.5))
dvisR_heatmap(res)
dvisR_heatmap(res, unlabeled = T)

## you could also run it in "debugging mode" if you want to skip the interactive parts 
##   (this is purely for debugging purposes, not for practical use)
set.seed(10)
res2 <- dvisR::dvisR_system(dat, 
                           system_options = system_options_default(minimum_instances_first_phase = 4,
                                                                   ntrials = 10, classifier = classifier_xgboost_closure(verbose = 0)),
                           debugging_inputs = list(round_inputs = list("5,6,9", "", "3,5,9,6",
                                                                       "y", "y", "n", "y", "y", "y", "n", "y", "n", "n")))
head(res2$df)

dvisR_embedding(res2)
dvisR_embedding(res2, embedding_method = embedding_umap)

