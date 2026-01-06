########### R script for cell type prediction. ###########
########### Created by Danwen Qian on 04-08-2025. ###########
########### Last modified by Danwen Qian on 04-08-2025. ###########

library(reticulate)

use_python("~/miniconda3/envs/scRNA/bin/python", required = TRUE)

library(Seurat)
#library(SeuratDisk)
library(ggplot2)
library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)


setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets")


combined<-readRDS("combined_filter.rds")
DSRCT<-subset(combined,study=="DSRCT_in_house")

DSRCT<-SplitObject(DSRCT,split.by="sample_ID")

for (i in 1:length(DSRCT)){
    sparse_matrix <- GetAssayData(DSRCT[[i]], assay = "RNA", slot = "counts")
    cell_predictions <- run_scATOMIC(sparse_matrix,mc.cores = 10)
    results_DSRCT <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = T, modify_results = T, mc.cores = 10, raw_counts = sparse_matrix, min_prop = 0.5, known_cancer_type = unique(DSRCT[[i]]$cancer_type))
    write.csv(results_DSRCT, paste0("annotation/DSRCT/",names(DSRCT)[i],"_CNV.csv"))
}

