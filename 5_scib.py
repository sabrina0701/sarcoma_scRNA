########### Python script for scib. ###########
########### Created by Danwen Qian on 4-09-2023. ###########
########### Last modified by Danwen Qian on 4-09-2023. ###########

import numpy as np
import pandas as pd
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import matplotlib.pyplot as plt
import scib



adata = sc.read_h5ad("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/immune_integrated_log.h5ad")

adata = sc.read_h5ad("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/Myeloids_SCT.h5ad")

adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
adata.obsm["CCA"] = adata.obsm["X_integrated.cca"]
adata.obsm["rPCA"] = adata.obsm["X_integrated.rpca"]
adata.obsm["Harmony"] = adata.obsm["X_harmony"]
adata.obsm["FastMNN"] = adata.obsm["X_integrated.mnn"]
adata.obsm["scVI"] = adata.obsm["X_integrated.scvi"]
 
bm = Benchmarker(
    adata,
    batch_key="study",
    label_key="scATOMIC_pred",
    embedding_obsm_keys=["Unintegrated", "CCA","rPCA","Harmony","FastMNN", "scVI"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    n_jobs=1,
)

bm.benchmark()

bm.plot_results_table(min_max_scale=False,save_dir="/mnt/beegfs01/scratch/d_qian/SingleCell/plots")


df = bm.get_results(min_max_scale=False)
df.to_csv("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/metrics_immune.csv")




adata = sc.read_h5ad("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/immune_integrated_SCT.h5ad")

adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
adata.obsm["CCA"] = adata.obsm["X_integrated.cca"]
adata.obsm["rPCA"] = adata.obsm["X_integrated.rpca"]
adata.obsm["Harmony"] = adata.obsm["X_harmony"]

 
bm = Benchmarker(
    adata,
    batch_key="study",
    label_key="scATOMIC_pred",
    embedding_obsm_keys=["Unintegrated", "CCA","rPCA","Harmony"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    n_jobs=1,
)
bm.benchmark()

bm.plot_results_table(min_max_scale=False,save_dir="/mnt/beegfs01/scratch/d_qian/SingleCell/plots")


df = bm.get_results(min_max_scale=False)
df.to_csv("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/metrics_Myeloid_SCT.csv")
