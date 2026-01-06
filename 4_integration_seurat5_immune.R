########### R script for integrate seurat object. ###########
########### Created by Danwen Qian on 10-10-2025. ###########
########### Last modified by Danwen Qian on 10-10-2025. ###########

setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
options(future.globals.maxSize = 5*1024^3)
options(Seurat.object.assay.version = "v5")

#Load the seurat object
blood<-readRDS("immune_cells.rds")

blood[["RNA"]] <- split(blood[["RNA"]], f = blood$study)
blood <- NormalizeData(blood)
blood <- FindVariableFeatures(blood)
blood <- ScaleData(blood)
blood <- RunPCA(blood)

blood <- FindNeighbors(blood, dims = 1:30, reduction = "pca")
blood <- FindClusters(blood, resolution = 1, cluster.name = "unintegrated_clusters")
blood <- RunUMAP(blood, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
head(blood@reductions[["umap.unintegrated"]]@cell.embeddings)
plot1<-DimPlot(blood, reduction = "umap.unintegrated", group.by = c("study", "scATOMIC_pred","unintegrated_clusters"))
ggsave("../plots/blood.umap.unintegrated.pdf", plot=plot1, height=10, width=40, dpi=600)

##five method
blood <- IntegrateLayers(
  object = blood, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

blood <- IntegrateLayers(
  object = blood, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

blood <- IntegrateLayers(
  object = blood, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

blood <- IntegrateLayers(
  object = blood, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

blood <- IntegrateLayers(
  object = blood, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/d_qian/miniconda3/envs/scRNA", verbose = FALSE
)

blood <- FindNeighbors(blood, reduction = "integrated.cca", dims = 1:30)
blood <- FindClusters(blood, resolution = 1, cluster.name = "cca_clusters")
blood <- RunUMAP(blood, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  blood,
  reduction = "umap.cca",
  group.by = c("study", "scATOMIC_pred", "cca_clusters"),
  combine = T
)
ggsave("../plots/immune.umap.cca.pdf", plot=p1, height=10, width=40, dpi=600)

blood <- FindNeighbors(blood, reduction = "integrated.rpca", dims = 1:30)
blood <- FindClusters(blood, resolution = 1, cluster.name = "rpca_clusters")
blood <- RunUMAP(blood, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p1 <- DimPlot(
  blood,
  reduction = "umap.rpca",
  group.by = c("study", "scATOMIC_pred", "rpca_clusters"),
  combine = T
)
ggsave("../plots/immune.umap.rpca.pdf", plot=p1, height=10, width=40, dpi=600)

blood <- FindNeighbors(blood, reduction = "harmony", dims = 1:30)
blood <- FindClusters(blood, resolution = 1, cluster.name = "harmony_clusters")
blood <- RunUMAP(blood, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p1 <- DimPlot(
  blood,
  reduction = "umap.harmony",
  group.by = c("study", "scATOMIC_pred", "harmony_clusters"),
  combine = T
)
ggsave("../plots/immune.umap.harmony.pdf", plot=p1, height=10, width=40, dpi=600)

blood <- FindNeighbors(blood, reduction = "integrated.mnn", dims = 1:30)
blood <- FindClusters(blood, resolution = 1, cluster.name = "mnn_clusters")
blood <- RunUMAP(blood, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
p1 <- DimPlot(
  blood,
  reduction = "umap.mnn",
  group.by = c("study", "scATOMIC_pred", "mnn_clusters"),
  combine = T
)
ggsave("../plots/immune.umap.mnn.pdf", plot=p1, height=10, width=40, dpi=600)


blood <- FindNeighbors(blood, reduction = "integrated.scvi", dims = 1:30)
blood <- FindClusters(blood, resolution = 1, cluster.name = "scvi_clusters")
blood <- RunUMAP(blood, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
p1 <- DimPlot(
  blood,
  reduction = "umap.scvi",
  group.by = c("study", "scATOMIC_pred", "scvi_clusters"),
  combine = T
)
ggsave("../plots/immune.umap.scvi.pdf", plot=p1, height=10, width=40, dpi=600)

blood <- JoinLayers(blood)
saveRDS(blood, "immune_integrated_log.rds")

blood[["RNA"]] <- as(object = blood[["RNA"]], Class = "Assay")
library(sceasy)

sceasy::convertFormat(blood, from="seurat", to="anndata",
                       outFile='immune_integrated_log.h5ad')




#SaveH5Seurat(blood, filename = "blood_integrated.h5Seurat")
#Convert("blood_integrated.h5Seurat", dest = "h5ad")

