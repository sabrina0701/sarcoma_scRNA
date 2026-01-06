########### R script for creating seurat object. ###########
########### Created by Danwen Qian###########

library(Seurat)
library(stringr)
library(dplyr)
library(scDblFinder)
library(BiocParallel)
library(Matrix)
library(biomaRt)



##Gruel_NatCommun_2024 need to exclude benign lipoma and normal adipose tissues
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Gruel_NatCommun_2024")
data <- Read10X(data.dir = "/mnt/beegfs01/scratch/d_qian/SingleCell/Gruel_NatCommun_2024/downloads")
data <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
data$study<-"Gruel_liposarcoma"
data$cancer_type<-"Liposarcoma"
data$patient_ID<-sub("\\..*$", "", rownames(data@meta.data))
data$sample_ID<-sub("^(([^.]*\\.[^.]*)).*", "\\1", rownames(data@meta.data))
data$platform<-"10X"
data$tissue<-"tumor"
sce<-as.SingleCellExperiment(data)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
data <- AddMetaData(
  object = data,
  metadata = scDblFinder
)
saveRDS(data,"/mnt/beegfs01/scratch/d_qian/SingleCell/Gruel_NatCommun_2024/Gruel_liposarcoma.rds")

##Jerby_NatMed_2021
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Jerby_NatMed_2021/downloads")
counts <- read.csv("GSE131309_RAW/GSM3770932_SyS.tumors10x_counts.csv.gz", row.names = 1, check.names = FALSE)
data <- CreateSeuratObject(counts = counts, project = "sarcoma", min.cells = 3, min.features = 200)
data$study<-"Jerby_synovial_sarcoma"
data$cancer_type<-"Synovial_sarcoma"
data$patient_ID<-sub("\\..*$", "", rownames(data@meta.data))
data$sample_ID<-data$patient_ID
data$platform<-"10X"
data$tissue<-"tumor"
sce<-as.SingleCellExperiment(data)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
data <- AddMetaData(
  object = data,
  metadata = scDblFinder
)
saveRDS(data,"Jerby_synovial_sarcoma.rds")


##Wu_SciAdv_2022
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Wu_SciAdv_2022/downloads")
dir<-dir("/mnt/beegfs01/scratch/d_qian/SingleCell/Wu_SciAdv_2022/downloads/GSE179033_RAW")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("GSE179033_RAW/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Wu_MPNST"
  scRNAlist[[i]]$cancer_type<-sub("-.*", "", dir[i])
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Wu_MPNST.rds")


##Zhou_NatCommun_2020
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Zhou_NatCommun_2020")
dir<-dir("/mnt/beegfs01/scratch/d_qian/SingleCell/Zhou_NatCommun_2020/downloads/GSE152048")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE152048/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Zhou_osteosarcoma"
  scRNAlist[[i]]$cancer_type<-"Osteosarcoma"
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Zhou_osteosarcoma.rds")

##Epithelioid_sarcoma_in_house
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Epithelioid_sarcoma")
dir<-dir("/mnt/beegfs01/scratch/d_qian/SingleCell/Epithelioid_sarcoma/SWISARC_scRNAseq/scRNA")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data<- readMM(paste0("SWISARC_scRNAseq/scRNA/",dir[i],"/matrix.mtx"))
  gene_names = readLines(paste0("SWISARC_scRNAseq/scRNA/",dir[i],"/genes.txt"))
  barcode_names = readLines(paste0("SWISARC_scRNAseq/scRNA/",dir[i],"/barcodes.txt"))
  data <- Matrix::t(data)
  rownames(data) = gene_names
  colnames(data) = barcode_names
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Epithelioid_sarcoma_in_house"
  scRNAlist[[i]]$cancer_type<-"Epithelioid_sarcoma"
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Epithelioid_sarcoma.rds")


##DSRCT
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/DSRCT")
dir<-dir("downloads/GGSE263521_RAW")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE263521_RAW/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"DSRCT_in_house"
  scRNAlist[[i]]$cancer_type<-"DSRCT"
  scRNAlist[[i]]$patient_ID<-sub("_.*", "", dir[i])
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)

load("220320_DSRCT_PEZ_JV.RData")
sarc <- UpdateSeuratObject(sarc)
counts<-sarc@assays[["RNA"]]@counts
colnames(counts) <- paste0("PEZ_JV_", colnames(counts))
data <- CreateSeuratObject(counts = counts, project = "sarcoma", min.cells = 3, min.features = 200)
data$study<-"DSRCT_in_house"
data$cancer_type<-"DSRCT"
data$patient_ID<-"PEZ_JV"
data$sample_ID<-"PEZ_JV"
data$platform<-"10X"
data$tissue<-"tumor"

scRNA<-merge(scRNA,data)
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"DSRCT_in_house.rds")


## Liu_FrontOncol_2021
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Liu_FrontOncol_2021")
dir<-dir("downloads/GSE162454_RAW")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE162454_RAW/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Liu_osteosarcoma"
  scRNAlist[[i]]$cancer_type<-"Osteosarcoma"
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Liu_osteosarcoma.rds")


##Taylor_Cancer_2025
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Taylor_Cancer_2025")
dir<-dir("downloads/GSE299023_RAW")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE299023_RAW/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Taylor_osteosarcoma"
  scRNAlist[[i]]$cancer_type<-"Osteosarcoma"
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Taylor_osteosarcoma.rds")


##Wei_NatCancer_2022
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Wei_NatCancer_2022")
dir<-dir("downloads/GSE195709_RAW/human")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE195709_RAW/human/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Wei_rhabdomyosarcoma"
  scRNAlist[[i]]$cancer_type<-"Rhabdomyosarcoma"
  scRNAlist[[i]]$patient_ID <-"Patient_8"
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
scRNA$patient_ID[scRNA$sample_ID=="29806"]<-"Patient_9"
scRNA$patient_ID[scRNA$sample_ID=="20082"]<-"Patient_10"
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Wei_rhabdomyosarcoma.rds")


##Patel_DevCell_2022
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Patel_DevCell_2022")
dir<-dir("downloads/GSE174376_RAW/human")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/GSE174376_RAW/human/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Patel_rhabdomyosarcoma"
  scRNAlist[[i]]$cancer_type<-"Rhabdomyosarcoma"
  scRNAlist[[i]]$patient_ID<-sub("_.*", "", dir[i])
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Patel_rhabdomyosarcoma.rds")

##Goodspeed_ClinCancerRes_2025
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Goodspeed_ClinCancerRes_2025")
counts <-counts <- read.table("downloads/GSE277083/GSE277083_counts_matrix_primary.txt.gz",header = TRUE,row.names = 1,sep = "\t")
data <- CreateSeuratObject(counts = counts, project = "sarcoma", min.cells = 3, min.features = 200)
data$study<-"Goodspeed_Ewing_sarcoma"
data$cancer_type<-"Ewing_sarcoma"
meta_info <- read.table("downloads/GSE277083/GSE277083_cell_metadata_primary.txt",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)
rownames(meta_info) <- meta_info$Cell
data$patient_ID<-meta_info[colnames(data), "Sample"]
data$sample_ID<-data$patient_ID
data$platform<-"10X"
data$tissue<-"tumor"
data$scDblFinder.sample<-data$sample_ID
data$scDblFinder.class<-meta_info[colnames(data), "scDblFinder.class"]
saveRDS(data,"Goodspeed_Ewing_sarcoma.rds")

##Lu_Oncogene_2024
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Lu_Oncogene_2024")
dir<-dir("downloads")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("downloads/",dir[i]))
  colnames(data) <- paste0(dir[i], "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Lu_UPS"
  scRNAlist[[i]]$cancer_type<-"UPS"
  scRNAlist[[i]]$patient_ID<-dir[i]
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)
sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Lu_UPS.rds")


##inhouse dataset
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/inhouse")
dir<-dir("raw_mtx")
data <- Read10X(data.dir ="raw_mtx")
scRNA <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
scRNA$Sample_Name <- sub("^((?:[^_]+_){4}[^_]+).*", "\\1", colnames(scRNA))
scRNA$cells<-colnames(scRNA)

map_table <- data.frame(
  Sample_Name = c(
    "OBS_000_2c8a49_01_10x3p",
    "OBS_000_c7ce57_01_10x3p",
    "OBS_000_29a6f5_01_10x3p",
    "OBS_000_662b8f_01_10x3p",
    "OBS_000_d08107_01_10x3p",
    "OBS_000_27ff1a_01_10x3p",
    "OBS_000_ac4122_01_10x3p",
    "OBS_000_1ba747_01_10x3p",
    "OBS_000_dd315e_01_10x3p",
    "OBS_000_1c25f9_01_10x3p",
    "OBS_009_b538b6_01_10x3p"
  ),
study = c(
    "Myxoid_liposarcoma_in_house",
    "Myxoid_liposarcoma_in_house",
    "Myxoid_liposarcoma_in_house",
    "Myxoid_liposarcoma_in_house",
    "Myxoid_liposarcoma_in_house",
    "Myxoid_liposarcoma_in_house",
    "Clear_cell_sarcoma_in_house",
    "Clear_cell_sarcoma_in_house",
    "Clear_cell_sarcoma_in_house",
    "Clear_cell_sarcoma_in_house",
    "Clear_cell_sarcoma_in_house"
   ),
  cancer_type = c(
    "Myxoid_liposarcoma",
    "Myxoid_liposarcoma",
    "Myxoid_liposarcoma",
    "Myxoid_liposarcoma",
    "Myxoid_liposarcoma",
    "Myxoid_liposarcoma",
    "Clear_cell_sarcoma",
    "Clear_cell_sarcoma",
    "Clear_cell_sarcoma",
    "Clear_cell_sarcoma",
    "Clear_cell_sarcoma"
  ),
  patient_ID=c(
    "1-TTP338647",
    "1-TTP338647",
    "1-TTP11H10043",
    "1-TTP11H10043",
    "1-TTP09H07630",
    "1-TTP09H07630",
    "4-TTP08H00769",
    "1-TTP12H10602",
    "1-TTP17H08805",
    "1-TTP18H08442",
    "1-TTP19H05075"
  ),
  sample_ID=c(
    "1-TTP338647",
    "1-TTP339297",
    "1-TTP16H04854",
    "1-TTP11H10043",
    "1-TTP09H07630",
    "1-TTP10H02019",
    "4-TTP08H00769",
    "1-TTP12H10602",
    "1-TTP17H08805",
    "1-TTP18H08442",
    "1-TTP19H05075"
))
scRNA@meta.data <- scRNA@meta.data %>%
  left_join(map_table, by = "Sample_Name")
rownames(scRNA@meta.data) <- scRNA$cells
scRNA@meta.data <- scRNA@meta.data [,-c(4,5)]
scRNA$platform<-"10X"
scRNA$tissue<-"tumor"

sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Myxoid_liposarcoma_CCS_in_house.rds")


##Leruste_CancerCell_2019
setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/Leruste_CancerCell_2019")
dir<-dir("matrix")
scRNAlist <- list()
for(i in 1:length(dir)){
  print(dir[i])
  data <- Read10X(data.dir = paste0("matrix/",dir[i]))
  colnames(data) <- paste0(sub("_.*", "", dir[i]), "_", colnames(data))
  scRNAlist[[i]] <- CreateSeuratObject(counts = data, project = "sarcoma", min.cells = 3, min.features = 200)
  scRNAlist[[i]]$study<-"Leruste_ATRT"
  scRNAlist[[i]]$cancer_type<-"ATRT"
  scRNAlist[[i]]$patient_ID<-sub("_.*", "", dir[i])
  scRNAlist[[i]]$sample_ID<-dir[i]
  scRNAlist[[i]]$platform<-"10X"
  scRNAlist[[i]]$tissue<-"tumor"
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- JoinLayers(scRNA)

a<-readRDS("Seurat/sc10x.rna.seurat_INI287_Myeloid_QC_done.RDS")
colnames(a) <- paste0("INI287_Myeloid_", colnames(a))
a <- GetAssayData(a, assay = "RNA", slot = "counts")
a<-CreateSeuratObject(counts = a, project = "sarcoma", min.cells = 3, min.features = 200)

b<-readRDS("Seurat/sc10x.rna.seurat_INI287_T_cells_QC_done.RDS")
colnames(b) <- paste0("INI287_T_", colnames(b))
b <- GetAssayData(b, assay = "RNA", slot = "counts")
b <-CreateSeuratObject(counts = b, project = "sarcoma", min.cells = 3, min.features = 200)

c<-readRDS("Seurat/sc10x.rna.seurat_INI287_Tumor_cells_QC_done.RDS")
colnames(c) <- paste0("INI287_Tumor_", colnames(c))
c <- GetAssayData(c, assay = "RNA", slot = "counts")
c <-CreateSeuratObject(counts = c, project = "sarcoma", min.cells = 3, min.features = 200)

INI287<-merge(a,list(b,c))

a<-readRDS("Seurat/sc10x.rna.seurat_INI256_myeloid_QC_done.RDS")
colnames(a) <- paste0("INI256_Myeloid_", colnames(a))
a <- GetAssayData(a, assay = "RNA", slot = "counts")
a<-CreateSeuratObject(counts = a, project = "sarcoma", min.cells = 3, min.features = 200)

b<-readRDS("Seurat/sc10x.rna.seurat_INI256_T_cells_QC_done.RDS")
colnames(b) <- paste0("INI256_T_", colnames(b))
b <- GetAssayData(b, assay = "RNA", slot = "counts")
b <-CreateSeuratObject(counts = b, project = "sarcoma", min.cells = 3, min.features = 200)

c<-readRDS("Seurat/sc10x.rna.seurat_INI256_Tumor_cells_QC_done.RDS")
colnames(c) <- paste0("INI256_Tumor_", colnames(c))
c <- GetAssayData(c, assay = "RNA", slot = "counts")
c <-CreateSeuratObject(counts = c, project = "sarcoma", min.cells = 3, min.features = 200)

INI256<-merge(a,list(b,c))

INI262<-readRDS("Seurat/sc10x.rna.seurat_INI262_tumor_live_cells_QC_done.RDS")
colnames(INI262) <- paste0("INI262_", colnames(INI262))
INI262 <- GetAssayData(INI262, assay = "RNA", slot = "counts")
INI262 <-CreateSeuratObject(counts = INI262, project = "sarcoma", min.cells = 3, min.features = 200)

data<-merge(INI287,list(INI256,INI262))
data <- JoinLayers(data)
data$study<-"Leruste_ATRT"
data$cancer_type<-"ATRT"
data$patient_ID<-sub("_.*", "", colnames(data))
data$sample_ID<-data$patient_ID
data$platform<-"10X"
data$tissue<-"tumor"

scRNA<-merge(scRNA,data)
scRNA <- JoinLayers(scRNA)

sce<-as.SingleCellExperiment(scRNA)
sce<-scDblFinder(sce,samples="sample_ID")
scDblFinder<-data.frame(sce@colData[,c("scDblFinder.sample","scDblFinder.class")])
scRNA <- AddMetaData(
  object = scRNA,
  metadata = scDblFinder
)
saveRDS(scRNA,"Leruste_ATRT.rds")
