########### R script for combine datasets. ###########
########### Created by Danwen Qian on 28-07-2025. ###########
########### Last modified by Danwen Qian on 28-07-2025. ###########

library(Seurat)
library(biomaRt)
library(dplyr)
library(ggplot2)


setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets")

DSRCT<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/DSRCT/DSRCT_in_house.rds")
DSRCT<-subset(DSRCT,scDblFinder.class  == "singlet")

Eps<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Epithelioid_sarcoma/Epithelioid_sarcoma.rds")
ensembl_ids <- rownames(Eps)
ensembl_ids_clean <- gsub("\\.\\d+$", "", ensembl_ids)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotation <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_clean,
  mart = mart
)
annotation <- annotation[annotation$hgnc_symbol != "" & !is.na(annotation$hgnc_symbol),]
annotation <- annotation[!duplicated(annotation$ensembl_gene_id), ]
gene_map <- data.frame(
  original_id = ensembl_ids,
  clean_id = ensembl_ids_clean,
  stringsAsFactors = FALSE
)
gene_map <- left_join(gene_map, annotation, by = c("clean_id" = "ensembl_gene_id"))
gene_map_filtered <- gene_map %>% filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
gene_map_filtered <- gene_map_filtered[!duplicated(gene_map_filtered$hgnc_symbol), ]
Eps <- Eps[gene_map_filtered$original_id, ]
rownames(Eps) <- gene_map_filtered$hgnc_symbol
gene<-read.table("gene_Eps.txt",head=TRUE)
rownames(Eps) <- gene$gene
saveRDS(Eps,"/mnt/beegfs01/scratch/d_qian/SingleCell/DSRCT/DSRCT_in_house_f.rds")

Eps<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Epithelioid_sarcoma/Epithelioid_sarcoma_f.rds")
Eps<-subset(Eps,scDblFinder.class  == "singlet")

Gruel<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Gruel_NatCommun_2024/Gruel_liposarcoma.rds")
Gruel <- subset(Gruel, subset = !(sample_ID %in% c("CAR.nl","DET.nl","GAL.nl","MAG.lipoma","MAT.nl", "VAR.nl")))
gene<-read.table("gene_Gruel.txt",head=TRUE)
rownames(Gruel) <- gene$gene
Gruel<-subset(Gruel,scDblFinder.class  == "singlet")

Jerby<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Jerby_NatMed_2021/Jerby_synovial_sarcoma.rds")
gene<-read.table("gene_Jerby.txt",head=TRUE)
rownames(Jerby) <- gene$gene
Jerby<-subset(Jerby,scDblFinder.class  == "singlet")


Wu<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Wu_SciAdv_2022/Wu_MPNST.rds")
gene<-read.table("gene_Wu.txt",head=TRUE)
rownames(Wu) <- gene$gene
Wu<-subset(Wu,scDblFinder.class  == "singlet")
Wu<-subset(Wu,cancer_type=="MPNST")


Zhou<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Zhou_NatCommun_2020/Zhou_osteosarcoma.rds")
gene<-read.table("gene_Zhou.txt",head=TRUE)
rownames(Zhou) <- gene$gene
Zhou<-subset(Zhou,scDblFinder.class  == "singlet")

Liu<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Liu_FrontOncol_2021/Liu_osteosarcoma.rds")
gene<-read.table("gene_Liu.txt",head=TRUE)
rownames(Liu) <- gene$gene
Liu<-subset(Liu,scDblFinder.class  == "singlet")

Wei<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Wei_NatCancer_2022/Wei_rhabdomyosarcoma.rds")
gene<-read.table("gene_Wei.txt",head=TRUE)
rownames(Wei) <- gene$gene
Wei<-subset(Wei,scDblFinder.class  == "singlet")

Patel<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Patel_DevCell_2022/Patel_rhabdomyosarcoma.rds")
gene<-read.table("gene_Patel.txt",head=TRUE)
rownames(Patel) <- gene$gene
Patel<-subset(Patel,scDblFinder.class  == "singlet")

Taylor<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Taylor_Cancer_2025/Taylor_osteosarcoma.rds")
gene<-read.table("gene_Taylor.txt",head=TRUE)
rownames(Taylor) <- gene$gene
Taylor<-subset(Taylor,scDblFinder.class  == "singlet")

Goodspeed<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Goodspeed_ClinCancerRes_2025/Goodspeed_Ewing_sarcoma.rds")
gene<-read.table("gene_Goodspeed.txt",head=TRUE)
rownames(Goodspeed) <- gene$gene
Goodspeed<-subset(Goodspeed,scDblFinder.class  == "singlet")

Lu<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Lu_Oncogene_2024/Lu_UPS.rds")
gene<-read.table("gene_Lu.txt",head=TRUE)
rownames(Lu) <- gene$gene
Lu<-subset(Lu,scDblFinder.class  == "singlet")

in_house<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/inhouse/Myxoid_liposarcoma_CCS_in_house.rds")
gene<-read.table("gene_inhouse.txt",head=TRUE)
rownames(in_house) <- gene$gene
in_house<-subset(in_house,scDblFinder.class  == "singlet")

Leruste<-readRDS("/mnt/beegfs01/scratch/d_qian/SingleCell/Leruste_CancerCell_2019/Leruste_ATRT.rds")
Leruste<-subset(Leruste,scDblFinder.class  == "singlet")

combined<-merge(DSRCT,y=list(Eps,Gruel,Jerby,Wu,Zhou,Liu,Wei,Patel,Taylor,Goodspeed,Lu,in_house))
combined <- JoinLayers(combined)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(combined,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'study', ncol = 3)

df <- FetchData(combined, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "study"))

library(ggplot2)
library(patchwork)  

VlnPlot_custom <- function(object, features, group.by, ncol = 3, pt.size = 0.1) {
  plots <- list()
  
  for (f in features) {
    df <- FetchData(object, vars = c(f, group.by))
    colnames(df) <- c("value", "group")
    
    p <- ggplot(df, aes(x = group, y = value, fill = group)) +
      geom_violin(trim = TRUE, scale = "width") +
      geom_jitter(width = 0.2, size = pt.size, alpha = 0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = "", y = f, title = f)
    
    plots[[f]] <- p
  }
  
  wrap_plots(plots, ncol = ncol)
}

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
plot1<-VlnPlot(combined, features = features, group.by = "study", ncol = 3)
ggsave("../plots/QC_metrics.png", plot=plot1, height=8, width=15, dpi=600)
# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'study')
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'study')
ggsave("../plots/FeatureScatter.pdf", plot=plot1+plot2, height=6, width=15, dpi=600)

clinical<-read.csv("Clinical_information.csv")
combined$cells_ID<-rownames(combined@meta.data)
combined@meta.data<-left_join(combined@meta.data,clinical[,4:14],by="sample_ID")
meta<-read.csv("metadata_public.csv")
combined@meta.data<-left_join(combined@meta.data,meta,by="cells_ID")
rownames(combined@meta.data)<-combined$cells_ID
combined@meta.data<-combined@meta.data[,-13]

saveRDS(combined,"combined.rds")

#filter the data based on the QC metrics
combined<- subset(combined,subset= nFeature_RNA > 200&nFeature_RNA <8000 &nCount_RNA>300&nCount_RNA<100000&percent.mt<30)

saveRDS(combined,"combined_filter.rds")


#Normalizing the data and draw a PCA plot
cat ("Normalizing\n")
combined <- NormalizeData(combined,normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)


cat ("Addcellcycle\n")
load("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/cycle.rda")
combined <- CellCycleScoring(combined,
                                g2m.features = g2m_genes,
                                s.features = s_genes)


cat ("RunPCA \n")
combined <- RunPCA(combined,features = VariableFeatures(object = combined))
print(combined[["pca"]], dims = 1:5, nfeatures = 5)
plot1<-DimPlot(combined, reduction = "pca", group.by="study")
ggsave("../plots/original_PCA.png", plot=plot1, height=10, width=10, dpi=600)
plot1<-DimPlot(combined, reduction = "pca", group.by= "Phase",split.by = "Phase")
ggsave("../plots/original_PCA_cellcycle.png", plot=plot1, height=10, width=10, dpi=600)

cat ("RunUMAP \n")
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
plot1<-DimPlot(combined, reduction = "umap", group.by = "study")
ggsave("../plots/original_umap.png", plot=plot1, height=10, width=10, dpi=600)



