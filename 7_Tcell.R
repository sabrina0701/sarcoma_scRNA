setwd("/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggsci)
library(dittoSeq)
library(ggbeeswarm)
library(ggpubr)
options(future.globals.maxSize = 39*1024^3)


#Tcells
Tcell<-subset(combined_integrated,clusters_coarse=="T/NK_cells")
a<-subset(combined_integrated,clusters_coarse=="Proliferating")
a<-subset(a,CD3D>0|CD3E>0)
Tcell<-merge(Tcell,a)
Tcell <- JoinLayers(Tcell)
Tcell <- subset(Tcell, study !="Jerby_synovial_sarcoma")


###select log CCA for the downstream analysis
Tcell<-readRDS("Tcell.rds")
Tcell[["RNA"]] <- split(Tcell[["RNA"]], f = Tcell$study)
Tcell <- NormalizeData(Tcell)
Tcell<- FindVariableFeatures(Tcell)
Tcell <- ScaleData(Tcell)
Tcell <- RunPCA(Tcell)

Tcell <- IntegrateLayers(
  object = Tcell, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

saveRDS(Tcell,"Tcell_integrated.rds")

Tcell <- FindNeighbors(Tcell, reduction = "integrated.cca", dims = 1:30)
Tcell <- FindClusters(Tcell, resolution = 1, cluster.name = "cca_clusters")

##select resolution
getSSEfromExp <- function (exp, maxPC) {
  sumres <- 0
  for (i in levels(Idents(exp))) {       # Clusters
    for (j in 1:maxPC) {                 # PC numbers
      m <- mean (exp@reductions$pca@cell.embeddings[Idents(exp) == i,j])  # get mean for cluster, pc
      res <- exp@reductions$pca@cell.embeddings[Idents(exp) == i,j] - m   # vector of residuals
      sumres <- sumres + sum(res ** 2)   # add sum of SSE for this dimension
    }
  }
  sumres
}
getSSEfromExp(Tcell,30)

SSE<-read.csv("/Users/qiandanwen/Documents/postdoc in france/scRNA_sarcoma/SSE_results_Tcell.csv")
plot1<-ggplot(data = SSE, mapping = aes(x = cluster, y = Variance)) + geom_line()
ggsave("../plots/ElbowPlot_T.pdf", plot=plot1, height=4, width=4, dpi=600)


Tcell<-readRDS("Tcell_integrated.rds")
Tcell<- FindNeighbors(Tcell, reduction = "integrated.cca", dims = 1:30)
Tcell <- FindClusters(Tcell, resolution = 0.4, cluster.name = "cca_clusters_0.4")

Tcell <- RunUMAP(Tcell, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

p1 <- DimPlot(
  Tcell,
  reduction = "umap.cca",
  group.by = c("cca_clusters_0.4"),
  combine = T,label=T
)
ggsave("../plots/Tcell_UNAP_0.4.pdf", plot=p1, height=4, width=4, dpi=600)

Tcell <- JoinLayers(Tcell)

Idents(object = Tcell) <- "cca_clusters_0.4"

markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file="/mnt/beegfs01/scratch/d_qian/SingleCell/markers/markers_cca_Tcell.csv")
top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,file="/mnt/beegfs01/scratch/d_qian/SingleCell/markers/top20_cca_Tcell.csv")

markers.to.plot <- c("CD3D","CD3E","CD8A","CD8B","CD4","GZMK","CCL4","CCL5","CX3CR1","FGFBP2","FCGR3A","PRF1","CD40LG","IL7R","TCF7", "MAL", "LTB", "CCR7","SELL","LEF1","XCL1","XCL2","KLRD1","KLRC1","FOXP3","CTLA4","ENTPD1","CCR10","CCR4","CXCL13","HAVCR2","LAG3","TRDC","TRGV9","TRGC1","SLC4A10","ZBTB16","MKI67","STMN1","TRAC","BACH2","IFIT1","IFIT3","ISG15","MT-ND2","MT-ND1","MT-ND3","CD14","CD163")
p<-DotPlot(Tcell, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("../plots/Tcell_Bubble_heatmap_before.pdf", p, height=6, width=15, dpi=600)


Tcell <- subset(Tcell, cca_clusters_0.4 != 17)

#Tcell <- RenameIdents(object = Tcell,
"0" = "CD8_GZMK+_Tem",
"1" = "CD8_CX3CR1+_Tte",
"2" = "CD4_activated_Tm",
"3" = "CD4_Tcm",
"4" = "CD4/CD8_Tn",
"5" = "XCL1+_NK",
"6" = "CD4_activated_Treg",
"7" = "CD4_tissue-homing_Treg",
"8" = "CD8_Ttex",
"9" = "innate_like_gamma_delta_T",
"10" = "CD4/CD8_Proliferating",
"11" = "CD4/CD8_activated_T",
"12" = "CD4/CD8_Tisg",
"13" = "FGFBP2+_NK",
"14" = "mito_high_T",
"15" = "CD4_Tcm",
"16" = "CD8_GZMK+_Tem",
"18" = "CD8_Tcm")

Tcell$Tcell_subsets<-Tcell@active.ident

mypal =pal_d3("category20c")(20)

umap.df <- FetchData(Tcell, vars = c("umapcca_1", "umapcca_2", "Tcell_subsets"))
umap.cols <- DiscretePalette(length(unique(umap.df$Tcell_subsets)), palette = "polychrome")
ggrepel.df <- umap.df %>% group_by(Tcell_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean)
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(Tcell_subsets = ""))

p<-ggplot(umap.df, aes(x = umapcca_1, y = umapcca_2)) +
geom_point(aes(fill = Tcell_subsets), size = 1,  pch = 21) +
geom_label_repel(umap.df %>% group_by(Tcell_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean), size = 2, mapping = aes(x = umapcca_1, y = umapcca_2, label = Tcell_subsets)) +guides(colour = guide_legend(override.aes = list(size=2))) + theme_classic()  + NoLegend()+theme(panel.border =  element_rect(colour = "black", fill = NA, size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(text = element_text(size = 10))+scale_fill_manual(values = mypal)

ggsave("../plots/Tcells_custom_UMAP.pdf", p, height=4, width=4, dpi=600)

mypal1<-c("#3182BDFF", "#E6550DFF","#756BB1FF", "#31A354FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF", "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF","#BDBDBDFF", "#C6DBEFFF", "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")
cell_number<-data.frame(table(Tcell$Tcell_subsets))
cell_number <- cell_number[order(cell_number$Freq, decreasing = TRUE), ]
piepercent<- round(100*cell_number$Freq/sum(cell_number$Freq), 1)
labels<-paste0(cell_number$Var1," ",piepercent,"%")
pdf(file = "../plots/pie_chart_cells_T.pdf")
pie(cell_number$Freq, labels=labels,col =mypal,cex=1)
dev.off()

## bubble plot
Idents(Tcell) <- factor(Idents(Tcell), levels = c("CD4/CD8_Tn","CD4_Tcm","CD4_activated_Tm","CD4_activated_Treg","CD4_tissue-homing_Treg","CD8_GZMK+_Tem","CD8_CX3CR1+_Tte","CD8_Ttex", "CD8_Tcm","CD4/CD8_Proliferating","CD4/CD8_activated_T", "CD4/CD8_Tisg","mito_high_T","XCL1+_NK","FGFBP2+_NK","innate_like_gamma_delta_T"))

markers.to.plot <- c("CD3D","CD3E","CD8A","CD8B","CD4","TCF7","MAL","LTB","CCR7","SELL","LEF1","CD40LG","IL7R","FOXP3","CTLA4","ENTPD1","CCR10","CCR4","GZMK", "CCL4","CCL5","CX3CR1","FGFBP2","FCGR3A","PRF1","CXCL13","HAVCR2","LAG3","MKI67","STMN1","TRAC","BACH2","IFIT1","IFIT3","ISG15","XCL1","XCL2","KLRD1","KLRC1","TRDC","TRGV9","TRGC1","SLC4A10","ZBTB16","MT-ND2","MT-ND1","MT-ND3")
p<-DotPlot(Tcell, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("../plots/Tcell_Bubble_heatmap.pdf", p, height=6, width=15, dpi=600)



##NK-like T cells (T_NK_like), comprising both TRDC⁺ γδ T cells and CD8⁺ T cells expressing an NK-like transcriptional program (XCL1, XCL2, KLRD1, KLRF1, KLRC1, FCER1G, TYROBP).”
#Activated Tregs (Treg_activated) expressing FOXP3, IL2RA, CTLA4, LAYN, TIGIT and ENTPD1; and Tissue-homing Tregs (Treg_CCR10) characterized by CCR10, CCR4, PI16, IL6R and HPGD expression.
#“Conventional activated T cells (T_conv_activated), comprising both CD4⁺ and CD8⁺ T cells and characterized by TRAC/TRBC, ICOS, CD69 and TNFAIP3 expression but lacking clear cytotoxic, regulatory, exhausted or γδ/NK-like signatures.”


##"CCR7","SELL" naive;
#"FOXP3", Treg;
#"GZMK", "CCL5", CD8_GZMK+_Tem;
#"CD40LG","IL7R", CD4_Tm;
#"FCGR3A","CX3CR1", FCGR3A_NK;
#"NCAM1","KLRC1", NCAM1_NK;
#"CXCR3","ZNF683", CD8_Trm;
#"KLRG1", CD8_Temra;
#"CCL20","CCR6","RORC" CD4_Th17;
#KLRG1,"RORA","CXCR6","KLRB1","SLC4A10",CD8_Tc17;
#"IFIT3","IFIT1", CD8_Tisg;
#AHR、CD83、TNFRSF4、TNFRSF18、TOX2 和 MAFF,ICOS, CD4 Tfh
#"IFNG",TNFRSF9 "CD8_Teff",
#"CD44", memory T cells

#Tcm ：IL7R、GPR183,Tem ：CD40LG、ICOS、TNFRSF25、TRADD
#MAIT:"SLC4A10","ZBTB16"

plot<-dittoBarPlot(Tcell, "Tcell_subsets", group.by = "cancer_type",var.labels.reorder=c(10,9,3,1,8,16,2,4,12,14,6,5,7,13,15,11))+scale_fill_manual(values = mypal1)
ggsave("../plots/Tcell_barplot_cancers.pdf", plot, height=6, width=7, dpi=600)

plot<-dittoBarPlot(Tcell, "Tcell_subsets", group.by = "study",var.labels.reorder=c(10,9,3,1,8,16,2,4,12,14,6,5,7,13,15,11))+scale_fill_manual(values = mypal1)
ggsave("../plots/Tcell_barplot_study.pdf", plot, height=6, width=7, dpi=600)

a <- subset(Tcell, subset = Tumour_site != "")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_site",var.labels.reorder=c(10,9,3,1,8,16,2,4,12,14,6,5,7,13,15,11))+scale_fill_manual(values = mypal1)
ggsave("../plots/Tcell_barplot_Tumour_site.pdf", plot, height=6, width=7, dpi=600)

a <- subset(Tcell, subset = Sample_type != "")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Sample_type",var.labels.reorder=c(10,9,3,1,8,16,2,4,12,14,6,5,7,13,15,11))+scale_fill_manual(values = mypal1)
ggsave("../plots/Tcell_barplot_Sample_type.pdf", plot, height=6, width=6, dpi=600)

a <- subset(Tcell, subset = Treatment != "")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Treatment",var.labels.reorder=c(10,9,3,1,8,16,2,4,12,14,6,5,7,13,15,11))+scale_fill_manual(values = mypal1)
ggsave("../plots/immune_barplot_Treatment.pdf", plot, height=6, width=6, dpi=600)

###barbox
proportion<-as.data.frame.array(prop.table(table(Tcell$sample_ID,Tcell$Tcell_subsets),margin=1))
proportion$sample_ID<-rownames(proportion)
tumor_type<-Tcell@meta.data[,c(4:7,13:20)]
tumor_type<-tumor_type %>% distinct(sample_ID, .keep_all = TRUE)
proportion<-merge(proportion,tumor_type,by="sample_ID")

plot<-ggplot(data=proportion,aes_string(y="CD8_Ttex",x="cancer_type",color="cancer_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/Tcell_boxplot_cancers_Ttex.pdf"), plot, height=5, width=12, dpi=600)


plot<-ggplot(data=proportion,aes_string(y="CD8_Ttex",x="study",color="study"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/Tcell_boxplot_study_Ttex.pdf"), plot, height=5, width=12, dpi=600)

a <- subset(proportion, subset = Tumour_site != "")
plot<-ggplot(data=a,aes_string(y="CD8_Ttex",x="Tumour_site",color="Tumour_site"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/Tcell_boxplot_Tumour_site_Ttex.pdf"), plot, height=5, width=12, dpi=600)

a <- subset(proportion, subset = Sample_type != "")
plot<-ggplot(data=a,aes_string(y="CD8_Ttex",x="Sample_type",color="Sample_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/Tcell_boxplot_Sample_type_Tcell.pdf"), plot, height=5, width=6, dpi=600)

a <- subset(proportion, subset = Treatment != "")
plot<-ggplot(data=a,aes_string(y="CD8_Ttex",x="Treatment",color="Treatment"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/Tcell_boxplot_Treatment_Tcell.pdf"), plot, height=5, width=8, dpi=600)


##pathway analysis
#df <- read.csv("/Users/qiandanwen/Downloads/plots/monocle3/pathway_markers.csv", stringsAsFactors = FALSE)
df <- read.csv("/mnt/beegfs01/scratch/d_qian/SingleCell/markers/pathway_markers.csv", stringsAsFactors = FALSE)

marker.list <- lapply(df, function(col) {
    genes <- col[!is.na(col) & col != ""]
    return(genes)
})

Tcell <-  AddModuleScore(Tcell,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")

for(i in 1:length(marker.list)){
    colnames(Tcell@meta.data)[colnames(Tcell@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}

Differentiation <- c("naive", "Activation.Effector.function", "Exhaustion")
Function <- c("TCR.Signaling", "Cytotoxicity", "Cytokine.Cytokine.receptor",
              "Chemokine.Chemokine.receptor", "Senescence", "Anergy",
              "NFKB.Signaling", "Stress.response", "MAPK.Signaling", "Adhesion",
              "IFN.Response","Treg.signature", "Costimulatory.molecules")
Metabolism <- c("Oxidative.phosphorylation", "Glycolysis",  "Lipid.metabolism")
Apoptosis <- c("Pro.apoptosis", "Anti.apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)

FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(Tcell$Tcell_subsets)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- levels(Idents(Tcell))
rownames(FunctionScoreMatrix) <- MarkerNameVector

for(ci in 1:ncol(FunctionScoreMatrix)){
    for(ri in 1:nrow(FunctionScoreMatrix)){
        FunctionVec <- as_tibble(Tcell@meta.data) %>% pull(MarkerNameVector[ri])
        fv <- mean(FunctionVec[Tcell$Tcell_subsets == levels(Idents(Tcell))[ci]])
        FunctionScoreMatrix[ri, ci] <- fv
    }
}

library(scales)
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))


my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
    colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
signatureType_row <- data.frame(Signature.type = c(
                                    rep("Differentiation", length(Differentiation)),
                                    rep("Function", length(Function)),
                                    rep("Metabolism", length(Metabolism)),
                                    rep("Apoptosis", length(Apoptosis))))

rownames(signatureType_row) <- MarkerNameVector

library(pheatmap)
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 16, 19),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 8,
         height = 6,
         filename = file.path("../plots/Tcell_FunctionScore_heatmap.pdf"))

