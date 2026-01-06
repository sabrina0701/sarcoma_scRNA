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

combined_integrated<-readRDS("immune_integrated_log.rds")

combined_integrated <- FindNeighbors(combined_integrated, reduction = "integrated.cca", dims = 1:30)
combined_integrated <- FindClusters(combined_integrated, resolution = 0.4, cluster.name = "cca_clusters_0.4")

combined_integrated <- RunUMAP(combined_integrated, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
combined_integrated,
  reduction = "umap.cca",
  group.by = c( "cca_clusters_0.4"),
  combine = T, label=T
)
ggsave("../plots/immune.umap.cca0.4.pdf", plot=p1, height=10, width=10, dpi=600)

Idents(object = combined_integrated) <- "cca_clusters_0.4"

markers.to.plot <- c("PTPRC","CD3D","NCAM1","XCL1","CD68","ITGAM","CD19","CD79A","IGHG1","JCHAIN","XCR1","CLEC9A","CCR7","LAMP3","CLEC10A","CD1C","LILRA4","CLEC4C","TPSB2","KIT","MKI67","STMN1")
p<-DotPlot(combined_integrated, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("../plots/Bubble_heatmap_before.pdf", p, height=6, width=12, dpi=600)

p1<-FeaturePlot(combined_integrated,features=c("PTPRC","CD3D","NCAM1","XCL1","CD68","ITGAM","CD19","CD79A","IGHG1","JCHAIN","TPSB2","KIT","MKI67","STMN1"),ncol=4,cols = c("lightgrey","coral2"),reduction = "umap.cca")

ggsave("../plots/FeaturePlot_cca_log.pdf", p1, height=15, width=20, dpi=600)

markers <- FindAllMarkers(combined_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file="/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/markers/immune_markers_cca.csv")
top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,file="/mnt/beegfs01/scratch/d_qian/SingleCell/combined_datasets/markers/innune_top20_cca.csv")


combined_integrated <- RenameIdents(object = combined_integrated,
"0" ="Myeloids",
"1" = "T/NK_cells",
"2" ="Myeloids",
"3" ="T/NK_cells",
"4" ="Myeloids",
"5" ="T/NK_cells",
"6" ="T/NK_cells",
"7" ="Proliferating",
"8" ="T/NK_cells",
"9" ="Myeloids",
"10"="B_cells/Plasmablasts",
"11"="Myeloids",
"12"="Myeloids",
"13"="T/NK_cells",
"14"="Mast_cells",
"15"="T/NK_cells",
"16"="Myeloids",
"17"="Myeloids",
"18"="Myeloids",
"19"="T/NK_cells"
)

combined_integrated$clusters_coarse<-combined_integrated@active.ident

mypal =pal_d3("category20c")(20)

umap.df <- FetchData(combined_integrated, vars = c("umapcca_1", "umapcca_2", "clusters_coarse"))
umap.cols <- DiscretePalette(length(unique(umap.df$clusters_coarse)), palette = "polychrome")
ggrepel.df <- umap.df %>% group_by(clusters_coarse) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean)
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(clusters_coarse = ""))
p<-ggplot(umap.df, aes(x = umapcca_1, y = umapcca_2)) +
geom_point(aes(fill = clusters_coarse), size = 1,  pch = 21)+
geom_label_repel(umap.df %>% group_by(clusters_coarse) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean), size = 4, mapping = aes(x = umapcca_1, y = umapcca_2, label = clusters_coarse)) +guides(colour = guide_legend(override.aes = list(size=2))) + theme_classic() + NoLegend() +theme(panel.border =  element_rect(colour = "black", fill = NA, size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(text = element_text(size = 15))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_custom_UMAP_log.pdf", p, height=6, width=6, dpi=600)

cell_number<-data.frame(table(combined_integrated$clusters_coarse))
cell_number <- cell_number[order(cell_number$Freq, decreasing = TRUE), ]
piepercent<- round(100*cell_number$Freq/sum(cell_number$Freq), 1)
labels<-paste0(cell_number$Var1," ",piepercent,"%")
pdf(file = "../plots/immune_pie_chart_cells_annotated.pdf")
pie(cell_number$Freq, labels=labels,col =mypal,cex=1)
dev.off()

Idents(combined_integrated) <- factor(Idents(combined_integrated), levels = c("Myeloids","T/NK_cells","B_cells/Plasmablasts","Mast_cells","Proliferating"))

markers.to.plot <- c("PTPRC","CD68","ITGAM","CD3D","NCAM1","XCL1","CD19","CD79A","IGHG1","JCHAIN","TPSB2","KIT","MKI67","STMN1")
p<-DotPlot(combined_integrated, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("../plots/immune_Bubble_heatmap.pdf", p, height=6, width=10, dpi=600)

p1<-FeaturePlot(combined_integrated,features=c("PTPRC","CD68","ITGAM","CD3D","NCAM1","XCL1","CD19","CD79A","IGHG1","JCHAIN","TPSB2","KIT","MKI67","STMN1"),ncol=4,cols = c("lightgrey","coral2"),reduction = "umap.cca")
ggsave("../plots/immune_FeaturePlot.pdf", p1, height=15, width=20, dpi=600)

plot<-dittoBarPlot(combined_integrated, "clusters_coarse", group.by = "cancer_type",var.labels.reorder=c(3,5,4,1,2))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_barplot_cancers.pdf", plot, height=4, width=6, dpi=600)

plot<-dittoBarPlot(combined_integrated, "clusters_coarse", group.by = "study",var.labels.reorder=c(3,5,4,1,2))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_barplot_study.pdf", plot, height=4, width=6, dpi=600)

a <- subset(combined_integrated, subset = Tumour_site != "")
plot<-dittoBarPlot(a, "clusters_coarse", group.by = "Tumour_site",var.labels.reorder=c(3,5,4,1,2))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_barplot_Tumour_site.pdf", plot, height=4, width=6, dpi=600)

a <- subset(combined_integrated, subset = Sample_type != "")
plot<-dittoBarPlot(a, "clusters_coarse", group.by = "Sample_type",var.labels.reorder=c(3,5,4,1,2))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_barplot_Sample_type.pdf", plot, height=4, width=6, dpi=600)

a <- subset(combined_integrated, subset = Treatment != "")
plot<-dittoBarPlot(a, "clusters_coarse", group.by = "Treatment",var.labels.reorder=c(3,5,4,1,2))+scale_fill_manual(values = mypal)
ggsave("../plots/immune_barplot_Treatment.pdf", plot, height=4, width=6, dpi=600)

###barbox
proportion<-as.data.frame.array(prop.table(table(combined_integrated$sample_ID,combined_integrated$clusters_coarse),margin=1))
proportion$sample_ID<-rownames(proportion)
tumor_type<-combined_integrated@meta.data[,c(4:7,13:20)]
tumor_type<-tumor_type %>% distinct(sample_ID, .keep_all = TRUE)
proportion<-merge(proportion,tumor_type,by="sample_ID")

plot<-ggplot(data=proportion,aes_string(y="Myeloids",x="cancer_type",color="cancer_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/immune_boxplot_cancers_Myeloids.pdf"), plot, height=5, width=12, dpi=600)


plot<-ggplot(data=proportion,aes_string(y="Myeloids",x="study",color="study"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/immune_boxplot_study_Myeloids.pdf"), plot, height=5, width=12, dpi=600)

a <- subset(proportion, subset = Tumour_site != "")
plot<-ggplot(data=a,aes_string(y="Myeloids",x="Tumour_site",color="Tumour_site"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/immune_boxplot_Tumour_site_Myeloids.pdf"), plot, height=5, width=12, dpi=600)

a <- subset(proportion, subset = Sample_type != "")
plot<-ggplot(data=a,aes_string(y="Myeloids",x="Sample_type",color="Sample_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/immune_boxplot_Sample_type_Myeloids.pdf"), plot, height=5, width=12, dpi=600)

a <- subset(proportion, subset = Treatment != "")
plot<-ggplot(data=a,aes_string(y="Myeloids",x="Treatment",color="Treatment"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust = 1),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("../plots/immune_boxplot_Treatment_Myeloids.pdf"), plot, height=5, width=12, dpi=600)


##SCT
combined_integrated<-readRDS("immune_integrated_SCT.rds")

combined_integrated <- FindNeighbors(combined_integrated, reduction = "integrated.rpca", dims = 1:30)
combined_integrated <- FindClusters(combined_integrated, resolution = 0.4, cluster.name = "rpca_clusters_0.4")

combined_integrated <- RunUMAP(combined_integrated, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p1 <- DimPlot(
combined_integrated,
  reduction = "umap.rpca",
  group.by = c( "rpca_clusters_0.4"),
  combine = T, label=T
)
ggsave("../plots/immune.umap.sct.rpca0.4.pdf", plot=p1, height=10, width=10, dpi=600)





#Tcells
Tcell<-subset(combined_integrated,clusters_coarse=="T/NK_cells")
a<-subset(combined_integrated,clusters_coarse=="Proliferating")
a<-subset(a,CD3D>0|CD3E>0)
Tcell<-merge(Tcell,a)
Tcell <- JoinLayers(Tcell)


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

Tcell <- JoinLayers(Tcell)

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

SSE<-read.csv("SSE_results_Tcell_all.csv")
plot1<-ggplot(data = SSE, mapping = aes(x = cluster, y = Variance)) + geom_line()
ggsave("plots/ElbowPlot_T_all.pdf", plot=plot1, height=4, width=4, dpi=600)


Tcell<-readRDS("Tcells_all.rds")
Tcell<- FindNeighbors(Tcell, reduction = "integrated.cca", dims = 1:30)
Tcell <- FindClusters(Tcell, resolution = 0.3, cluster.name = "cca_clusters_0.3")

Tcell <- RunUMAP(Tcell, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

p1 <- DimPlot(
  Tcell,
  reduction = "umap.cca",
  group.by = c("cca_clusters_0.3"),
  combine = T
)


Idents(object = Tcell) <- "cca_clusters_0.3"

markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file="/SAN/colcc/NexTGen/WP6/scRNA/markers/markers_cca_Tcell.csv")
top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,file="/SAN/colcc/NexTGen/WP6/scRNA/markers/top20_cca_Tcell.csv")

cat("gene set enrichment analysis")
library(msigdbr)
library(fgsea)
library(tibble)
library(pheatmap)

m_df <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")## 定义KEGG基因集
m_df <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")## 定义GO基因集
m_df <- msigdbr(species = "Homo sapiens", category = "C8")##C7/C8

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

selected_names <- names(fgsea_sets)[grep("_T_CELL|_NK_", names(fgsea_sets))]
filtered_list <- fgsea_sets[selected_names]

fgseaRes<-data.frame()
for(i in c(0:15)){
cluster.genes<- markers %>% dplyr::filter(cluster == i) %>%arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC)%>% distinct(gene, .keep_all = TRUE) #基因按logFC排序
ranks<- deframe(cluster.genes)
fgseaRes1<- fgsea(filtered_list, stats = ranks, nperm = 1000) #运行fgsea
fgseaRes1$cluster<-i
fgseaRes<-bind_rows(fgseaRes, fgseaRes1)
}

pathway <- fgseaRes %>%
    select(pathway, NES, cluster) %>%
    pivot_wider(names_from = cluster, values_from = NES,names_prefix = "Cluster_")
pathway<-data.frame(pathway)
pathway[is.na(pathway)] <- 0
rownames(pathway)<-pathway$pathway
pathway<-pathway[,-1]

p <- pheatmap(pathway,
cluster_rows = T,
cluster_cols = T,
show_rownames = T,
show_colnames = T,
fontsize = 8)
ggsave("plots/Tcell_C8_heatmap.pdf", p, height=12, width=8, dpi=600)


FeaturePlot(Tcell,features=c("CCR7","LEF1","SELL"),ncol=4,cols = c("lightgrey","coral2"),reduction = "umap.cca")

p1<-VlnPlot(Tcell, features=c("CCR7","SELL","FOXP3","GZMK","IFNG","CCL5","CD40LG","IL7R","FCGR3A","PRF1","NCAM1","XCL1","CXCR3","ZNF683","KLRG1","CCL20","CCR6","RORA","SLC4A10","KLRB1","IFIT3","IFIT1","ICOS","TOX2","MAFF","MT-ATP6","HSPA1A"),group.by = "cca_clusters_0.3",ncol=2)
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
#Temra细胞通常通过以下几个维度表现出显著的基因特征： 细胞毒性基因（如 GZMB, PRF1, GZMK）趋化因子（如 CCL5, CX3CR1） NK细胞标志基因（如 KLRD1, KLRF1）衰竭相关基因（如 LAG3, TIM-3）,FGFBP2/CX3CR1/PRF1
#TRGC1, TRDC, NCR3, KLRG1, 和 CXCR6 等基因的细胞更可能是 γδ T细胞
#Tcm ：IL7R、GPR183,Tem ：CD40LG、ICOS、TNFRSF25、TRADD
#MAIT:"SLC4A10","ZBTB16"

Tcell <- RenameIdents(object = Tcell,
"0" ="CD4/CD8_Tn",
"1" = "CD8_GZMK+_Tem",
"2" = "CD4_Tcm",
"3" = "CD8_Trm",
"4" = "CD8_Trm",
"5" = "NCAM1_NK",
"6" = "CD8_Temra",
"7" = "dying_cells",
"8" = "CD4_Th17",
"9" = "CD8_Tc17",
"10"= "CD8_Tisg",
"11"="CD4_Treg",
"12"= "CD4_Tfh",
"13"= "unknown"
)

Tcell <- RenameIdents(object = Tcell,
"0" ="CD8_GZMK+_Tem",
"1" = "CD4/CD8_Tn",
"2" = "CD4_Tm",
"3" = "CD8_Trm",
"4" = "FCGR3A+_NK",
"5" = "NCAM1+_NK",
"6" = "CD8_TNFRSF9+_Tem",
"7" = "CD4/CD8_Tn",
"8" = "CD4_Th17",
"9" = "Stressed_T",
"10"= "CD8_KLRG1+_Tem",
"11"="CD8_Tc17",
"12"= "CD8_Tisg",
"13"= "CD4_Treg",
"14"=
"15"=
)

Tcell <- RenameIdents(object = Tcell,
"0" = "CD8_GZMK+_Tem",
"1" = "CD4/CD8_Tn",
"2" = "CD4_Tem",
"3" = "CD8_KLRG1+_Tem",
"4" = "CD8_Trm",
"5" = "Stressed_T",
"6" = "NK",
"7" = "CD8_Temra",
"8" = "CD8_TNFRSF9+_Tem",
"9" = "CD4_Th17-like",
"10"= "CD4/CD8_Proliferating",
"11"= "MAIT/gamma_delta_T",
"12"= "CD8_Tisg",
"13"= "CD4_Treg",
"14"= "CD4_Tcm",
"15"= "gamma_delta_T"
)

Tcell$Tcell_subsets<-Tcell@active.ident

mypal =pal_d3("category20c")(20)

umap.df <- FetchData(Tcell, vars = c("umapcca_1", "umapcca_2", "Tcell_subsets"))
umap.cols <- DiscretePalette(length(unique(umap.df$Tcell_subsets)), palette = "polychrome")
ggrepel.df <- umap.df %>% group_by(Tcell_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean)
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(Tcell_subsets = ""))

p<-ggplot(umap.df, aes(x = umapcca_1, y = umapcca_2)) +
geom_point(aes(fill = Tcell_subsets), size = 1.5,  pch = 21) +
geom_label_repel(umap.df %>% group_by(Tcell_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean), size = 2, mapping = aes(x = umapcca_1, y = umapcca_2, label = Tcell_subsets)) +guides(colour = guide_legend(override.aes = list(size=2))) + theme_classic()  + NoLegend()+theme(panel.border =  element_rect(colour = "black", fill = NA, size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(text = element_text(size = 10))+scale_fill_manual(values = mypal)

ggsave("plots/custom_UMAP_Tcells.pdf", p, height=4, width=4, dpi=600)

cell_number<-data.frame(table(Tcell$Tcell_subsets))
cell_number <- cell_number[order(cell_number$Freq, decreasing = TRUE), ]
piepercent<- round(100*cell_number$Freq/sum(cell_number$Freq), 1)
labels<-paste0(cell_number$Var1," ",piepercent,"%")
pdf(file = "plots/pie_chart_cells_T.pdf")
pie(cell_number$Freq, labels=labels,col =mypal,cex=1)
dev.off()

##Figure1e bubble plot
Idents(Tcell) <- factor(Idents(Tcell), levels = c("CD4/CD8_Tn","CD4_Tem","CD4_Tcm", "CD4_Th17-like","CD4_Treg","CD8_GZMK+_Tem", "CD8_KLRG1+_Tem", "CD8_TNFRSF9+_Tem", "CD8_Trm","CD8_Temra","CD8_Tisg","MAIT/gamma_delta_T","Stressed_T","CD4/CD8_Proliferating","NK","gamma_delta_T"))


markers.to.plot <- c("CD3E","CD8A","CD4","CCR7","SELL","CD44","CD40LG","ICOS","IL7R","GPR183","CCR6","CCL20","IL17A","RORC","FOXP3","GZMK","CCL4L2","CCL5","KLRG1","GNLY","TNFRSF9","ZNF683","XCL1","FGFBP2","CX3CR1","IFIT3","IFIT1","SLC4A10","ZBTB16","TRGC1", "TRDC", "KLRD1", "KLRF1", "KLRC1","HSPA6","DNAJA4","MKI67","STMN1")
p<-DotPlot(Tcell, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("plots/Bubble_heatmap_T.pdf", p, height=6, width=12, dpi=600)

markers.to.plot <- c("ENTPD1","ITGAE","TNFRSF9","CXCL13","PDCD1","HAVCR2","LAG3","CTLA4","TIGIT","GZMB","IFNG","TNF")
p<-DotPlot(Tcell, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("plots/Bubble_tumor_reactive_T.pdf", p, height=6, width=8, dpi=600)

proportion<-as.data.frame.array(prop.table(table(Tcell$sample_ID,Tcell$Tcell_subsets),margin=1))
proportion$sample_ID<-rownames(proportion)
tumor_type<-Tcell@meta.data[,c(4:7,9,13:15,19)]
tumor_type<-tumor_type %>% distinct(sample_ID, .keep_all = TRUE)
proportion<-merge(proportion,tumor_type,by="sample_ID")

plot<-dittoBarPlot(Tcell, "Tcell_subsets", group.by = "cancer_type",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_cancers_Tcells.pdf", plot, height=5, width=6, dpi=600)
plot<-ggplot(data=proportion,aes_string(y="`CD4/CD8_Tn`",x="cancer_type",color="cancer_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_cancers_Tcells_naive.pdf"), plot, height=5, width=12, dpi=600)

plot<-dittoBarPlot(Tcell, "Tcell_subsets", group.by = "study",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))+scale_fill_manual(values = mypal)
ggsave("plots/barplot_study_Tcells.pdf", plot, height=5, width=6, dpi=600)
plot<-ggplot(data=proportion,aes_string(y="`CD4/CD8_Tn`",x="study",color="study"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_study_Tcells_naive.pdf"), plot, height=5, width=12, dpi=600)

plot<-dittoBarPlot(Tcell, "Tcell_subsets", group.by = "Time_point",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_Time_point_Tcell.pdf", plot, height=5, width=6, dpi=600)
plot<-ggplot(data=proportion %>% filter(!is.na(Time_point)), aes_string(y="`CD4/CD8_Tn`",x="Time_point",color="Time_point"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Time_point_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell,Tumour_site %in% c("primary","metastatic"))
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_site",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_Tumour_site_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,Tumour_site %in% c("primary","metastatic"))
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Tumour_site",color="Tumour_site"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_site_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell, Treatment %in% c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Treatment",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))+scale_x_discrete(limits = c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
ggsave("plots/barplot_Treatment_Tcell.pdf", plot, height=6, width=6, dpi=600)
a<-subset(proportion,Treatment %in% c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
a$Treatment <- factor(a$Treatment, levels = c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))  # 按照你想要的顺序替换各组名称
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Treatment",color="Treatment"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Treatment_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)
plot<-ggplot(data=a, aes_string(y="`CD8_GZMK+_Tem`",x="Treatment",color="Treatment"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Treatment_Tcells_GZMK.pdf"), plot, height=5, width=6, dpi=600)


a<-subset(Tcell,cancer_type=="Ependymoma")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(4,1,2,5,8,6,7,10,12,3,11,9))
ggsave("plots/barplot_Tumour_subgroup_ENPs_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Ependymoma")
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_ENPs_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)


a<-subset(Tcell,cancer_type=="Hepatoblastoma")
plot<-dittoBarPlot(a,"Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_Tumour_subgroup_HB_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Hepatoblastoma")
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_HB_Tcells_navie.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell,cancer_type=="Medulloblastoma")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(4,1,2,5,8,6,7,10,12,3,11,9))
ggsave("plots/barplot_Tumour_subgroup_Medulloblastoma_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Medulloblastoma")
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Medulloblastoma_Tcells_navie.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell,cancer_type=="Neuroblastoma")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_Tumour_subgroup_NB_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Neuroblastoma")
plot<-ggplot(data=a, aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Neuroblastoma_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell,cancer_type=="Rhabdomyosarcoma")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(6,2,1,3,4,7,8,11,12,9,10,14,16,5,15,13))
ggsave("plots/barplot_Tumour_subgroup_RMS_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Rhabdomyosarcoma")
plot<-ggplot(data=a %>% filter(!is.na(Tumour_subgroup)), aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Rhabdomyosarcoma_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Tcell,cancer_type=="Wilms_tumor")
plot<-dittoBarPlot(a, "Tcell_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(4,1,2,5,8,6,7,10,12,3,11,9))
ggsave("plots/barplot_Tumour_subgroup_Wilms_tumor_Tcell.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Wilms_tumor")
plot<-ggplot(data=a %>% filter(!is.na(Tumour_subgroup)), aes_string(y="`CD4/CD8_Tn`",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Wilms_tumor_Tcells_naive.pdf"), plot, height=5, width=6, dpi=600)

##pathway analysis
df <- read.csv("pathway_markers.csv", stringsAsFactors = FALSE)
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
Metabolism <- c("Oxidative.phosphorylation", "Glycolysis", "Fatty.acid.metabolism","OXPHOS", "Lipid.metabolism")
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
         gaps_row = c(3, 16, 21),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 8,
         height = 6,
         filename = file.path("plots/Tcell_FunctionScore_heatmap.pdf"))



Myeloid<-subset(combined_integrated,clusters_coarse_2=="Myeloids")
a<-subset(combined_integrated,clusters_coarse_2=="Proliferating")
Myeloid<-merge(Myeloid,a)
Myeloid <- JoinLayers(Myeloid)

Myeloid[["RNA"]] <- split(Myeloid[["RNA"]], f = Myeloid$study)
Myeloid <- NormalizeData(Myeloid)
Myeloid<- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid)

Myeloid  <- IntegrateLayers(
  object = Myeloid , method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

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
getSSEfromExp(Myeloid,30)

SSE<-read.csv("SSE_results_Myeloids.csv")
plot1<-ggplot(data = SSE, mapping = aes(x = cluster, y = Variance)) + geom_line()
ggsave("plots/ElbowPlot_Myeloid.pdf", plot=plot1, height=4, width=4, dpi=600)

Myeloid <- JoinLayers(Myeloid)


Myeloid <- RunUMAP(Myeloid, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

Myeloid <- readRDS("Myeloids.rds")
Myeloid<- FindNeighbors(Myeloid, reduction = "integrated.cca", dims = 1:30)
Myeloid <- FindClusters(Myeloid, resolution = 0.4, cluster.name = "cca_clusters_0.4")


p1 <- DimPlot(
  Myeloid,
  reduction = "umap.cca",
  group.by = c("cca_clusters_0.4"),
  combine = T
)

Idents(object = Myeloid) <- "cca_clusters_0.4"

markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file="/SAN/colcc/NexTGen/WP6/scRNA/markers/markers_cca_Myeloids.csv")
top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,file="/SAN/colcc/NexTGen/WP6/scRNA/markers/top20_cca_Myeloids.csv")

p1<-VlnPlot(Myeloid, features=c("ARG1","MRC1","CD274","CX3CR1",),group.by = "cca_clusters_0.3",ncol=2)
##CCL2/3/4/5/20, CCL3L1, CCL3L3, CCL4L2, CCL4L4, CXCL1/2/3/5/8, G0S2, IL1B, IL1RN, IL6, INHBA, KLF2/6, NEDD9, PMAIP1, S100A8/A9, SPP1, inflam, "IL10","TGFB1" immunsuppressive
##"APOE","APOC1","FABP5",LA TAM,8,13
##ARG1、MRC1、CD274 和 CX3CR1 Reg Mac;
#"ISG15","CXCL10", IFN_TAM; 7
#"CLEC10A","CD1C",cDC2; 6
#"LILRA4","CLEC4C", pDC; 9
#""XCR1","CLEC9A", cDC1/CCR7+ DC; 11
#"TPSB2","KIT", mast cells; 12
#"CD14","S100A8","CDA",CD14 mono; 1
#"OCSTAMP","ACP5",osteoclasts; 4
#"VEGFA","SPP1","CEBPB","SLC2A1", angio TAM, 5,14 ;
#CCL18, CCL23, CD52, FABP4, FBP1, LGALS3, MACRO, MCEMP1, MRC1, MSR1, PPARG, RBP4,RTM-TAM, 2,3;
##ECM-remodeling COL3A1,COL5A2,SERPINH1


Myeloid <- RenameIdents(object = Myeloid,
"0" ="Inflam_TAMs",
"1" = "CD14_mono",
"2" = "RTM_TAMs",
"3" = "RTM_TAMs",
"4" = "osteoclasts",
"5" = "Angio_TAMs",
"6" = "cDC2",
"7" = "IFN_TAMs",
"8" = "LA_TAMs",
"9" = "pDC",
"10"= "Inflam_TAMs",
"11"="cDC1/CCR7+_DC",
"12"= "Mast_cells",
"13"= "LA_TAMs",
"14"= "Angio_TAMs"
)


Myeloid <- RenameIdents(object = Myeloid,
                        "0" = "Inflam_TAMs",
                        "1" = "LA_TAMs",
                        "2" = "Reg_TAMs",
                        "3" = "Angio_TAMs",
                        "4" = "CD14_mono",
                        "5" = "Osteoclasts",
                        "6" = "TAM/T_doublets",
                        "7" = "cDC2",
                        "8" = "Stressed_TAMs",
                        "9" = "Proliferating_TAMs",
                        "10"= "IFN_TAMs",
                        "11"= "Inflam_TAMs",
                        "12"= "pDC",
                        "13"= "ECM-remodeling_TAMs",
                        "14"= "cDC1/CCR7+_DC",
                        "15"= "Mast_cells",
                        "16"= "RTM_TAMs"
)

Myeloid$Myeloid_subsets<-Myeloid@active.ident

mypal =pal_d3("category20c")(20)

umap.df <- FetchData(Myeloid, vars = c("umapcca_1", "umapcca_2", "Myeloid_subsets"))
umap.cols <- DiscretePalette(length(unique(umap.df$Myeloid_subsets)), palette = "polychrome")
ggrepel.df <- umap.df %>% group_by(Myeloid_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean)
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(Myeloid_subsets = ""))

p<-ggplot(umap.df, aes(x = umapcca_1, y = umapcca_2)) +
geom_point(aes(fill = Myeloid_subsets), size = 1.5,  pch = 21) +
geom_label_repel(umap.df %>% group_by(Myeloid_subsets) %>% summarise_at(c("umapcca_1", "umapcca_2"), mean), size = 2, mapping = aes(x = umapcca_1, y = umapcca_2, label = Myeloid_subsets)) +guides(colour = guide_legend(override.aes = list(size=2))) + theme_classic()  + NoLegend()+theme(panel.border =  element_rect(colour = "black", fill = NA, size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(text = element_text(size = 10))+scale_fill_manual(values = mypal)

ggsave("plots/custom_UMAP_Myeloid.pdf", p, height=4, width=5, dpi=600)


##FOLR2和TREM2是LA-TAMs的核心marker
##CD209, CCL18, SLC40A1 Reg_TAMs
##SLC2A1,SPP1,VCAN, Angio_TAMs
##Figure1e bubble plot
Idents(Myeloid) <- factor(Idents(Myeloid), levels = c("IFN_TAMs","Inflam_TAMs","LA_TAMs","Reg_TAMs","Angio_TAMs","Stressed_TAMs","Proliferating_TAMs","ECM-remodeling_TAMs","RTM_TAMs","CD14_mono","Osteoclasts","cDC2","pDC","cDC1/CCR7+_DC","Mast_cells","TAM/T_doublets"))


markers.to.plot <- c("ISG15","CXCL9","CXCL10","CCL4L2","CCL3L1","CCL3","IL1B","APOE","FOLR2","TREM2","CD209","CCL18","SLC40A1","SLC2A1","SPP1","VCAN","HSPA6","DNAJA4","MKI67","STMN1","COL3A1","COL5A2","SERPINH1","MARCO","MRC1","FABP4","CD14","FCN1","S100A8","S100A9","OCSTAMP","DCSTAMP","CLEC10A","CD1C","LILRA4","CLEC4C","XCR1","CLEC9A","CCR7","LAMP3","TPSB2","MS4A2","KIT","CD3D","CD3E")

p<-DotPlot(Tcell, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
    RotatedAxis()
ggsave("plots/Bubble_heatmap_Myeloid.pdf", p, height=6, width=14, dpi=600)

cell_number<-data.frame(table(Myeloid$Myeloid_subsets))
cell_number <- cell_number[order(cell_number$Freq, decreasing = TRUE), ]
piepercent<- round(100*cell_number$Freq/sum(cell_number$Freq), 1)
labels<-paste0(cell_number$Var1," ",piepercent,"%")
pdf(file = "plots/pie_chart_cells_Myeloid.pdf")
pie(cell_number$Freq, labels=labels,col =mypal,cex=1)
dev.off()

proportion<-as.data.frame.array(prop.table(table(Myeloid$sample_ID,Myeloid$Myeloid_subsets),margin=1))
proportion$sample_ID<-rownames(proportion)
tumor_type<-Myeloid@meta.data[,c(5:7,9,13:15,19)]
tumor_type<-tumor_type %>% distinct(sample_ID, .keep_all = TRUE)
proportion<-merge(proportion,tumor_type,by="sample_ID")

plot<-dittoBarPlot(Myeloid, "Myeloid_subsets", group.by = "cancer_type",var.labels.reorder=c(6,7,8,13,1,15,12,5,14,2,10,4,11,3,9,16))
ggsave("plots/barplot_cancers_Myeloid.pdf", plot, height=5, width=6, dpi=600)
plot<-ggplot(data=proportion,aes_string(y="Inflam_TAMs",x="cancer_type",color="cancer_type"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_cancers_Myeloids_Inflam.pdf"), plot, height=5, width=12, dpi=600)

a<-subset(Myeloid,Tumour_site %in% c("primary","metastatic"))
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_site",var.labels.reorder=c(6,7,8,13,1,15,12,5,14,2,10,4,11,3,9,16))
ggsave("plots/barplot_Tumour_site_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,Tumour_site %in% c("primary","metastatic"))
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Tumour_site",color="Tumour_site"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_site_Myeloids_Inflam.pdf"), plot, height=5, width=6, dpi=600)

plot<-dittoBarPlot(Myeloid, "Myeloid_subsets", group.by = "Time_point",var.labels.reorder=c(6,7,8,13,1,15,12,5,14,2,10,4,11,3,9,16))
ggsave("plots/barplot_Time_point_Myeloid.pdf", plot, height=5, width=6, dpi=600)
plot<-ggplot(data=proportion  %>% filter(!is.na(Time_point)), aes_string(y="Inflam_TAMs",x="Time_point",color="Time_point"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Time_point_Myeloids_Inflam.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Myeloid, Treatment %in% c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Treatment",var.labels.reorder=c(6,7,8,13,1,14,12,5,2,10,4,11,3,9,15))+scale_x_discrete(limits = c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
ggsave("plots/barplot_Treatment_Myeloid.pdf", plot, height=6, width=6, dpi=600)
a<-subset(proportion,Treatment %in% c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))
a$Treatment <- factor(a$Treatment, levels = c("naive","Surgery","Chemotherapy","Surgery+Chemotherapy","Surgery+Radiotherapy"))  # 按照你想要的顺序替换各组名称
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Treatment",color="Treatment"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Treatment_Inflam_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)


a<-subset(Myeloid,cancer_type=="Ependymoma")
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(5,6,7,12,1,13,11,2,9,4,10,3,8,14))
ggsave("plots/barplot_Tumour_subgroup_ENPs_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Ependymoma")
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_ENPs_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)


a<-subset(Myeloid,cancer_type=="Hepatoblastoma")
plot<-dittoBarPlot(a,"Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(5,6,7,12,1,13,11,2,9,4,10,3,8,14))
ggsave("plots/barplot_Tumour_subgroup_HB_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Hepatoblastoma")
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Hepatoblastoma_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Myeloid,cancer_type=="Medulloblastoma")
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(5,6,7,11,1,12,10,2,4,9,3,8,13))
ggsave("plots/barplot_Tumour_subgroup_Medulloblastoma_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Medulloblastoma")
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Medulloblastoma_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Myeloid,cancer_type=="Neuroblastoma")
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(6,7,8,13,1,14,12,5,2,10,4,11,3,9,15))
ggsave("plots/barplot_Tumour_subgroup_NB_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion,cancer_type=="Neuroblastoma")
plot<-ggplot(data=a, aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Neuroblastoma_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Myeloid,cancer_type=="Rhabdomyosarcoma")
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(6,7,8,13,1,14,12,5,2,10,4,11,3,9,15))
ggsave("plots/barplot_Tumour_subgroup_RMS_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion ,cancer_type=="Rhabdomyosarcoma")
plot<-ggplot(data=a %>% filter(!is.na(Tumour_subgroup)), aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Rhabdomyosarcoma_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)

a<-subset(Myeloid,cancer_type=="Wilms_tumor")
plot<-dittoBarPlot(a, "Myeloid_subsets", group.by = "Tumour_subgroup",var.labels.reorder=c(5,6,7,11,1,12,10,2,4,9,3,8,13))
ggsave("plots/barplot_Tumour_subgroup_Wilms_tumor_Myeloid.pdf", plot, height=5, width=6, dpi=600)
a<-subset(proportion ,cancer_type=="Wilms_tumor")
plot<-ggplot(data=a %>% filter(!is.na(Tumour_subgroup)), aes_string(y="Inflam_TAMs",x="Tumour_subgroup",color="Tumour_subgroup"))+theme_classic()+geom_boxplot(outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ggtitle("")+xlab("")
ggsave(paste0("plots/boxplot_Tumour_subgroup_Wilms_tumor_TAMs_Inflam.pdf"), plot, height=5, width=6, dpi=600)




##stacked barplot


Idents(object = tumor) <- "cancer_type"

markers.to.plot <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-DMA","HLA-DMB")
p<-DotPlot(tumor, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
RotatedAxis()
ggsave("plots/HLA.pdf", p, height=6, width=10, dpi=600)
