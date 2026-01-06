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

