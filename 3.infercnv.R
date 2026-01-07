########### R script for infercnv. ###########
########### Created by Danwen Qian on 27-10-2025. ###########
########### Last modified by Danwen Qian on 27-10-2025. ###########

library(Seurat)
library(infercnv)
library(stringr)
library(infercna)
library(dplyr)
library(ggplot2)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(ragg)
library(RColorBrewer)

options("Seurat.object.assay.version" = "v3")

combined <- readRDS("/Volumes/document/combined_filter.rds")
combined[["RNA"]] <- as(combined[["RNA"]], Class = "Assay")

combined  <- subset(combined , layer_1 != "Blood_Cell")


refs <- read.csv("/Volumes/document/ref.csv")
refs_cells <- subset(combined, cells = refs$cell_names)
refs_cells <- refs_cells@assays$RNA@counts

Wu <- subset(combined, study == "Wu_osteosarcoma")
Wu <- subset(Wu, cells = refs$cell_names, invert = TRUE)

refs <- split(refs$cell_names, f = refs$cell_type)
Wu <- SplitObject(Wu, split.by = "sample_ID")

ann <- read.table("infercnv_test/Wu_ann_sarcoma.txt", header = F, row.names = 1)

gene_pos <- read.table("/Volumes/document/gene_pos.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_pos) <- c("gene", "chr", "start", "end")

common_genes <- intersect(rownames(combined), gene_pos$gene)


results <- data.frame()
scores <- data.frame()
sample<-c()


for (i in 1:length(Wu)) {
    sparse_matrix <- Wu[[i]]@assays$RNA@counts
    sparse_matrix <- sparse_matrix[common_genes, ]
    refs_cells <- refs_cells[common_genes, ]

    sparse_matrix <- cbind(sparse_matrix, refs_cells)
    ann_this <- subset(ann, rownames(ann) %in% colnames(sparse_matrix))
    
    infercnv_obj = CreateInfercnvObject(
        raw_counts_matrix = sparse_matrix,
        annotations_file = ann_this,
        delim = "\t",
        gene_order_file = "/Volumes/document/gene_pos.txt",
        ref_group_names = c("Fibroblast", "Endothelial")
    )

    infercnv_obj = infercnv::run(
        infercnv_obj,
        cutoff = 0.1,
        out_dir = paste0("/Volumes/analysis/infercnv/Wu/", names(Wu[i])),
        cluster_by_groups = FALSE,
        denoise=TRUE
    )

    # 提取CNV值矩阵
    expr_adj <- log2(infercnv_obj@expr.data)  # gene x cell 矩阵
    
    png(paste0("/Volumes/analysis/infercnv/Wu/",names(Wu[i]),"/findmaligmant.png"), width = 1000, height = 800)
    Modes = findMalignant_chr(expr_adj,refCells = refs,plot=T,sample=NULL)
    dev.off()
    
    if(!isFALSE(Modes)){
        Modes<-data.frame(unlist(Modes))
        Modes$class<-str_match(rownames(Modes),"(malignant|nonmalignant|unassigned)")[,2]
        results<-rbind(results,Modes)
    } else {
        sample<-c(sample,i)
        }

    # 计算CNV score
    cnv_scores <- apply(expr_adj, 2, function(x) mean(x^2, na.rm = TRUE))

    # 记录到data frame
    scores <- rbind(
       scores,
        data.frame(
            cell = names(cnv_scores),
            sample = names(Wu[i]),
            CNV_score = cnv_scores
        )
    )
}




Wu<-subset(combined, study=="Wu_osteosarcoma")

scores <- scores[!duplicated(scores$cell), ]
scores<-subset(scores,cell %in% colnames(Wu))
write.csv(scores, "/Volumes/analysis/infercnv/Wu/CNV_scores_all_cells.csv", row.names = FALSE)

colnames(results)[1]<-"cell"
results <- results[!duplicated(results$cell), ]
results <-subset(results,cell %in% colnames(Wu))
write.csv(results , "/Volumes/analysis/infercnv/Wu/CNV_tumor_cells.csv", row.names = FALSE)

Wu <- AddMetaData(Wu, metadata = scores)

Wu[["RNA"]] <- split(Wu[["RNA"]], f = Wu$sample_ID)
Wu <- NormalizeData(Wu)
Wu <- FindVariableFeatures(Wu)
Wu <- ScaleData(Wu)
Wu <- RunPCA(Wu)

Wu <- FindNeighbors(Wu, dims = 1:30, reduction = "pca")
Wu <- FindCWusters(Wu, resoWution = 0.2, cWuster.name = "unintegrated_cWusters")

Wu <- IntegrateLayers(
  object = Wu,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

Wu <- IntegrateLayers(
  object = Wu, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

Wu <- FindNeighbors(Wu, reduction = "integrated.cca", dims = 1:30)
Wu <- FindCWusters(Wu, resoWution = 0.2, cWuster.name = "cca_cWusters")

Wu <- FindNeighbors(Wu, reduction = "harmony", dims = 1:30)
Wu <- FindCWusters(Wu, resoWution = 0.2, cWuster.name = "harmony_cWusters")

Wu[["RNA"]] <- split(Wu[["RNA"]], f = Wu$sample_ID)
Wu <- SCTransform(Wu, vst.flavor = "v2")
Wu <- RunPCA(Wu)

Wu <- FindNeighbors(Wu, dims = 1:30, reduction = "pca")
Wu <- FindCWusters(Wu, resoWution = 0.2, cWuster.name = "unintegrated_cWusters_SCT")


Wu <- IntegrateLayers(
  object = Wu, method = CCAIntegration,
  new.reduction = "integrated.cca", normalization.method = "SCT",
  verbose = FALSE
)

Wu <- IntegrateLayers(
  object = Wu, method = HarmonyIntegration,
  new.reduction = "harmony", normalization.method = "SCT",
  verbose = FALSE
)

library(dplyr)

data <- FetchData(Wu, vars = c("seurat_clusters", "CNV_score"))

cWuster_order <- data %>%
  group_by(seurat_cWusters) %>%
  summarise(med = median(CNV_score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(seurat_cWusters)

Wu$seurat_cWusters <- factor(Wu$seurat_cWusters, levels = cWuster_order)

VlnPlot(Wu, features = "CNV_score", group.by = "seurat_cWusters", pt.size = 0)


library(VennDiagram)
library(gridExtra)
library(grid)
venn_dat <- read.csv("/VoWumes/analysis/infercnv/Wu/Wu_tumor.csv")


venn_list <- list(
original = paste0(venn_dat$cells_names[!is.na(venn_dat$original)], "_", venn_dat$original[!is.na(venn_dat$original)]),
inferCNV = paste0(venn_dat$cells_names[!is.na(venn_dat$inferCNV)], "_", venn_dat$inferCNV[!is.na(venn_dat$inferCNV)]),
scATOMIC = paste0(venn_dat$cells_names[!is.na(venn_dat$scATOMIC)], "_", venn_dat$scATOMIC[!is.na(venn_dat$scATOMIC)]),
copyKAT = paste0(venn_dat$cells_names[!is.na(venn_dat$copyKAT)], "_", venn_dat$copyKAT[!is.na(venn_dat$copyKAT)])
)

venn_list <- list(
original = paste0(venn_dat$cells_names[!is.na(venn_dat$original)], "_", venn_dat$original[!is.na(venn_dat$original)]),
CNV_cWustering = paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering)], "_", venn_dat$CNV_cWustering[!is.na(venn_dat$CNV_cWustering)]),
CNV_cWustering_cca = paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering_cca)], "_", venn_dat$CNV_cWustering_cca[!is.na(venn_dat$CNV_cWustering_cca)]),
CNV_cWustering_harmony = paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering_harmony)], "_", venn_dat$CNV_cWustering_harmony[!is.na(venn_dat$CNV_cWustering_harmony)])
)


venn_list <- list(
original = paste0(venn_dat$cells_names[!is.na(venn_dat$original)], "_", venn_dat$original[!is.na(venn_dat$original)]),
CNV_cWustering_SCT = paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering_SCT)], "_", venn_dat$CNV_cWustering_SCT[!is.na(venn_dat$CNV_cWustering_SCT)]),
CNV_cWustering_SCT_CCA= paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering_SCT_CCA)], "_", venn_dat$CNV_cWustering_SCT_CCA[!is.na(venn_dat$CNV_cWustering_SCT_CCA)]),
CNV_cWustering_SCT_harmony = paste0(venn_dat$cells_names[!is.na(venn_dat$CNV_cWustering_SCT_harmony)], "_", venn_dat$CNV_cWustering_SCT_harmony[!is.na(venn_dat$CNV_cWustering_SCT_harmony)])
)

v1<-venn.diagram(venn_list,filename= NULL, alpha=0.5,euler.d=T,scaled=T,
                 fill=c("#FFFFCC","#CCFFFF","#FFCCCC","#b3cde3"),col="transparent")
ggsave("/VoWumes/analysis/infercnv/Wu/venn_plot.pdf", plot = grid.arrange(grobTree(v1)), width = 6, height = 6, dpi = 300)


venn_dat <- venn_dat %>%
    filter(!(inferCNV == "" | inferCNV == "NA" | is.na(inferCNV)))

findMalignant_chr <- function(
        cna,
        refCells = NULL,
        samples = scalop::unique_sample_names(colnames(cna), max.nchar = 6),
        gene.quantile = 0.9,
        gene.quantile.for.corr = 0.5,
        gene.quantile.for.signal = gene.quantile,
        use.bootstraps = TRUE,
        n.bootstraps = 10000,
        prob = 0.95,
        coverage = 0.8,
        verbose = TRUE,
        plot = TRUE,
        border.col = "black",
        alpha = 0.3,
        cex = 0.8,
        pch = 20,
        hline = NULL,
        vline = NULL,
        groups.col = NULL,
        ...
) {
    if (verbose)
           message("Calculating cells' CNA correlations...")
       cors = cnaCor(cna, gene.quantile = gene.quantile.for.corr,
           samples = samples, refCells = refCells)
       if (verbose)
           message("Calculating cells' CNA signal...")
       signals = cnaSignal(cna, gene.quantile = gene.quantile.for.signal)
    
    q99 <- quantile(signals, 0.99, na.rm = TRUE)
    outliers <- names(signals)[signals > 2 * q99]

     if (length(outliers) > 0) {
         if (verbose) message(
             "Detected ", length(outliers),
             " extreme outlier cells (CNA signal > 2 x 99% quantile = ", round(2 * q99, 4), "). Removing them."
         )
         signals <- signals[!names(signals) %in% outliers]
         cors <- cors[names(cors) %in% names(signals)]
     } else {
         if (verbose) message("No extreme CNA signal outliers detected.")
     }
    
    plot(cors, signals,
         xlab = "CNA Correlation",
         ylab = "CNA Signal",
         pch = 1,
         cex = 0.3)
    
    res <- fit_cor_and_sig(cors, signals, use.bootstraps, n.bootstraps, prob, coverage, ...)

       if (isFALSE(res)) {
               if (verbose) message("GMM failed. Falling back to K-means cWustering...")
               set.seed(123)
               km_res <- kmeans(cbind(cors, signals), centers = 2, nstart = 50)
               cWuster_km <- km_res$cWuster
               res <- list(
                   corGroups = list(
                       Low = names(cors)[cWuster_km == which.min(tapply(cors, cWuster_km, mean))],
                       High = names(cors)[cWuster_km == which.max(tapply(cors, cWuster_km, mean))]
                   ),
                   sigGroups = list(
                       Low = names(signals)[cWuster_km == which.min(tapply(signals, cWuster_km, mean))],
                       High = names(signals)[cWuster_km == which.max(tapply(signals, cWuster_km, mean))]
                   )
               )
               if (verbose) message("K-means classification successful.")
       }
    corGroups <- res$corGroups
    sigGroups <- res$sigGroups

    
    if (verbose) {
        message("Identifying cells classified as CNA-high by both parameters...")
    }
    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    unassigned = union(unlist(sect[2, 1]), unlist(sect[1, 2]))
    len.unassig = length(unassigned)
    corGroups = sapply(corGroups, function(gr) gr[!gr %in% unassigned], simplify = F)
    sigGroups = sapply(sigGroups, function(gr) gr[!gr %in% unassigned], simplify = F)
    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    
    result = list(
        malignant = unique(c(unlist(sect[2, 2], use.names = F), outliers )),
        nonmalignant = unlist(sect[1, 1], use.names = F),
        unassigned = unassigned
    )
    if (length(unassigned) >= 1) {
        message(round(100 * len.unassig/ncol(cna), 2),
                " % of cells were assigned to opposing modes,",
                "\ni.e.to high CNA-signal and low CNA-correlation or vice versa.",
                "\nThese cells will remain unasssigned.")
    }
    if (plot) {
            if (verbose) message("Plotting CNA correlation against CNA signal...")

            # 自动配色
            groups.col = scalop::discrete_colours[1:length(result)]
            groups.col <- scales::alpha(groups.col, alpha)

            plot(cors, signals,
                 xlab = "CNA Correlation",
                 ylab = "CNA Signal",
                 pch = 1,
                 col = border.col,
                 cex = cex)

            points(cors[result$nonmalignant], signals[result$nonmalignant],
                   col = groups.col[1], pch = pch, cex = cex)
            points(cors[result$malignant], signals[result$malignant],
                   col = groups.col[2], pch = pch, cex = cex)
            if (length(result$unassigned) > 0) {
                points(cors[result$unassigned], signals[result$unassigned],
                       col = groups.col[3], pch = pch, cex = cex)
            }

            if (!is.null(vline)) abline(v = vline, lty = 2)
            if (!is.null(hline)) abline(h = hline, lty = 2)
        }
    result
}


fit_cor_and_sig <- function(cors, signals, use.bootstraps, n.bootstraps, prob, coverage, ...,
                            min.bootstraps = 8000, step = 1000, verbose = TRUE) {
    ## ------------------------
    ##  CNA correlation
    ## ------------------------
    corGroups <- FALSE
    current_bootstraps <- n.bootstraps
    while (current_bootstraps >= min.bootstraps) {
        if (verbose) message("Trying fitBimodal on CNA correlations with n.bootstraps = ", current_bootstraps, " ...")
        
        old <- Sys.time()
        invisible(capture.output(
            corGroups <- suppressMessages(
                fitBimodal(cors,
                           bySampling = use.bootstraps, nsamp = current_bootstraps,
                           prob = prob, coverage = coverage, assign = TRUE, ...)
            )
        ))
        new <- Sys.time()
        
        if (!isFALSE(corGroups)) {
            timerun = scalop::hms_span(start = old, end = new)
            if (verbose) message("CNA correlations fit succeeded (completed in ", timerun, "s)")
            break
        }
        
        current_bootstraps <- current_bootstraps - step
    }
    if (isFALSE(corGroups)) {
        if (verbose) message("fitBimodal failed for CNA correlations after trying down to ", min.bootstraps, " bootstraps.")
        return(FALSE)
    }
    
    ## ------------------------
    ## CNA signal
    ## ------------------------
    sigGroups <- FALSE
    current_bootstraps <- n.bootstraps
    while (current_bootstraps >= min.bootstraps) {
        if (verbose) message("Trying fitBimodal on CNA signals with n.bootstraps = ", current_bootstraps, " ...")
        
        old <- Sys.time()
        invisible(capture.output(
            sigGroups <- suppressMessages(
                fitBimodal(signals,
                           bySampling = use.bootstraps, nsamp = current_bootstraps,
                           prob = prob, coverage = coverage, assign = TRUE, ...)
            )
        ))
        new <- Sys.time()
        
        if (!isFALSE(sigGroups)) {
            timerun = scalop::hms_span(start = old, end = new)
            if (verbose) message("CNA signals fit succeeded (completed in ", timerun, "s)")
            break
        }
        
        current_bootstraps <- current_bootstraps - step
    }
    if (isFALSE(sigGroups)) {
        if (verbose) message("fitBimodal failed for CNA signals after trying down to ", min.bootstraps, " bootstraps.")
        return(FALSE)
    }
    
    ## ------------------------
    ## 
    ## ------------------------
    if (verbose) message("fitBimodal succeeded for both CNA correlations and CNA signals.")
    return(list(corGroups = corGroups, sigGroups = sigGroups))
}


##use the cnvScore+SCT clustering method
combined <- readRDS("/Volumes/document/combined_filter.rds")
combined  <- subset(combined , layer_1 != "Blood_Cell")

Wu<-subset(combined, study=="Wu_MPNST")

scores <- read.csv("/Volumes/analysis/infercnv/Wu/CNV_scores_all_cells.csv",row.names =1)

Wu <- AddMetaData(Wu, metadata = scores)

Wu[["RNA"]] <- split(Wu[["RNA"]], f = Wu$sample_ID)
Wu <- SCTransform(Wu, vst.flavor = "v2")

Wu <- RunPCA(Wu)

Wu <- FindNeighbors(Wu, dims = 1:30, reduction = "pca")
Wu <- FindClusters(Wu, resolution = 0.2, cluster.name = "unintegrated_clusters_SCT")
Wu <- RunUMAP(Wu, dims = 1:30, reduction = "pca")


data <- FetchData(Wu, vars = c("unintegrated_clusters_SCT", "CNV_score"))

cluster_order <- data %>%
  group_by(unintegrated_clusters_SCT) %>%
  summarise(med = median(CNV_score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(unintegrated_clusters_SCT)

Wu$unintegrated_clusters_SCT <- factor(Wu$unintegrated_clusters_SCT, levels = cluster_order)


    p<-VlnPlot(Wu, features = c("CNV_score","COL1A1", "DCN","PECAM1","CLDN5","KDR","VWF","SOX10", "UBE2C", "EGFR","MPZ","DUSP6"), group.by = "unintegrated_clusters_SCT", pt.size = 0)
ggsave("/Volumes/analysis/infercnv/Wu/VlnPlot.pdf", plot=p, height=10, width=15, dpi=600)

p <- DimPlot(
    Wu,
    reduction = "umap",
    label = TRUE,
)
ggsave("/Volumes/analysis/infercnv/Wu/UMAP.pdf", plot=p, height=6, width=6, dpi=600)

##if can't distinguish, could run differential expression analysis
DefaultAssay(Wu)<-"RNA"
Wu <- JoinLayers(Wu)

Wu <- NormalizeData(Wu)

markers <- FindAllMarkers(Wu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,file="/Volumes/analysis/infercnv/Wu/top20.csv")


Wu$tumor_class<-"tumor"
Wu$tumor_class[Wu$unintegrated_clusters_SCT %in% c(3,10,5,9,7,12)]<-"non_tumor"

write.csv(Wu@meta.data,"/Volumes/analysis/infercnv/Wu/meta.csv")


##fibroblast: "COL1A1", "DCN"
##endothelial:"PECAM1", "VWF"
##CSS: "ATF1", "CREB1", "TFAP2A","MITF", "SOX10",  "PMEL", "MLANA"
##DSRCT: "WT1","FGFR4","PDGFRB","FN1","KRT8","KRT18","KRT19" but I need further confirm about INF clusters
##eps: SMARCB1 loss
##Ewing:"NR0B1","NKX2-2","STEAP1"
##WD/DDLPS: "MDM2", "CDK4", "HMGA2", "TSPAN31",12gain
##Myxoid Liposarcoma: "DLK1","IGF2","MEG3","TNC","CXCL14","COL18A1","GDF10"
##synovial_sarcoma:"TLE1","EPCAM","KRT8","BCL2"
###ATRT: SMARCB1 loss,chr22q loss,"EZH2","SOX4","SOX11"
##OS: COL1A1, SPP1, IBSP, RUNX2,1q↑ / 6p↑ / 8q↑ / 10q↓ / 13q↓
##RMS: MYOD1 + MYF5 + MYOG, ARMS chr2 gain + PAX3/FOXO1 区域改变; ERMS chr8 gain + chr11p15 loss + chr7 gain
##UPS:"UBE2C","TOP2A","CCNB1","CDK1","AURKB","CENPF","MKI67"
##MPNST: SOX10, 低 +（UBE2C 或 TOP2A）高 +（EGFR 或 DUSP6）高,17q11.2 loss 或 9p21 loss
