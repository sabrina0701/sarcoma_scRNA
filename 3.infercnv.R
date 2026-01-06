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



# 保存结果

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

# 提取数据
data <- FetchData(Wu, vars = c("seurat_clusters", "CNV_score"))

# 计算每个cWuster的中位数
cWuster_order <- data %>%
  group_by(seurat_cWusters) %>%
  summarise(med = median(CNV_score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(seurat_cWusters)

# 按中位数顺序重新设置因子
Wu$seurat_cWusters <- factor(Wu$seurat_cWusters, levels = cWuster_order)

# 重新画图
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

    # -------- 修改部分开始 --------
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
       # -------- 修改部分结束 --------
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

            # 给不同分组上色
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
    ## 第一阶段：先跑 CNA correlation
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
    ## 第二阶段：再跑 CNA signal
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
    ## 如果两个都成功，返回
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


# 提取数据
data <- FetchData(Wu, vars = c("unintegrated_clusters_SCT", "CNV_score"))

# 计算每个cluster的中位数
cluster_order <- data %>%
  group_by(unintegrated_clusters_SCT) %>%
  summarise(med = median(CNV_score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(unintegrated_clusters_SCT)

# 按中位数顺序重新设置因子
Wu$unintegrated_clusters_SCT <- factor(Wu$unintegrated_clusters_SCT, levels = cluster_order)

# 重新画图
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


##for DSRCT
p <- DotPlot(
DSRCT,
  features = list(
    "tumor_core" = c("WT1","FGFR4","IGF2","PDGFRB","FN1"), # Wu 融合/核心程序（WT1若掉测也看其余几个）
    "epi_mix" = c("KRT8","KRT18","KRT19","EPCAM")  ,# 肿瘤常见的上皮混合程序
    "Proliferation"  = c("MKI67","CCNB1","CDK1","TOP2A","UBE2C") ,# 增殖（肿瘤常高）
    "Endothelial"    = c("PECAM1","VWF","KDR","CLDN5","ESAM","PLVAP"),   # 内皮
    "Fibroblast"       = c("COL1A1","COL1A2","DCN","LUM","COL3A1","THBS2"),# 成纤维/CAF
    "Pericyte / SMC"         = c("ACTA2","TAGLN","MYH11","RGS5")       # 周细胞/平滑肌
  )
) +
  RotatedAxis() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text  = element_text(face = "bold")
  )
ggsave("/Volumes/analysis/infercnv/Wu/dotplot_Wu_states.pdf",
       plot = p, height = 7, width = 12, dpi = 600)

##for Wu
p <- DotPlot(
  Eps,
  features = list(
    "tumor_core_INI1_loss" = c("SMARCB1" ,"EZH2", "EED", "SUZ12","IGF2BP3" ),
    "tumor_epithelial" = c("KRT8","KRT18","KRT19","EPCAM","KRT7","MUC1" )  ,# 肿瘤常见的上皮混合程序
#  Classic-type Wu 的特征性支持（非所有病例都有，但出现即强支持）
    "tumor_classic_support" = c("EGFR",  "MET", "CD34"),       # 经典型 ES (~50–70%) 阳性,
    "Proliferation"  = c("MKI67","CCNB1","CDK1","TOP2A","UBE2C") ,# 增殖（肿瘤常高）
    "Endothelial"    = c("PECAM1","VWF","KDR","CLDN5","ESAM","PLVAP"),   # 内皮
    "Fibroblast"       = c("COL1A1","COL1A2","DCN","LUM","COL3A1","THBS2"),# 成纤维/CAF
    "Pericyte / SMC"         = c("ACTA2","TAGLN","MYH11","RGS5")       # 周细胞/平滑肌
  )
) +
  RotatedAxis() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text  = element_text(face = "bold")
  )

ggsave("/Volumes/analysis/infercnv/Wu/dotplot_Wu_states.pdf",
       plot = p, height = 7, width = 12, dpi = 600)

##Ewing
p <- DotPlot(
  Wu,
  features = list(
    "tumor_core" = c( "NR0B1",        # = DAX1，EWS-FLI1直接靶点，EwS高
                      "NKX2-2",       # 经典EwS转录因子靶基因（诊断常用）
                      "STEAP1",       # EwS高表达的膜蛋白，受EWS-FLI1调控
                      "GLI1"    ),
    "IGF axis / receptor" = c("IGF1", "IGF1R"), # EwS依赖IGF/IGF1R轴，常上调
    "Surface / diagnostic aid" = c( "CD99",         # = MIC2；蛋白层面最经典，但RNA可变
                                    "NCAM1" ),        # 可变（辅助）),       # 经典型 ES (~50–70%) 阳性,
    "Proliferation"  = c("MKI67","CCNB1","CDK1","TOP2A","UBE2C") ,# 增殖（肿瘤常高）
    "Endothelial"    = c("PECAM1","VWF","KDR","CLDN5","ESAM","PLVAP"),   # 内皮
    "Fibroblast"       = c("COL1A1","COL1A2","DCN","LUM","COL3A1","THBS2"),# 成纤维/CAF
    "Pericyte / SMC"         = c("ACTA2","TAGLN","MYH11","RGS5"),       # 周细胞/平滑肌
    "Osteoblast_lineage" = c("ALPL","RUNX2","BGLAP"),
    "Schwann/Glial" = c("SOX10","MPZ","PLP1"),
    "Neuron_like" = c("TUBB3","STMN2","DCX","SNAP25")
  )
) +
  RotatedAxis() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text  = element_text(face = "bold")
  )

ggsave("/Volumes/analysis/infercnv/Wu/dotplot_states.pdf",
       plot = p, height = 7, width = 12, dpi = 600)


# 自动扫描目录下的所有 run.final.infercnv_obj.rds 或 run.final.infercnv_obj
base_dir <- "/Volumes/analysis/infercnv/Wu"
samples <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
path  <- file.path(base_dir, samples, "run.final.infercnv_obj")

# 读取 & 提取矩阵 ---------------------------------------------------
extract_expr_matrix <- function(path, sample_id) {
  message("Reading ", sample_id, " ...")
  obj <- readRDS(path)
  
  mat <- as.matrix(obj@expr.data)
  
  # 返回矩阵与gene_order
  go <- obj@gene_order
  if (!"gene" %in% colnames(go)) go$gene <- rownames(go)
  list(expr = mat, gene_order = go)
  
  mat <- mat[rownames(mat) %in% go$gene, , drop = FALSE]
  go  <- go[match(rownames(mat), go$gene), ]
  go$chr <- gsub("^chr", "", go$chr, ignore.case = TRUE)
  go <- go[!go$chr %in% c("M","MT","m","mt"), ]
  mat <- mat[go$gene, , drop = FALSE]

  list(expr = mat, gene_order = go[, c("gene","chr","start","stop")])
}

lst <- imap(path, extract_expr_matrix)

# 统一基因集合与顺序 ----------------------------------------------
common_genes <- Reduce(intersect, lapply(lst, function(x) rownames(x$expr)))
go0 <- lst[[1]]$gene_order
go0 <- go0[go0$gene %in% common_genes, ]
go0 <- go0[order(factor(go0$chr, levels = c(as.character(1:22),"X","Y")),
                 as.numeric(go0$start)), ]
genes_ordered <- go0$gene

# 合并所有样本矩阵
expr_list <- lapply(lst, function(x) x$expr[genes_ordered, , drop = FALSE])
gene_cell_all <- do.call(cbind, expr_list)
rownames(gene_cell_all) <- genes_ordered
message("Merged matrix: ", nrow(gene_cell_all), " genes × ", ncol(gene_cell_all), " cells")

# 细胞对齐到 Seurat 对象 ------------------------------------------
common_cells <- intersect(colnames(Wu), colnames(gene_cell_all))
gene_cell_all <- gene_cell_all[, common_cells, drop = FALSE]
clusters <- Idents(Wu)[colnames(gene_cell_all)]

# 每群下采样最多500个细胞 ----------------------------------------
set.seed(42)
by_cl <- split(names(clusters), clusters)
sel_cells <- unlist(lapply(by_cl, function(v) {
  k <- min(500, length(v))
  if (length(v) > k) sample(v, k) else v
}), use.names = FALSE)

gene_cell_ds <- gene_cell_all[, sel_cells, drop = FALSE]
clusters_ds <- clusters[sel_cells]

# =========================
# 5) 构造 inferCNV 主图方向的矩阵：
#    行=细胞（按cluster排序与分块），列=基因（按chr→start固定）
# =========================
# 列 = 基因：基因顺序已经按 go0 固定
# 行 = 细胞：按 cluster 排序并切块
ord_rows <- order(clusters_ds)
mat <- t(gene_cell_ds)[ord_rows, , drop = FALSE]   # 行=细胞，列=基因
cl_sorted <- clusters_ds[ord_rows]

# =========================
# 6) inferCNV风格配色：中心 x.center，1%~99% 对称裁剪
#    若你的矩阵是比例刻度（二倍体≈1），用 x.center=1；如是log刻度，用0
# =========================
x.center <- 1
vals <- as.vector(mat); vals <- vals[is.finite(vals)]
vals_neq <- vals[vals != x.center]
if (length(vals_neq) < 10) {
  delta <- max(0.05, 0.5 * sd(vals, na.rm = TRUE))
} else {
  qs <- quantile(vals_neq, probs = c(0.01, 0.99), na.rm = TRUE)
  delta <- max(abs(x.center - qs[1]), abs(qs[2] - x.center))
}
low_th  <- x.center - delta
high_th <- x.center + delta

mat_clip <- mat
mat_clip[mat_clip < low_th]  <- low_th
mat_clip[mat_clip > high_th] <- high_th

col_fun <- colorRamp2(c(low_th, x.center, high_th), c("darkblue","white","darkred"))

# =========================
# 7) 画图（推荐用 ComplexHeatmap 的栅格输出，文件小）
#    行按 cluster 分块显示（相当于 inferCNV 的分组）
# =========================
row_split <- factor(cl_sorted, levels = unique(cl_sorted))
col_split <- factor(go0$chr[match(colnames(mat_clip), go0$gene)],
                    levels = c(as.character(1:22), "X", "Y"))
# 2) 列分块：按染色体（确保顺序为 1..22, X, Y）
outfile <- file.path(base_dir, "inferCNV_merged_downsampled.png")
ragg::agg_png(outfile, width = 2000, height = 3000, res = 220)
ht <- Heatmap(
  mat_clip,
  name = "CNV",
  col = col_fun,
  cluster_rows = FALSE,            # 行不聚类（我们已按cluster排好）
  cluster_columns = FALSE,         # 列不聚类（我们已按chr→start排好）
  row_split = row_split,           # 行按cluster分块
  column_split = col_split,        # 列按染色体分块 ✅ 关键修正
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,               # 栅格化，文件小
  raster_device = "png",
  raster_quality = 2,
  border = NA,
  # 分块间距（可调）
  row_gap = unit(1, "mm"),
  column_gap = unit(1, "mm"),
  # 染色体分块标题样式（默认会显示 1,2,...,X,Y）
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  # 可选：行分块标题（cluster编号）
  row_title_gp = gpar(fontsize = 9, fontface = "bold")
)
draw(ht)
dev.off()

##correlation
# 1) 为每个 cluster 先得到基因层面的平均（沿用上面 cluster_profiles）
## 先生成每个 cluster 的基因层面平均型谱：genes × clusters
clusters_lvls <- levels(row_split)  # 与主图一致的cluster顺序
cluster_profiles <- sapply(clusters_lvls, function(cl) {
  colMeans(mat_clip[cl_sorted == cl, , drop = FALSE], na.rm = TRUE)  # 长度=基因数
})
# 补上行名=基因，列名=cluster
rownames(cluster_profiles) <- colnames(mat_clip)
colnames(cluster_profiles) <- clusters_lvls


# 2) 构造每条染色体的聚合（中位数或均值都行）
chr_levels <- c(as.character(1:22), "X", "Y")
gene_chr  <- go0$chr[match(rownames(cluster_profiles), go0$gene)]
gene_chr  <- factor(gene_chr, levels = chr_levels)

# 聚合函数：按染色体对行聚合（对每个 cluster 列取每条 chr 的中位数）
agg_by_chr <- apply(cluster_profiles, 2, function(v) {
  tapply(v, gene_chr, median, na.rm = TRUE)
})
# 结果维度：染色体 × cluster

# 3) 染色体维度上的 cluster 相关性
cor_mat_chr <- cor(agg_by_chr, method = "spearman", use = "pairwise.complete.obs")

# 4) 画热图
ht_cor_chr <- ComplexHeatmap::Heatmap(
  cor_mat_chr,
  name = "r (chromosomes)",
  col = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  rect_gp = grid::gpar(col = NA)
)

outfile_cor2 <- file.path(base_dir, "cluster_CNV_correlation_byChr.png")
ragg::agg_png(outfile_cor2, width = 1600, height = 1400, res = 220)
draw(ht_cor_chr)
dev.off()


rho <- cor(cluster_profiles, method = "spearman", use = "pairwise.complete.obs")

# 4) 画热图
ht_cor_chr <- ComplexHeatmap::Heatmap(
  rho,
  name = "r (chromosomes)",
  col = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  rect_gp = grid::gpar(col = NA)
)

outfile_cor2 <- file.path(base_dir, "cluster_CNV_correlation.png")
ragg::agg_png(outfile_cor2, width = 1600, height = 1400, res = 220)
draw(ht_cor_chr)
dev.off()


cells<-read.csv("/Volumes/analysis/infercnv/tumor_cell.csv")
tumor<-subset(combined,cells = cells$cell_name)
saveRDS(tumor,"/Volumes/analysis/infercnv/tumor.rds")
