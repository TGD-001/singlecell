library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(DropletUtils)
library(tidyr)
library(clustree)
library(presto)
library(Matrix)
library("loupeR")
#参数设置
dim_range <- 1:20
resolution_arg<- 0.2
percentmt<- 5
features_numbers <- 2000
#调色盘
palette1 <- brewer.pal(9, "Set1")
palette2 <- brewer.pal(8, "Set2")
palette3 <- brewer.pal(9, "Set3")
custom_colors <- c(palette1, palette2, palette3)

#输出文件夹
outfolder <- "/histor/wangxh/tangguodong/singlecell/data/GS_combin/test_GS/G/"

#qc过滤及去双胞流程
config <- read.table("config.csv", sep = "\t", header = FALSE)
sample <- config[, 1]
sample_list <- list()
sample_matrix_list <- list()
for (i in 1:length(sample)) {
  # 修正路径拼接：用paste0正确连接所有部分
  data.file <- paste0(getwd(), "/cellbender/", sample[i], "/output_filtered_seurat.h5")
  # 读取.h5文件（确保Seurat包已加载）
  counts <- Read10X_h5(filename = data.file, use.names = TRUE)
  # 创建Seurat对象并预处理
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample[i], min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  pdf(paste0(outfolder, sample[i], "proQC_vinplot.pdf"), width = 20, height = 10)
  p<-VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  dev.off()
  
  print(paste0(sample[i],"现有：", dim(seurat_obj)[1], "基因 x", dim(seurat_obj)[2], "细胞"))
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < percentmt)
  print(paste0(sample[i],"质控后有：", dim(seurat_obj)[1], "基因 x", dim(seurat_obj)[2], "细胞"))

  pdf(paste0(outfolder, sample[i], "postQC_vinplot.pdf"), width = 20, height = 10)
  p<-VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)  
  dev.off()
   
  print(head(seurat_obj@meta.data))
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = features_numbers)
  
  pdf(paste0(outfolder, sample[i], "VariableFeatures.pdf"), width = 20, height = 10) 
  top10 <- head(VariableFeatures(seurat_obj), 10)
  plot1 <- VariableFeaturePlot(seurat_obj)  
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  dev.off()

  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
   
  pdf(paste0(outfolder, sample[i], "ElbowPlot.pdf"), width = 20, height = 10)
  p<-ElbowPlot(seurat_obj, ndims = 30)
  print(p)
  dev.off()

  seurat_obj <- FindNeighbors(seurat_obj, dims = dim_range)
  
  #聚类树
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) {
  seurat_obj <- FindClusters(seurat_obj, graph.name = "RNA_snn", resolution = res, algorithm = 1)
  }
  apply(seurat_obj@meta.data[,grep("RNA_snn_res",colnames(seurat_obj@meta.data))],2,table)
  p_tree <- clustree(seurat_obj@meta.data,  prefix = "RNA_snn_res.")
  pdf(paste0(outfolder,sample[i], "resolution.pdf"), width = 12, height = 10)
  print(p_tree)
  dev.off()

  seurat_obj <- FindClusters(seurat_obj, resolution = resolution_arg)
  seurat_obj <- RunUMAP(seurat_obj, dims = dim_range)
  seurat_obj <- RunTSNE(seurat_obj, dims = dim_range)
    
  pdf(paste0(outfolder, sample[i], "umap_tsne.pdf"), width = 20, height = 10)
  p1<-DimPlot(seurat_obj, reduction = "umap" ,label = TRUE)
  p2<-DimPlot(seurat_obj, reduction = "tsne",label = TRUE)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
  #DoubletFinder去除双胞
  sweep.res.list_seurat_obj <- paramSweep(seurat_obj, PCs = dim_range, sct = FALSE)
  sweep.stats_seurat_obj <- summarizeSweep(sweep.res.list_seurat_obj, GT = FALSE)
  bcmvn_seurat_obj <- find.pK(sweep.stats_seurat_obj)
  mpK<-as.numeric(as.vector(bcmvn_seurat_obj$pK[which.max(bcmvn_seurat_obj$BCmetric)]))
  homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)
  doublerate <- ncol(seurat_obj)*4e-06
  print(paste("细胞总数为:",ncol(seurat_obj)))
  print(paste("估计双胞率为:",doublerate))
  nExp_poi <- round(doublerate *nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seurat_obj <- doubletFinder(seurat_obj, PCs = dim_range, pN = 0.25, pK = mpK,nExp = nExp_poi.adj)
  df_col <- grep("^DF.classifications_", colnames(seurat_obj@meta.data), value = TRUE)
  
  doublet_table <- table(seurat_obj[[df_col]])
  print("双胞检测结果:")
  print(doublet_table)
  doublet_rate <- doublet_table["Doublet"] / sum(doublet_table)
  print(paste("双胞比例:", round(doublet_rate * 100, 2), "%"))

  pdf(paste0(outfolder, sample[i], "doublefinder.pdf"), width = 20, height = 15)
  p<-DimPlot(seurat_obj, reduction = "umap", group.by = df_col,
        cols = c("Singlet" = "blue", "Doublet" = "red"))
  print(p)
  dev.off()
  
  seurat_obj_filtered <- seurat_obj[, seurat_obj@meta.data[, df_col] == "Singlet"]
  
  #保存单个对象
  mtx <- as.matrix(seurat_obj_filtered@assays$RNA@layers$counts)
  rownames(mtx) <- rownames(seurat_obj_filtered@assays$RNA@features)
  colnames(mtx) <- rownames(seurat_obj_filtered@assays$RNA@cells)
  mtx <- Matrix(mtx, sparse = TRUE)
  DropletUtils:::write10xCounts(paste0(outfolder, sample[i], "singlet_matrix"), mtx, version="3")
  
  singlet_seurat <- CreateSeuratObject(counts = mtx,project = sample[i])
  singlet_seurat[["percent.mt"]] <- PercentageFeatureSet(singlet_seurat, pattern = "^MT-") 

  #保存对象，为数据整合做准备
  sample_matrix_list[[sample[i]]] <- singlet_seurat 
  
  singlet_seurat <- NormalizeData(singlet_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  singlet_seurat <- FindVariableFeatures(singlet_seurat, selection.method = "vst", nfeatures = features_numbers)
  singlet_seurat <- ScaleData(singlet_seurat, features = rownames(singlet_seurat))
  singlet_seurat <- RunPCA(singlet_seurat, features = VariableFeatures(singlet_seurat))
  
  #PCA选择的量化指标
  pc_variance <- singlet_seurat@reductions$pca@stdev^2 / sum(singlet_seurat@reductions$pca@stdev^2) * 100
  cumulative_variance <- cumsum(pc_variance)
  threshold <- 80
  pc_num <- which(cumulative_variance >= threshold)[1]
  cat("当累计方差解释比例达到", threshold, "%时，应选择前", pc_num, "个主成分。\n")
  print(data.frame(PC = 1:length(cumulative_variance), 累计方差解释比例 = cumulative_variance))
    
  #肘形图
  pdf(paste0(outfolder,"clean_", sample[i], "_ElbowPlot.pdf"), width = 20, height = 10)
  p<-ElbowPlot(seurat_obj, ndims = 30)
  print(p)
  dev.off()

  singlet_seurat <- FindNeighbors(singlet_seurat, dims = dim_range)
  
  #聚类树
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) {
  singlet_seurat <- FindClusters(singlet_seurat, graph.name = "RNA_snn", resolution = res, algorithm = 1)
  }
  apply(singlet_seurat@meta.data[,grep("RNA_snn_res",colnames(singlet_seurat@meta.data))],2,table)
  p_tree <- clustree(singlet_seurat@meta.data,  prefix = "RNA_snn_res.")
  pdf(paste0(outfolder,"clean_",sample[i], "_resolution.pdf"), width = 12, height = 10)
  print(p_tree)
  dev.off()


  singlet_seurat <- FindClusters(singlet_seurat, resolution = resolution_arg)
  singlet_seurat <- RunUMAP(singlet_seurat, dims = dim_range)
  singlet_seurat <- RunTSNE(singlet_seurat, dims = dim_range)
  
  pdf(paste0(outfolder,"clean_", sample[i], "_UMI_percent.mt.pdf"), width = 20, height = 10)
  P <- FeaturePlot(object = singlet_seurat, features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 3, pt.size = 0.5)
  print(P)
  dev.off()

  pdf(paste0(outfolder,"clean_", sample[i], "_umap_tsne.pdf"), width = 20, height = 10)
  p1<-DimPlot(singlet_seurat, reduction = "umap", label = TRUE)
  p2<-DimPlot(singlet_seurat, reduction = "tsne", label = TRUE)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  #创造一个loupe文件
  create_loupe_from_seurat(singlet_seurat,output_dir=outfolder,output_name=paste0(sample[i],".cloupe"))
  sample_list[[sample[i]]] <- singlet_seurat
}

save(sample_list, file = paste0(outfolder, "sample_list.RData"))
save(sample_matrix_list, file = paste0(outfolder, "sample_matrix_list.RData"))
