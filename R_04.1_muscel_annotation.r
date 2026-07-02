library(Seurat)
library(ggplot2)

outfolder <-"./R/R_06_subcell/"
seurat_obj <- readRDS("./R/R_06_subcell/muscle_newSCT.rds")
# 1. 获取当前的聚类结果（比如分辨率resolution = 0.2时的结果）
cluster_ids <- seurat_obj@meta.data$seurat_clusters
# 2. 创建一个命名字符向量，将每个聚类编号映射到具体的细胞类型
cluster_to_celltype <- c(
  "0" = "Nep2+", 
  "1" = "flr+", 
  "3" = "Hdc+", 
  "4" = "LG13254+", 
  "5" = "lz+",
  "6" = "sesB+"
)
# 将聚类编号映射为细胞类型，并添加到元数据中
seurat_obj@meta.data$sub_celltype <- cluster_to_celltype[as.character(cluster_ids)]
# 将新创建的细胞类型列设置为默认标识
Idents(seurat_obj) <- seurat_obj@meta.data$sub_celltype
# 验证结果
print(table(seurat_obj@meta.data$sub_celltype))
#画图展示
pdf(paste0(outfolder, "Dimplot_muscle_annot_umap.pdf"), width = 20, height = 15)
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) +ggtitle("UMAP muscle annotation")+
  theme(plot.title = element_text(size = 20),  # 调整图标题字体大小
        axis.title = element_text(size = 16),  # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12), # 图例文字大小
        axis.text = element_text(size = 14))   # 调整坐标轴刻度标签字体
print(p1)
dev.off()

saveRDS(seurat_obj, file = paste0(outfolder, "muscle_annotation.rds"))
