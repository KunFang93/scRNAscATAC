library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(writexl)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(ggrepel)

############# Author: Kun Fang; Date: 07/06/23 ################################
# Load data
NTs.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/NTs.h5seurat')
PTs.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PTs.h5seurat')
RTs.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs.h5seurat')

# add cell type annotation
NTs.seu <- RenameIdents(object = NTs.seu,
                        "0" = "Luminal cells", 
                        "1" = "Basal cells",
                        "2" = "Fibroblasts", "3" = "Endothelial cells",
                        "4" = "Myeloid cells", "5" = "Luminal cells",
                        "6" = "Luminal Progenitor",
                        "7" = "Luminal Progenitor",
                        "8" = "Basal cells", 
                        "9" = "Luminal Progenitor",
                        "10" = "Luminal cells")
NTs.seu$cell.annot <- Idents(NTs.seu)

PTs.seu <- RenameIdents(object = PTs.seu,
                        "0" = "Breast cancer cells", "1" = "Myeloid cells",
                        "2" = "Basal cells", "3" = "Adipocyte-like cells",
                        "4" = "Breast cancer cells", 
                        "5" = "Breast cancer cells",
                        "6" = "Luminal cells","7" = "Breast cancer cells",
                        "8" = "Breast cancer cells", 
                        "9" = "Breast cancer cells",
                        "10" = "Neutrophils","11"="Breast cancer cells",
                        "12" = "Breast cancer cells", 
                        "13" = "Fibroblast-like cells",
                        "14" = "Breast cancer cells")
PTs.seu$cell.annot <- Idents(PTs.seu)

RTs.seu <- RenameIdents(object = RTs.seu,
                        "0" = "Breast cancer cells", "1" = "Myeloid cells",
                        "2" = "Breast cancer cells", "3" = "Luminal cells",
                        "4" = "Stromal cells", 
                        "5" = "Stromal cells",
                        "6" = "Breast cancer cells",
                        "7" = "Stromal cells",
                        "8" = "Breast cancer cells", "9" = "Stromal cells",
                        "10" = "Luminal cells","11"="Breast cancer cells",
                        "12" = "NK cells")
RTs.seu$cell.annot <- Idents(RTs.seu)

# read module list
# ER signaling, HER2 signaling, proliferation module,tumor invasion/metastasis
# immune response, angiogenesis, apoptosis phenotypes,resistance 
module_list <- c("ESR1","ERBB2","AURKA","PLAU","STAT1","VEGF")
module_genes_list <- lapply(X=module_list,FUN = function(x){
  geneslist <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/Modules.xlsx',
                          sheet = x, col_names="genes")$genes
  # get the common gene bewteen dataset and modules
  geneslist <- intersect(geneslist, rownames(TTs.cancer.seu))
})
names(module_genes_list) <- module_list

# Fig.3A
for (i in module_list){
  print(i)
  TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                   features = list(module_genes_list[[i]]),
                                   name=paste0(i,"_module"))
}
for (i in module_list[1:6]){
  current_features<-paste0(i,"_module1")
  print(FeaturePlot(TTs.cancer.seu, features=current_features, repel = TRUE,
                    min.cutoff = "q30", max.cutoff = "q99", pt.size=0.5, 
                    order=TRUE) + theme(legend.position = "right") +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
}


RTs_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(0,1,2,5,6),]$barcodes
PTs_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(3,4,7,8,9,10,12),]$barcodes

RTsvsPTs_markers <- FindMarkers(TTs.cancer.seu,
                                ident.1 = RTs_cells,
                                ident.2 = PTs_cells,
                                assay='SCT',
                                only.pos = TRUE,
                                logfc.threshold = 0.25)
RTsvsPTs_markers<-RTsvsPTs_markers[order(RTsvsPTs_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]

sheetlist <- list("RTs.vs.PTs"=RTsvsPTs_markers)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs_vs_PTs.xlsx')
SaveH5Seurat(TTs.cancer.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223')

# TTs.cancer.seu <-LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223.h5seurat')
# analysis based on module
all_meta <- TTs.cancer.seu@meta.data
modules_df <- all_meta[,c("ESR1_module1", "ERBB2_module1","AURKA_module1",
                          "STAT1_module1","VEGF_module1","CASP3_module1")]
modules_df_means <- aggregate(modules_df, by = list(all_meta$seurat_clusters), FUN = mean)
rownames(modules_df_means) <- paste0("cluster",modules_df_means$Group.1)
modules_df_means$Group.1 <- NULL
modules_df_means_z <- apply(modules_df_means, 2, scale)
rownames(modules_df_means_z) <- rownames(modules_df_means)
cluster_df2<- data.frame("clusters"=rownames(modules_df_means_z),
                         "sample"=c("RT","RT","RT","PT","PT","RT","RT",
                                    "PT","PT","PT","PT","PT&RT","PT"))

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=13)
names(color_list) <- cluster_df2$clusters
ra1 <- rowAnnotation(cluster=cluster_df2$clusters,
                     width = unit(20, "mm"), ## annotation bar width
                     col = list(cluster=color_list),
                     show_annotation_name = F,
                     gap=unit(1, "mm")
)

ra2 <- rowAnnotation(sample=cluster_df2$sample,
                     width = unit(20, "mm"), ## annotation bar width
                     col = list(sample=c("PT"="bisque2","RT"="bisque4","PT&RT"="bisque3")),
                     show_annotation_name = F,
                     gap=unit(1, "mm")
)
# Fig.3B
Heatmap(modules_df_means_z,
        left_annotation = ra2,
        right_annotation = ra1,
        show_row_names = T, 
        cluster_columns = F,
        col = rev(brewer.pal(n = 11, name = "RdBu")))
draw(ht_list, merge_legend = TRUE)
# Fig.3C
modules_pca <- prcomp(modules_df)
cluster_means <- aggregate(modules_pca$x, by = list(all_meta$seurat_clusters), FUN = mean)
cluster_means$shape <- c("0","0","0","1","2","3","4","5","6","6","6","7","8")
ggplot(cluster_means, aes(x=PC1, y=PC2, color=Group.1, shape=shape)) +
  geom_point(size=8)+theme(text = element_text(face="bold",size=16))+
  scale_shape_manual(values=seq(0,15))+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# cell state dimplot
TTs.cancer.seu <- RenameIdents(object = TTs.cancer.seu,
                               "0" = "RT_CS1", "1" = "RT_CS1",
                               "2" = "RT_CS1", "3" = "PT_CS2",
                               "4" = "PT_CS4", 
                               "5" = "RT_CS2",
                               "6" = "RT_CS3","7" = "PT_CS3",
                               "8" = "PT_CS1", 
                               "9" = "PT_CS1",
                               "10" = "PT_CS1","11"="PRT_CS",
                               "12" = "PT_CS5")
TTs.cancer.seu$cell.state <- Idents(TTs.cancer.seu)
# Fig.3D
DimPlot(TTs.cancer.seu,label.size = 8,label = T,repel=T) + NoLegend()

# Fig.3E
# RT_CS1/2/3 and RPT_CS module
PTs_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(3,4,7,8,9,10,12),]$barcodes
RT_CS1_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(0,1,2),]$barcodes
RT_CS2_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(5),]$barcodes
RT_CS3_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(6),]$barcodes
RPT_CS_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(11),]$barcodes
RTs_cells <- TTs.cancer.seu@meta.data[TTs.cancer.seu$seurat_clusters %in% c(0,1,2,5,6),]$barcodes
# RT_CS1
RT_CS1_markers <- FindMarkers(TTs.cancer.seu,
                              ident.1 = RT_CS1_cells,
                              ident.2 = PTs_cells,
                              assay='SCT',
                              logfc.threshold = 0.25)

EnhancedVolcano(RT_CS1_markers , 
                rownames(RT_CS1_markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
RT_CS1_markers$gene <- rownames(RT_CS1_markers)

RT_CS1_markers<-RT_CS1_markers[order(RT_CS1_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS1_module <- rbind(head(RT_CS1_markers,100)) # highest values

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(rownames(RT_CS1_module)),
                                 name="RT_CS1_module")

FeaturePlot(TTs.cancer.seu, features="RT_CS1_module1", repel = TRUE,
            min.cutoff = "q30", max.cutoff = "q96", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#RT_CS2
RT_CS2_markers <- FindMarkers(TTs.cancer.seu,
                              ident.1 = RT_CS2_cells,
                              ident.2 = PTs_cells,
                              assay='SCT',
                              logfc.threshold = 0.25)

EnhancedVolcano(RT_CS2_markers , 
                rownames(RT_CS2_markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
RT_CS2_markers$gene <- rownames(RT_CS2_markers)

RT_CS2_markers<-RT_CS2_markers[order(RT_CS2_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS2_module <- rbind(head(RT_CS2_markers,100)) # highest values

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(rownames(RT_CS2_module)),
                                 name="RT_CS2_module")

FeaturePlot(TTs.cancer.seu, features="RT_CS2_module1", repel = TRUE,
            min.cutoff = "q30", max.cutoff = "q98", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#RT_CS3
RT_CS3_markers <- FindMarkers(TTs.cancer.seu,
                              ident.1 = RT_CS3_cells,
                              ident.2 = PTs_cells,
                              assay='SCT',
                              logfc.threshold = 0.25)

EnhancedVolcano(RT_CS3_markers , 
                rownames(RT_CS3_markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
RT_CS3_markers$gene <- rownames(RT_CS3_markers)

RT_CS3_markers<-RT_CS3_markers[order(RT_CS3_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS3_module <- rbind(head(RT_CS3_markers,100)) # highest values

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(rownames(RT_CS3_module)),
                                 name="RT_CS3_module")

FeaturePlot(TTs.cancer.seu, features="RT_CS3_module1", repel = TRUE,
            min.cutoff = "q30", max.cutoff = "q99", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#RPT_CS
RPT_CS_markers <- FindMarkers(TTs.cancer.seu,
                              ident.1 = RPT_CS_cells,
                              ident.2 = c(PTs_cells, RTs_cells),
                              assay='SCT',
                              logfc.threshold = 0.25)

EnhancedVolcano(RPT_CS_markers , 
                rownames(RPT_CS_markers),
                x ="avg_log2FC", 
                y ="p_val_adj")
RPT_CS_markers$gene <- rownames(RPT_CS_markers)

RPT_CS_markers<-RPT_CS_markers[order(RPT_CS_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RPT_CS_module <- rbind(head(RPT_CS_markers,100)) # highest values

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(rownames(RPT_CS_module)),
                                 name="RPT_CS_module")

FeaturePlot(TTs.cancer.seu, features="RPT_CS_module1", repel = TRUE,
            min.cutoff = "q30", max.cutoff = "q99", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


sheetlist <- list("RT_CS1_all"=RT_CS1_markers, "RT_CS2_all"=RT_CS2_markers,
                  "RT_CS3_all"=RT_CS3_markers, "RPT_CS_all"=RPT_CS_markers)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_marker_all.xlsx')

sheetlist <- list("RT_CS1_module"=RT_CS1_module, "RT_CS2_module"=RT_CS2_module,
                  "RT_CS3_module"=RT_CS3_module, "RPT_CS_module"=RPT_CS_module)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_module_all.xlsx')

rt_cs1 <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_module_all.xlsx', sheet='RT_CS1_module')
rt_cs2 <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_module_all.xlsx', sheet='RT_CS2_module')
rt_cs3 <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_module_all.xlsx', sheet='RT_CS3_module')

x <- list(
  RT_CS1 = rt_cs1$gene, 
  RT_CS2 = rt_cs2$gene, 
  RT_CS3 = rt_cs3$gene
)

p1<- ggVennDiagram(x, label_alpha = 0, edge_size=0,label_size = 10) 
p1 + ggplot2::scale_fill_gradient(low="white",high = "red")

getVennOverlap <- function(lsvenn) {
  ItemsList <- gplots::venn(lsvenn, show.plot = FALSE)
  print(lengths(attributes(ItemsList)$intersections))
  return(attributes(ItemsList)$intersections)
}

x_intersections <- getVennOverlap(x)
# core signatures
list_df <- list('RT_CS1'=rt_cs1[rt_cs1$gene %in% x_intersections$`RT_CS1:RT_CS2:RT_CS3`,c('avg_log2FC','gene')],
                'RT_CS2'=rt_cs2[rt_cs2$gene %in% x_intersections$`RT_CS1:RT_CS2:RT_CS3`,c('avg_log2FC','gene')],
                'RT_CS3'=rt_cs3[rt_cs3$gene %in% x_intersections$`RT_CS1:RT_CS2:RT_CS3`,c('avg_log2FC','gene')])
cores_signature <- Reduce(function(x, y) merge(x,y,by='gene'),list_df)
colnames(cores_signature) <- c('gene', 'avg_log2FC.RT_CS1', 'avg_log2FC.RT_CS2',
                               'avg_log2FC.RT_CS3')
write.csv(cores_signature,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/core_signatures.csv',row.names = F)
