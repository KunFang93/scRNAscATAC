library(Seurat)
library(ggplot2)
library(writexl)
library(SeuratDisk)
library(sctransform)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(foreach)
library(doParallel)
library(readxl)
library(ggpatt)

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

DimPlot(NTs.seu, label = TRUE, repel = TRUE, label.size = 6)

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

DimPlot(PTs.seu, label = TRUE, repel = TRUE, label.size = 6)

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

DimPlot(RTs.seu, label = TRUE, repel = TRUE, label.size = 6)

# Fig.2A calculate the composition barplot
my_levels <- c("Breast cancer cells",
               "Basal cells","Luminal cells",
               "Stromal cells","Adipocyte-like cells","Fibroblast-like cells",
               "Myeloid cells","NK cells","Neutrophils")
all_meta<- rbind(PTs.seu@meta.data[,c('orig.ident','cell.annot')],
                 RTs.seu@meta.data[,c('orig.ident','cell.annot')])
all_meta$cell.annot <- factor(all_meta$cell.annot, levels= my_levels)
ggplot(all_meta, aes(x=orig.ident, fill=cell.annot)) + 
  geom_bar()+theme(text = element_text(face="bold",size=16))+
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# Fig.2B plot correlation map for TTs
PTs_mat <- PTs.seu@assays$integrated@scale.data
RTs_mat <- RTs.seu@assays$integrated@scale.data
TTs_mat <- cbind(PTs_mat,RTs_mat)

bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
  library(ff, quietly = TRUE)
  require(doMC)
  if(ncore=="all"){
    ncore = parallel::detectCores()
    registerDoMC(cores = ncore)
  } else{
    registerDoMC(cores = ncore)
  }

  NCOL <- ncol(x)

  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))

  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }

  gc()
  return(corMAT)
}

# remove few cells in order to set nblocks
TTs_cor_mat <- bigcorPar(TTs_mat[1:22583,1:18650],nblocks = 10)
# TTs_cor_mat <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cor_mat.rds')
TTs_cor_df <- as.ffdf(TTs_cor_mat)[,]
saveRDS(TTs_cor_df,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cor_mat.rds')

# TTs_cor_mat  <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cor_mat.rds')
TTs_cor_mat <- as.matrix(TTs_cor_mat)

# other information
cluster_PTs <- data.frame(cell_id = colnames(PTs_mat))
cluster_PTs$sample <- separate(data = cluster_PTs, col = cell_id, into = c("left", "right"), sep = "_")$left
PTs.seu$barcodes <- colnames(PTs.seu)
cluster_PTs$celltype <- plyr::mapvalues(
  x = cluster_PTs$cell_id, 
  from = PTs.seu$barcodes, 
  to = PTs.seu$cell.annot
)
cluster_PTs$celltype <- ifelse(cluster_PTs$celltype %in% c(1), "Cancer Cells", 
                               "Other Cells")

cluster_RTs <- data.frame(cell_id = colnames(RTs_mat))
cluster_RTs$sample <- separate(data = cluster_RTs, col = cell_id, into = c("left", "right"), sep = "_")$left
RTs.seu$barcodes <- colnames(RTs.seu)
cluster_RTs$celltype <- plyr::mapvalues(
  x = cluster_RTs$cell_id, 
  from = RTs.seu$barcodes, 
  to = RTs.seu$cell.annot
)
cluster_RTs$celltype <- ifelse(cluster_RTs$celltype %in% c(1), "Cancer Cells", 
                               "Other Cells")

cluster_all <- rbind(cluster_PTs,cluster_RTs)
cluster_all <- head(cluster_all,18650)

colnames(TTs_cor_mat) <- cluster_all$cell_id[1:18650]
rownames(TTs_cor_mat) <- cluster_all$cell_id[1:18650]

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=12)

# since the heatmap size is too big, random select 2000 cells from each to visualze
set.seed(34)
PTs_Cancers <- cluster_PTs[cluster_PTs$celltype=="Cancer Cells",]
PTs_Cancers.sample <- PTs_Cancers[sample(nrow(PTs_Cancers),1200),]
PTs_Others <- cluster_PTs[cluster_PTs$celltype=="Other Cells",]
PTs_Others.sample <- PTs_Others[sample(nrow(PTs_Others),900),]
RTs_Cancers <- cluster_RTs[cluster_RTs$celltype=="Cancer Cells",]
RTs_Cancers.sample <- RTs_Cancers[sample(nrow(RTs_Cancers),1000),]
RTs_Others <- cluster_RTs[cluster_RTs$celltype=="Other Cells",]
RTs_Others.sample <- RTs_Others[sample(nrow(RTs_Others),1000),]

cluster_PTs_sampled <- rbind(PTs_Cancers.sample,PTs_Others.sample)
cluster_PTs_sampled <- cluster_PTs_sampled[ order(as.numeric(row.names(cluster_PTs_sampled))), ]

cluster_RTs_sampled <- rbind(RTs_Cancers.sample,RTs_Others.sample)
cluster_RTs_sampled <- cluster_RTs_sampled[ order(as.numeric(row.names(cluster_RTs_sampled))), ]

cluster_all_sampled <- rbind(cluster_PTs_sampled,cluster_RTs_sampled)

plotCorr <- function(cluster_df,cluster_org,corr_origin_mat){
  samples <- unique(cluster_df$sample)
  cluster_sub_df <- cluster_org[cluster_org$cell_id %in% cluster_df$cell_id,]
  ha <- HeatmapAnnotation(sample=cluster_sub_df$sample,
                          col = list("sample" = setNames(color_list[1:length(samples)],samples)),
                          show_annotation_name = F,
                          annotation_name_gp = gpar(fontsize = 24, fontface="bold"),
                          annotation_legend_param = list(title_gp = gpar(fontsize = 16, fontface="bold"), labels_gp = gpar(fontsize = 16),grid_height=unit(18,"mm"), grid_width=unit(20,"mm")),
                          gap=unit(1, "mm")
                          
  )
  
  ra <- rowAnnotation(celltype=cluster_sub_df$celltype,
                      #use_raster = TRUE,
                      #text = row_anno_text(LETTERS[1:length(groups)], gp = gpar(fontsize = 20)),
                      width = unit(20, "mm"), ## annotation bar width
                      col = list(celltype=c("Cancer Cells"="gray11","Other Cells"="gray70")),
                      #gp = gpar(fontsize = 1),
                      show_annotation_name = F,
                      #annotation_name_gp = gpar(fontsize = 1),
                      annotation_legend_param = list(title_gp = gpar(fontsize = 16, fontface="bold"), labels_gp = gpar(fontsize = 16),grid_height=unit(18,"mm"), grid_width=unit(20,"mm")),
                      gap=unit(1, "mm")
  )
  
  cor_mat <- corr_origin_mat[as.numeric(rownames(cluster_sub_df)),
                             as.numeric(rownames(cluster_sub_df))]
  print(Heatmap(cor_mat,
                name="corr_mat",
                top_annotation = ha,
                left_annotation = ra,
                cluster_rows=FALSE, 
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface="bold"), labels_gp = gpar(fontsize = 16),grid_height=unit(20,"mm"), grid_width=unit(22,"mm")),
                show_row_names = F,
                show_column_names = F,
                use_raster = F
                # use_raster=TRUE, raster_quality=1, raster_by_magick=TRUE
  ))
}

plotCorr(cluster_all_sampled,cluster_all,TTs_cor_mat)
# only cancer cells
clusters_cancer <- rbind(PTs_Cancers.sample,RTs_Cancers.sample)
clusters_other <- rbind(PTs_Others.sample,RTs_Others.sample)

plotCorr(clusters_cancer,cluster_all,TTs_cor_mat)
plotCorr(clusters_other,cluster_all,TTs_cor_mat)

# Create cell to cell correlation data.frame with tidyr::pivot_longer function, that will give the correlation score 
# of any given "cell of origin" with any "other cell" 
TTs_cor_df <- as.data.frame(TTs_cor_mat)
rownames(TTs_cor_df) <- cluster_all$cell_id
colnames(TTs_cor_df) <- cluster_all$cell_id
TTs_cor_df$cell_of_origin <- rownames(TTs_cor_df)
TTs_cor_df <- tidyr::pivot_longer(TTs_cor_df, cols = seq_len(ncol(TTs_cor_df)-1),
                                  names_to="other_cell",values_to="correlation")
# Remove self correlations (e.g. cell_1 with cell_1), as it is always 1
TTs_cor_df <- TTs_cor_df[-which(TTs_cor_df$cell_of_origin == TTs_cor_df$other_cell),]
# Add cluster information (cluster of the cell of origin & cluster of the other cell)
TTs_cor_df$sample_origin <- cluster_all$sample[match(TTs_cor_df$cell_of_origin,cluster_all$cell_id)]
TTs_cor_df$sample_other <- cluster_all$sample[match(TTs_cor_df$other_cell,cluster_all$cell_id)]
TTs_cor_df$ct_origin <- cluster_all$celltype[match(TTs_cor_df$cell_of_origin,cluster_all$cell_id)]
TTs_cor_df$ct_other <- cluster_all$celltype[match(TTs_cor_df$other_cell,cluster_all$cell_id)]

intra_corr <- TTs_cor_df[TTs_cor_df$sample_origin==TTs_cor_df$sample_other,]
intra_corr1 <- intra_corr[intra_corr$ct_origin==intra_corr$ct_other,]
# Fig.2C Plot all and cancer cell boxplot
intra_corr2 <- intra_corr1[intra_corr1$ct_origin=="Cancer Cells",]
intra_corr2$type <- "Cancer Cells"
intra_corr1$type <- "All Cells"
intra_corr_all <- rbind(intra_corr1,intra_corr2)

p = ggplot(data=intra_corr_all)

p + stat_boxplot(aes(x=sample_origin,y=correlation,fill=type),
                 geom = "errorbar") + 
  geom_boxplot(aes(x=sample_origin,y=correlation,fill=type,pattern = type),
               color="black",
               outlier.shape = NA) + ylim(-0.3,0.6) + theme_classic()
all_meta<- rbind(PTs.seu@meta.data[,c('orig.ident','cell.annot')],
                 RTs.seu@meta.data[,c('orig.ident','cell.annot')])
all_meta$barcodes <- rownames(all_meta)

# load integrated TTs
TTs.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs.h5seurat')
# add cell annotaion
TTs.seu$barcodes <- colnames(TTs.seu)
TTs.seu$cell.annot <- plyr::mapvalues(
  x = TTs.seu$barcodes, 
  from = all_meta$barcodes, 
  to = as.character(all_meta$cell.annot)
)
# extract cancers cells
TTs.cancer.seu <- subset(TTs.seu, subset = cell.annot %in% c("Breast cancer cells"))
TTs.cancer.seu <- RunPCA(TTs.cancer.seu, npcs = 30, verbose = FALSE)
TTs.cancer.seu <- RunUMAP(TTs.cancer.seu, reduction = "pca", dims = 1:30)
TTs.cancer.seu <- FindNeighbors(TTs.cancer.seu, reduction = "pca", dims = 1:30)
TTs.cancer.seu <- FindClusters(TTs.cancer.seu, resolution = 0.9, algorithm = 2)

# Fig.2D
DimPlot(TTs.cancer.seu, label = TRUE, repel = TRUE, 
        label.size = 6, group.by = c('orig.ident','seurat_clusters'))

# find all markers
# Find markers
TTs.cancer.seu <- PrepSCTFindMarkers(TTs.cancer.seu)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
cluster_markers <- FindAllMarkers(TTs.cancer.seu,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  logfc.threshold = 0.25)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

DoHeatmap(TTs.cancer.seu, features = top_markers$gene) + theme(axis.text.y = element_text(size = 6))
TTs_cancer_list <- list("All"=cluster_markers,"Top20"=top_markers)
write_xlsx(TTs_cancer_list,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_cluster_markers.xlsx')

vis.markers2 <- c('FXYD3','MUCL1',
                  'NEAT1',
                  'S100A10',
                  "ESR1","PGR",'BCL2',"ERBB2",
                  'PRKD1',
                  'BRIP1',
                  'MKI67','AURKA','BIRC5', 'CCNB1',
                  'BRCA1'
)
# Fig.2E
VlnPlot(TTs.cancer.seu, features = vis.markers2,
        stack=T,pt.size=0,
        flip = T,
        add.noise = T)+#横纵轴不标记任何东西
  theme(axis.text.y = element_blank(), #不显示坐标刻度
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 14),
        legend.position = 'none')

# Fig.2F samples composition
ggplot(TTs.cancer.seu@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + 
  geom_bar(position = "fill")+theme(text = element_text(face="bold",size=16))+
  theme(panel.border = element_blank(), 
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
