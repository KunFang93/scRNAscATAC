library(Seurat)
library(tidyverse)
library(scuttle)
library(scDblFinder)
library(ggplot2)
library(writexl)
library(SeuratDisk)
library(celldex)
library(sctransform)
library(readxl)

############# Author: Kun Fang; Date: 07/06/23 ################################
# Load scRNA-seq
PT1 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PT1')
PT2 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PT2')
RT3 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT3')
RT4 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT4')

NT7 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/NT7')
NT8 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/NT8')
PT5 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PT5')
RT6 <- Read10X(data.dir = '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT6')

PT1.seu <- CreateSeuratObject(counts = PT1, project = "PT1", min.cells = 3, min.features = 200)
PT2.seu <- CreateSeuratObject(counts = PT2, project = "PT2", min.cells = 3, min.features = 200)
PT3.seu <- CreateSeuratObject(counts = PT5, project = "PT3", min.cells = 3, min.features = 200)
RT1.seu <- CreateSeuratObject(counts = RT3, project = "RT1", min.cells = 3, min.features = 200)
RT2.seu <- CreateSeuratObject(counts = RT4, project = "RT2", min.cells = 3, min.features = 200)
RT3.seu <- CreateSeuratObject(counts = RT6, project = "RT3", min.cells = 3, min.features = 200)
NT1.seu <- CreateSeuratObject(counts = NT7, project = "NT1", min.cells = 3, min.features = 200)
NT2.seu <- CreateSeuratObject(counts = NT8, project = "NT2", min.cells = 3, min.features = 200)

rm(PT1)
rm(PT2)
rm(RT3)
rm(RT4)
rm(NT7)
rm(NT8)
rm(PT5)
rm(RT6)

seurat_list <- c('PT1' = PT1.seu, 'PT2' = PT2.seu, 'PT3' = PT3.seu,
                 'RT1' = RT1.seu, 'RT2' = RT2.seu, 'RT3' = RT3.seu,
                 'NT1' = NT1.seu, 'NT2' = NT2.seu)

finddoublet <- function(seurat_obj){
  set.seed(34) 
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder(sce_obj)
  seu_new <- as.Seurat(sce_obj)
  return(seu_new)
}

nMAD <- function(x,nmads=3){
  xm <- median(x)
  md <- median(abs(x-xm))
  mads <- xm+nmads*md
  return(mads)
}

## profiling information for QC and perform QC
nmad=5
count = 1
sample.nfeaure.cut <- c()
sample.ncount.cut <- c()
sample.mt.cut <- c()
seurat_list_qc <- c()
for (obj in seurat_list){
  print(names(seurat_list)[count])
  # add percent.mt for qc
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # visual nfeature ncount percent.mt
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  # add doublet info
  obj = finddoublet(obj)
  # nMAD cut
  nfeature.upcut <- ceiling(nMAD(obj$nFeature_RNA,nmad))
  ncount.upcut <- ceiling(nMAD(obj$nCount_RNA,nmad))
  permt.upcut <- ceiling(nMAD(obj$percent.mt,nmad))
  sample.nfeaure.cut <- c(sample.nfeaure.cut, nfeature.upcut)
  sample.ncount.cut <- c(sample.ncount.cut, ncount.upcut)
  sample.mt.cut <- c(sample.mt.cut,permt.upcut)
  # perform QC
  obj.filt <- subset(obj, subset = nFeature_RNA <= nfeature.upcut & nFeature_RNA > nfeature.upcut/20)
  obj.filt <- subset(obj.filt, subset = nCount_RNA <= ncount.upcut & nCount_RNA > ncount.upcut/20)
  obj.filt <- subset(obj.filt, subset = percent.mt <= min(permt.upcut,25))
  obj.filt <- subset(obj.filt, subset = scDblFinder.class %in% c('singlet'))
  Idents(obj.filt) <- names(seurat_list)[count]
  # visualize again
  print(VlnPlot(obj.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  seurat_list_qc <- c(seurat_list_qc,obj.filt)
  count = count + 1
}
names(sample.mt.cut) <- names(seurat_list)
names(sample.ncount.cut) <- names(seurat_list)
names(sample.nfeaure.cut) <- names(seurat_list)
names(seurat_list_qc) <- names(seurat_list)

############################ Plot Fig.S1 For QC ################################
samples <- names(seurat_list)
samples.ncount.cut.summary <- data.frame(orig.ident=c(samples,samples), 
                                         types=c(rep('up.cut',length(samples)),
                                                 rep('down.cut',length(samples))),
                                         values=c(sample.ncount.cut,sample.ncount.cut/20))
samples.nfeature.cut.summary <- data.frame(orig.ident=c(samples,samples), 
                                           types=c(rep('up.cut',length(samples)),
                                                   rep('down.cut',length(samples))),
                                           values=c(sample.nfeature.cut,sample.nfeature.cut/20))
# mt don't need lower cut
samples.mt.cut.summary <- data.frame(orig.ident=c(samples), 
                                     types=c(rep('up.cut',length(samples))),
                                     values=c(sample.mt.cut))


metadata <- merge.seu@meta.data
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)

metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.nfeature.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)

metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.mt.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)

metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  geom_smooth(se=TRUE,level=0.9) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F)+
  geom_hline(data = samples.nfeature.cut.summary, 
             aes(yintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol = 4)

merge.qc.seu <- merge(x=seurat_list_qc[[1]], y=seurat_list_qc[2:length(seurat_list_qc)])

metadata.qc <- merge.qc.seu@meta.data
metadata.qc %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)

metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.nfeature.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)

metadata.qc %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha=.5) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(data = samples.mt.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)

metadata.qc %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  geom_smooth(se=TRUE,level=0.9) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(data = samples.ncount.cut.summary, 
             aes(xintercept = values), 
             linetype = "dashed",show.legend = F)+
  geom_hline(data = samples.nfeature.cut.summary, 
             aes(yintercept = values), 
             linetype = "dashed",show.legend = F) +
  facet_wrap(~orig.ident,ncol=4)

################## Main Figs ###############################
########### NTs ###########
NTs.seu <- merge(seurat_list_qc$NT1, y = c(seurat_list_qc$NT2), 
                 add.cell.ids = c("NT1", "NT2"), 
                 project = "scRNA_NTs")

# Run the standard workflow for visualization and clustering
NTs.seu.sct <- SCTransform(NTs.seu, method = "glmGamPoi", vars.to.regress = "percent.mt",
                           verbose = FALSE)
NTs.seu <- RunPCA(NTs.seu, npcs = 30, verbose = FALSE)
NTs.seu <- RunUMAP(NTs.seu, reduction = "pca", dims = 1:30)
NTs.seu <- RunTSNE(NTs.seu, reduction = "pca", dims = 1:30)
NTs.seu <- FindNeighbors(NTs.seu, reduction = "pca", dims = 1:30)
NTs.seu <- FindClusters(NTs.seu, resolution = 0.9, algorithm = 2)

# check the cluster location
DimPlot(NTs.seu, label=TRUE, repel = T)
# Find markers
NTs.seu <- PrepSCTFindMarkers(NTs.seu)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
NTs_markers <- FindAllMarkers(object = NTs.seu, 
                              assay='SCT',
                              only.pos = TRUE,
                              logfc.threshold = 0.25)

NTs_top <- NTs_markers %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

sheetlist <- list("allPosMarkers"=NTs_markers,"Top20Markers"=NTs_top)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/NTs_markers.xlsx')

# plot stackviolin plot
NTs.marker <- c("KRT15","KRT16","KRT17", # Basal
                "ESRP1","ELF3","RARRES1",       # Luminal
                "TOP2A", "CDK1","MKI67","CENPF", #Luminal Progenitor
                "COL4A6","COL4A5",
                "S100A10", "ID1",
                "CXCL1", "CXCL8", 
                "CD24")

VlnPlot(NTs.seu, features = NTs.marker,
        stack=T,pt.size=0,
        flip = T,
        add.noise = T) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 14),
        legend.position = 'none')

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
SaveH5Seurat(NTs.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/NTs',overwrite = TRUE)

NTs_top <- NTs_markers %>%
  group_by(cluster) %>%
  top_n(n = 10,
        wt = avg_log2FC)

DoHeatmap(NTs.seu, features = NTs_top$gene,label=FALSE) + theme(text = element_text(face="bold",size=8))

######### TTs ##########
# process TTs
seu.TTs.part1 <- merge(seurat_list_qc$PT1, y = c(seurat_list_qc$PT2, 
                                                 seurat_list_qc$RT1,
                                                 seurat_list_qc$RT2), 
                       add.cell.ids = c("PT1", "PT2", "RT1", "RT2"), 
                       project = "scRNA_09302020")

seu.TTs.part2 <- merge(seurat_list_qc$PT3, y = c(seurat_list_qc$RT3), 
                       add.cell.ids = c("PT3", "RT3"), 
                       project = "scRNA_01112021")


seu.batchs.list <- c("seu.TTs.batch1"=seu.TTs.part1, "seu.TTs.batch2"=seu.TTs.part2) 
seurat.qc.sct.list <- lapply(X = seu.batchs.list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seurat.qc.sct.list, nfeatures = 3000)
seurat.qc.sct.list <- PrepSCTIntegration(object.list = seurat.qc.sct.list, anchor.features = features)

batchs.anchors <- FindIntegrationAnchors(object.list = seurat.qc.sct.list, normalization.method = "SCT",
                                         anchor.features = features)

common_genes <- Reduce(intersect, list(rownames(seu.TTs.part1),rownames(seu.TTs.part2)))

TTs.combined.sct <- IntegrateData(anchorset = batchs.anchors, 
                                  normalization.method = "SCT",
                                  features.to.integrate = common_genes)

TTs.combined.sct <- RunPCA(TTs.combined.sct, npcs = 30, verbose = FALSE)
TTs.combined.sct <- RunTSNE(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- RunUMAP(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- FindNeighbors(TTs.combined.sct, reduction = "pca", dims = 1:30)
TTs.combined.sct <- FindClusters(TTs.combined.sct, resolution = 0.9)
SaveH5Seurat(TTs.combined.sct,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs',overwrite = TRUE)

#### PTs ####
PTs.combined.sct <- subset(TTs.combined.sct, subset = ident %in% c('PT1','PT2','PT3'))
DefaultAssay(PTs.combined.sct) <- "integrated"

PTs.combined.sct <- RunPCA(PTs.combined.sct, npcs = 30, verbose = FALSE)
PTs.combined.sct <- RunUMAP(PTs.combined.sct, reduction = "pca", dims = 1:30)
PTs.combined.sct <- FindNeighbors(PTs.combined.sct, reduction = "pca", dims = 1:30)
PTs.combined.sct <- FindClusters(PTs.combined.sct, resolution = 0.9, algorithm = 2)

# check the cluster location
DimPlot(PTs.combined.sct, label = TRUE, repel = TRUE,label.size = 9)

# Find markers
PTs.combined.sct <- PrepSCTFindMarkers(PTs.combined.sct)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
PTs_markers <- FindAllMarkers(object = PTs.combined.sct, 
                              assay='SCT',
                              only.pos = TRUE,
                              logfc.threshold = 0.25)

PTs_top <- PTs_markers %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

sheetlist <- list("allPosMarkers"=PTs_markers,"Top20Markers"=PTs_top)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PTs_markers.xlsx')
SaveH5Seurat(PTs.combined.sct,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PTs',overwrite = TRUE)

# PTs plot stackviolin plot
PTs.marker <- c("BRCA1","ESR1","BCL2","PGR","CCNB1","MKI67","AURKA",
                "BAG1", "KRT8","BRIP1", "TUBB", "TUBA1B", "KRT17", "CAV1",
                "IFI27","IGFBP4","TFF1","FTL","FGF13", "SLC3A2")

DotPlot(PTs.seu, features = PTs.marker, group.by = 'seurat_clusters', 
        cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

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

PTs_top <- PTs_top %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(PTs.seu, features = PTs_top$gene, label=FALSE) + theme(text = element_text(face="bold",size=8))

#### RTs ####
RTs.combined.sct <- subset(TTs.combined.sct, subset = ident %in% c('RT1','RT2','RT3'))
DefaultAssay(RTs.combined.sct) <- "integrated"

RTs.combined.sct <- RunPCA(RTs.combined.sct, npcs = 30, verbose = FALSE)
RTs.combined.sct <- RunUMAP(RTs.combined.sct, reduction = "pca", dims = 1:30)
RTs.combined.sct <- FindNeighbors(RTs.combined.sct, reduction = "pca", dims = 1:30)
RTs.combined.sct <- FindClusters(RTs.combined.sct, resolution = 0.9, algorithm = 2)

# check the cluster location
DimPlot(RTs.combined.sct, label = TRUE, repel = TRUE,label.size = 9)

# Find markers
RTs.combined.sct <- PrepSCTFindMarkers(RTs.combined.sct)
# Find markers for every cluster compared to all remaining cells, report only the positive ones
RTs_markers <- FindAllMarkers(object = RTs.combined.sct, 
                              assay='SCT',
                              only.pos = TRUE,
                              logfc.threshold = 0.25)

RTs_top <- RTs_markers %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

sheetlist <- list("allPosMarkers"=RTs_markers,"Top20Markers"=RTs_top)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs_markers.xlsx')
SaveH5Seurat(RTs.combined.sct,'/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs',overwrite = TRUE)

RTs.marker <- c("MUCL1","KRT19","XBP1","NEAT1","ESR1","CCNB1","MKI67","AURKA",
                "TUBB", "TUBA1B","UBE2T", "PCNA", "MSH6", "AGR2", "PDIA4",
                "KRT80","SLC3A2","CD83","NFKBIA","NFKBIZ")


DotPlot(RTs.seu, features = RTs.marker, group.by = 'seurat_clusters', 
        cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


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

RTs_top <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs_markers.xlsx',sheet = "Top20Markers")
RTs_top <- RTs_top %>%
  group_by(cluster) %>%
  top_n(n = 8,
        wt = avg_log2FC)

DoHeatmap(RTs.seu, features = RTs_top$gene, label=FALSE) +theme(text = element_text(face="bold",size=8))
