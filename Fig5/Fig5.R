library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(SeuratDisk)
library(tidyverse)
library(plyr)
library(ggVennDiagram)
library(stringr)
library(ComplexHeatmap)
library(reshape2)
library(readxl)
library(writexl)
set.seed(1234)

############# Author: Kun Fang; Date: 07/07/23 ################################
############# Part 1: create ATAC Seurat obj #################
find_common_peaks <- function(peakfiles, width.upcut=10000, width.lowcut=20){
  peaks.list <- lapply(peakfiles, FUN = function(files) {
    read.table(files, col.names = c("chr", "start", "end"))
  })
  peaks.gr.list <- lapply(peaks.list, FUN = function(files) {
    makeGRangesFromDataFrame(files)
  })
  myGRangesList<-GRangesList(peaks.gr.list)   
  combined.peaks <- reduce(unlist(myGRangesList))
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < width.upcut & peakwidths > width.lowcut]
  return(combined.peaks)
}

create_merge_atac_seurat <- function(sc.csv, frag.tsv, combined.peaks, cells.cut=500){
  md <- read.table(
    file = sc.csv,
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  # perform an initial filtering of low count cells
  md <- md[md$passed_filters > cells.cut, ]
  
  frags <- CreateFragmentObject(
    path = frag.tsv,
    cells = rownames(md)
  )
  
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    cells = rownames(md)
  )
  
  sample_assay <- CreateChromatinAssay(counts, fragments = frags)
  sample.seu <- CreateSeuratObject(sample_assay, assay = "ATAC", meta.data=md)
  return(sample.seu)
} 

NTs.combined.peaks <- find_common_peaks(
  c('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT7/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT8/peaks.bed')
)

NT7.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT7/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT7/fragments.tsv.gz",
                                         NTs.combined.peaks)
NT8.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT8/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NT8/fragments.tsv.gz",
                                         NTs.combined.peaks)

# add information to identify dataset of origin
NT7.atac.seu$dataset <- 'NT1'
NT8.atac.seu$dataset <- 'NT2'

# merge all datasets, adding a cell ID to make sure cell names are unique
NTs.combined.seu <- merge(
  x = NT7.atac.seu,
  y = list(NT8.atac.seu),
  add.cell.ids = c("NT1","NT2")
)

NTs.combined.seu <- RunTFIDF(NTs.combined.seu)
NTs.combined.seu <- FindTopFeatures(NTs.combined.seu, min.cutoff = 20)
NTs.combined.seu <- RunSVD(NTs.combined.seu)
NTs.combined.seu <- RunUMAP(NTs.combined.seu, dims = 1:50, reduction = 'lsi')

saveRDS(NTs.combined.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/NTs.atac.rds')

TTs.combined.peaks <- find_common_peaks(
  c('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT1/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT2/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT5/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT3/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT4/peaks.bed',
    '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT6/peaks.bed')
)

PT1.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT1/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT1/fragments.tsv.gz",
                                         TTs.combined.peaks)
PT2.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT2/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT2/fragments.tsv.gz",
                                         TTs.combined.peaks)

PT5.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT5/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT5/fragments.tsv.gz",
                                         TTs.combined.peaks)
RT3.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT3/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT3/fragments.tsv.gz",
                                         TTs.combined.peaks)
RT4.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT4/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT4/fragments.tsv.gz",
                                         TTs.combined.peaks)
RT6.atac.seu <- create_merge_atac_seurat("/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT6/singlecell.csv",
                                         "/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT6/fragments.tsv.gz",
                                         TTs.combined.peaks)

# add information to identify dataset of origin
PT1.atac.seu$dataset <- 'PT1'
PT2.atac.seu$dataset <- 'PT2'
PT5.atac.seu$dataset <- 'PT3'
RT3.atac.seu$dataset <- 'RT1'
RT4.atac.seu$dataset <- 'RT2'
RT6.atac.seu$dataset <- 'RT3'

# regionstate bugs not accept h5seurat
saveRDS(PT1.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT1.atac.rds')
saveRDS(PT2.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT2.atac.rds')
saveRDS(PT5.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT5.atac.rds')
saveRDS(RT3.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT3.atac.rds')
saveRDS(RT4.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT4.atac.rds')
saveRDS(RT6.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT6.atac.rds')

############# Part 2: combine ATAC Seurat obj #################
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# annotations <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/annotation.rds')

TTs.batch1.seu <- merge(
  x = PT1.atac.seu,
  y = list(PT2.atac.seu,RT3.atac.seu,RT4.atac.seu),
  add.cell.ids = c("PT1","PT2","RT1","RT2")
)

TTs.batch2.seu <- merge(
  x = PT5.atac.seu,
  y = list(RT6.atac.seu),
  add.cell.ids = c("PT3","RT3")
)

profile_atac_qc <- function(sample.seu, annotations){
  # add the gene information to the object
  Annotation(sample.seu) <- annotations
  
  # compute nucleosome signal score per cell
  sample.seu <- NucleosomeSignal(object = sample.seu)
  
  # compute TSS enrichment score per cell
  sample.seu <- TSSEnrichment(object = sample.seu, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  sample.seu$pct_reads_in_peaks <- sample.seu$peak_region_fragments / sample.seu$passed_filters * 100
  sample.seu$blacklist_ratio <- sample.seu$blacklist_region_fragments / sample.seu$peak_region_fragments
  
  sample.seu$high.tss <- ifelse(sample.seu$TSS.enrichment > 2, 'High', 'Low')
  print(TSSPlot(sample.seu, group.by = 'high.tss') + NoLegend())
  
  sample.seu$nucleosome_group <- ifelse(sample.seu$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  print(FragmentHistogram(object = sample.seu, group.by = 'nucleosome_group'))
  
  print(VlnPlot(
    object = sample.seu,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  ))
  return(sample.seu)
}

TTs.batch1.seu <- profile_atac_qc(TTs.batch1.seu,annotations)
TTs.batch2.seu <- profile_atac_qc(TTs.batch2.seu,annotations)

TTs.batch1.seu <- subset(
  x = TTs.batch1.seu,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

TTs.batch2.seu <- subset(
  x = TTs.batch2.seu,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

TTs.batch1.seu <- FindTopFeatures(TTs.batch1.seu, min.cutoff = 10)
TTs.batch1.seu <- RunTFIDF(TTs.batch1.seu)
TTs.batch1.seu <- RunSVD(TTs.batch1.seu)
TTs.batch1.seu$dataset <- "batch1"

TTs.batch2.seu <- FindTopFeatures(TTs.batch2.seu, min.cutoff = 10)
TTs.batch2.seu <- RunTFIDF(TTs.batch2.seu)
TTs.batch2.seu <- RunSVD(TTs.batch2.seu)
TTs.batch2.seu$dataset <- "batch2"

# simple merge
TTs.atac.combined <- merge(TTs.batch1.seu, TTs.batch2.seu)

# process the combined dataset
TTs.atac.combined <- FindTopFeatures(TTs.atac.combined, min.cutoff = 10)
TTs.atac.combined <- RunTFIDF(TTs.atac.combined)
TTs.atac.combined <- RunSVD(TTs.atac.combined)
TTs.atac.combined <- RunUMAP(TTs.atac.combined, reduction = "lsi", dims = 2:30)
TTs.atac.combined$sample <- str_split_i(rownames(TTs.atac.combined@meta.data),"_",1)
TTs.atac.combined$dataset <- NULL
p1 <- DimPlot(TTs.atac.combined, group.by = "sample")

# anchor integration
peaks.use <- Reduce(intersect, list(rownames(TTs.batch1.seu), 
                                    rownames(TTs.batch2.seu)))

integration.anchors <- FindIntegrationAnchors(
  object.list = list(TTs.batch1.seu, TTs.batch2.seu),
  anchor.features = peaks.use,
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = TTs.atac.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated$sample <- str_split_i(rownames(integrated@meta.data),"_",1)
integrated$dataset <- NULL

p2 <- DimPlot(integrated, group.by = "sample")

(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
DimPlot(integrated, group.by = "sample")

TTs.atac.seu <- integrated
TTs.atac.seu <- RunTFIDF(TTs.atac.seu)
TTs.atac.seu <- FindTopFeatures(TTs.atac.seu, min.cutoff = 'q0')
TTs.atac.seu <- RunSVD(TTs.atac.seu)

DepthCor(TTs.atac.seu)

TTs.atac.seu <- RunUMAP(object = TTs.atac.seu, reduction = 'integrated_lsi', reduction.name = 'umap.atac',dims = 2:30)
TTs.atac.seu <- FindNeighbors(object = TTs.atac.seu, reduction = 'lsi', dims = 2:30)
TTs.atac.seu <- FindClusters(object = TTs.atac.seu, verbose = FALSE, algorithm = 3)
DimPlot(object = TTs.atac.seu, label =TRUE,group.by = c('sample','seurat_clusters'),repel=T) 

# add the gene activity matrix to the Seurat object as a new assay and normalize it
TTs.gene.activities <- GeneActivity(TTs.atac.seu)
TTs.atac.seu[['ACTIVITY']] <- CreateAssayObject(counts = TTs.gene.activities)
DefaultAssay(TTs.atac.seu) <- "ACTIVITY"
TTs.atac.seu <- NormalizeData(
  object = TTs.atac.seu,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(TTs.atac.seu$nCount_ACTIVITY)
)
TTs.atac.seu <- ScaleData(TTs.atac.seu, features = rownames(TTs.atac.seu))
saveRDS(TTs.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.combined.rds')

############# Part 3: integrate ATAC and RNA Seurat obj and get pure cell states atac seu.obj #################
PTs.rna.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/PTs.h5seurat')
RTs.rna.seu <- LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RTs.h5seurat')
PTs.rna.seu <- RenameIdents(object = PTs.rna.seu,
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
PTs.rna.seu$cell.annot <- Idents(PTs.rna.seu)

RTs.rna.seu <- RenameIdents(object = RTs.rna.seu,
                            "0" = "Breast cancer cells", "1" = "Myeloid cells",
                            "2" = "Breast cancer cells", "3" = "Luminal cells",
                            "4" = "Stromal cells", 
                            "5" = "Stromal cells",
                            "6" = "Breast cancer cells",
                            "7" = "Stromal cells",
                            "8" = "Breast cancer cells", "9" = "Stromal cells",
                            "10" = "Luminal cells","11"="Breast cancer cells",
                            "12" = "NK cells")
RTs.rna.seu$cell.annot <- Idents(RTs.rna.seu)

# bulid atac.seu with transferred labels
PTs.atac.seu <- subset(TTs.atac.seu, subset = sample %in% c("PT1","PT2","PT3"))
RTs.atac.seu <- subset(TTs.atac.seu, subset = sample %in% c("RT1","RT2","RT3"))

PTs.atac.seu <- RunTFIDF(PTs.atac.seu)
PTs.atac.seu <- FindTopFeatures(PTs.atac.seu, min.cutoff = 'q0')
PTs.atac.seu <- RunSVD(PTs.atac.seu)
PTs.atac.seu <- RunUMAP(object = PTs.atac.seu, reduction = 'lsi', reduction.name = 'umap.atac',dims = 2:30)

RTs.atac.seu <- RunTFIDF(RTs.atac.seu)
RTs.atac.seu <- FindTopFeatures(RTs.atac.seu, min.cutoff = 'q0')
RTs.atac.seu <- RunSVD(RTs.atac.seu)
RTs.atac.seu <- RunUMAP(object = RTs.atac.seu, reduction = 'lsi', reduction.name = 'umap.atac',dims = 2:30)

# dimplot before integration
p1 <- DimPlot(PTs.rna.seu, group.by = "cell.annot", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(PTs.atac.seu, group.by = "sample", label = TRUE,repel=TRUE) + NoLegend() + ggtitle("ATAC")
p1 + p2

p1 <- DimPlot(RTs.rna.seu, group.by = "cell.annot", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(RTs.atac.seu, group.by = "sample", label = TRUE,repel=TRUE) + NoLegend() + ggtitle("ATAC")
p1 + p2

# SCT version Identify anchors
PTs.transfer.anchors.sct <- FindTransferAnchors(
  reference = PTs.rna.seu,
  query = PTs.atac.seu,
  features = VariableFeatures(object = PTs.rna.seu),
  reference.assay = "integrated", query.assay = "ACTIVITY",
  reduction = 'cca',
  normalization.method ='SCT'
)

PTs.predicted.labels.sct <- TransferData(
  anchorset = PTs.transfer.anchors.sct,
  refdata = PTs.rna.seu$seurat_clusters,
  weight.reduction = PTs.atac.seu[['lsi']],
  dims = 2:30
)

PTs.atac.seu <- AddMetaData(object = PTs.atac.seu, metadata = PTs.predicted.labels.sct)

plot1 <- DimPlot(
  object = PTs.rna.seu,
  group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = PTs.atac.seu,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot1 + plot2
saveRDS(PTs.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PTs.atac.rds')

RTs.transfer.anchors.sct <- FindTransferAnchors(
  reference = RTs.rna.seu,
  query = RTs.atac.seu,
  features = VariableFeatures(object = RTs.rna.seu),
  reference.assay = "integrated", query.assay = "ACTIVITY",
  reduction = 'cca',
  normalization.method ='SCT'
)

RTs.predicted.labels.sct <- TransferData(
  anchorset = RTs.transfer.anchors.sct,
  refdata = RTs.rna.seu$seurat_clusters,
  weight.reduction = RTs.atac.seu[['lsi']],
  dims = 2:30
)

RTs.atac.seu <- AddMetaData(object = RTs.atac.seu, metadata = RTs.predicted.labels.sct)

plot1 <- DimPlot(
  object = RTs.rna.seu,
  group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = RTs.atac.seu,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot1 + plot2
saveRDS(RTs.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTs.atac.rds')

# co-embedding only for visualization purpose 
rna_atac_coemb <- function(rna.seu,atac.seu){
  genes.use <- VariableFeatures(rna.seu)
  refdata <- GetAssayData(rna.seu, assay = "integrated", slot = "data")[genes.use, ]
  # Identify anchors
  transfer.anchors <- FindTransferAnchors(reference = rna.seu, 
                                          query = atac.seu, 
                                          features = VariableFeatures(object = rna.seu),
                                          reference.assay = "integrated", 
                                          query.assay = "ACTIVITY", 
                                          reduction = "cca",
                                          normalization.method ='SCT')
  # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
  # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
  imputation <- TransferData(anchorset = transfer.anchors, 
                             refdata = refdata, 
                             weight.reduction = atac.seu[["lsi"]],
                             dims = 2:30)
  atac.seu[["integrated"]] <- imputation
  
  # add predicted.id
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = rna.seu$seurat_clusters,
    weight.reduction = atac.seu[['lsi']],
    dims = 2:30
  )
  atac.seu$predicted.id <- predicted.labels$predicted.id
  rna.seu$Libtype <- 'RNA'
  atac.seu$Libtype <- 'ATAC'
  coembed <- merge(x = rna.seu, y = atac.seu)
  
  # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
  # datasets
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  return(coembed)
}
# PTs
PTs.coembed.seu <- rna_atac_coemb(PTs.rna.seu,PTs.atac.seu)
PTs.coembed.seu$seurat_clusters <- ifelse(is.na(PTs.coembed.seu$predicted.id), 
                                          PTs.coembed.seu$integrated_snn_res.0.9, 
                                          PTs.coembed.seu$predicted.id)
Idents(PTs.coembed.seu) <- PTs.coembed.seu$seurat_clusters
PTs.coembed.seu <- RenameIdents(object = PTs.coembed.seu,
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
PTs.coembed.seu$cell.annot <- Idents(PTs.coembed.seu)
DimPlot(PTs.coembed.seu, group.by = c("Libtype", "cell.annot"), label=T, label.size = 5)
# RTs
RTs.coembed.seu <- rna_atac_coemb(RTs.rna.seu,RTs.atac.seu)
RTs.coembed.seu$seurat_clusters <- ifelse(is.na(RTs.coembed.seu$predicted.id), 
                                          RTs.coembed.seu$integrated_snn_res.0.9, 
                                          RTs.coembed.seu$predicted.id)
Idents(RTs.coembed.seu) <- RTs.coembed.seu$seurat_clusters
RTs.coembed.seu <- RenameIdents(object = RTs.coembed.seu,
                                "0" = "Breast cancer cells", "1" = "Myeloid cells",
                                "2" = "Breast cancer cells", "3" = "Luminal cells",
                                "4" = "Stromal cells", 
                                "5" = "Stromal cells",
                                "6" = "Breast cancer cells",
                                "7" = "Stromal cells",
                                "8" = "Breast cancer cells", "9" = "Stromal cells",
                                "10" = "Luminal cells","11"="Breast cancer cells",
                                "12" = "NK cells")
RTs.coembed.seu$cell.annot <- Idents(RTs.coembed.seu)
DimPlot(RTs.coembed.seu, group.by = c("Libtype", "cell.annot"), label=T, label.size = 5)

saveRDS(PTs.coembed.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PTs.coembed.rds')
saveRDS(RTs.coembed.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTs.coembed.rds')

# only select cancer
# select only cancer cells atac
PTs.atac.cancer.seu <- subset(PTs.atac.seu, subset = predicted.id %in% c(0,4,5,7,8,9,11,12,14))
RTs.atac.cancer.seu <- subset(RTs.atac.seu, subset = predicted.id %in% c(0,2,6,8,11))
PTs.atac.cancer.cells <- colnames(PTs.atac.cancer.seu)
RTs.atac.cancer.cells <- colnames(RTs.atac.cancer.seu)
rm(PTs.atac.cancer.seu)
rm(RTs.atac.cancer.seu)
TTs.atac.seu$barcodes <- colnames(TTs.atac.seu)
TTs.atac.cancer.seu <- subset(TTs.atac.seu, subset = barcodes %in% c(PTs.atac.cancer.cells, RTs.atac.cancer.cells))

# load processed rna object from Fig3.R
TTs.rna.cancer.seu <-LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223.h5seurat')
# Fig.4A for TTs co-embeding visualization 
TTs.transfer.anchors.sct <- FindTransferAnchors(
  reference = TTs.rna.cancer.seu,
  query = TTs.atac.cancer.seu,
  features = VariableFeatures(object = TTs.rna.cancer.seu),
  reference.assay = "integrated", query.assay = "ACTIVITY",
  reduction = 'cca',
  normalization.method ='SCT'
)

TTs.predicted.labels.sct <- TransferData(
  anchorset = TTs.transfer.anchors.sct,
  refdata = TTs.rna.cancer.seu$seurat_clusters,
  weight.reduction = TTs.atac.cancer.seu[['lsi']],
  dims = 2:30
)
TTs.atac.cancer.seu$predicted.id <- TTs.predicted.labels.sct$predicted.id

genes.use <- VariableFeatures(TTs.rna.cancer.seu)
refdata <- GetAssayData(TTs.rna.cancer.seu, assay = "integrated", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = TTs.transfer.anchors.sct, 
                           refdata = refdata, 
                           weight.reduction = TTs.atac.cancer.seu[["lsi"]],
                           dims = 2:30)

TTs.atac.cancer.seu[["integrated"]] <- imputation

TTs.rna.cancer.seu$Libtype <- 'RNA'
TTs.atac.cancer.seu$Libtype <- 'ATAC'
TTs.cancer.coembed.seu <- merge(x = TTs.rna.cancer.seu, y = TTs.atac.cancer.seu)

TTs.cancer.coembed.seu <- ScaleData(TTs.cancer.coembed.seu, features = genes.use, do.scale = FALSE)
TTs.cancer.coembed.seu <- RunPCA(TTs.cancer.coembed.seu, features = genes.use, verbose = FALSE)
TTs.cancer.coembed.seu <- RunUMAP(TTs.cancer.coembed.seu, dims = 1:30)
DimPlot(TTs.cancer.coembed.seu, group.by = c("Libtype", "seurat_clusters"), label=T, label.size = 5)

TTs.cancer.coembed.seu$seurat_clusters <- ifelse(is.na(TTs.cancer.coembed.seu$predicted.id), 
                                                 TTs.cancer.coembed.seu$integrated_snn_res.0.9, 
                                                 TTs.cancer.coembed.seu$predicted.id)

DimPlot(TTs.cancer.coembed.seu, group.by = c("Libtype", "seurat_clusters"), label=T, label.size = 5)
saveRDS(TTs.cancer.coembed.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.coembed.rds')

# TTs.cancer.coembed.seu only has cancer cells
Idents(TTs.cancer.coembed.seu) <- TTs.cancer.coembed.seu$seurat_clusters
TTs.cancer.coembed.seu <- RenameIdents(object = TTs.cancer.coembed.seu,
                                       "0" = "RT_CS1", "1" = "RT_CS1",
                                       "2" = "RT_CS1", "3" = "PT_CS2",
                                       "4" = "PT_CS4", 
                                       "5" = "RT_CS2",
                                       "6" = "RT_CS3","7" = "PT_CS3",
                                       "8" = "PT_CS1", 
                                       "9" = "PT_CS1",
                                       "10" = "PT_CS1","11"="PRT_CS1",
                                       "12" = "PT_CS5")
TTs.cancer.coembed.seu$cell.state <- Idents(TTs.cancer.coembed.seu)
DimPlot(TTs.cancer.coembed.seu, group.by = c("Libtype", "cell.state"),  
        order = rev(c("PT_CS1","PT_CS2","PT_CS3","PT_CS4","PT_CS5","PRT_CS1","RT_CS1","RT_CS2",
                      "RT_CS3")),label.size = 5,label=T,repel=T) + NoLegend()

TTs.cancer.atac.seu <- subset(TTs.cancer.coembed.seu, subset = Libtype %in% c("ATAC"))
TTs.cancer.atac.metadata <- TTs.cancer.atac.seu@meta.data
rm(TTs.cancer.atac.seu)
# get pure atac seu.obj and assign cell states to it
TTs.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.combined.rds')
TTs.atac.seu$barcodes <- rownames(TTs.atac.seu@meta.data)
TTs.cancer.atac.seu <- subset(TTs.atac.seu,subset = barcodes %in% TTs.cancer.atac.metadata$barcodes)
TTs.cancer.atac.seu$cell.state  <- plyr::mapvalues(
  x = TTs.cancer.atac.seu$barcodes, 
  from = as.character(TTs.cancer.atac.metadata$barcodes), 
  to = as.character(TTs.cancer.atac.metadata$cell.state)
)
DefaultAssay(TTs.cancer.atac.seu) <- "ATAC"

# Peak calling
TTs.peaks <- CallPeaks(
  object = TTs.cancer.atac.seu,
  group.by = "cell.state",
  macs2.path = '/hpc/apps/miniconda3/4.9.2/envs/macs2-2.2.6/bin/macs2'
)
saveRDS(TTs.peaks,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.peaks.rds')
# TTs.peaks <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.peaks.rds')
saveRDS(TTs.cancer.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.0615.rds')

############# Part 4: identify DAs #################
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)

# TTs.cancer.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.0615.rds')
# TTs.peaks <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.peaks.rds')
# TTs.rna.cancer.seu <-LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223.h5seurat')
TTs.rna.cancer.seu <- RenameIdents(object = TTs.rna.cancer.seu,
                                   "0" = "RT_CS1", "1" = "RT_CS1",
                                   "2" = "RT_CS1", "3" = "PT_CS2",
                                   "4" = "PT_CS4", 
                                   "5" = "RT_CS2",
                                   "6" = "RT_CS3","7" = "PT_CS3",
                                   "8" = "PT_CS1", 
                                   "9" = "PT_CS1",
                                   "10" = "PT_CS1","11"="PRT_CS1",
                                   "12" = "PT_CS5")

Idents(TTs.rna.cancer.seu) <- factor(x = Idents(TTs.rna.cancer.seu), 
                                     levels = c("PT_CS1","PT_CS2","PT_CS3",
                                                "PT_CS4","PT_CS5","PRT_CS1",
                                                "RT_CS1","RT_CS2","RT_CS3"))

# add peaks layer
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
TTs.peaks <- keepStandardChromosomes(TTs.peaks, pruning.mode = "coarse")
TTs.peaks <- subsetByOverlaps(x = TTs.peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(TTs.cancer.atac.seu),
  features = TTs.peaks,
  cells = colnames(TTs.cancer.atac.seu)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
PT1_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT1/fragments.tsv.gz")
PT2_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT2/fragments.tsv.gz")
PT3_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/PT5/fragments.tsv.gz")
RT1_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT3/fragments.tsv.gz")
RT2_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT4/fragments.tsv.gz")
RT3_fragpath <- CreateFragmentObject(path ="/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RT6/fragments.tsv.gz")
TTs.cancer.atac.seu[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = list(PT1_fragpath, PT2_fragpath, PT3_fragpath,
                   RT1_fragpath, RT2_fragpath, RT3_fragpath),
  annotation = annotation
)
saveRDS(TTs.cancer.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.0620.rds')
saveRDS(TTs.peaks,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.filt.peaks.rds')

pwm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = 9606, all_versions = FALSE)
)

# get genes from pwm
gene_symbols <- unlist(lapply(pwm, function(x) x@name))
gene_symbols_split <- unlist(strsplit(gene_symbols, "::"))

# add motif information
TTs.cancer.atac.seu <- AddMotifs(TTs.cancer.atac.seu, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

# find differential peaks
TTs.cancer.atac.meta <- TTs.cancer.atac.seu@meta.data
rt_cs1_bc <- TTs.cancer.atac.meta[TTs.cancer.atac.meta$cell.state=='RT_CS1',]$barcodes
rt_cs2_bc <- TTs.cancer.atac.meta[TTs.cancer.atac.meta$cell.state=='RT_CS2',]$barcodes
rt_cs3_bc <- TTs.cancer.atac.meta[TTs.cancer.atac.meta$cell.state=='RT_CS3',]$barcodes
pt_bc <- TTs.cancer.atac.meta[TTs.cancer.atac.meta$cell.state %in% c('PT_CS1',
                                                                     'PT_CS2',
                                                                     'PT_CS3',
                                                                     'PT_CS4',
                                                                     'PT_CS5'),]$barcodes
rt_cs1_dpeaks <- FindMarkers(
  object = TTs.cancer.atac.seu,
  ident.1 = rt_cs1_bc,
  ident.2 = pt_bc,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

rt_cs2_dpeaks <- FindMarkers(
  object = TTs.cancer.atac.seu,
  ident.1 = rt_cs2_bc,
  ident.2 = pt_bc,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

rt_cs3_dpeaks <- FindMarkers(
  object = TTs.cancer.atac.seu,
  ident.1 = rt_cs3_bc,
  ident.2 = pt_bc,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

rt_cs1_dpeaks$region <- rownames(rt_cs1_dpeaks)
rt_cs2_dpeaks$region <- rownames(rt_cs2_dpeaks)
rt_cs3_dpeaks$region <- rownames(rt_cs3_dpeaks)

sheetlist <- list("RT_CS1"=rt_cs1_dpeaks,"RT_CS2"=rt_cs2_dpeaks,"RT_CS3"=rt_cs3_dpeaks)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/CS_dpeaks_macs2.xlsx')
saveRDS(TTs.cancer.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.addmotif.0620.rds')


############# Part 5: assign genomic location for peaks and find DA, DEA #################
# this part is done by Fig4.py

############# Part 6: Plot panels in Fig.4 #################
# TTs.cancer.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.addmotif.0620.rds')
Idents(TTs.cancer.atac.seu) <- factor(x = Idents(TTs.cancer.atac.seu), 
                                      levels = c("PT_CS1","PT_CS2","PT_CS3",
                                                 "PT_CS4","PT_CS5","PRT_CS1",
                                                 "RT_CS1","RT_CS2","RT_CS3"))
TTs.peaks <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.atac.filt.peaks.rds')
TTs.peaks.annot <- read.csv('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs_peaks_filt_annot_df.csv',row.names = 1)

# Load dp_all from part5
rtcs1_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all.xlsx',sheet = "RT_CS1")
rtcs1_dp_all$...1<-NULL
rtcs2_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all.xlsx',sheet = "RT_CS2")
rtcs2_dp_all$...1<-NULL
rtcs3_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all.xlsx',sheet = "RT_CS3")
rtcs3_dp_all$...1<-NULL

# Load dp_de_all from part5
rtcs1_dp_de_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_all.xlsx',sheet = "RT_CS1")
rtcs1_dp_de_all$...1<-NULL
rtcs2_dp_de_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_all.xlsx',sheet = "RT_CS2")
rtcs2_dp_de_all$...1<-NULL
rtcs3_dp_de_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_all.xlsx',sheet = "RT_CS3")
rtcs3_dp_de_all$...1<-NULL

# Load dp_de_module from part5
rtcs1_dp_de_mod <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_module.xlsx',sheet = "RT_CS1")
rtcs1_dp_de_mod$...1<-NULL
rtcs2_dp_de_mod <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_module.xlsx',sheet = "RT_CS2")
rtcs2_dp_de_mod$...1<-NULL
rtcs3_dp_de_mod <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_De_module.xlsx',sheet = "RT_CS3")
rtcs3_dp_de_mod$...1<-NULL

cs_upsetplot <- function(rtcs1_df,rtcs2_df,rtcs3_df,labels){
  # dp_de upset
  distal_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Distal',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Distal',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Distal',]$region
  )
  distal_m = make_comb_mat(distal_lt)
  proximal_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Proximal',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Proximal',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Proximal',]$region
  )
  proximal_m = make_comb_mat(proximal_lt)
  promoter_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Promoter',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Promoter',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Promoter',]$region
  )
  promoter_m = make_comb_mat(promoter_lt)
  top_ha = HeatmapAnnotation(
    "Distal" = anno_barplot(comb_size(distal_m), 
                            gp = gpar(fill = "chocolate"), height = unit(3, "cm"),
                            axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    "Proximal" = anno_barplot(comb_size(proximal_m), 
                              gp = gpar(fill = "darkorange"), height = unit(3, "cm"),
                              axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    "Promoter" = anno_barplot(comb_size(promoter_m), 
                              gp = gpar(fill = "orange"), height = unit(3, "cm"),
                              axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0,
    annotation_name_gp= gpar(fontsize = 16,fontface="bold"))
  # the same for using m2 or m3
  ss = c("RT_CS1"=length(rtcs1_df$region),
         "RT_CS2"=length(rtcs2_df$region),
         "RT_CS3"=length(rtcs3_df$region))
  UpSet(distal_m, 
        top_annotation = top_ha,
        set_order = c("RT_CS1","RT_CS2","RT_CS3"),
        left_annotation = rowAnnotation(
          "Total Counts" = anno_barplot(-ss, 
                                        baseline = 0,
                                        axis_param = list(
                                          at = -1 * labels,
                                          labels = labels,
                                          labels_rot = 0),
                                        border = FALSE, 
                                        gp = gpar(fill = c("cornflowerblue","darkorchid1","deeppink")), 
                                        width = unit(4, "cm")
                                        
          ),
          set_name = anno_text(c("RT_CS1","RT_CS2","RT_CS3"), 
                               location = 0.5, 
                               just = "center",
                               gp = gpar(fontsize = 14,
                                         fontface="bold"),
                               width = max_text_width(set_name(distal_m)) + unit(4, "mm")),
          annotation_name_gp= gpar(fontsize = 14, fontface="bold")
        ), 
        right_annotation = NULL,
        show_row_names = FALSE
  )
}

cs_upsetplot(rtcs1_dp_all,rtcs2_dp_all,rtcs3_dp_all,c(0,5000,10000,15000))
cs_upsetplot(rtcs1_dp_de_all,rtcs2_dp_de_all,rtcs3_dp_de_all,c(0,1000,2000,15000))

# choosing background peaks
open.peaks <- AccessiblePeaks(TTs.cancer.atac.seu, idents = c("PT_CS1","PT_CS2","PT_CS3","PT_CS4","PT_CS5"))
# match the overall GC content in the peak set
meta.feature <- GetAssayData(TTs.cancer.atac.seu, assay = "peaks", slot = "meta.features")

findEnrichMotif <- function(rt_cs_dpeaks, lc.cut, meta.feature){
  rt_cs_dpeaks_upp <- rt_cs_dpeaks[rt_cs_dpeaks$avg_log2FC>=lc.cut,]$region
  rt_cs_dpeaks_dwp <- rt_cs_dpeaks[rt_cs_dpeaks$avg_log2FC<=-lc.cut,]$region
  if (length(grep("\\.", rt_cs_dpeaks_upp))>=1){
    # remove extension https://github.com/stuart-lab/signac/issues/1109
    rt_cs_dpeaks_upp <- rt_cs_dpeaks_upp[-grep("\\.", rt_cs_dpeaks_upp)]
  } else {
    rt_cs_dpeaks_upp <- rt_cs_dpeaks_upp
  }
  if (length(grep("\\.", rt_cs_dpeaks_dwp))>=1){
    # remove extension https://github.com/stuart-lab/signac/issues/1109
    rt_cs_dpeaks_dwp <- rt_cs_dpeaks_dwp[-grep("\\.", rt_cs_dpeaks_dwp)]
  } else {
    rt_cs_dpeaks_dwp <- rt_cs_dpeaks_dwp
  }
  rt_cs_peaks_upp.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[rt_cs_dpeaks_upp, ],
    n = 50000
  )
  rt_cs_peaks_dwp.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[rt_cs_dpeaks_dwp, ],
    n = 50000
  )
  rt_cs_upp_enriched.motifs <- FindMotifs(
    object = TTs.cancer.atac.seu,
    features = rt_cs_dpeaks_upp,
    background = rt_cs_peaks_upp.matched
  )
  rt_cs_dwp_enriched.motifs <- FindMotifs(
    object = TTs.cancer.atac.seu,
    features = rt_cs_dpeaks_dwp,
    background = rt_cs_peaks_dwp.matched
  )
  # only keep pval < 0.05, fold.enrichment > 1.5
  rt_cs_upp_enriched.motifs <- rt_cs_upp_enriched.motifs[(rt_cs_upp_enriched.motifs$p.adjust<0.05)
                                                         &(rt_cs_upp_enriched.motifs$fold.enrichment>1.5),]
  rt_cs_dwp_enriched.motifs <- rt_cs_dwp_enriched.motifs[(rt_cs_dwp_enriched.motifs$p.adjust<0.05)
                                                         &(rt_cs_dwp_enriched.motifs$fold.enrichment>1.5),]
  rt_cs_upp_enriched.motifs$type <- "UpA"
  rt_cs_dwp_enriched.motifs$type <- "DownA"
  # ordered by fold enrichment
  rt_cs_upp_enriched.motifs <- rt_cs_upp_enriched.motifs[order(rt_cs_upp_enriched.motifs$fold.enrichment,decreasing = T),]
  rt_cs_dwp_enriched.motifs <- rt_cs_dwp_enriched.motifs[order(rt_cs_dwp_enriched.motifs$fold.enrichment,decreasing = T),]
  rt_cs_enriched.motifs <- rbind(rt_cs_upp_enriched.motifs, rt_cs_dwp_enriched.motifs)
  return(rt_cs_enriched.motifs)
}

plotMotif <- function(TTs.cancer.atac.seu,enriche.region.motifs,title){
  p1<-MotifPlot(
    object = TTs.cancer.atac.seu,
    motifs = c(enriche.region.motifs[(enriche.region.motifs$type=='UpA')
                                     &(enriche.region.motifs$percent.observed>30),]$motif[1:5]),
    ncol=5
  )
  p2<-MotifPlot(
    object = TTs.cancer.atac.seu,
    motifs = c(enriche.region.motifs[(enriche.region.motifs$type=='DownA')
                                     &(enriche.region.motifs$percent.observed>30),]$motif[1:5]),
    ncol=5
  )
  print(p1 + ggtitle(title) + p2 + plot_layout(nrow = 2))
}

plotvenn_motif <- function(rtcs1_df,rtcs2_df,rtcs3_df,region){
  rt_cs1_enriche.region.motifs <- findEnrichMotif(rtcs1_df[rtcs1_df$genomeLoc_annot==region,],1,meta.feature)
  rt_cs2_enriche.region.motifs <- findEnrichMotif(rtcs2_df[rtcs2_df$genomeLoc_annot==region,],1,meta.feature)
  rt_cs3_enriche.region.motifs <- findEnrichMotif(rtcs3_df[rtcs3_df$genomeLoc_annot==region,],1,meta.feature)
  
  dp_motif_up_region <- list(
    RT_CS1 = rt_cs1_enriche.region.motifs[rt_cs1_enriche.region.motifs$type=="UpA",]$motif,
    RT_CS2 = rt_cs2_enriche.region.motifs[rt_cs2_enriche.region.motifs$type=="UpA",]$motif,
    RT_CS3 = rt_cs3_enriche.region.motifs[rt_cs3_enriche.region.motifs$type=="UpA",]$motif
  )
  
  dp_motif_dw_region <- list(
    RT_CS1 = rt_cs1_enriche.region.motifs[rt_cs1_enriche.region.motifs$type=="DownA",]$motif,
    RT_CS2 = rt_cs2_enriche.region.motifs[rt_cs1_enriche.region.motifs$type=="DownA",]$motif,
    RT_CS3 = rt_cs3_enriche.region.motifs[rt_cs1_enriche.region.motifs$type=="DownA",]$motif
  )
  
  p1<- ggVennDiagram(dp_motif_up_region, label_alpha = 0, edge_size=0,label_size = 10)
  p2<- ggVennDiagram(dp_motif_dw_region, label_alpha = 0, edge_size=0,label_size = 10)
  print(p1 + ggplot2::scale_fill_gradient(low="white",high = "red") + NoLegend() + ggtitle("UpA") +
          p2 + ggplot2::scale_fill_gradient(low="white",high = "red") + NoLegend() + ggtitle("DownA"))
  
  
  plotMotif(TTs.cancer.atac.seu,rt_cs1_enriche.region.motifs,"RT_CS1")
  plotMotif(TTs.cancer.atac.seu,rt_cs2_enriche.region.motifs,"RT_CS2")
  plotMotif(TTs.cancer.atac.seu,rt_cs3_enriche.region.motifs,"RT_CS3")
  motif_list <- list("RT_CS1"=rt_cs1_enriche.region.motifs,
                     "RT_CS2"=rt_cs2_enriche.region.motifs,
                     "RT_CS3"=rt_cs3_enriche.region.motifs)
  return(motif_list)
  
}
# all dp
dp_distal_motif_list <- plotvenn_motif(rtcs1_dp_all,rtcs2_dp_all,rtcs3_dp_all,'Distal')
dp_proximal_motif_list <- plotvenn_motif(rtcs1_dp_all,rtcs2_dp_all,rtcs3_dp_all,'Proximal')
dp_promoter_motif_list <- plotvenn_motif(rtcs1_dp_all,rtcs2_dp_all,rtcs3_dp_all,'Promoter')

write_xlsx(dp_distal_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_distal_motifs.xlsx')
write_xlsx(dp_proximal_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_proximal_motifs.xlsx')
write_xlsx(dp_promoter_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_promoter_motifs.xlsx')

# dp de
dp_de_distal_motif_list <- plotvenn_motif(rtcs1_dp_de_all,rtcs2_dp_de_all,rtcs3_dp_de_all,'Distal')
dp_de_proximal_motif_list <- plotvenn_motif(rtcs1_dp_de_all,rtcs2_dp_de_all,rtcs3_dp_de_all,'Proximal')
dp_de_promoter_motif_list <- plotvenn_motif(rtcs1_dp_de_all,rtcs2_dp_de_all,rtcs3_dp_de_all,'Promoter')
write_xlsx(dp_de_distal_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_de_distal_motifs.xlsx')
write_xlsx(dp_de_proximal_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_de_proximal_motifs.xlsx')
write_xlsx(dp_de_promoter_motif_list, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/dp_de_promoter_motifs.xlsx')

# motif heatmap 
## footprint
pwm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = 9606, all_versions = FALSE)
)

# get genes from pwm
gene_symbols <- unlist(lapply(pwm, function(x) x@name))
motif_name <- unlist(lapply(pwm, function(x) x@ID))
motif_gene <- data.frame(motif_name,gene_symbols)

dp_de_distal_merged_upa <- Reduce(function(x, y) full_join(x, y, by = "motif"), 
                                  lapply(dp_de_distal_motif_list, function(df) df[df$type=="UpA", c("motif", "fold.enrichment")])
)
rownames(dp_de_distal_merged_upa) <- dp_de_distal_merged_upa$motif
dp_de_distal_merged_upa$motif<- NULL
colnames(dp_de_distal_merged_upa) <- c("RT_CS1","RT_CS2","RT_CS3")
dp_de_distal_merged_upa[is.na(dp_de_distal_merged_upa)] <- 0

# dp_de_distal_merged_upa_z <- as.data.frame(t(scale(t(dp_de_distal_merged_upa))))
set.seed(34)
built_ht_c <- function(z_map,ht,cluster_name,motif_gene){
  ht_c <- as.data.frame(t(t(z_map[row_order(ht)[[cluster_name]],])))
  ht_c$sum <- rowSums(ht_c)
  ht_c$gene<- mapvalues(x = rownames(ht_c), from=motif_gene$motif_name,to=motif_gene$gene_symbols)
  ht_c <- ht_c[order(ht_c$sum,decreasing = T),]
  return(ht_c)
}
ddd_ht <- Heatmap(dp_de_distal_merged_upa, km=11, show_row_names = F, cluster_columns = F)
# saveRDS(ddd_ht, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/ddd_ht_062523.rds')
ddd_ht <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/ddd_ht_062523.rds')
ddd_ht <- draw(ddd_ht)
ddd_ht_c1 <- built_ht_c(dp_de_distal_merged_upa,ddd_ht,"1",motif_gene)
ddd_ht_c2 <- built_ht_c(dp_de_distal_merged_upa,ddd_ht,"2",motif_gene)
ddd_ht_c4 <- built_ht_c(dp_de_distal_merged_upa,ddd_ht,"4",motif_gene)
ddd_ht_c10 <- built_ht_c(dp_de_distal_merged_upa,ddd_ht,"10",motif_gene)
ddd_ht_c11 <- built_ht_c(dp_de_distal_merged_upa,ddd_ht,"11",motif_gene)

ddd_ht_c_list <- list("C1"=ddd_ht_c1,"C2"=ddd_ht_c2,"C4"=ddd_ht_c4,"C10"=ddd_ht_c10,
                      "C11"=ddd_ht_c11)
write_xlsx(ddd_ht_c_list,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/ddd_ht_c_list.xlsx')

MotifPlot(
  object = subset(TTs.cancer.atac.seu,subset=cell.state %in% c("RT_CS1")),
  motifs = c("HES5","ZNF85","KLF3","TFAP2C","PAX1"),
  ncol=5
)

MotifPlot(
  object = subset(TTs.cancer.atac.seu,subset=cell.state %in% c("RT_CS2")),
  motifs = c("ZBTB7B","MYBL1","TFDP1","ZBTB7C","ESR1"),
  ncol=5
)

MotifPlot(
  object = subset(TTs.cancer.atac.seu,subset=cell.state %in% c("RT_CS3")),
  motifs = c("SMAD5","HOXB2","HOXA2","NFKB1","NR2F1"),
  ncol=5
)

MotifPlot(
  object = subset(TTs.cancer.atac.seu,subset=cell.state %in% c("RT_CS1","RT_CS2","RT_CS3")),
  motifs = c("GLIS1","FOSB::JUNB","FOSL1::JUN","SOX12","HIF1A"),
  ncol=5
)

MotifPlot(
  object = subset(TTs.cancer.atac.seu,subset=cell.state %in% c("RT_CS1","RT_CS2","RT_CS3")),
  motifs = c("EGR3","FOSL2::JUNB","FOS","FOSL2::JUND","JUNB"),
  ncol=5
)

# add chromatin accessiblity score
rtcs1_dp_de_mod <- rtcs1_dp_de_mod[(rtcs1_dp_de_mod$sync_score>0),]
rtcs2_dp_de_mod <- rtcs2_dp_de_mod[(rtcs2_dp_de_mod$sync_score>0),]
rtcs3_dp_de_mod <- rtcs3_dp_de_mod[(rtcs3_dp_de_mod$sync_score>0),]

rtcs_intersection <- getVennOverlap(list("RT_CS1"=rtcs1_dp_de_mod$region,
                                         "RT_CS2"=rtcs2_dp_de_mod$region,
                                         "RT_CS3"=rtcs3_dp_de_mod$region))

rtcs1_spe_mod <- rtcs1_dp_de_mod[rtcs1_dp_de_mod$region %in% rtcs_intersection$RT_CS1,]
rtcs2_spe_mod <- rtcs2_dp_de_mod[rtcs2_dp_de_mod$region %in% rtcs_intersection$RT_CS2,]
rtcs3_spe_mod <- rtcs3_dp_de_mod[rtcs3_dp_de_mod$region %in% rtcs_intersection$RT_CS3,]
# add module scores to metadata
TTs.atac.modscore <- AddChromatinModule(object = TTs.cancer.atac.seu, 
                                        features = list("RT_CS1_module"=rtcs_intersection$RT_CS1,
                                                        "RT_CS2_module"=rtcs_intersection$RT_CS2,
                                                        "RT_CS3_module"=rtcs_intersection$RT_CS3
                                        ), 
                                        genome = BSgenome.Hsapiens.UCSC.hg38, 
                                        assay = "peaks")
metadata<- TTs.atac.modscore@meta.data
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=9)
p1 <- ggplot(TTs.atac.modscore@meta.data, 
             aes(x=cell.state, y=RT.CS1.module, fill=cell.state)) +
  geom_violin() +   
  geom_boxplot(width=0.2,outlier.shape = NA) +
  geom_hline(yintercept=mean(metadata[metadata$cell.state=="RT_CS1",]$RT.CS1.module), 
             linetype="dashed", color = color_list[[7]]) +
  
  theme_classic() + theme(legend.position="none",
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text = element_text(face="bold",size=14))
p2 <- ggplot(TTs.atac.modscore@meta.data, 
             aes(x=cell.state, y=RT.CS2.module, fill=cell.state)) +
  geom_violin() +   
  geom_boxplot(width=0.2,outlier.shape = NA) +
  geom_hline(yintercept=mean(metadata[metadata$cell.state=="RT_CS2",]$RT.CS2.module), 
             linetype="dashed", color = color_list[[8]]) +
  
  theme_classic() + theme(legend.position="none",
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text = element_text(face="bold",size=14))

p3 <- ggplot(TTs.atac.modscore@meta.data, 
             aes(x=cell.state, y=RT.CS3.module, fill=cell.state)) +
  geom_violin() +    
  geom_boxplot(width=0.2,outlier.shape = NA) + 
  geom_hline(yintercept=mean(metadata[metadata$cell.state=="RT_CS3",]$RT.CS3.module), 
             linetype="dashed", color = color_list[[9]]) +
  theme_classic() + theme(legend.position="none",
                          axis.text = element_text(face="bold",size=14))
p1+p2+p3+plot_layout(nrow=3)

TTs.rna.cancer.seu <-LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223.h5seurat')
TTs.rna.cancer.seu <- RenameIdents(object = TTs.rna.cancer.seu,
                                   "0" = "RT_CS1", "1" = "RT_CS1",
                                   "2" = "RT_CS1", "3" = "PT_CS2",
                                   "4" = "PT_CS4", 
                                   "5" = "RT_CS2",
                                   "6" = "RT_CS3","7" = "PT_CS3",
                                   "8" = "PT_CS1", 
                                   "9" = "PT_CS1",
                                   "10" = "PT_CS1","11"="PRT_CS1",
                                   "12" = "PT_CS5")

Idents(TTs.rna.cancer.seu) <- factor(x = Idents(TTs.rna.cancer.seu), 
                                     levels = c("PT_CS1","PT_CS2","PT_CS3",
                                                "PT_CS4","PT_CS5","PRT_CS1",
                                                "RT_CS1","RT_CS2","RT_CS3"))
plotall <- function(genesymbol,region,atac.seu,rna.seu,extend.up,extend.dw){
  cov_plot<-CoveragePlot(
    object = atac.seu,
    region = region,
    assay = "ATAC",
    annotation = FALSE,
    peaks = FALSE,
    extend.upstream = extend.up,
    extend.downstream = extend.dw
  )
  
  annot_plot <- AnnotationPlot(
    object = atac.seu,
    region = region,
    extend.upstream = extend.up,
    extend.downstream = extend.dw
  )
  
  peak_plot <- PeakPlot(
    object = atac.seu,
    region = region,
    extend.upstream = extend.up,
    extend.downstream = extend.dw
  )
  
  expr_plot <- ExpressionPlot(
    object = rna.seu,
    features = genesymbol,
    assay = "SCT"
  )
  
  print(CombineTracks(
    plotlist = list(cov_plot, peak_plot, annot_plot),
    expression.plot = expr_plot,
    heights = c(10, 1, 2),
    widths = c(10, 1)
  ))
}

rtcs1_mod_select <- rtcs1_dp_de_mod[rtcs1_dp_de_mod$region %in% rtcs_intersection$RT_CS1,]
rtcs2_mod_select <- rtcs2_dp_de_mod[rtcs2_dp_de_mod$region %in% rtcs_intersection$RT_CS2,]
rtcs3_mod_select <- rtcs3_dp_de_mod[rtcs3_dp_de_mod$region %in% rtcs_intersection$RT_CS3,]

plotall("BMP7","chr20-57168753-57266641",TTs.cancer.atac.seu,
        TTs.rna.cancer.seu,200000,20000)

plotall("ALCAM","chr3-105366909-105576900",TTs.cancer.atac.seu,
        TTs.rna.cancer.seu,200000,20000)

CoveragePlot(
  object = TTs.cancer.atac.seu,
  region = "chr20-57124811-57125206",
  assay = "ATAC",
  extend.upstream = 1000,
  extend.downstream = 1000
)

CoveragePlot(
  object = TTs.cancer.atac.seu,
  region = "chr3-105386817-105387849",
  assay = "ATAC",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# TTs.cancer.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.motif.0615.rds')
# pwm <- getMatrixSet(
#   x = JASPAR2022,
#   opts = list(species = 9606, all_versions = FALSE)
# )
# gene_symbols <- unlist(lapply(pwm, function(x) x@name))
# motif_name <- unlist(lapply(pwm, function(x) x@ID))
# motif_gene <- data.frame(motif_name,gene_symbols)
# plot footprint
TTs.cancer.atac.seu <- Footprint(
  object = TTs.cancer.atac.seu,
  motif.name = c("MA0476.1","MA1951.1", #FOS
                 "MA1135.1","MA1136.1", # FOSB::JUNB
                 "MA0821.2", # HES5
                 "MA0524.2", "MA0815.1", "MA0814.2", #TFAP2C
                 "MA0112.3", # ESR1
                 "MA1557.1", # SMAD5
                 "MA0017.2", "MA1537.1", "MA1538.1"  # NR2F1
  ),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(TTs.cancer.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.motif.footprint.0615.rds')

Idents(TTs.cancer.atac.seu) <- factor(x = Idents(TTs.cancer.atac.seu), 
                                      levels = c("PT_CS1","PT_CS2","PT_CS3",
                                                 "PT_CS4","PT_CS5","PRT_CS1",
                                                 "RT_CS1","RT_CS2","RT_CS3"))

PlotFootprint(TTs.cancer.atac.seu, features = c("MA0476.1"),
              label=F)
PlotFootprint(TTs.cancer.atac.seu, features = c("MA1135.1"),
              label=F)


