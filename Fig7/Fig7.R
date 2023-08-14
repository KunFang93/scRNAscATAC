library(InstaPrism)
library(SeuratDisk)
library(Seurat)
library(readxl)
library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)

PTs.seu <- LoadH5Seurat('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/PTs.h5seurat')
RTs.seu <- LoadH5Seurat('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/RTs.h5seurat')
core_signature <- read.csv('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig3/core_signatures.csv')

# assign cell types
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

# assign cell states
PTs.seu$barcodes <- rownames(PTs.seu@meta.data)
RTs.seu$barcodes <- rownames(RTs.seu@meta.data)
cell.states <- read.csv('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/annotation/barcodes_cs.csv',row.names = 1)
PTs.seu$cell.states<- ifelse(PTs.seu$cell.annot=='Breast cancer cells',
                             plyr::mapvalues(
                                  x = PTs.seu$barcodes, 
                                  from = as.character(cell.states$barcodes), 
                                  to = as.character(cell.states$cell.state)
                                ),
                             as.character(PTs.seu$cell.annot)
                             )
RTs.seu$cell.states<- ifelse(RTs.seu$cell.annot=='Breast cancer cells',
                             plyr::mapvalues(
                               x = RTs.seu$barcodes, 
                               from = as.character(cell.states$barcodes), 
                               to = as.character(cell.states$cell.state)
                             ),
                             as.character(RTs.seu$cell.annot)
)                          

# seurat to sc
PTs_exprs <- as.matrix(PTs.seu[['RNA']]@counts) %>% apply(2,function(x)((x/sum(x))*1e+09)) 
RTs_exprs <- as.matrix(RTs.seu[['RNA']]@counts) %>% apply(2,function(x)((x/sum(x))*1e+09)) 

PT_refPhi_obj = refPrepare(PTs_exprs, 
                           PTs.seu$cell.states, 
                           PTs.seu$cell.states)
RT_refPhi_obj = refPrepare(RTs_exprs, 
                           RTs.seu$cell.states, 
                           RTs.seu$cell.states)
# bulk data
construct_exprs <- function(exprs_cohort, pdata_cohort,PT_refPhi_obj, RT_refPhi_obj){
  # instaprism
  print('PT')
  print("prism")
  PT.InstaPrism.res = InstaPrism(input_type = 'refPhi', 
                                 bulk_Expr = exprs_cohort, 
                                 refPhi = PT_refPhi_obj,
                                 outlier.cut=0.01,
                                 n.iter=1000)
  print("reconstruct")
  PTCSs_list <- list()
  PTCSs <- c('PT_CS1','PT_CS2','PT_CS3','PT_CS4','PT_CS5')
  for (i in seq_along(PTCSs)){
    PTCSs_list[[i]] <-reconstruct_Z_ct_initial(InstaPrism_obj = PT.InstaPrism.res,
                                               cell.type.of.interest = PTCSs[i])
  }
  
  # Calculate the sum of all data frames
  PTCS_expr_sum <- Reduce("+", PTCSs_list)
  # Calculate the average
  PTCS_expr_avg <- PTCS_expr_sum / length(PTCSs_list)
  
  print("RTs")
  print("prism")
  RT.InstaPrism.res = InstaPrism(input_type = 'refPhi', 
                                 bulk_Expr = exprs_cohort, 
                                 refPhi = RT_refPhi_obj,
                                 outlier.cut=0.05,
                                 n.iter=1000)
  print('reconstruct')
  RTCSs_list <- list()
  RTCSs <- c('RT_CS1','RT_CS2','RT_CS3')
  for (i in seq_along(RTCSs)){
    RTCSs_list[[i]] <- reconstruct_Z_ct_initial(InstaPrism_obj = RT.InstaPrism.res,
                                                        cell.type.of.interest = RTCSs[i])
  }
  # Calculate the sum of all data frames
  RTCS_expr_sum <- Reduce("+", RTCSs_list)
  # Calculate the average
  RTCS_expr_avg <- RTCS_expr_sum / length(RTCSs_list)
  
  return(list('PT'=PTCS_expr_avg,"RT"=RTCS_expr_avg,'PT.res'= PT.InstaPrism.res,
              'RT.res'=RT.InstaPrism.res))
}

plotheatmap <- function(pdata_cohort,PTCS_expr_avg,RTCS_expr_avg,sortmean=F){
  relapse_samples <- pdata_cohort[pdata_cohort$rfs==1,]$geo_accession
  relapsefree_samples <- pdata_cohort[pdata_cohort$rfs==0,]$geo_accession
  PTCS_expr_avg_csig <- PTCS_expr_avg[rownames(PTCS_expr_avg) %in% core_signature$gene,relapsefree_samples]
  RTCS_expr_avg_csig <- RTCS_expr_avg[rownames(RTCS_expr_avg) %in% core_signature$gene,relapse_samples]
  PTCS_expr_avg_csig <- as.data.frame(PTCS_expr_avg_csig)
  RTCS_expr_avg_csig <- as.data.frame(RTCS_expr_avg_csig)
  PTCS_expr_avg_csig$gene <- rownames(PTCS_expr_avg_csig)
  RTCS_expr_avg_csig$gene <- rownames(RTCS_expr_avg_csig)
  RTPTCS1 <- merge(RTCS_expr_avg_csig,PTCS_expr_avg_csig,by="gene",all=TRUE)
  rownames(RTPTCS1) <- RTPTCS1$gene
  RTPTCS1$gene <- NULL

  # find the smallest non-NA value in the data frame
  RTPTCS1[, -1] <- t(apply(RTPTCS1[, -1], 1, function(x) replace(x, is.na(x), min(x, na.rm = TRUE))))

  # RTPTCS1 <- cbind(RTCS_expr_avg_csig,PTCS_expr_avg_csig)
  RTPTCS1_zscore <- t(scale(t(RTPTCS1)))
  RTPTCS1_zscore <- RTPTCS1_zscore[,c(relapse_samples,relapsefree_samples)]
  # plot
  relapse_mean = rowMeans(RTPTCS1_zscore[,relapse_samples])
  relapsefree_mean = rowMeans(RTPTCS1_zscore[,relapsefree_samples])
  mean_mat = cbind(relapse_mean,relapsefree_mean)

  if (sortmean){
    mean_mat = mean_mat[order(mean_mat[,"relapse_mean"]-mean_mat[,"relapsefree_mean"],decreasing = T),]
    RTPTCS1_zscore <- RTPTCS1_zscore[rownames(mean_mat),]
  }
  ta <- HeatmapAnnotation(sum = anno_points(colSums(RTPTCS1_zscore)))
  ht_list = Heatmap(RTPTCS1_zscore, name = "mat",
                    column_split = c(rep("Relapse", length(relapse_samples)),
                                     rep("Relapse-Free", length(relapsefree_samples))),
                    cluster_rows = F,
                    cluster_columns = F,
                    show_column_names = F,
                    show_row_names = T,
                    column_gap = unit(5, "mm"),
                    top_annotation = ta,
                    row_names_gp = gpar(fontsize = 6),
                    col = colorRamp2(c(-4, 0, 4), c("green", "white", "red"))
  ) +
    Heatmap(mean_mat, name = "relapse mean", row_names_gp = gpar(fontsize = 6),
            cluster_columns = F, width = unit(15, "mm"))
  draw(ht_list)
}


# ChanrionGeneExpTamBca
exprs_cohort1 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/ChanrionGeneExpTamBca/ChanrionGeneExpTamBca.xlsx',
                            sheet = "Expression_matrix") 
pdata_cohort1 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/ChanrionGeneExpTamBca/ChanrionGeneExpTamBca.xlsx',
                            sheet = "Process_Pdata")

exprs_cohort1_combined <- as.data.table(exprs_cohort1)[, lapply(.SD, mean, na.rm = TRUE), by = Gene_Symbol]
exprs_cohort1_combined <- as.data.frame(exprs_cohort1_combined)
rownames(exprs_cohort1_combined) <- exprs_cohort1_combined[,c('Gene_Symbol')]
exprs_cohort1_combined$Gene_Symbol <- NULL
exprs_cohort1_combined <- as.matrix(exprs_cohort1_combined)

cohort1_avg <- construct_exprs(exprs_cohort1_combined, pdata_cohort1,
                               PT_refPhi_obj, RT_refPhi_obj)
plotheatmap(pdata_cohort1,cohort1_avg$PT,cohort1_avg$RT,sortmean = T)

build_combined_exprs <- function(exprs_cohort,genecol="symbol"){
  # Check if any column other than 'symbol' is non-numeric and convert them to numeric
  exprs_cohort <- as.data.table(exprs_cohort)
  numeric_columns <- setdiff(names(exprs_cohort), genecol)
  for(col in numeric_columns) {
    if(!is.numeric(exprs_cohort[[col]])) {
      warning(paste("Converting non-numeric column:", col))
      exprs_cohort[[col]] <- as.numeric(exprs_cohort[[col]])
    }
  }
  
  # Group by 'symbol' and calculate the geometric mean of numeric columns, handling NAs and negative values
  exprs_cohort_combined <- exprs_cohort[, lapply(.SD, mean, na.rm = TRUE), .SDcols = numeric_columns, by = genecol]
  exprs_cohort_combined <- as.data.frame(exprs_cohort_combined)
  rownames(exprs_cohort_combined) <- exprs_cohort_combined[[genecol]]
  exprs_cohort_combined[[genecol]] <- NULL
  exprs_cohort_combined <- as.matrix(exprs_cohort_combined)
  return(exprs_cohort_combined)
}


# CreightonSignatureClinical
exprs_cohort2 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/CreightonSignatureClinical/CreightonSignatureClinical.xlsx',
                            sheet = "Expression_matrix") 
pdata_cohort2 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/CreightonSignatureClinical/CreightonSignatureClinical.xlsx',
                            sheet = "Process_Pdata")

exprs_cohort2_combined <- build_combined_exprs(exprs_cohort2,"symbol")
cohort2_avg <- construct_exprs(exprs_cohort2_combined, pdata_cohort2,
                               PT_refPhi_obj, RT_refPhi_obj)
plotheatmap(pdata_cohort2,cohort2_avg$PT,cohort2_avg$RT,sortmean = T)


# LoiGeneProfilingBCaTam
exprs_cohort3 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/LoiGeneProfilingBCaTam/LoiGeneProfilingBCaTam.xlsx',
                            sheet = "Expression_matrix") 
pdata_cohort3 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/LoiGeneProfilingBCaTam/LoiGeneProfilingBCaTam.xlsx',
                            sheet = "Process_Pdata")
exprs_cohort3_combined <- build_combined_exprs(exprs_cohort3,"symbol")
pdata_cohort3 <- pdata_cohort3 %>% drop_na(rfs)
pdata_cohort3$geo_accession <- pdata_cohort3$samplename
cohort3_avg <- construct_exprs(exprs_cohort3_combined, pdata_cohort3,
                               PT_refPhi_obj, RT_refPhi_obj)
plotheatmap(pdata_cohort3,cohort3_avg$PT,cohort3_avg$RT,sortmean = T)

# NagallaInteractionImmunitySubtypeBRCA
exprs_cohort4 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/NagallaInteractionImmunitySubtypeBRCA/NagallaInteractionImmunitySubtypeBRCA.xlsx',
                            sheet = "Expression_matrix") 
pdata_cohort4 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/NagallaInteractionImmunitySubtypeBRCA/NagallaInteractionImmunitySubtypeBRCA.xlsx',
                            sheet = "Process_Pdata")
exprs_cohort4_combined <- build_combined_exprs(exprs_cohort5,"Gene Symbol")
pdata_cohort4 <- pdata_cohort4 %>% drop_na(dfs)
pdata_cohort4$rfs <- pdata_cohort4$dfs
cohort4_avg <- construct_exprs(exprs_cohort4_combined, pdata_cohort4,
                               PT_refPhi_obj, RT_refPhi_obj)
plotheatmap(pdata_cohort4,cohort4_avg$PT,cohort4_avg$RT,sortmean = T)


# ChinAberrationBRcAPathophysiologies
exprs_cohort5 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/ChinAberrationBRcAPathophysiologies/ChinAberrationBRcAPathophysiologies.xlsx',
                            sheet = "Expression_matrix") 
pdata_cohort5 <- read_excel('/Users/kfang/Documents/lab/Jin_lab/BRCA_cohort/ChinAberrationBRcAPathophysiologies/ChinAberrationBRcAPathophysiologies.xlsx',
                            sheet = "Process_Pdata")
exprs_cohort5<-exprs_cohort5[ , colSums(is.na(exprs_cohort5)) == 0]
exprs_cohort5_combined <- build_combined_exprs(exprs_cohort5,"Gene Symbol")
pdata_cohort5 <- pdata_cohort5 %>% drop_na(rfs)

pdata_cohort5 <- pdata_cohort5[pdata_cohort5$geo_accession %in% colnames(exprs_cohort5),]
cohort5_avg <- construct_exprs(exprs_cohort5_combined, pdata_cohort5,
                               PT_refPhi_obj, RT_refPhi_obj)
plotheatmap(pdata_cohort5,cohort5_avg$PT,cohort5_avg$RT,sortmean = T)

