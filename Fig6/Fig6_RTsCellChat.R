library(optparse)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(purrr)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ComplexHeatmap)
source("styling.R")

plot_contribution_celltype <- function(cell_type = "NB", signif = 0.05,prob_cut=0.01,ypos=1.68) {
  # prepare data
  vis_data <- 
    sig_data %>%
    filter(
      pathway_name %in% cellchat@netP$pathways, 
      source == {{cell_type}}, 
      target != {{cell_type}},
      pval <= signif
    ) %>%
    group_by(interaction = interaction_name_2) %>%
    summarise(
      prob.norm = sum(prob.norm),
      count = n(),
      pathway_name = first(pathway_name)
    ) %>%
    ungroup() %>%
    mutate(prob.avg = prob.norm / count) %>% 
    mutate(interaction = fct_reorder(interaction, prob.norm))
    
  #print(vis_data)
  # export source data
  #vis_data_cut <- vis_data[vis_data$prob.norm>=prob_cut, ]
  #vis_data_cut %>% 
  #vis_data %>% 
  #  select(interaction, pathway_name, total_outgoing_comm_score = prob.norm) %>% 
  #  arrange(desc(total_outgoing_comm_score)) %>% 
  #  save_table("source_data/figure_3b", "Figure 3b")

  vis_data_cut <- vis_data %>% 
    filter(prob.norm>=prob_cut) %>%
    arrange(desc(prob.norm))


  # make plot
  ggplot(vis_data_cut, aes(interaction, prob.norm)) +
    #geom_col(fill = "#ffb03e", width = 0.85) +
    geom_col(fill = "#3e8dff", width = 0.85) +
    geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +
    geom_text(
      aes(y = ypos, label = pathway_name),
      size = BASE_TEXT_SIZE_MM,
      fontface = "italic"
    ) +
    xlab("Ligand-receptor pair") +
    scale_y_continuous(
      paste0("Total outgoing communication score from ",cell_type),
      expand = expansion(0)
    ) +
    coord_flip() +
    theme_nb(grid = FALSE) +
    theme(
      axis.title.x = element_text(hjust = 1),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        size = BASE_LINE_SIZE
      )
    )
}

plot_n_interactions <- function(abbreviations) {
  # prepare data
  #rc_order <- names(CELL_TYPE_ABBREVIATIONS)
  rc_order <- names(abbreviations)
  mat <- cellchat@net$count[rc_order, rc_order]
  total_target <- colSums(mat)
  total_source <- rowSums(mat)
  bar_ylim <- c(0, max(total_source, total_target))
  
  # export source data
  mat %>% 
    as_tibble(rownames = "source") %>% 
    save_table("source_data/figure_3a", "Figure 3a")
  
  # make plot
  Heatmap(
    mat,
    col = RColorBrewer::brewer.pal(9, "PuBu"),
#RColorBrewer::brewer.pal(9, "YlOrBr"),
    name = "Number of\ninteractions",
    
    cluster_rows = FALSE, 
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    row_names_side = "left",
    row_title = "Source (ligand)", 
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    cluster_columns = FALSE, 
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = "Target (receptor)",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    width = unit(20, "mm"),
    height = unit(20, "mm"),
    
    top_annotation = HeatmapAnnotation(
      count_bar = anno_barplot(
        total_target,
        ylim = bar_ylim,
        border = FALSE, 
        axis = FALSE,
        #gp = gpar(fill = "gray70", col = "gray70")
        gp = gpar(fill = "#9fc5e8", col = "#9fc5e8")
      ),
      count_text = anno_text(
        total_target, 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ), 
      simple_anno_size_adjust = TRUE,
      show_annotation_name = FALSE
    ), 
    
    right_annotation = rowAnnotation(
      count_text = anno_text(
        total_source, 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ), 
      count_bar = anno_barplot(
        total_source,
        ylim = bar_ylim,
        border = FALSE, 
        axis = FALSE,
        #gp = gpar(fill = "gray70", col = "gray70")
        gp = gpar(fill = "#9fc5e8", col = "#9fc5e8")
      ),
      simple_anno_size_adjust = TRUE,
      show_annotation_name = FALSE
    ), 
    
    heatmap_legend_param = list(
      at = c(min(mat), max(mat)),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      title_position = "topleft", 
      border = NA, 
      legend_height = unit(20, "mm"), 
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      grid_width = unit(2, "mm")
    )
  )
}

option_list <- list(
  make_option(c("-f", "--flag"), type = "character", default = FALSE,
              action = "store", help = "Input the bedgraph file !"
  )
)

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(2, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(2, "pt"),
  TITLE_PADDING = unit(3, "mm")
)


opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))
dir="../"
readfile<-paste0(dir,opt$flag,".h5seurat")
RTsDict<-c(
           "0"="Breast cancer cells",
           "1"="Myeloid",
           "2"="Breast cancer cells",
           "3"="Luminal",
           "4"="Stromal",
           "5"="Stromal",
           "6"="Breast cancer cells",
           "7"="Stromal",
           "8"="Breast cancer cells",
           "9"="Stromal",
           "10"="Luminal",
           "11"="Breast cancer cells",
           "12"="NK"
)
#RTsDict<-c(
#           "0"="Breast cancer cells",
#           "1"="Myeloid cells",
#           "2"="Breast cancer cells",
#           "3"="Luminal cells",
#           "4"="Stromal cells",
#           "5"="Stromal cells",
#           "6"="Breast cancer cells",
#           "7"="Stromal cells",
#           "8"="Breast cancer cells",
#           "9"="Stromal cells",
#           "10"="Luminal cells",
#           "11"="Breast cancer cells",
#           "12"="NK cells"
#)

CDict<-c(
           "0"="C0",
           "1"="C1",
           "2"="C2",
           "3"="C3",
           "4"="C4",
           "5"="C5",
           "6"="C6",
           "7"="C7",
           "8"="C8",
           "9"="C9",
           "10"="C10",
           "11"="C11",
           "12"="C12",
           "13"="C13",
           "14"="C14"
)

PTsDict<-c(
           "0"="Breast cancer cells",
           "1"="Myeloid",
           "2"="Basal",
           "3"="Adipocyte",
           "4"="Breast cancer cells",
           "5"="Breast cancer cells",
           "6"="Luminal",
           "7"="Breast cancer cells",
           "8"="Breast cancer cells",
           "9"="Breast cancer cells",
           "10"="Neutrophils",
           "11"="Breast cancer cells",
           "12"="Breast cancer cells",
           "13"="Fibroblast",
           "14"="Breast cancer cells"
)

#PTsDict<-c(
#           "0"="Breast cancer cells",
#           "1"="Myeloid cells",
#           "2"="Basal cells",
#           "3"="Adipocyte-like cells",
#           "4"="Breast cancer cells",
#           "5"="Breast cancer cells",
#           "6"="Luminal cells",
#           "7"="Breast cancer cells",
#           "8"="Breast cancer cells",
#           "9"="Breast cancer cells",
#           "10"="Neutrophils",
#           "11"="Breast cancer cells",
#           "12"="Breast cancer cells",
#           "13"="Fibroblast-like cells",
#           "14"="Breast cancer cells"
#)

if(0)
{
ExprMatrix<-LoadH5Seurat(readfile,assays = "RNA")
############################################################
new_cluster_ids<-NULL
if(opt$flag=="PTs"){
new_cluster_ids <- levels(ExprMatrix) %>% map_chr(\(x) PTsDict[x])
}else if(opt$flag=="RTs"){
new_cluster_ids <- levels(ExprMatrix) %>% map_chr(\(x) RTsDict[x])
}else{
new_cluster_ids <- levels(ExprMatrix) %>% map_chr(\(x) CDict[x])
}
############################################################
names(new_cluster_ids) <- levels(ExprMatrix)
ExprMatrix <- RenameIdents(ExprMatrix, new_cluster_ids)

BCaCellStates<-read.csv(paste0(dir,"cellstates_barcodes.csv"))#"../cellstates_barcodes.csv")
newIdent<-map2_chr(.x=ExprMatrix@active.ident,.y=names(ExprMatrix@active.ident),
    .f=function(x,y){
       if(x=="Breast cancer cells")
       {
           #paste0("BCa_",BCaCellStates[BCaCellStates$barcode==y,]$cell_states)
           #paste0(BCaCellStates[BCaCellStates$barcode==y,]$cell_states,"      ")
           paste0(BCaCellStates[BCaCellStates$barcode==y,]$cell_states)
       } else {
           as.character(x)
       }
    })

Idents(ExprMatrix)<-factor(newIdent)
#ExprMatrix<-subset(ExprMatrix, idents = Idents(ExprMatrix))
#ExprMatrix<-subset(ExprMatrix, idents = setdiff(Idents(ExprMatrix),c("PRT_CS1      ","RT_CS1      ","RT_CS2      ")))
if(opt$flag=="PTs")
{
    ExprMatrix<-subset(ExprMatrix, idents = setdiff(Idents(ExprMatrix),c("PRT_CS1","RT_CS1","RT_CS2")))
}
if(opt$flag=="RTs")
{
    ExprMatrix<-subset(ExprMatrix, idents = setdiff(Idents(ExprMatrix),c("PRT_CS1")))
}
#ExprMatrix<-subset(ExprMatrix, idents = Idents(ExprMatrix))
cellchat <- createCellChat(object = ExprMatrix, group.by = "ident", assay = "RNA")
print(Idents(ExprMatrix))

#ExprMatrix@meta.data$annotated_clusters<-factor(newIdent)
#data.input <- ExprMatrix@assays$RNA$data
#meta<-data.frame(labels=ExprMatrix@meta.data$annotated_clusters)
#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#cellchat <- addMeta(cellchat, meta = meta)
#cellchat <- setIdent(cellchat, ident.use = "labels")

#CellChatDB <- CellChatDB.human
#CellChatDB.use <- subsetDB(CellChatDB, search = "Immune Signaling")
#CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.human #.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,seed.use=1,population.size = TRUE) #needs some time
#cellchat <- computeCommunProb(cellchat,population.size = TRUE) #needs some time
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

groupSize <- as.numeric(table(cellchat@idents))

#df.net <- subsetCommunication(cellchat)
#write.csv(df.net, paste0("df-",opt$flag,".csv"), row.names=FALSE)

df_net <-
  cellchat %>% 
  subsetCommunication(thresh = NA) %>% 
  as_tibble() %>% 
  group_by(interaction_name_2) %>% 
  mutate(prob.norm = prob / max(prob)) %>% 
  ungroup()
cellchat %>% saveRDS(paste0("cellchat_",opt$flag,".rds"))
df_net %>% saveRDS(paste0("signaling_",opt$flag,".rds"))
#saveRDS(cellchat, file = )#"CommunProb_cellchat_RTs.rds")
#saveRDS(df.net, file = )#"CommunProb_cellchat_RTs.rds")
saveRDS(ExprMatrix, file = paste0("ExprMatrix_",opt$flag,".rds"))
sig_data <- readRDS(paste0("signaling_",opt$flag,".rds"))

}else{

future::plan("multicore", workers = 4)
cellchat<-readRDS(paste0("cellchat_",opt$flag,".rds"))
sig_data <- readRDS(paste0("signaling_",opt$flag,".rds"))
write.csv(sig_data, paste0("df-",opt$flag,".csv"), row.names=FALSE)
groupSize <- as.numeric(table(cellchat@idents))

}

pathways.show.all <- cellchat@netP$pathways
TypeOfCells<-levels(cellchat@idents)

sink(file=paste0("./Pathway_",opt$flag,".txt"))
print(pathways.show.all)
sink(file=NULL)
sink(file=paste0("./CellTypes_",opt$flag,".txt"))
print(TypeOfCells)
sink(file=NULL)

PTCancerStates=c("PT_CS1","PT_CS2","PT_CS3","PT_CS4","PT_CS5")
RTCancerStates=c("RT_CS1","RT_CS2","RT_CS3")
if(opt$flag=="PTs") {CancerStates=PTCancerStates}
if(opt$flag=="RTs") {CancerStates=RTCancerStates}
#CancerStates=(stop("Flag should be PTs or RTs") if (opt$flag!="PTs") and (opt$flag!="RTs") else (PTCancerStates if opt$flag=="PTs" else RTCancerStates))
#print(vertex.receiver)
vertex.receiver=match(CancerStates,TypeOfCells)


if(opt$flag=="PTs"){
    (p <- plot_n_interactions(CELL_TYPE_ABBREVIATIONS_PTs))
}
if(opt$flag=="RTs"){
    (p <- plot_n_interactions(CELL_TYPE_ABBREVIATIONS_RTs))
}
ggsave_publication(paste0("heatmap_",opt$flag), plot = p, width = 8, height = 7,type="png")

if(opt$flag=="PTs"){
    prob_cut=c(2.6,2.6,0.6,1.7,0.3)
    ypos=c(1.68,2.0,1.0,1.68,0.4)
}
if(opt$flag=="RTs"){
    prob_cut=c(1.68,0.7,0.7)
    ypos=c(1.45,0.7,0.7)
}
for (i in 1:length(CancerStates)) {
    plot_contribution_celltype(cell_type = CancerStates[i],signif=5e-8,prob_cut=prob_cut[i],ypos=ypos[i])
    ggsave_publication(paste0(CancerStates[i]), width = 5.5, height = 5,type="png")
}

if(1)
{
#netVisual_heatmap(cellchat, signaling = NULL, color.heatmap = "Reds",font.size=20,font.size.title=24)
for (i in 1:length(pathways.show.all)) {
  #options(repr.plot.width = 8, repr.plot.height =16)
  g2<-plotGeneExpression(cellchat,signaling=pathways.show.all[i], enriched.only = TRUE,angle.x=0)
  g2<-g2+plot_annotation(paste(pathways.show.all[i],"signaling pathway"),theme=theme(plot.title=element_text(hjust=0.5,size=24)))
  #g2<- g2 + ggtitle(paste("Expression of enriched genes in ",pathways.show.all[i],"signaling"))
  ggsave(filename=paste0(pathways.show.all[i], "_GeneExpression.pdf"), plot=g2, width = 12, height = 10, units = 'in', dpi = 300)

  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  #g <- netVisual_bubble(cellchat, sources.use = 5, targets.use=c(1:7), remove.isolate = FALSE,signaling=pathways.show.all[i],font.size=15,font.size.title=16,vjust.x=0,hjust.x=0)
  #ggsave(filename=paste0(pathways.show.all[i], "_Bubble.pdf"), plot=g, width = 9, height = 6, units = 'in', dpi = 300)

  {

  netVisual(cellchat, signaling = pathways.show.all[i],vertex.receiver = vertex.receiver, layout = "hierarchy",out.format="svg",vertex.weight=NULL,vertex.label.cex = 0.8,title.space = 6,label.dist=4,space.v=1.5,space.h = 1.6)

  }

  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 14,font.size.title=16)
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  #gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  #ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 9, height = 6, units = 'in', dpi = 300)
}

}


