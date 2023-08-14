library(CellChat)
library(writexl)

cellchat.rt <- readRDS('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig5/cellchat_RTs.rds')
cellchat.pt <- readRDS('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig5/cellchat_PTs.rds')
TypeOfCells<-levels(cellchat.rt@idents)
vertex.receiver=match(c("RT_CS1","RT_CS2","RT_CS3"),TypeOfCells)
netVisual(cellchat.rt, signaling = "MK",
                     vertex.receiver = vertex.receiver, 
                     layout = "hierarchy",
                     vertex.weight=NULL,
                     out.format="pdf",
                     vertex.label.cex = 0.8,
                     title.space = 6,
                     label.dist=4,
                     space.v=1.5,
                     space.h = 1.6)


rankplot <- function(cellchat.obj, sources, targets, measure, color){
  ranked <- rankNet(cellchat.obj,mode = 'single',
                    measure = measure,
                    sources.use = sources, 
                    targets.use = targets,
                    return.data = T)
  ranked_gg <- ggplot_build(ranked$gg.obj)
  ranked_gg$data[[1]]$fill <- color
  plot(ggplot_gtable(ranked_gg))
  return(ranked)
}
# rt overall
# source
rtcs_source <- rankplot(cellchat.rt,c('RT_CS1','RT_CS2','RT_CS3'),
                    c("Myeloid","Luminal","Stromal","NK"),"weight","#d661a4")
# target
rtcs_target <- rankplot(cellchat.rt,c("Myeloid","Luminal","Stromal","NK"),
                    c('RT_CS1','RT_CS2','RT_CS3'),"weight","#638b88")


# rt_cs1
rtcs1_source <- rankplot(cellchat.rt,c('RT_CS1'),
                        c("Myeloid","Luminal","Stromal","NK"),"weight","#d661a4")
rtcs1_target <- rankplot(cellchat.rt,c("Myeloid","Luminal","Stromal","NK"),
                        c('RT_CS1'),"weight","#638b88")

# rt_cs2
rtcs2_source <- rankplot(cellchat.rt,c('RT_CS2'),
                         c("Myeloid","Luminal","Stromal","NK"),"count","#d661a4")
rtcs2_target <- rankplot(cellchat.rt,c("Myeloid","Luminal","Stromal","NK"),
                         c('RT_CS2'),"count","#638b88")

# rt_cs3
rtcs3_source <- rankplot(cellchat.rt,c('RT_CS3'),
                         c("Myeloid","Luminal","Stromal","NK"),"weight","#d661a4")
rtcs3_target <- rankplot(cellchat.rt,c("Myeloid","Luminal","Stromal","NK"),
                         c('RT_CS3'),"weight","#638b88")
# pt overall
ptcs_source <- rankplot(cellchat.pt,c('PT_CS1','PT_CS2','PT_CS3','PT_CS4','PT_CS5'),
         c("Myeloid","Luminal","Basal","Adipocyte","Neutrophils","Fibroblast"),"#d661a4")
ptcs_target <- rankplot(cellchat.pt,c("Myeloid","Luminal","Basal","Adipocyte","Neutrophils","Fibroblast"),
         c('PT_CS1','PT_CS2','PT_CS3','PT_CS4','PT_CS5'),"#638b88")

sheet_list <- list("RTCS_Source"=rtcs_source$signaling.contribution,
                "RTCS_Target"=rtcs_target$signaling.contribution,
                "PTCS_Source"=ptcs_source$signaling.contribution,
                "PTCS_Target"=ptcs_target$signaling.contribution)
write_xlsx(sheet_list,'/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig5/Pathway_contribution_rank.xlsx')

#as source
netVisual_aggregate(cellchat.rt,signaling = "LAMININ",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5) 

netVisual_aggregate(cellchat.rt,signaling = "COLLAGEN",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "FN1",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "GRN",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "WNT",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "EGF",
                    sources.use = c("RT_CS1","RT_CS2","RT_CS3"), 
                    targets.use = c("Myeloid","Luminal","Stromal","NK"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)
# as target
netVisual_aggregate(cellchat.rt,signaling = "BMP",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "EGF",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "LAMININ",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "OCLN",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)

netVisual_aggregate(cellchat.rt,signaling = "PSAP",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)
netVisual_aggregate(cellchat.rt,signaling = "CDH1",
                    sources.use = c("Myeloid","Luminal","Stromal","NK"), 
                    targets.use = c("RT_CS1","RT_CS2","RT_CS3"),
                    vertex.label.cex=2,edge.label.cex=2,arrow.size=1.5)
