library(CellChat)
library(patchwork)
library(Matrix)
filtered.counts=readMM("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/ROSMAPMM.mtx")
cells=read.csv('/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv')
genes=read.csv('/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv')
colnames(filtered.counts)=1:ncol(filtered.counts)
rownames(filtered.counts)=genes$x
filtered.counts=filtered.counts[!duplicated(rownames(filtered.counts)),]
filtered.counts=normalizeData(filtered.counts)
cellchat <- createCellChat(object = filtered.counts, meta = cells, group.by = "broad_class")
cellchat <- addMeta(cellchat, meta = cells)
cellchat <- setIdent(cellchat, ident.use = "broad_class")

# Part I
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human 
pdf("databasecategories_newdt.pdf")
showDatabaseCategory(CellChatDB)
dev.off()
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
print(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 

# Part II
options(future.globals.maxSize= 891289600)
cellchat1 <- computeCommunProb(cellchat)  
cellchat <- filterCommunication(cellchat1, min.cells = 20)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"net_lr2.csv")
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat,slot.name="netP")
write.csv(df.netP,"net_pathway2.csv")

levels(cellchat@idents)
pdf("new_multiLR_mic.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =5, targets.use = c(1,2,3,4,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,3,4,6,7,8,9), targets.use = 5, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_Ast.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =1, targets.use = c(2,3,4,5,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(2,3,4,5,6,7,8,9), targets.use = 1, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_End.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =2, targets.use = c(1,3,4,5,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3,4,5,6,7,8,9), targets.use = 2, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_Ex.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =3, targets.use = c(1,2,4,5,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,4,5,6,7,8,9), targets.use = 3, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_In.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,3,5,6,7,8,9), targets.use = 4, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_None.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =6, targets.use = c(1,2,3,4,5,7,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,3,4,5,7,8,9), targets.use = 6, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_Oli.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =7, targets.use = c(1,2,3,4,5,6,8,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,3,4,5,6,8,9), targets.use = 7, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_Opc.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =8, targets.use = c(1,2,3,4,5,6,7,9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,2,3,4,5,6,7,9), targets.use =8, remove.isolate = FALSE)
dev.off()

pdf("new_multiLR_Per.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use =9, targets.use = c(1,2,3,4,5,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =9, targets.use = c(1,2,3,4,5,6,7,8), remove.isolate = FALSE)
dev.off()

cellchat@netP$pathways
head(cellchat@LR$LRsig)

# cell_types=c('Astr', 'Endo', 'Ex', 'Inh', 'Micr', 'None','Olig', 'OPC', 'PerI')
# cell_types_numbers=c(1:length(cell_types))

# for(i in cell_types_numbers){
#     cell_type=cell_types[i]
#     from=cell_types_numbers[cell_types_numbers!=i]
#     to=i
#     pdf(paste0('net_vis_bubble_full_to_',cell_type,'.pdf'),height=10,width=10)
#     netVisual_bubble(cellchat, sources.use = from, targets.use = to, remove.isolate = FALSE)
#     dev.off()
#     from=i
#     to=cell_types_numbers[cell_types_numbers!=i]
#     pdf(paste0('net_vis_bubble_full_from_',cell_type,'.pdf'),height=10,width=10)
#     netVisual_bubble(cellchat, sources.use = from, targets.use = to, remove.isolate = FALSE)
#     dev.off()
# }



pdf("new_interactions.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("new_interactions2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("new_aggregatedcommu_singlec.pdf",height=12)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat@netP$pathways
pathways.show <- c("PSAP") 
vertex.receiver = seq(1,4)
pdf("new_PSAP_hierarchy.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()
pdf("new_PSAP_circle.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("new_PSAP_chord.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
pdf("new_PSAP_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("new_PSAP_LR_contribution.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pairLR.PSAP <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.PSAP[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
pdf("new_PSAP_hierarchy2.pdf")
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
pdf("new_PSAP_circle2.pdf")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()
