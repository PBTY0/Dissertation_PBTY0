library(CellChat)
library(patchwork)
library(Matrix)
library(NMF)
library(ggalluvial)

# Part I
filtered.counts <- readMM("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_count_matrix.mtx")
rownames(filtered.counts) <- readLines("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_column_metadata.txt")
filtered.counts=normalizeData(filtered.counts)


cellchat <- createCellChat(object = filtered.counts, meta = filtered.colMetadata, group.by = "broad.cell.type")
cellchat <- addMeta(cellchat, meta = filtered.colMetadata)
cellchat <- setIdent(cellchat, ident.use = "broad.cell.type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human


CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)
cellchat1 <- identifyOverExpressedGenes(cellchat)
cellchat2 <- identifyOverExpressedInteractions(cellchat1)
cellchat3 <- projectData(cellchat2, PPI.human) 
cellchat_prob <- computeCommunProb(cellchat3) 
cellchat_prob <- filterCommunication(cellchat_prob, min.cells = 10)
df.net <- subsetCommunication(cellchat_prob)
write.csv(df.net,"net_lr.csv")
cellchat_prob <- computeCommunProbPathway(cellchat_prob)
df.netP <- subsetCommunication(cellchat_prob,slot.name="netP")
write.csv(df.netP,"net_pathway.csv")


# library(mindr)
# (out <- capture.output(str(cellchat_prob)))
# out2 <- paste(out, collapse="\n")
# p <- mm(gsub("\\.\\.@","# ",gsub("\\.\\. ","#",out2)),type ="text")
# saveTXT("data.pdf",p)

#
levels(cellchat_prob@idents)

pdf("multiLR_mic.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =5, targets.use = c(1,2,3,4,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,3,4,6,7,8), targets.use = 5, remove.isolate = FALSE)
dev.off()

pdf("multiLR_Ast.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =1, targets.use = c(2,3,4,5,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(2,3,4,5,6,7,8), targets.use = 1, remove.isolate = FALSE)
dev.off()

pdf("multiLR_End.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =2, targets.use = c(1,3,4,5,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,3,4,5,6,7,8), targets.use = 2, remove.isolate = FALSE)
dev.off()

pdf("multiLR_Ex.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =3, targets.use = c(1,2,4,5,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,4,5,6,7,8), targets.use = 3, remove.isolate = FALSE)
dev.off()

pdf("multiLR_In.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =4, targets.use = c(1,2,3,5,6,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,3,5,6,7,8), targets.use = 4, remove.isolate = FALSE)
dev.off()

pdf("multiLR_Oli.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =6, targets.use = c(1,2,3,4,5,7,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,3,4,5,7,8), targets.use = 6, remove.isolate = FALSE)
dev.off()

pdf("multiLR_Opc.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =7, targets.use = c(1,2,3,4,5,6,8), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,3,4,5,6,8), targets.use =7, remove.isolate = FALSE)
dev.off()

pdf("multiLR_Per.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat_prob, sources.use =8, targets.use = c(1,2,3,4,5,6,7), remove.isolate = FALSE)
netVisual_bubble(cellchat_prob, sources.use =c(1,2,3,4,5,6,7), targets.use =8, remove.isolate = FALSE)
dev.off()

cellchat_prob@netP$pathways
head(cellchat_prob@LR$LRsig)
#


cell_types=c('Ast', 'End', 'Ex', 'In', 'Mic', 'Oli', 'Opc', 'Per')
cell_types_numbers=c(1:length(cell_types))

# for(i in cell_types_numbers){
#     cell_type=cell_types[i]
#     from=cell_types_numbers[cell_types_numbers!=i]
#     to=i
#     pdf(paste0('net_vis_bubble_full_to_',cell_type,'.pdf'),height=10,width=10)
#     netVisual_bubble(cellchat_prob, sources.use = from, targets.use = to, remove.isolate = FALSE)
#     dev.off()
#     from=i
#     to=cell_types_numbers[cell_types_numbers!=i]
#     pdf(paste0('net_vis_bubble_full_from_',cell_type,'.pdf'),height=10,width=10)
#     netVisual_bubble(cellchat_prob, sources.use = from, targets.use = to, remove.isolate = FALSE)
#     dev.off()
# }


# netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30,thresh = 0.001)


# for(direction in c('incoming','outgoing')){

#     pdf(paste0(direction,'galluvial.pdf'))
#     selectK(cellchat_prob, pattern=direction)
#     dev.off()

#     pdf(paste0(direction,'identify_patterns.pdf'))
#     nPatterns = 3
#     cellchat <- identifyCommunicationPatterns(cellchat_prob, pattern =direction, k = nPatterns)
#     dev.off()

#     pdf(paste0(direction,'river.pdf'))
#     netAnalysis_river(cellchat_prob, pattern=direction)
#     dev.off()


# }


pdf("aggregatedcommu.pdf")
netVisual_circle(cellchat_prob@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("aggregatedcommu2.pdf")
netVisual_circle(cellchat_prob@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("aggregatedcommu_singlec.pdf",height=12)
mat <- cellchat_prob@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat_prob@netP$pathways
pathways.show <- c("PSAP") 
vertex.receiver = seq(1,4)
pdf("PSAP_hierarchy.pdf")
netVisual_aggregate(cellchat_prob, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()
pdf("PSAP_circle.pdf")
netVisual_aggregate(cellchat_prob, signaling = pathways.show, layout = "circle")
dev.off()
pdf("PSAP_chord.pdf")
netVisual_aggregate(cellchat_prob, signaling = pathways.show, layout = "chord")
dev.off()
pdf("PSAP_heatmap.pdf")
netVisual_heatmap(cellchat_prob, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("PSAP_LR_contribution.pdf")
netAnalysis_contribution(cellchat_prob, signaling = pathways.show)
dev.off()

pairLR.PSAP <- extractEnrichedLR(cellchat_prob, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.PSAP[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
pdf("PSAP_hierarchy2.pdf")
netVisual_individual(cellchat_prob, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
pdf("PSAP_circle2.pdf")
netVisual_individual(cellchat_prob, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()

# p <- plotGeneExpression(cellchat_prob, signaling = "NRG")
# ggsave("NRG_GeneExpression_vln.pdf",p, width=8, height=8)
# p <- plotGeneExpression(cellchat_prob, signaling = "NRG",type="dot")
# ggsave("NRG_GeneExpression_dot.pdf",p, width=8, height=6)
# error:Error in new(Class = "seurat", raw.data = raw.data, is.expr = is.expr,: argument "raw.data" is missing, with no default

# ht1 <- netAnalysis_signalingRole_heatmap(cellchat_prob, pattern = "outgoing")
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat_prob, pattern = "incoming")
# p <- ht1 + ht2
# ggsave("outin_heatmap.pdf",p, width=6, height=4)