library(CellChat)
library(patchwork)
library(Matrix)
​
# For AD data
filtered.counts <- readMM("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_count_matrix.mtx")
rownames(filtered.counts) <- readLines("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/filtered_column_metadata.txt",stringsAsFactors=F)
translate_cell_types=data.frame(
  control=c('astrocytes', 'endothelial', 'fetal_quiescent', 'fetal_replicating', 'hybrid', 'microglia', 'neurons','neurons', 'oligodendrocytes', 'OPC', NA),
  ad=c('Ast',        'End',          NA,                NA,                  NA,      'Mic',       'Ex',     'In',      'Oli',              'Opc','Per'),
  stringsAsFactors=F)
​
translate_cell_types_full=na.omit(translate_cell_types)
​
for(cell_type in na.omit(translate_cell_types_full$ad)){
  print(cell_type)
  filtered.colMetadata$broad.cell.type[filtered.colMetadata$broad.cell.type==cell_type]=translate_cell_types_full$control[translate_cell_types_full$ad==cell_type]
}

#colnames(filtered.counts)=filtered.colMetadata$broad.cell.type
filtered.counts=filtered.counts[,(filtered.colMetadata$broad.cell.type %in% translate_cell_types_full$control)]
filtered.counts=normalizeData(filtered.counts)
filtered.colMetadata=filtered.colMetadata[(filtered.colMetadata$broad.cell.type %in% translate_cell_types_full$control),]

​

cellchat <- createCellChat(object = filtered.counts, meta = filtered.colMetadata, group.by = "broad.cell.type")
cellchat <- addMeta(cellchat, meta = filtered.colMetadata)
cellchat <- setIdent(cellchat, ident.use = "broad.cell.type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
​
​
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)
cellchat1 <- identifyOverExpressedGenes(cellchat)
cellchat2 <- identifyOverExpressedInteractions(cellchat1)
cellchat3 <- projectData(cellchat2, PPI.human) 
cellchat_prob <- computeCommunProb(cellchat3) 
cellchat_prob <- filterCommunication(cellchat_prob, min.cells = 10)
cellchat_prob <- computeCommunProbPathway(cellchat_prob)
cellchat_prob <- aggregateNet(cellchat_prob)
cellchat_prob <- netAnalysis_computeCentrality(cellchat_prob)
cellchat_case <-cellchat_prob
saveRDS(cellchat_case, "cellchat_case.rds")
​
# For healthy data
filtered.healthy=read.csv('/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/normal/2015_dataset/matrix.csv',stringsAsFactors=F)
genes=filtered.healthy[,1]
filtered.healthy=as.matrix(filtered.healthy[,-1])
rownames(filtered.healthy)=genes
print(dim(filtered.healthy))
metasamp=read.table('/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/normal/2015_dataset/GSM_cell_types',header=T,stringsAsFactors=F)
print(dim(metasamp))
filtered.healthy=filtered.healthy[,(metasamp$cell_type %in% translate_cell_types_full$control)]
filtered.healthy=normalizeData(filtered.healthy)
metasamp=metasamp[(metasamp$cell_type %in% translate_cell_types_full$control),]
colnames(metasamp)[2]='broad.cell.type'
#colnames(filtered.healthy)=metasamp$broad.cell.type
filtered.healthy=as.matrix(filtered.healthy)
cellchat <- createCellChat(object = filtered.healthy, meta = metasamp, group.by = "broad.cell.type")
cellchat<- addMeta(cellchat, meta = metasamp)
cellchat<- setIdent(cellchat, ident.use = "broad.cell.type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
​
​
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat<- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)
cellchat1 <- identifyOverExpressedGenes(cellchat)
cellchat2 <- identifyOverExpressedInteractions(cellchat1)
cellchat3 <- projectData(cellchat2, PPI.human) 
cellchat_prob <- computeCommunProb(cellchat3) 
cellchat_prob <- filterCommunication(cellchat_prob, min.cells = 10)
cellchat_prob <- computeCommunProbPathway(cellchat_prob)
cellchat_prob <- aggregateNet(cellchat_prob)
cellchat_prob <- netAnalysis_computeCentrality(cellchat_prob)
cellchat_control <-cellchat_prob
saveRDS(cellchat_control, "cellchat_control.rds")
​
​
#merge two datasets
cellchat_control=readRDS('cellchat_control.rds')
cellchat_case=readRDS('cellchat_case.rds')
levels(cellchat_case@idents)
levels(cellchat_control@idents)
colnames(cellchat_case@data)=filtered.colMetadata$broad.cell.type
metasamp=read.table('/lustre/projects/DTAdb/resources/epigenetics/scRNAseq_microglia/normal/2015_dataset/GSM_cell_types',header=T,stringsAsFactors=F)
metasamp=metasamp[(metasamp$cell_type %in% translate_cell_types_full$control),]
colnames(metasamp)[2]='broad.cell.type'
colnames(cellchat_control@data)=metasamp$broad.cell.type
object.list <- list(control=cellchat_control, case=cellchat_case)
source('modified.R')
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("overview_number_strength.png",p, width=6, height=4)

png("Diff_number_strength_count.png")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
dev.off()

png("Diff_number_strength_weight.png")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

png("Diff_number_strength_number_heatmap1.png")
netVisual_heatmap(cellchat)
dev.off()

png("Diff_number_strength_number_heatmap2.png")
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  pdf(paste0('counts',i,'.pdf'))
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  dev.off()
}

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("compare_pathway_strength.png",p, width=10, height=6)

# png("pathways_similarity.png")
# rankSimilarity(cellchat, type = "functional")
# dev.off()

# library(ComplexHeatmap)
# pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# gg1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
# gg2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
# p <- gg1 + gg2
# ggsave("compare_signal_pattern_all.png",p, width=10, height=6)

gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2,3,4,5,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6), targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_ast_prob.png",p, width=12, height=10)

gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4,5,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6), targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_endo_prob.png",p, width=12, height=10)

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4,5,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2,4,5,6), targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_mic_prob.png",p, width=12, height=10)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,3,5,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,5,6), targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_neuron_prob.png",p, width=12, height=10)

gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2,3,4,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,6), targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_olig_prob.png",p, width=12, height=10)

gg1 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1,2,3,4,5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in case", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5), targets.use = 6,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in case", angle.x = 45, remove.isolate = T)
p <- gg1 + gg2
ggsave("compare_LR_OPC_prob.png",p, width=12, height=10)

pos.dataset = "case"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "case",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets ="control",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(2,3,4,5,6), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[,"interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use =1, targets.use =c(2,3,4,5,6), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
p <- gg1 + gg2
ggsave("compare_LR_ast_gene.png",p, width=12, height=10)

# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# png("signalgp_similarity.png")
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# dev.off()

pathways.show <- c("PSAP") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  pdf(paste0('PSAP',i,'.pdf'))
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off()
}
