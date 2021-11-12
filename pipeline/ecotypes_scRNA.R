suppressPackageStartupMessages({
library(cluster)
library(ggplot2)
library(viridis)
library(reshape2)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("discovery_scRNA_CRC", "Cell_type_specific_genes", "0.05") 
args = commandArgs(T)   
dataset = args[1]
fractions = args[2]
p_val_cutoff = as.numeric(as.character(args[3]))

if(is.na(p_val_cutoff))
{
	p_val_cutoff = 1
}

key_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
states_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery")
output_dir = file.path("../EcoTyper", dataset, fractions, "Ecotypes", "discovery")
dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path(key_dir, "rank_data.txt"))

all_mapping = NULL
all_classes = NULL
all_classes_filt = NULL
	
for(cell_type in key[,1])
{	
	#print(cell_type) 
	n_states = key[key[,1] == cell_type,2]

	mapping_path = file.path(states_dir, cell_type, n_states, 'mapping_to_initial_states.txt')
	classes_path = file.path(states_dir, cell_type, n_states, 'initial_state_assignment.txt')
	classes_filt_path = file.path(states_dir, cell_type, n_states, 'state_assignment.txt')

	mapping = read.delim(mapping_path)
	mapping = mapping[,c("State", "InitialState")] 
	mapping$CellType = cell_type
	all_mapping = rbind(all_mapping, mapping)

	classes = read.delim(classes_path)
	clinical = read_clinical(classes$ID, dataset = dataset)
	classes$Sample = clinical$Sample
	
	classes = as.data.frame(table(classes$Sample, classes$State))
	colnames(classes) = c("ID", "State", "Freq")
	splits = split(classes, classes$ID)
	classes = do.call(rbind, lapply(splits, function(spl)
	{
		spl$Frac = spl$Freq / sum(spl$Freq)
		spl$Max= ifelse((!is.na(max(spl$Freq))) & (max(spl$Freq) > 0) & (max(spl$Freq) == spl$Freq), 1, 0)
		spl
		}))
	classes$CellType = cell_type
	all_classes = rbind(all_classes, classes)
}

all_mapping$InitialID = paste(all_mapping$CellType, all_mapping$InitialState, sep = "_")
all_mapping$ID = paste(all_mapping$CellType, all_mapping$State, sep = "_")
write.table(all_mapping, file.path(output_dir, "mapping_all_states.txt"), sep = "\t", row.names = F)

casted = dcast(all_classes, ID ~ CellType + State, value.var = "Max")
clusters = t(casted[,-1])
clusters[is.na(clusters)] = 0

clusters = clusters[match(all_mapping$InitialID, rownames(clusters)),]
rownames(clusters) = all_mapping$ID

colnames(clusters) = casted[,1]
write.table(clusters, file.path(output_dir, "binary_classification_all_states.txt"), sep = "\t")

jaccard = matrix(NA, nrow(clusters), nrow(clusters))
for(i in 1:(nrow(clusters)))
{
	for(j in (i):(nrow(clusters)))
	{
		int <- sum(clusters[i,] & clusters[j,])
		idx <- int / (sum(clusters[i,]) + sum(clusters[j,]) - int)
		
		p = 1 - phyper(int, sum(clusters[i,]), ncol(clusters) - sum(clusters[i,]), sum(clusters[j,]))
		if(is.na(p) | (p >= p_val_cutoff)){			
			idx = 0
		}		
		jaccard[i, j] = jaccard[j, i] = idx
	}
}
jaccard[is.na(jaccard)] = 0

rownames(jaccard) = colnames(jaccard) = rownames(clusters)
write.table(jaccard, file.path(output_dir, "jaccard_matrix.txt"), sep = "\t")

hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{	
	silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
		sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
		tmp <<- summary(sil)
		tmp$avg.width
	}))

	g2 <- ggplot(silh, aes(x = K, y = Silhouette)) + 
		geom_point() +		
		geom_line() + 
		geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") + 
		theme_bw() +
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 1) 		
		
	pdf(file.path(output_dir, paste0("nclusters_", name, ".pdf")), width = 4, height = 4)
	plot(g2)
	tmp = dev.off()
	png(file.path(output_dir, paste0("nclusters_", name, ".png")), width = 4, height = 4, res = 300, units = "in")
	plot(g2)
	tmp = dev.off()
	silh[which.max(silh$Silhouette), 1]	
}
hc = hclust(as.dist(1-jaccard), method = "average")
n_clust = choose_clusters(jaccard, "jaccard", range = 2:(nrow(jaccard) - 1))
clust =  hclusCut(jaccard, n_clust)$cluster

sil = silhouette(clust, as.dist(1-jaccard))
avg_silhouette = summary(sil)
write.table(avg_silhouette$avg.width, file.path(output_dir, "silhouette_initial.txt"), sep = "\t", row.names = F)

top_ann = as.data.frame(t(sapply(rownames(jaccard), function(x) {
	s= strsplit(x, "_")[[1]]
	c(paste0(s[-length(s)],collapse = "_"), s[length(s)])
})))

colnames(top_ann) = c("CellType","State")
top_ann$InitialEcotype = as.factor(sprintf("IE%02d", clust))
write.table(top_ann, file.path(output_dir, "initial_ecotypes.txt"), sep = "\t")

top_ann$ID = rownames(top_ann)
top_ann = top_ann[order(top_ann$InitialEcotype),]
write.table(top_ann, file.path(output_dir, "ecotypes.txt"), sep = "\t", row.names = F)

jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]

top_ann$"Cell type" = top_ann$CellType
diag(jaccard) = 1
pdf(file.path(output_dir, "initial_jaccard_matrix.pdf"), width = 7, height = 7, family = "Helvetica")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("InitialEcotype", "Cell type", "State"), 
 legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
color_palette = c("gray", viridis(4)), raster_quality = 5,
color_range = c(0, 0.1, 0.2, 0.3))
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	

ord = top_ann$InitialEcotype
dup = (which(!duplicated(ord)) - 1)
fract = dup / nrow(top_ann)
width =  c(fract[-1], 1) - fract
decorate_heatmap_body("ht1", {
    #grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 2))
})
tmp = dev.off()

tb = table(top_ann$InitialEcotype)
tb = tb[tb > 2]

top_ann = top_ann[top_ann$InitialEcotype %in% names(tb),]

nm = unique(top_ann$InitialEcotype)
mapping = sprintf("E%d", 1:length(nm))
names(mapping) = nm

top_ann$Ecotype = mapping[as.character(top_ann$InitialEcotype)]
top_ann$Ecotype = ecotype_to_factor(top_ann$Ecotype)
top_ann = top_ann[order(top_ann$Ecotype),]
write.table(top_ann, file.path(output_dir, "ecotypes.txt"), sep = "\t", row.names = F)

jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]

sil <- silhouette(as.numeric(as.character(gsub("E", "", as.character(top_ann$Ecotype)))), as.dist(1-jaccard))
avg_silhouette <<- summary(sil)
write.table(avg_silhouette$avg.width, file.path(output_dir, "silhouette.txt"), sep = "\t", row.names = F)

top_ann$"Cell type" = top_ann$CellType

pdf(file.path(output_dir, "jaccard_matrix.pdf"), width = 7, height = 7, family = "Helvetica")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("Ecotype", "Cell type", "State"), 
 legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
color_palette = c("gray", viridis(4)), raster_quality = 5,
color_range = c(0, 0.1, 0.2, 0.3))
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	

ord = top_ann$Ecotype
dup = (which(!duplicated(ord)) - 1)
fract = dup / nrow(top_ann)
width =  c(fract[-1], 1) - fract
decorate_heatmap_body("ht1", {
    grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 2))
})
tmp = dev.off()

png(file.path(output_dir, "jaccard_matrix.png"), width = 5, height = 7, res = 300, units = "in", family = "Helvetica")
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	
decorate_heatmap_body("ht1", {
    grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 2))
})
tmp = dev.off()
