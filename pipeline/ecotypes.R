suppressPackageStartupMessages({
library(cluster)
library(ggplot2)
library(viridis)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("MyDiscovery", "Carcinoma")
args = c("MyDiscovery", "Carcinoma")
args = commandArgs(T)  
dataset = args[1]
fractions = args[2]

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
	colnames(classes) = c('ID', cell_type)
	if(is.null(all_classes))
	{
		all_classes = classes
	}else{
		all_classes = merge(all_classes, classes, by = 'ID')
	}
}

all_mapping$InitialID = paste(all_mapping$CellType, all_mapping$InitialState, sep = "_")
all_mapping$ID = paste(all_mapping$CellType, all_mapping$State, sep = "_")
write.table(all_mapping, file.path(output_dir, "mapping_all_states.txt"), sep = "\t", row.names = F)

#clinical = read_clinical((all_classes_filt[,1]))
#all_classes_filt = all_classes_filt[clinical$Tissue == "Tumor",]

rownames(all_classes) = all_classes[,1]
all_classes = all_classes[,-1]
all_classes = t(all_classes) 

clusters = do.call(rbind, lapply(1:nrow(all_classes), function(i){
	vals = sort(unique(as.vector(all_classes[i,])))
	tmp = do.call(rbind, lapply(vals, function(x){
		all_classes[i,] == x
	}))
	rownames(tmp) = paste0(rownames(all_classes)[i], '_', vals)
	tmp
}))

clusters = clusters[match(all_mapping$InitialID, rownames(clusters)),]
rownames(clusters) = all_mapping$ID
write.table(clusters, file.path(output_dir, "binary_classification_all_states.txt"), sep = "\t")

jaccard = matrix(NA, nrow(clusters), nrow(clusters))
for(i in 1:(nrow(clusters)))
{
	for(j in (i):(nrow(clusters)))
	{
		int <- sum(clusters[i,] & clusters[j,])
		idx <- int / (sum(clusters[i,]) + sum(clusters[j,]) - int)
		p = 1 - phyper(int, sum(clusters[i,]), ncol(all_classes) - sum(clusters[i,]), sum(clusters[j,]))
		if(p > 0.01){
			idx = 0
		}
		jaccard[i, j] = jaccard[j, i] = idx
	}
}
rownames(jaccard) = colnames(jaccard) = rownames(clusters)
write.table(jaccard, file.path(output_dir, "jaccard_matrix.txt"), sep = "\t")

hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{
	#gap = clusGap(jaccard, hclusCut, 50, B = 100) 
	silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
		sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
		tmp <<- summary(sil)
		tmp$avg.width
	}))

	g2 <- ggplot(silh, aes(x = K, y = Silhouette)) + 
		geom_point() +
		#ylab("Sparseness (coef)") + 
		geom_line() + 
		geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") + 
		theme_bw() +
		theme(panel.grid = element_blank()) + 
		theme(aspect.ratio = 1) 
		#facet_wrap(CellType, ncol = 4) + 
		
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
#n_clust = 5
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

tb = table(top_ann$InitialEcotype)
tb = tb[tb > 2]

top_ann = top_ann[top_ann$InitialEcotype %in% names(tb),]

top_ann$ID = rownames(top_ann)
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
pdf(file.path(output_dir, "jaccard_matrix.pdf"), width = 8, height = 7, family = "Helvetica")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("Ecotype", "Cell type", "State"), 
 legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
color_palette = c("white", viridis(4)), 
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

png(file.path(output_dir, "jaccard_matrix.png"), width = 8, height = 7, res = 300, units = "in", family = "Helvetica")
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	
decorate_heatmap_body("ht1", {
    grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 2))
})
tmp = dev.off()
