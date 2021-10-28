suppressPackageStartupMessages({
library(cluster)
library(ggplot2)
library(viridis)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("Lung", "TR12", "Ecotype", "Histology", "Tissue")
args = c("Lung", "TR12", "Ecotype", "Histology", "Tissue")
args = commandArgs(T) 
dataset = args[1]
fractions = args[2]
fraction_processing = args[3]
top_cols = args[4:length(args)]


if(is.na(top_cols[1]))
{
	top_cols = c("Ecotype")
}
top_cols = unique(c(top_cols, "Ecotype"))

key_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
states_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery")
output_dir = file.path("../EcoTyper", dataset, fractions, "Ecotypes", "discovery")
dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path(key_dir, "rank_data.txt"))
ecotypes = read.delim(file.path(output_dir, 'ecotypes.txt'), sep = '\t')
mapping = read.delim(file.path(output_dir, 'mapping_all_states.txt'), sep = '\t')
ecotypes$Ecotype = ecotype_to_factor(ecotypes$Ecotype)
ecotypes = ecotypes[order(ecotypes$Ecotype),]

all_H = NULL
all_classes_filt = NULL
for(cell_type in key[,1])
{	
	#print(cell_type) 
	n_clusters = key[key[,1] == cell_type,2]
	classes = read.delim(file.path(states_dir, cell_type, n_clusters, 'state_abundances.txt'))
	rownames(classes) = paste0(cell_type, "_", rownames(classes))
	all_H = rbind(all_H, classes)

	classes = read.delim(file.path(states_dir, cell_type, n_clusters, 'state_assignment.txt'))
	clusters = ecotypes[ecotypes$CellType == cell_type,] 
	classes = classes[classes$State %in% clusters$State,]
	classes = classes[,c("ID", "State")]
	colnames(classes) = c('ID', cell_type)

	if(is.null(all_classes_filt))
	{
		all_classes_filt = classes
	}else{
		all_classes_filt = merge(all_classes_filt, classes, by = 'ID', all = T)
	}
} 

all_H = all_H[match(ecotypes$ID, rownames(all_H)),]
write.table(all_H, file.path(output_dir, "combined_state_abundances.txt"), sep = "\t")

H = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
	clst <<- clst
	#print(clst)
	inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
	apply(all_H[rownames(all_H) %in% inc,,drop = F], 2, mean)
}))
rownames(H) = levels(ecotypes$Ecotype)
#write.table(H, file.path(output_dir, "ecotype_abundance.txt"), sep = "\t")
H = apply(H, 2, function(x) x / sum(x))
write.table(H, file.path(output_dir, "ecotype_abundance.txt"), sep = "\t")

p_vals = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
	clst <<- clst
	#print(clst)
	inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
	
	apply(all_H, 2, function(x){
		x <<- x
		err <<- F
		p <<- 1
		tryCatch({
			p <<- t.test(x[rownames(all_H) %in% inc], x[!(rownames(all_H) %in% inc)])$p.value
		}, error = function(x) err <<- T)
		p
	})
	
})) 
rownames(p_vals) = levels(ecotypes$Ecotype)
write.table(p_vals, file.path(output_dir, "assignment_p_vals.txt"), sep = "\t")

assignment = as.data.frame(apply(H, 2, function(x) rownames(H)[which.max(x)]))

clinical = data.frame(ID = rownames(assignment), MaxEcotype = assignment[,1])
clinical$AssignmentP = sapply(1:ncol(H), function(i) {
	p_vals[which.max(H[,i]), i] 
	})
clinical$AssignmentQ = p.adjust(clinical$AssignmentP, method = "BH")
clinical$AssignedToEcotypeStates = clinical$ID %in% all_classes_filt$ID
clinical$Ecotype = ifelse((clinical$AssignmentQ < 0.25) & clinical$AssignedToEcotypeStates, as.character(clinical$MaxEcotype), "Unassigned")

tmp = read_clinical(clinical$ID, dataset = dataset)
to_rem = colnames(tmp)[colnames(tmp) %in% colnames(clinical)]
to_rem = to_rem[to_rem != "ID"]
tmp = tmp[,!colnames(tmp) %in% to_rem]
clinical = merge(clinical, tmp, by = "ID", all.x = T)

clinical = clinical[order(clinical$MaxEcotype),]
H = H[,match(clinical$ID, colnames(H))]
all_H = all_H[,match(clinical$ID, colnames(all_H))]
all_H = all_H[match(ecotypes$ID, rownames(all_H)),]

write.table(clinical, file.path(output_dir, "initial_ecotype_assignment.txt"), sep = "\t")

rownames(clinical) = clinical$ID
rownames(ecotypes) = ecotypes$ID

h <- heatmap_simple(all_H, top_annotation = clinical, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	column_split = ifelse(clinical$Ecotype == "Unassigned", "Unassigned", "Assigned"),
	width = unit(7, "in"), height = unit(4, "in"),
	legend_name = "State abundance",
	color_range = seq(0, quantile(as.matrix(all_H), .9), length.out = 8), color_palette = c("white", viridis(8)), raster_quality = 5)

pdf(file.path(output_dir, "heatmap_all_samples.pdf"), width = 12, height = 7)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
tmp = dev.off()

clinical_filt = clinical[clinical$Ecotype != "Unassigned",]
write.table(clinical_filt, file.path(output_dir, "ecotype_assignment.txt"), sep = "\t")

clinical_filt$Ecotype = ecotype_to_factor(clinical_filt$Ecotype)
clinical_filt = clinical_filt[order(clinical_filt$Ecotype),]

small_H = as.matrix(all_H[,match(clinical_filt$ID, colnames(all_H))])
h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	width = unit(5, "in"), height = unit(3, "in"),
	legend_name = "State abundance",
	color_range = seq(0, quantile(as.matrix(all_H), .9), length.out = 8), color_palette = c("white", viridis(8)), raster_quality = 20)

pdf(file.path(output_dir, "heatmap_assigned_samples_viridis.pdf"), width = 8, height = 6)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
#draw(h1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
suppressWarnings({
ord = ecotypes[row_order(h),]
})
dup = (which(!duplicated(ord$Ecotype)) - 1)
fract = dup / nrow(ecotypes)
height =  c(fract[-1], 1) - fract

suppressWarnings({
ord = clinical_filt[column_order(h),]
})
dup = (which(!duplicated(ord$Ecotype)) - 1)
fract2 = dup / nrow(clinical_filt)
width =  c(fract2[-1], 1) - fract2
decorate_heatmap_body("hmap", {
    grid.rect(unit(fract2, "native"), unit(1-fract, "native"), width = unit(width, "native"), height = unit(height, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(output_dir, "heatmap_assigned_samples_viridis.png"), width = 8, height = 6, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
    grid.rect(unit(fract2, "native"), unit(1-fract, "native"), width = unit(width, "native"), height = unit(height, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	width = unit(5, "in"), height = unit(3, "in"), 
	legend_name = "State abundance",
	color_range = c(seq(0, quantile(small_H, .8), length.out = 8)), color_palette = c("black", brewer.pal(8, "YlGnBu")), raster_quality = 20)

pdf(file.path(output_dir, "heatmap_assigned_samples_YlGnBu.pdf"), width = 8, height = 6)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
    grid.rect(unit(fract2, "native"), unit(1-fract, "native"), width = unit(width, "native"), height = unit(height, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(output_dir, "heatmap_assigned_samples_YlGnBu.png"), width = 8, height = 6, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
    grid.rect(unit(fract2, "native"), unit(1-fract, "native"), width = unit(width, "native"), height = unit(height, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()


