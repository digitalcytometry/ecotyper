suppressPackageStartupMessages({
library(NMF)
library(ComplexHeatmap)
source("lib/misc.R")
source("lib/heatmaps.R")
source("lib/scRNA.R")
source("lib/ecotyper.R")
})

args = c("Carcinoma", "Carcinoma_Fractions", "Epithelial.cells", "8", "single_cell_test_data", "FALSE")
args = commandArgs(T) 
dataset = args[1]
fractions = args[2]
cell_type = args[3]
n_states = args[4]
scRNA_dataset = args[5]
calculate_significance_flag = as.logical(args[6])
top_cols = args[7:length(args)]

cat(paste0("Running cell state recovery on cell type: ", cell_type, "...\n"))

n_bootstraps = 10

if(is.na(top_cols[1]))
{
	top_cols = c("State")
}

states_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type, n_states) 
output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "recovery", "scRNA", scRNA_dataset, cell_type, n_states) 
dir.create(output_dir, recursive = T, showWarning = F) 

W = read.delim(file.path(states_dir, "W.txt"))
mapping = read.delim(file.path(states_dir, "mapping_to_initial_states.txt"))
initial_gene_info = read.delim(file.path(states_dir, "initial_gene_info.txt"))

tmp = read_scRNA(scRNA_dataset, cell_type, log2 = T) 
if(is.null(tmp))
{
	quit()
}
data = tmp[[1]]
scRNA_clinical = tmp[[2]]

data = data[match(initial_gene_info$Gene, rownames(data)),]
rownames(data) = initial_gene_info$Gene
data[is.na(data)] = 0

if(ncol(data) > 2500)
{
	set.seed(1234)
	data = data[,sample(1:ncol(data), 2500)]
	scRNA_clinical = scRNA_clinical[match(colnames(data), scRNA_clinical$ID),]
}

unscaled_data = data

write.table(data, file.path(output_dir, "expression_matrix.txt"), sep = "\t")
write.table(scRNA_clinical, file.path(output_dir, "annotation.txt"), sep = "\t")

tmp_cp = colnames(data)
data = t(apply(data, 1, scale))
colnames(data) = tmp_cp
data[is.na(data)] = 0

write.table(data, file.path(output_dir, "expression_matrix_scaled.txt"), sep = "\t")

nmf_reference <- NMFpredict(W, data) 
save(nmf_reference, file = file.path(file.path(output_dir, "nmf.RData")))

H = nmf_reference@fit@H
rownames(H) = colnames(W)
write.table(H, file.path(output_dir, "H.txt"), sep = "\t")

H <- apply(H, 2, function(x) (x/sum(x)))
write.table(H, file.path(output_dir, "initial_state_abundances.txt"), sep = "\t")

assignment = data.frame(ID = colnames(H), State = apply(H, 2, function(x) rownames(H)[which.max(x)]))
assignment = assignment[order(assignment$State),]
write.table(assignment, file.path(output_dir, "initial_state_assignment.txt"), sep = "\t", row.names= F)
initial_assignment = assignment

H = H[match(mapping$InitialState, rownames(H)),]
rownames(H) = mapping$State
write.table(H, file.path(output_dir, "state_abundances.txt"), sep = "\t")

assignment$InitialState = assignment$State
assignment$State = mapping[match(assignment$InitialState, mapping$InitialState), "State"]
assignment = assignment[!is.na(assignment$State),]
assignment = assignment[order(assignment$State),]
write.table(assignment, file.path(output_dir, "state_assignment.txt"), sep = "\t", row.names= F)

#this section plots the heatmap

sum_squares <- function(hmcs)
{	
	#hmcs[is.na(hmcs)] <- 0
	for(i in 1:10){  
		#print(i)
	    hmcs <- t(apply(hmcs, 1, function(x) x/sqrt(sum(x^2, na.rm = T))))
	    #hmcs[is.na(hmcs)] <- 0
	    hmcs <- hmcs - apply(hmcs,1,mean, na.rm = T)
	    hmcs[which(hmcs<0)] <- 0
	}

	hmcs <- hmcs - apply(hmcs,1,mean, na.rm = T)

	thres <- .025
	hmcs[(hmcs>thres)] <- thres
	hmcs[(hmcs<= -thres)] <- -thres

	#sliding window
	if(ncol(hmcs) > 100)
	{
		w <- 1
		if(ncol(hmcs) > 250)
		{
			w <- 2
		}		
		if(ncol(hmcs) > 500)
		{
			w <- 3
		}
	    hmcs[(hmcs<= 0)] <- 0	    
	    for(i in 1:nrow(hmcs))
	    {
			d <- sapply(c(1:(ncol(hmcs)-w)), function(j) mean(hmcs[i,c(j:(j+w))], na.rm = T))			
			if(i==1) smooth <- d
			else{smooth <- rbind(smooth,d)}
	    }
	    smooth = cbind(smooth, hmcs[,(ncol(hmcs) - w + 1):ncol(hmcs)])
	    rownames(smooth) = rownames(hmcs)	    
	    colnames(smooth) = colnames(hmcs)[1:(ncol(hmcs))]
	    hmcs <- smooth
	    hmcs <- (hmcs - sapply(1:nrow(hmcs),function(i) mean(hmcs[i,], na.rm = T)))
	    thres <- .005
	    hmcs[which(hmcs>thres)] <- thres
	    hmcs[which(hmcs<= -thres)] <- -thres
	}	
	hmcs <- hmcs - apply(hmcs,1,mean, na.rm = T)
	hmcs * 1000
}

suppressWarnings({
discovery_data = fread(file.path(states_dir, "heatmap_data.txt"), data.table = F)
})
rownames(discovery_data) = discovery_data[,1]
discovery_data = discovery_data[,-1]
discovery_annotation = read.delim(file.path(states_dir, "heatmap_top_ann.txt"))
gene_info = read.delim(file.path(states_dir, "gene_info.txt"))

assignment = merge(assignment, scRNA_clinical, by = "ID", all.x = T)
assignment = assignment[order(assignment$State),]
rownames(assignment) = assignment$ID

scRNA_plot_data = unscaled_data
scRNA_plot_data = scRNA_plot_data[match(rownames(discovery_data), rownames(scRNA_plot_data)),match(assignment$ID, colnames(scRNA_plot_data))]
rownames(scRNA_plot_data) = rownames(discovery_data)
smoothed_data = sum_squares(scRNA_plot_data)

color_palette = c("black", "#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4", "#005af2", "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4", "#0049e2", "#0046e0", "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6", "#003ed4", "#003dd3", "#003cd1", "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9", "#0036c7", "#0034c5", "#0033c4", "#0131ba", "#0c2eae", "#132ca1", "#172995", "#192688", "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e", "#191743", "#171539", "#15122e", "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a", "#26290e", "#30360f", "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500", "#8ba400", "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000", "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00", "#d5fc00", "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00", "#e0ff00", "#e2ff00", "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02", "#eeff07", "#f0ff0b", "#f2ff0e", "#f4ff12", "#f5ff14", "#f7ff17")

suppressWarnings({
p <- heatmap_simple(discovery_data, top_annotation = discovery_annotation, dataset = dataset,
	top_columns = unique(c(top_cols[top_cols %in% colnames(discovery_annotation)], "State")), 
	column_title = "Discovery dataset\n",
	name = "hmap1",
scale_rows = T, width = unit(3, "in"), height = unit(3, "in"),
raster_quality = 7, legend_name = "Relative expression",
color_palette = color_palette, show_annotation_name = F,
color_range = 2 * seq(-1, 1, length.out = length(color_palette) - 1))

q <- heatmap_simple(smoothed_data, top_annotation = assignment, dataset = dataset,
	top_columns = unique(c(top_cols[top_cols %in% colnames(assignment)], "State")),
	column_title = paste0("Single cells assigned to cell states\nn = ", nrow(assignment), " (", format(nrow(assignment) * 100 / ncol(H), digits = 2), "%)"),
	name = "hmap2",
scale_rows = F, scale_columns = F, width = unit(3, "in"), height = unit(3, "in"),
raster_quality = 10, legend_name = "Relative expression",
color_palette = color_palette,
color_range = 5 * (seq(-1, 1, length.out = length(color_palette) - 1)))
})

pdf(file.path(output_dir, "state_assignment_heatmap.pdf"), width = 9, height = 6)
suppressWarnings({
	draw(top_markers_left(scRNA_plot_data, 5, gene_info, "State") + p + q, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = F)
})

gene_info$State = as.factor(as.character(gene_info$State))
discovery_annotation$State = as.factor(as.character(discovery_annotation$State))
assignment$State = factor(as.character(assignment$State), levels = levels(discovery_annotation$State))

rect = rectangle_annotation_coordinates(gene_info$State, discovery_annotation$State)
rect2 = rectangle_annotation_coordinates(gene_info$State, assignment$State)
 
decorate_heatmap_body("hmap1", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 1))
})
decorate_heatmap_body("hmap2", {
    grid.rect(x = unit(rect2$x, "native"), y = unit(rect2$y, "native"), width = unit(rect2$w, "native"), height = unit(rect2$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 1))
})

tmp = dev.off()

png(file.path(output_dir, "state_assignment_heatmap.png"), width = 9, height = 6, units = "in", res = 300)
suppressWarnings({
	draw(top_markers_left(scRNA_plot_data, 5, gene_info, "State") + p + q , heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = F)
})
decorate_heatmap_body("hmap1", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 1))
})
decorate_heatmap_body("hmap2", {
    grid.rect(x = unit(rect2$x, "native"), y = unit(rect2$y, "native"), width = unit(rect2$w, "native"), height = unit(rect2$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 1))
})

tmp = dev.off()  

####

initial_assignment = read.delim(file.path(output_dir, "initial_state_assignment.txt"))
initial_gene_info = read.delim(file.path(states_dir, "initial_gene_info.txt"))
data = data[match(initial_gene_info$Gene, rownames(data)), match(initial_assignment$ID, colnames(data))]
rownames(data) = initial_gene_info$Gene
data[is.na(data)] = 0

if(calculate_significance_flag)
{
	calculate_recovery_significance(data, initial_gene_info, initial_assignment, output_dir, n_bootstraps = 30)
}

