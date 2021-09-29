suppressPackageStartupMessages({
library(NMF)
library(data.table)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("Lung", "TR12", "two_tiered", "Fibroblasts", "5", "State", "Histology", "Tissue")
args = commandArgs(T) 
dataset = args[1]
fractions = args[2]
cell_type = args[3]
n_states = args[4]
top_cols = args[5:length(args)]

if(is.na(top_cols[1]))
{
	top_cols = c("State")
}
top_cols = unique(c(top_cols, "State"))

csx_dir = file.path("../CIBERSORTx/hires", dataset, fractions)
exp_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type)
input_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type, n_states) 
output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type, n_states) 
dir.create(output_dir, recursive = T, showWarning = F)

raw_data = read.delim(file.path(exp_dir, "expression_top_genes_scaled.txt"))

if(!file.exists(file.path(input_dir, "estim.RData")))
{
	stop(paste("NMF output not found in '",file.path(input_dir, "estim.RData"), "'. Please make sure that cell state discovery ran successfully!"))	
}
#print(n_states) 

csx_data = read.delim(file.path(exp_dir, "expression_full_matrix_log2.txt"))
#rownames(csx_data) = csx_data[,1]
#csx_data = csx_data[,-1]
csx_data = csx_data[match(rownames(raw_data), rownames(csx_data)),]
gene_filt = apply(csx_data, 1, function(x) rev(sort(table(unlist(x))))[1] > 0.5 * ncol(csx_data))
gene_filt = gene_filt[!gene_filt]

gene_info = read.delim(file.path(output_dir, "initial_gene_info.txt"))
gene_info = gene_info[rownames(gene_info) %in% names(gene_filt),]

classes = read.delim(file.path(output_dir, "initial_state_assignment.txt"))
W = read.delim(file.path(output_dir, "W.txt"))
agg = apply(W, 2, function(x) sum(x[grepl("_neg", rownames(W))]) / sum(x[grepl("_pos", rownames(W))]))
agg = data.frame(InitialState = names(agg), PosRatio = agg)

tb = as.data.frame(table(gene_info$State))
colnames(tb) = c("State", "NumberMarkers")
tb = merge(tb, agg, by.x = "State", by.y = "InitialState", all.y = T)
tb$NumberMarkers = ifelse(is.na(tb$NumberMarkers), 0, tb$NumberMarkers)
write.table(tb, file.path(output_dir, "first_filter_data.txt"), sep = "\t", row.names = F)

to_keep = as.character(tb[tb$NumberMarkers > 10, "State"])
if(length(to_keep) < 2)
{
	cat(paste0("Warning: There would be less than two states left for '", cell_type, "', after filtering states with less than 10 marker genes! This indicates that the Cophenetic coefficient cutoff used at step 4 (choosing the number of cell states) might be set too high. You could try to lower it and re-run EcoTyper from step 4. Alternatively, the discovery dataset might be too small, or the cell type sepcific profiles contain too few genes. Not applying the filter for this cell type at the moment.\n"))
	to_keep = as.character(tb[,"State"])
}

to_keep2 = as.character(tb[(tb$PosRatio < 1), "State"])
if(length(intersect(to_keep2, to_keep)) < 2)
{
	cat(paste0("Warning: There would be less than two states left for '", cell_type, "', after applying the AFI filter! This indicates that the Cophenetic coefficient cutoff used at step 4 (choosing the number of cell states) might be set too high. You could try to lower it and re-run EcoTyper from step 4. Alternatively, the discovery dataset might be too small, or the cell type sepcific profiles contain too few genes. Not applying the filter for this cell type at the moment.\n"))	
}else{
	to_keep = intersect(to_keep2, to_keep)
}

classes = classes[classes$State %in% to_keep,]
gene_info = gene_info[gene_info$State %in% to_keep,]

data = raw_data
data = data[,match(classes$ID, colnames(data))]
data = data[match(gene_info$Gene, rownames(data)),]
classes$State = as.factor(as.character(classes$State))

agg_data = collapse_samples(data, classes$State)
h = hclust(dist(t(agg_data)))
mapping = sprintf("S%02d", 1:ncol(agg_data))
names(mapping) = colnames(agg_data)[h$ord]

tmp = as.data.frame(mapping)
tmp = data.frame(State = tmp[,1], InitialState = rownames(tmp))
write.table(tmp, file.path(output_dir, "mapping_to_initial_states.txt"), sep = "\t", row.names = F)

H = read.delim(file.path(output_dir, "initial_state_abundances.txt"))
H = H[match(tmp$InitialState, rownames(H)),]
rownames(H) = tmp$State
write.table(H, file.path(output_dir, "state_abundances.txt"), sep = "\t")

classes$InitialState = as.factor(as.character(classes$State))
classes$State = as.factor(as.character(mapping[as.character(classes$State)]))
classes = classes[order(classes$State),]
write.table(classes, file.path(output_dir, "state_assignment.txt"), sep = "\t")

data = raw_data
data = data[,match(classes$ID, colnames(data))]

clinical = read_clinical(colnames(data), dataset = dataset, dataset_type = "discovery")
top_ann = merge(clinical, classes, by = "ID", all.y = T)
rownames(top_ann) = top_ann$ID 

top_ann = top_ann[match(colnames(data), top_ann$ID),]
write.table(top_ann, file.path(output_dir, "heatmap_top_ann.txt"), sep = "\t")

gene_info$InitialState = gene_info$State
gene_info$State = as.factor(mapping[as.character(gene_info$State)])
gene_info = gene_info[order(gene_info$State),]
data = data[match(gene_info$Gene, rownames(data)),]

write.table(gene_info, file.path(output_dir, "gene_info.txt"), sep = "\t")
write.table(data, file.path(output_dir, "heatmap_data.txt"), sep = "\t")

color_palette = c("black", "#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4", "#005af2", "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4", "#0049e2", "#0046e0", "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6", "#003ed4", "#003dd3", "#003cd1", "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9", "#0036c7", "#0034c5", "#0033c4", "#0131ba", "#0c2eae", "#132ca1", "#172995", "#192688", "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e", "#191743", "#171539", "#15122e", "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a", "#26290e", "#30360f", "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500", "#8ba400", "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000", "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00", "#d5fc00", "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00", "#e0ff00", "#e2ff00", "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02", "#eeff07", "#f0ff0b", "#f2ff0e", "#f4ff12", "#f5ff14", "#f7ff17")
p <- heatmap_simple(data, top_annotation = top_ann, top_columns = top_cols, column_title = lookup_celltype(cell_type),
scale_rows = T, width = unit(4, "in"), height = unit(3, "in"),
raster_quality = 10, legend_name = "Relative expression",
color_palette = color_palette, name = "hmap1",
color_range = 1.5 * seq(-1, 1, length.out = length(color_palette) - 1))

pdf(file.path(output_dir, "state_assignment_heatmap.pdf"), width = 7, height = 6)
draw(top_markers_left(data, 5, gene_info, "State") + p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
rect = rectangle_annotation_coordinates(gene_info$State, top_ann$State)
decorate_heatmap_body("hmap1", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 1))
})
tmp = dev.off()

png(file.path(output_dir, "state_assignment_heatmap.png"), width = 7, height = 6, units = "in", res = 200)
draw(top_markers_left(data, 5, gene_info, "State") + p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
decorate_heatmap_body("hmap1", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 1))
})
tmp = dev.off()
