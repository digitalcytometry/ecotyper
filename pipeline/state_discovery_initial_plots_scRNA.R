suppressPackageStartupMessages({
library(NMF)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("Lung", "TR12", "two_tiered", "Fibroblasts", "2", "State", "Histology", "Tissue")
args = commandArgs(T) 
dataset_type = args[1]
dataset = args[2]
fractions = args[3]
cell_type = args[4]
n_states = args[5]
top_cols = args[6:length(args)]
if(is.na(top_cols[1]))
{
	top_cols = c("State")
}

exp_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type)
input_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type, n_states) 
output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type, n_states) 
dir.create(output_dir, recursive = T, showWarning = F)

raw_data = read.delim(file.path(exp_dir, "expression_top_genes_scaled.txt"))

if(!file.exists(file.path(input_dir, "estim.RData")))
{
	stop(paste("NMF output not found in '",file.path(input_dir, "estim.RData"), "'. Please make sure that cell state discovery ran successfully!"))
}

load(file.path(input_dir, "estim.RData"))

W = estim.r@fit@W
rownames(W) = c(paste0(rownames(raw_data), "__pos"),paste0(rownames(raw_data), "__neg"))
colnames(W) = sprintf("IS%02d", 1:ncol(W))
write.table(W, file.path(output_dir, "W.txt"), sep = "\t")

H = estim.r@fit@H
rownames(H) = sprintf("IS%02d", 1:nrow(H))
write.table(H, file.path(output_dir, "H.txt"), sep = "\t")
H = apply(H, 2, function(x) x / sum(x))
write.table(H, file.path(output_dir, "initial_state_abundances.txt"), sep = "\t")

classes = as.data.frame(predict(estim.r))
classes = data.frame(ID = rownames(classes), State = sprintf("IS%02d", classes[,1]))
classes = classes[order(classes$State),]
write.table(classes, file.path(output_dir, "initial_state_assignment.txt"), sep = "\t", row.names= F)

data = raw_data
data = data[,match(classes$ID, colnames(data))]

clinical = read_clinical(colnames(data), dataset = dataset, dataset_type = "discovery")
top_ann = merge(clinical, classes, by = "ID", all.y = T)
rownames(top_ann) = top_ann$ID

top_ann = top_ann[match(colnames(data), top_ann$ID),]

data = data[match(colnames(data), rownames(data)),]

#write.table(gene_info, file.path(output_dir, "initial_gene_info.txt"), sep = "\t")
#write.table(data, file.path(output_dir, "heatmap_data.txt"), sep = "\t")

color_palette = c("black", "#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4", "#005af2", "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4", "#0049e2", "#0046e0", "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6", "#003ed4", "#003dd3", "#003cd1", "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9", "#0036c7", "#0034c5", "#0033c4", "#0131ba", "#0c2eae", "#132ca1", "#172995", "#192688", "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e", "#191743", "#171539", "#15122e", "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a", "#26290e", "#30360f", "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500", "#8ba400", "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000", "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00", "#d5fc00", "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00", "#e0ff00", "#e2ff00", "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02", "#eeff07", "#f0ff0b", "#f2ff0e", "#f4ff12", "#f5ff14", "#f7ff17")
p <- heatmap_simple(data, top_annotation = top_ann, top_columns = top_cols, column_title = lookup_celltype(cell_type),
scale_rows = T, width = unit(4, "in"), height = unit(4, "in"),
raster_quality = 10, legend_name = "Relative expression",
color_palette = color_palette,
color_range = 1.5 * seq(-1, 1, length.out = length(color_palette) - 1))

pdf(file.path(output_dir, "initial_assignment_heatmap.pdf"), width = 7, height = 6)
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
tmp = dev.off()
png(file.path(output_dir, "initial_assignment_heatmap.png"), width = 7, height = 6, units = "in", res = 200)
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
tmp = dev.off()

