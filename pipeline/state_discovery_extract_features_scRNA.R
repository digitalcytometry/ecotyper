suppressPackageStartupMessages({
library(data.table)
source("lib/misc.R")
})

args = c("scRNA_CRC_Park", "scRNA_specific_genes", "CD4.T.cells", 7) 
args = commandArgs(T)
dataset = args[1] 
fractions = args[2]
cell_type = args[3]
n_states = args[4] 

dataset_type = "discovery"
 
input_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery_cross_cor", cell_type) 

log_data = read.delim(file.path(input_dir, "expression_full_matrix_log2.txt"))
scaled_data = read.delim(file.path(input_dir, "expression_top_genes_scaled_filt.txt"))
log_data = log_data[match(rownames(scaled_data), rownames(log_data)),]

if(nrow(scaled_data) < 50)
{
	stop(paste0("Only ", nrow(scaled_data), " genes were imputed for cell type: ", cell_type, ". At least 50 are required. Skipping this cell type!"))
}

output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type) 
dir.create(output_dir, recursive = T, showWarning = F)

annotation = read.delim(file.path(input_dir, n_states, "initial_state_assignment.txt"))
annotation = annotation[match(colnames(scaled_data), annotation$ID),]

FCs = get_FC_for_clusters(scaled_data, annotation$State)
splits = split(FCs, FCs$State)
FCs = do.call(rbind, lapply(splits, function(spl) {
	spl$Rank = rank(-spl$MaxFC)
	spl
}))
FCs = FCs[order(FCs$Rank, -FCs$MaxFC),]
write.table(FCs, file.path(output_dir, "full_FCs.txt"), sep = "\t", row.names = T)
FCs = FCs[FCs$MaxFC > 0,]
FCs = FCs[1:min(nrow(FCs), 1000),]
write.table(FCs, file.path(output_dir, "top_FCs.txt"), sep = "\t", row.names = T)

write.table(log_data, file.path(output_dir, "expression_full_matrix_log2.txt"), sep = "\t", row.names = T)
log_data = log_data[match(FCs$Gene, rownames(log_data)),]
write.table(log_data, file.path(output_dir, "expression_top_genes_log2.txt"), sep = "\t", row.names = T)

data = scaled_data[match(FCs$Gene, rownames(scaled_data)),]
write.table(data, file.path(output_dir, "expression_top_genes_scaled.txt"), sep = "\t", row.names = T)
write.table(scaled_data, file.path(output_dir, "expression_full_matrix_scaled.txt"), sep = "\t", row.names = T)
