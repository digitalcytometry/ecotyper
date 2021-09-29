suppressPackageStartupMessages({
library(doParallel)
library(NMF)
})

args = c("discovery", "scRNA_CRC_Park", "scRNA_specific_genes", "CD4.T.cells", "7", "2")
args = commandArgs(T) 
dataset_type = args[1]
dataset = args[2]
fractions = args[3]
cell_type = args[4]
n_clusters = as.integer(as.character(args[5]))
restart = as.integer(as.character(args[6]))

input_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type) 

if(!file.exists(file.path(input_dir, "expression_top_genes_scaled.txt")))
{
	stop(paste0("No input data for cell type: ", cell_type))
}

output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type, n_clusters, "restarts", restart) 
dir.create(output_dir, recursive = T, showWarning = F)

raw_data = read.delim(file.path(input_dir, "expression_top_genes_scaled.txt"))
data = posneg(as.matrix(raw_data))

cat(paste0("Running NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, ")...\n"))
seed = 1234 + restart
estim.r <- nmf(data, n_clusters, nrun = 1, method = "brunet", seed = seed, .opt='P1')
save(estim.r, file = file.path(output_dir, "estim.RData")) 


