suppressPackageStartupMessages({
library(NMF)
})

args = c("discovery", "MyDiscovery", "Carcinoma", "Epithelial.cells", "20", "5")
args = commandArgs(T) 

dataset_type = args[1]
dataset = args[2]
fractions = args[3]
cell_type = args[4]
max_n_clusters = as.integer(as.character(args[5]))
max_restarts = as.integer(as.character(args[6]))

base_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type, cell_type) 
if(!file.exists(file.path(base_dir, "expression_top_genes_scaled.txt")))
{
	stop(paste0("No input data for cell type: ", cell_type))
}

cat(paste0("Combining NMF restarts for '", cell_type, "'...\n"))
for(n_clusters in 2:max_n_clusters)
{
	#print(n_clusters)
	output_dir = file.path(base_dir, n_clusters) 

	all_coef = NULL
	used_reps = list()
	l = list()
	for(restart in 1:max_restarts)
	{
		input_dir = file.path(base_dir, n_clusters, "restarts", restart) 
		
		if(!file.exists(file.path(input_dir, "estim.RData")))
		{
			cat(paste0("Warning: NMF output file '", file.path(input_dir, "estim.RData"), "' is missing! Skipping it...\n"))
			next
		}
		used_reps = c(used_reps, restart)
		load(file.path(input_dir, "estim.RData"))
		l[[length(l) + 1]] = estim.r
		estim.r = NULL 
	}
	if(length(l) == 0)
	{
		next
	}
	write.table(t(as.data.frame(used_reps)), file.path(output_dir, "used.txt"), sep = "\t")
	
	suppressWarnings({
	estim.r = NMF:::NMFfitX(l, .merge = T)	
	})	
	save(estim.r, file = file.path(output_dir, "estim.RData"))
	
	coph = cophcor(estim.r)
	sis = summary(silhouette(estim.r, what = 'samples'))
	sig = summary(silhouette(estim.r, what = 'features'))
	sichc = summary(silhouette(estim.r, what = 'chc'))
	sicons = summary(silhouette(estim.r, what = 'consensus'))
	di = dispersion(estim.r) 
	all_coef = rbind(all_coef, data.frame(Cophenetic = coph, Silhouette_s = sis[4][[1]], 
		Silhouette_g = sig[4][[1]],
		Silhouette_chc = sichc[4][[1]],
		Silhouette_con = sicons[4][[1]],
		Dispersion = di, n_clusters = n_clusters))
	write.table(all_coef, file.path(output_dir, paste0("rank_data.txt")), sep = "\t", row.names = F)
}

