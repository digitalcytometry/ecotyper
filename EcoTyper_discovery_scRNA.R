suppressPackageStartupMessages({
library(config)
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/misc.R") 
source("pipeline/lib/multithreading.R")
})

parser <- ArgumentParser(add_help = F)

arguments = parser$add_argument_group('Arguments')

arguments$add_argument("-c", "--config", type = "character", metavar="<PATH>", 
    help="Path to the config files [required].")
arguments$add_argument("-h", "--help", action='store_true', help="Print help message.")

args <- parser$parse_args()
#print(args)
if(args$h || is.null(args$config))
{
	parser$print_help()
	quit()
}

config_file = abspath(args$config)

config <- config::get(file = config_file)
check_discovery_configuration(config)

discovery = config$Input$"Discovery dataset name"
discovery_type = config$Input$"Expression type"
scale_column = config$Input$"Annotation file column to scale by"
additional_columns = config$Input$"Annotation file column(s) to plot"

final_output = config$"Output"$"Output folder"

n_threads = config$"Pipeline settings"$"Number of threads"

nmf_restarts = config$"Pipeline settings"$"Number of NMF restarts"
max_clusters = config$"Pipeline settings"$"Maximum number of states per cell type"
cohpenetic_cutoff = config$"Pipeline settings"$"Cophenetic coefficient cutoff"
skip_steps = config$"Pipeline settings"$"Pipeline steps to skip"

suppressWarnings({
	final_output = abspath(final_output)	
})

#Starting EcoTyper
setwd("pipeline")
start = Sys.time()

if(config$"Pipeline settings"$"Filter non cell type specific genes")
{
	fractions = "Cell_type_specific_genes"
}else{
	fractions = "All_genes"
}

if(!1 %in% skip_steps & config$"Pipeline settings"$"Filter non cell type specific genes")
{
	cat("\nStep 1 (extract cell type specific genes)...\n")

	annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
	for(cell_type in cell_types)
	{
		print(cell_type)
		PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes.R", discovery, fractions, cell_type, scale_column))	
	}
	RunJobQueue()	
	cat("Step 1 (extract cell type specific genes) finished successfully!\n")
	
}else{
	cat("Skipping step 1 (extract cell type specific genes)...\n")
}

if(!2 %in% skip_steps)
{
	cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")

	annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

	for(cell_type in cell_types)
	{
		filter_genes = (fractions == "Cell_type_specific_genes")
		PushToJobQueue(paste("Rscript state_discovery_scRNA_distances.R", discovery, fractions, cell_type, filter_genes, scale_column))	
	}	
	RunJobQueue() 
	
	cat("Step 2 (cell state discovery on correrlation matrices): Running NMF (Warning: This step might take a long time!)...\n")

	for(cell_type in cell_types)
	{		
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}
		for(n_clusters in 2:max_clusters)
		{
			for(restart in 1:nmf_restarts)
			{
				if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, n_clusters, "restarts", restart, "estim.RData")))
				{
					PushToJobQueue(paste("Rscript state_discovery_NMF.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, restart))
				}else{					
					cat(paste0("Warning: Skipping NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, "), as the output file '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData"), "' already exists!\n"))
				}
			} 
		}			
	} 
	RunJobQueue()
		
	cat("Step 2 (cell state discovery on correrlation matrices): Aggregating NMF results...\n")
	for(cell_type in cell_types)
	{	
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}				
		PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery_cross_cor", discovery, fractions, cell_type, max_clusters, nmf_restarts))
	} 
	RunJobQueue()
	cat("Step 2 (cell state discovery on correrlation matrices) finished successfully!\n")
}else{
	cat("Skipping step 2 (cell state discovery on correrlation matrices)...\n")
}	

if(!3 %in% skip_steps)
{
	cat("\nStep 3 (choosing the number of cell states)...\n")
	PushToJobQueue(paste("Rscript state_discovery_rank_selection.R", "discovery_cross_cor", discovery, fractions, max_clusters, cohpenetic_cutoff))
	RunJobQueue()
	cat("Step 3 (choosing the number of cell states) finished successfully!\n")
}else{
	cat("Skipping step 3 (choosing the number of cell states)...\n")
}

if(!4 %in% skip_steps)
{
	cat("\nStep 4 (extracting cell state information)...\n")

	system(paste("cp -f ", config_file, file.path("../EcoTyper", discovery, "config_used.yml")))

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Extracting cell states information for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_initial_plots_scRNA.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 4 (extracting cell state information) finished successfully!\n")
}else{
	cat("\nSkipping step 4 (extracting cell state information)...\n")
}

if(!5 %in% skip_steps)
{
	cat("\nStep 5 (cell state re-discovery in expression matrices)...\n")

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Extracting marker genes for cell states defined in:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_extract_features_scRNA.R", discovery, fractions, cell_type, n_clusters)) 		 
	}	
	RunJobQueue()
	
	cat("\nStep 5 (cell state re-discovery in expression matrices): Running NMF on expression matrix...\n")

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}
		
		n_clusters = key[key[,1] == cell_type, 2]
		for(restart in 1:nmf_restarts)
		{
			if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData")))
			{
				PushToJobQueue(paste("Rscript state_discovery_NMF.R", "discovery", discovery, fractions, cell_type, n_clusters, restart))
			}else{					
				cat(paste0("Warning: Skipping NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, "), as the output file '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData"), "' already exists!\n"))
			}
		} 
	} 
	RunJobQueue()
		
	cat("Step 5 (cell state re-discovery in expression matrices): Aggregating NMF results...\n")
	for(cell_type in key[,1])
	{	
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}				
		PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery", discovery, fractions, cell_type, max_clusters, nmf_restarts))
	} 
	RunJobQueue()
	cat("Step 5 (cell state re-discovery in expression matrices) finished successfully!\n")
}else{
	cat("Skipping step 5 (cell state re-discovery in expression matrices)...\n")
}	

if(!6 %in% skip_steps) 
{
	cat("\nStep 6 (extracting information for re-discovered cell states)...\n")

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Extracting cell states information for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_initial_plots.R", "discovery", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 6 (extracting information for re-discovered cell states) finished successfully!\n")
}else{
	cat("\nSkipping step 6 (extracting information for re-discovered cell states)...\n")
}

if(!7 %in% skip_steps)
{
	cat("\nStep 7 (cell state QC filter)...\n")

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Filtering low-quality cell states for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_first_filter_scRNA.R", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 7 (cell state QC filter) finished successfully!\n")
}else{
	cat("\nSkipping step 7 (cell state QC filter)...\n")
}

if(!8 %in% skip_steps)
{
	cat("\nStep 8 (ecotype discovery)...\n")
	PushToJobQueue(paste("Rscript ecotypes_scRNA.R", discovery, fractions)) 
	RunJobQueue()
	PushToJobQueue(paste("Rscript ecotypes_assign_samples_scRNA.R", discovery, fractions, "State",paste(additional_columns, collapse = " "))) 
	cat("Step 8 (ecotype discovery) finished successfully!\n")
	RunJobQueue()
}else{
	cat("Skipping step 8 (ecotype discovery)...\n")
}

cat("\nCopying EcoTyper results to the output folder!\n")

if(file.exists(final_output) && length(list.files(final_output)) > 0)
{
	old_results_folder = paste0(final_output, format(Sys.time(), " %a %b %d %X %Y"))
	dir.create(old_results_folder, recursive = T, showWarnings = F)
	warning(paste0("The output folder contains files from a previous run. Moving those files to: '", old_results_folder, "'"))	
	system(paste0("mv -f ", final_output, "/* '", old_results_folder, "'"))
}

dir.create(final_output, recursive = T, showWarnings = F)

system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"), final_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.pdf"), final_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.png"), final_output))

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{		
	n_clusters = key[key[,1] == cell_type, 2]	
	ct_output = file.path(final_output, cell_type)
	dir.create(ct_output, recursive = T, showWarnings = F)
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "gene_info.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_abundances.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.pdf"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.png"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_data.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_top_ann.txt"), ct_output))	
}	

ct_output = file.path(final_output, "Ecotypes")
dir.create(ct_output, recursive = T, showWarnings = F)
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotypes.txt"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_assignment.txt"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_abundance.txt"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.pdf"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.png"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.pdf"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.png"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.png"), ct_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.pdf"), ct_output))

end = Sys.time()
cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 0), "\n"))

