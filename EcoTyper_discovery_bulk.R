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
CSx_username = config$Input$"CIBERSORTx username"
CSx_token = config$Input$"CIBERSORTx token"
scale_column = config$Input$"Annotation file column to scale by"
additional_columns = config$Input$"Annotation file column(s) to plot"
fractions = config$Input$"Cell type fractions"

final_output = config$"Output"$"Output folder"

n_threads = config$"Pipeline settings"$"Number of threads"
nmf_restarts = config$"Pipeline settings"$"Number of NMF restarts"
max_clusters = config$"Pipeline settings"$"Maximum number of states per cell type"
cohpenetic_cutoff = config$"Pipeline settings"$"Cophenetic coefficient cutoff"
skip_steps = config$"Pipeline settings"$"Pipeline steps to skip"
CSx_singularity_path_fractions = config$"Pipeline settings"$"CIBERSORTx fractions Singularity path"
CSx_singularity_path_hires = config$"Pipeline settings"$"CIBERSORTx hires Singularity path"
min_states = config$"Pipeline settings"$"Minimum number of states in ecotypes"

suppressWarnings({
	final_output = abspath(final_output)	
	if(!fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
	{
		fractions_path = abspath(fractions)
	}
})

#Starting EcoTyper
setwd("pipeline")
start = Sys.time()

if(!fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
{
	fractions = "Custom"
}

if(!1 %in% skip_steps)
{
	cat("\nStep 1 (cell type fraction estimation): Running CIBERSORTxFractions on the discovery dataset...\n")
	if(fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
	{
		
		if(discovery_type == "RNA-seq")
		{
			PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "TR4", "no_batch", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))	
			RunJobQueue()
			PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "LM22", "B_mode", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))
			RunJobQueue() 
			PushToJobQueue(paste("Rscript csx_fractions_two_tiered.R", "discovery", discovery, "TR4", "no_batch", "LM22", "B_mode", fractions))		
			RunJobQueue()
		}else{
			if(discovery_type == "Affymetrix")
			{
				PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "TR4", "B_mode", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))
				RunJobQueue()
				PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "LM22", "no_batch", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))
				RunJobQueue()
				PushToJobQueue(paste("Rscript csx_fractions_two_tiered.R", "discovery", discovery, "TR4", "B_mode", "LM22", "no_batch", fractions))
				RunJobQueue()
			}else{
				PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "LM22", "B_mode", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))
				RunJobQueue()
				PushToJobQueue(paste("Rscript csx_fractions.R", "discovery", discovery, "TR4", "B_mode", CSx_username, CSx_token, paste0("'", CSx_singularity_path_fractions, "'")))			
				RunJobQueue()
				PushToJobQueue(paste("Rscript csx_fractions_two_tiered.R", "discovery", discovery, "TR4", "B_mode", "LM22", "B_mode", fractions))
				RunJobQueue()
			}
		}
		
	}else{
		cat("Step 1 (cell type fraction estimation): Loading user-provided cell type fractions...\n")
				
		dir.create(file.path("../CIBERSORTx/fractions/discovery", discovery, fractions), recursive = T, showWarnings = F)
		system(paste("cp -f", paste0("'", fractions_path, "'"), file.path("../CIBERSORTx/fractions/discovery", discovery, fractions, "CIBERSORTx_Adjusted.txt")))	
	}
	cat("Step 1 (cell type fraction estimation) finished successfully!\n")
}else{
	cat("Skipping step 1 (cell type fraction estimation)...\n")
}

if(!2 %in% skip_steps)
{
	cat("\nStep 2 (cell type expression purification): Running CIBERSORTxHiRes...\n")
	PushToJobQueue(paste("Rscript csx_hires_scheduler.R", "discovery", discovery, fractions, CSx_username, CSx_token, n_threads, paste0("'", CSx_singularity_path_hires, "'")))
	RunJobQueue()

	cat("Step 2 (cell type expression purification): Aggregating CIBERSORTxHiRes results...\n")
	PushToJobQueue(paste("Rscript csx_hires_aggregate_worker_results.R", "discovery", discovery, fractions))
	RunJobQueue()
	
	cat("Step 2 (cell type expression purification) finished successfully!\n")
}else{
	cat("Skipping step 2 (cell type expression purification)...\n")
}

if(!3 %in% skip_steps)
{
	cat("\nStep 3 (cell state discovery): Preparing the NMF input...\n")

	classes_path = file.path(file.path("../CIBERSORTx/hires", discovery, fractions, "classes.txt"))
	classes = read.delim(classes_path)
	for(cell_type in colnames(classes))
	{
		PushToJobQueue(paste("Rscript state_discovery_extract_features.R", discovery, fractions, cell_type, scale_column))
	}
	RunJobQueue() 
	
	cat("Step 3 (cell state discovery): Running NMF (Warning: This step might take a long time!)...\n")
	classes_path = file.path(file.path("../CIBERSORTx/hires", discovery, fractions, "classes.txt"))
	classes = read.delim(classes_path) 
	for(cell_type in colnames(classes))
	{		
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}
		for(n_clusters in 2:max_clusters)
		{
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
	} 
	RunJobQueue()
	
	classes_path = file.path(file.path("../CIBERSORTx/hires", discovery, fractions, "classes.txt"))
	classes = read.delim(classes_path)
	cat("Step 3 (cell state discovery): Aggregating NMF results...\n")
	for(cell_type in colnames(classes))
	{	
		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
		{			
			next
		}				
		PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery", discovery, fractions, cell_type, max_clusters, nmf_restarts))
	} 
	RunJobQueue()
	cat("Step 3 (cell state discovery) finished successfully!\n")
}else{
	cat("Skipping step 3 (cell state discovery)...\n")
}	

if(!4 %in% skip_steps)
{
	cat("\nStep 4 (choosing the number of cell states)...\n")
	PushToJobQueue(paste("Rscript state_discovery_rank_selection.R", "discovery", discovery, fractions, max_clusters, cohpenetic_cutoff))
	RunJobQueue()
	cat("Step 4 (choosing the number of cell states) finished successfully!\n")
}else{
	cat("Skipping step 4 (choosing the number of cell states)...\n")
}

if(!5 %in% skip_steps)
{
	cat("\nStep 5 (extracting cell state information)...\n")

	system(paste("cp -f ", paste0("'", config_file, "'"), file.path("../EcoTyper", discovery, "config_used.yml")))

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Extracting cell states information for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_initial_plots.R", "discovery", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 5 (extracting cell state information) finished successfully!\n")
}else{
	cat("\nSkipping step 5 (extracting cell state information)...\n")
}

if(!6 %in% skip_steps)
{
	cat("\nStep 6 (cell state QC filter)...\n")

	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Filtering low-quality cell states for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_first_filter.R", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 6 (cell state QC filter) finished successfully!\n")
}else{
	cat("\nSkipping step 6 (cell state QC filter)...\n")
}

if(!7 %in% skip_steps)
{
	cat("\nStep 7 (advanced cell state QC filter). Warning: This filter is not recommended, unless your discovery dataset is composed of samples that can be confounded by serious biological or technical differences (e.g. it contains multiple tumor types).\n")

	PushToJobQueue(paste("Rscript state_discovery_calculate_dropout_score.R", discovery, fractions))
	RunJobQueue()
	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
	for(cell_type in key[,1])
	{	
		cat(paste("Filtering low-quality cell states for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_second_filter.R", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 7 (advanced cell state QC filter) finished successfully!\n")
}else{
	cat("\nSkipping step 7 (advanced cell state QC filter)...\n")
}

if(!8 %in% skip_steps)
{
	cat("\nStep 8 (ecotype discovery)...\n")
	PushToJobQueue(paste("Rscript ecotypes.R", discovery, fractions, min_states)) 
	RunJobQueue()
	PushToJobQueue(paste("Rscript ecotypes_assign_samples.R", discovery, fractions, "Ecotype", paste(additional_columns, collapse = " "))) 
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
	system(paste0("mv -f ", paste0("'", final_output, "'"), "/* '", old_results_folder, "'"))
}

dir.create(final_output, recursive = T, showWarnings = F)

system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"), paste0("'", final_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.pdf"), paste0("'", final_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.png"), paste0("'", final_output, "'")))

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{		
	n_clusters = key[key[,1] == cell_type, 2]	
	ct_output = file.path(final_output, cell_type)
	dir.create(ct_output, recursive = T, showWarnings = F)
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "gene_info.txt"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_abundances.txt"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment.txt"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.pdf"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.png"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_data.txt"), paste0("'", ct_output, "'")))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_top_ann.txt"), paste0("'", ct_output, "'")))
}	

ct_output = file.path(final_output, "Ecotypes")
dir.create(ct_output, recursive = T, showWarnings = F)
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotypes.txt"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_assignment.txt"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_abundance.txt"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.pdf"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.png"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.pdf"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.png"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.png"), paste0("'", ct_output, "'")))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.pdf"), paste0("'", ct_output, "'")))

end = Sys.time()
cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", paste0("'", final_output, "'"), "'.\nRun time: ", format(end - start, digits = 1), "\n"))

