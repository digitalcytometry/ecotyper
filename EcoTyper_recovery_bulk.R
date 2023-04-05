suppressPackageStartupMessages({
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/multithreading.R")
})

parser <- ArgumentParser(add_help = F)

arguments = parser$add_argument_group('Arguments')

arguments$add_argument("-d", "--discovery", type="character", default="Carcinoma", 
    help=paste0("The name of the discovery dataset used to define cell states and ecotypes. Accepted values: ",
    "'Carcinoma' will recover the cell states and ecotypes defined across carcinomas, as described in the EcoTyper carcinoma paper, ",
    "'Lymphoma' will recover the cell states and ecotypes defined in diffuse large B cell lymphoma (DLBCL), as described in the EcoTyper lymphoma paper, ",
    "'<MyDiscovery>' the value used in the field 'Discovery dataset name' of the config file used for running EcoTyper discovery ('EcoTyper_discovery.R') script. ",
    "[default: '%(default)s']"),
    metavar="<character>") 
arguments$add_argument("-m", "--matrix", type = "character", metavar="<PATH>",
    help="Path to a tab-delimited file containing the input bulk tissue expression matrix, with gene names on the first column and sample ids as column names [required].")
arguments$add_argument("-a", "--annotation", type = "character",  metavar="<PATH>", default="NULL",
    help="Path to a tab-delimited annotation file containing the annotation of samples in the input matrix. This file has to contain in column 'ID' the same ids used as column names in the input matrix, and any number of additional columns. The additional columns can be plotted as color bars in the output heatmaps. [default: '%(default)s']")
arguments$add_argument("-c", "--columns",  type = "character",  metavar="<character>", default="NULL", 
    help="A comma-spearated list of column names from the annotation file to be plotted as color bar in the output heatmaps. [default: '%(default)s']")
arguments$add_argument("-t", "--threads",  type = "integer",  metavar="<integer>", default=10, 
    help="Number of threads. [default: '%(default)s']")
arguments$add_argument("-o", "--output", type = "character",  metavar="<PATH>", default="RecoveryOutput",
    help="Output directory path. [default: '%(default)s']")
arguments$add_argument("-h", "--help", action='store_true', help="Print help message.")

args <- parser$parse_args()
#print(args)
if(is.null(args$matrix))
{
	parser$print_help()
	quit()
}

input_mat = normalizePath(args$matrix)
discovery = args$discovery
annotation_path = normalizePath(args$annotation)
columns = args$columns
n_threads = args$threads

if(!file.exists(input_mat))
{
	stop("Error: Path to the input expression matrix does not exist!")
}

if(! discovery %in% c("Carcinoma", "Lymphoma"))
{
	config_file = file.path("EcoTyper", discovery, "config_used.yml")	
	if(!file.exists(config_file))
	{
		stop("Error: Cannot read the config file used for the discovery of cell states and ecotypes. It should be in the following path: '", config_file, "'. Please make sure that the '--discovery (-d)' argument provided is correct!")
	}
	config <- config::get(file = config_file)
	
	fractions = config$"Input"$"Cell type fractions" 
	if(is.null(fractions))
	{
		if(config$"Pipeline settings"$"Filter genes" == "cell type specific")
		{
			fractions = "Cell_type_specific_genes"
		}else{
			if(config$"Pipeline settings"$"Filter genes" == "no filter")
			{
				fractions = "All_genes"
			}else{
				n_genes = as.integer(as.numeric(config$"Pipeline settings"$"Filter genes"))				
				fractions = paste0("Top_", n_genes)
			}
		}
	}else{
		if(!fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
		{
			fractions = "Custom"
		}
	}
}else{
	fractions = paste0(discovery, "_Fractions") 
}

if(annotation_path != "NULL" && columns != "NULL")
{
	if(!file.exists(annotation_path))
	{
		stop("Error: Path to the input annotation file does not exist!")
	}else{
		annotation = read.delim(annotation_path)
		#print(head(annotation))
		if(!"ID" %in% colnames(annotation))
		{
			stop("Error: The annotation file provided, does not contain the column 'ID'.")
		}

		additional_columns = strsplit(columns, ",")[[1]]
		
		if(!all(additional_columns %in% colnames(annotation)))
		{
			stop(paste0("The following columns are missing from the annotation file provided: ", "'", 
				paste(additional_columns[!additional_columns %in% colnames(annotation)], collapse = "'"), "'"))			
		} 
	}
}else{
	additional_columns = c()
}

recovery = gsub(".tsv$", "", gsub(".txt$", "", basename(input_mat)))

input_dir = file.path("datasets", "bulk", recovery)
dir.create(input_dir, recursive = T, showWarning = F)

dir.create(file.path(args$output, recovery), recursive = T, showWarnings = F)
final_output = normalizePath(file.path(args$output, recovery))

system(paste0("ln -sf '", input_mat, "' '", file.path(input_dir, "data.txt"), "'"))
system(paste0("ln -sf '", annotation_path, "' '", file.path(input_dir, "annotation.txt"), "'"))

start = Sys.time()
cur_dir = getwd()
setwd("pipeline")

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{
    n_states = key[key[,1] == cell_type, 2]
    PushToJobQueue((paste("Rscript state_recovery_bulk.R", discovery, fractions, cell_type, n_states, recovery, "FALSE", paste(additional_columns, collapse = " ")))) 
}   
RunJobQueue()

PushToJobQueue((paste("Rscript ecotypes_recovery.R", discovery, fractions, recovery, paste(additional_columns, collapse = " ")))) 
RunJobQueue()

cat("\nCopying EcoTyper results to the output folder!\n")

if(file.exists(final_output) && length(list.files(final_output)) > 0)
{
	old_results_folder = paste0(final_output, format(Sys.time(), " %a %b %d %X %Y"))
	dir.create(old_results_folder, recursive = T, showWarnings = F)
	warning(paste0("The output folder contains files from a previous run. Moving those files to: '", old_results_folder, "'"))	
	system(paste0("mv -f '", final_output, "'/* '", old_results_folder, "'"))

}

dir.create(final_output, recursive = T, showWarnings = F)

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{		
	n_clusters = key[key[,1] == cell_type, 2]	
	ct_output = file.path(final_output, cell_type)
	dir.create(ct_output, recursive = T, showWarnings = F)

	system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_abundances.txt"), "' '", ct_output, "'"))
    system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment.txt"), "' '", ct_output, "'"))
    system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment_heatmap.pdf"), "' '", ct_output, "'")) 
    system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment_heatmap.png"), "' '", ct_output, "'")) 
    system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "heatmap_data.txt"), "' '", ct_output, "'"))    
    system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "heatmap_top_ann.txt"), "' '", ct_output, "'"))    
	
}	

ct_output = file.path(final_output, "Ecotypes")
dir.create(ct_output, recursive = T, showWarnings = F)
system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_assignment.txt"), "' '", ct_output, "'"))
system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_abundance.txt"), "' '", ct_output, "'"))
system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.pdf"), "' '", ct_output, "'"))
system(paste0("cp -f '", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.png"), "' '", ct_output, "'"))

end = Sys.time()
cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 1), "\n"))

