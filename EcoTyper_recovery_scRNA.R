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
    help="Path to a tab-delimited file containing the input scRNA-seq expression matrix, with gene names on the first column and cell ids as column names [required].")
arguments$add_argument("-a", "--annotation", type = "character",  metavar="<PATH>",
    help="Path to a tab-delimited annotation file containing the annotation of cells in the input matrix. This file should contain at least two columns, 'ID' with the same values as the columns of the expression matrix, and 'CellType' (case sensitive) which contains the cell type for each cell. These values are limited to the set of cell types analyzed in the discovery dataset. If the argument '-d' is set to 'Carcinoma', then the accepted values for column 'CellType' are: 'B.cells', 'CD4.T.cells', 'CD8.T.cells', 'Dendritic.cells', 'Endothelial.cells', 'Epithelial.cells', 'Fibroblasts', 'Mast.cells', 'Monocytes.and.Macrophages', 'NK.cells', 'PCs' and 'PMNs'. If the argument '-d' is set to 'Lymphoma', then the accepted values for column 'CellType' are: 'B.cells', 'Plasma.cells', 'T.cells.CD8', 'T.cells.CD4', 'T.cells.follicular.helper', 'Tregs', 'NK.cells', 'Monocytes.and.Macrophages', 'Dendritic.cells', 'Mast.cells', 'Neutrophils', 'Fibroblasts', 'Endothelial.cells'. All other values will be ignored for these two cases. Additionally, this file can contain any number of columns, that can be used for plotting color bars in the output heatmaps (see argument '-c'). [required]")
arguments$add_argument("-c", "--columns",  type = "character",  metavar="<character>", default="NULL", 
    help="A comma-spearated list of column names from the annotation file to be plotted as color bar in the output heatmaps. [default: '%(default)s']")
arguments$add_argument("-z", "--z-score",  type = "character",  metavar="<bool>", default="FALSE", 
    help="A flag indicating whether the significance quantification procedure should be run. Note that this procedure might be slow, as the NMF model is applied 30 times on the same dataset. [default: '%(default)s']")
arguments$add_argument("-s", "--subsample",  type = "integer",  metavar="<integer>", default=-1, 
    help="An integer specifying the number of cells each cell type will be downsampled to. For values <50, no downsampling will be performed. [default: '%(default)s' (no downsampling)]")
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
z_flag = as.logical(args$z)
n_threads = as.numeric(args$threads)
subsample = as.integer(args$subsample)

if(!file.exists(input_mat))
{
	stop("Error: Path to the input expression matrix does not exist!")
}

if(!discovery %in% c("Carcinoma", "Lymphoma"))
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
		if(config$"Pipeline settings"$"Filter non cell type specific genes")
		{
			fractions = "Cell_type_specific_genes"
		}else{
			fractions = "All_genes"
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

if(!file.exists(annotation_path))
{
	stop("Error: Path to the input annotation file does not exist!")
}else{
	annotation = read.delim(annotation_path)	
	if(!all(c("ID", "CellType") %in% colnames(annotation)))
	{
		stop("Error: The annotation file provided, does not contain the columns 'ID' and 'CellType.")
	}

	if("Sample" %in% colnames(annotation))
	{		
		ecotypes = read.delim(file.path("EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotypes.txt"))
		disc_cell_types = table(ecotypes$CellType)
		rec_cell_types= table(annotation$CellType)
		if(length(rec_cell_types) * 2 >= length(disc_cell_types))
		{
			cat("The annotation file contains column 'Sample', and more than half of the cell types are present in the recovery dataset. Will perform ecotype recovery.\n")
			recover_ecotypes = T
		}else{
			cat("The annotation file contains column 'Sample', but less than half of the cell types are present in the recovery dataset. Will NOT perform ecotype recovery.\n")
			recover_ecotypes = F
		}

	}else{
		recover_ecotypes = F
		cat("The annotation file does not contain column 'Sample'. Will NOT perform ecotype recovery.\n")
	}

	if(columns != "NULL")
	{
		additional_columns = strsplit(columns, ",")[[1]]
	
		if(!all(additional_columns %in% colnames(annotation)))
		{
			stop(paste0("The following columns are missing from the annotation file provided: ", "'", 
				paste(additional_columns[!additional_columns %in% colnames(annotation)], collapse = "'"), "'"))			
		} 
	}else{
		additional_columns = c()
	}
}

recovery = gsub(".tsv$", "", gsub(".txt$", "", basename(input_mat)))

input_dir = file.path("datasets", "scRNA", recovery)
dir.create(input_dir, recursive = T, showWarning = F)

dir.create(file.path(args$output, recovery), recursive = T, showWarnings = F)
final_output = normalizePath(file.path(args$output, recovery))

PushToJobQueue(paste0("ln -sf ", input_mat, " ", file.path(input_dir, "data.txt")))
RunJobQueue()
PushToJobQueue(paste0("ln -sf ", annotation_path, " ", file.path(input_dir, "annotation.txt")))
RunJobQueue()

start = Sys.time()
cur_dir = getwd()
setwd("pipeline")

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
key = key[key[,1] %in% annotation$CellType,]
if(nrow(key) == 0)
{
	stop("Error: The column 'CellType' of the annotation file contains invalid values!")
}
for(cell_type in key[,1])
{
    n_states = key[key[,1] == cell_type, 2]
    PushToJobQueue((paste("Rscript state_recovery_scRNA.R", discovery, fractions, cell_type, n_states, recovery, z_flag, subsample, paste(additional_columns, collapse = " ")))) 
}   
RunJobQueue()

if(recover_ecotypes)
{
	PushToJobQueue((paste("Rscript ecotypes_recovery_scRNA.R", discovery, fractions, recovery))) 
	RunJobQueue()
}

cat("\nCopying EcoTyper results to the output folder...\n")

if(file.exists(final_output) && length(list.files(final_output)) > 0)
{
	old_results_folder = paste0(final_output, format(Sys.time(), " %a %b %d %X %Y"))
	dir.create(old_results_folder, recursive = T, showWarnings = F)
	warning(paste0("The output folder contains files from a previous run. Moving those files to: '", old_results_folder, "'"))	
	system(paste0("mv -f ", final_output, "/* '", old_results_folder, "'"))
}

dir.create(final_output, recursive = T, showWarnings = F)

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
key = key[key[,1] %in% annotation$CellType,] 
for(cell_type in key[,1])
{		
	n_clusters = key[key[,1] == cell_type, 2]	
	ct_output = file.path(final_output, cell_type)
	dir.create(ct_output, recursive = T, showWarnings = F)
		
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", "scRNA", recovery, cell_type, n_clusters, "state_assignment.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", "scRNA", recovery, cell_type, n_clusters, "state_assignment_heatmap.pdf"), ct_output))		
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", "scRNA", recovery, cell_type, n_clusters, "state_assignment_heatmap.png"), ct_output))		
	if(z_flag)
	{
		system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", "scRNA", recovery, cell_type, n_clusters, "recovery_z_scores.txt"), ct_output))	
	}
}	

if(recover_ecotypes)
{
	ct_output = file.path(final_output, "Ecotypes")
	dir.create(ct_output, recursive = T, showWarnings = F)
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_assignment.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_abundance.txt"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.pdf"), ct_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.png"), ct_output))
}

end = Sys.time()
cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 1), "\n"))

