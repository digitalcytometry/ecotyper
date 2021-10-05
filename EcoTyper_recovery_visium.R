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

if(args$h || is.null(args$config))
{
	parser$print_help()
	quit()
}

config_file = abspath(args$config)

config <- config::get(file = config_file)

discovery = config$Input$"Discovery dataset name"
recovery = config$Input$"Recovery dataset name"
input_path = config$Input$"Input Visium directory"
fractions = config$Input$"Recovery cell type fractions"
coo = config$Input$"Malignant cell of origin"
CSx_username = config$Input$"CIBERSORTx username"
CSx_token = config$Input$"CIBERSORTx token"
n_threads = config$"Pipeline settings"$"Number of threads"

final_output = config$"Output"$"Output folder"

suppressWarnings({
	input_path = abspath(input_path)	
	final_output = abspath(file.path(final_output, recovery))	
})

suppressWarnings({
	fractions_path = abspath(fractions)	
	})

if(!discovery %in% c("Carcinoma", "Lymphoma"))
{
	discovery_config_file = file.path("EcoTyper", discovery, "config_used.yml")	
	if(!file.exists(discovery_config_file))
	{
		stop("Error: Cannot read the configuration file used for the discovery of cell states and ecotypes. It should be in the following path: '", config_file, "'. Please make sure that the '--discovery (-d)' argument provided is correct!")
	}
	config <- config::get(file = discovery_config_file)
	discovery_fractions = config$"Input"$"Cell type fractions"
	if(!is.null(discovery_fractions) && discovery_fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
	{
		fractions = discovery_fractions
	}else{
		if(!is.null(config$"Input"$"Filter non cell type specific genes"))
		{
			if(config$"Input"$"Filter non cell type specific genes")
			{
				fractions = "Cell_type_specific_genes"
			}else{
				fractions = "All_genes"
			}
		}else{
			fractions = "Custom"	
		}
	}
}else{
	fractions = paste0(discovery, "_Fractions") 
}

#Starting EcoTyper
setwd("pipeline")
start = Sys.time()

cat("\nLoading visium data...\n")
PushToJobQueue(paste("Rscript spatial_load_visium_data.R", recovery, input_path))
RunJobQueue()

if(fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions") && !file.exists(fractions_path))
{
	cat("\nRunning CIBERSORTxFractions on the visium dataset...\n")
	PushToJobQueue(paste("Rscript csx_fractions.R", "visium", recovery, "LM22", "B_mode", CSx_username, CSx_token, TRUE))
	RunJobQueue()
	PushToJobQueue(paste("Rscript csx_fractions.R", "visium", recovery, "TR4", "B_mode", CSx_username, CSx_token, TRUE))			
	RunJobQueue()
	PushToJobQueue(paste("Rscript csx_fractions_two_tiered.R", "visium", recovery, "TR4", "B_mode", "LM22", "B_mode", fractions))
	RunJobQueue()
	coo = "Epithelial.cells"
	if(fractions %in% "Lymphoma_Fractions")
	{
		coo = "B.cells"
	}
}else{
	cat("\nLoading user-provided cell type fractions...\n")
		
	dir.create(file.path("../CIBERSORTx/fractions/visium", recovery, fractions), recursive = T, showWarnings = F)
	PushToJobQueue(paste("cp -f", fractions_path, file.path("../CIBERSORTx/fractions/visium", recovery, fractions, "CIBERSORTx_Adjusted.txt")))
	RunJobQueue()		
}

cat("\nRunning cell state recovery on the visium dataset...\n")
key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{
    n_states = key[key[,1] == cell_type, 2]
    PushToJobQueue((paste("Rscript state_recovery_visium.R", discovery, fractions, cell_type, n_states, recovery, "FALSE"))) 
}   
RunJobQueue()

cat("\nCalculating cell state abundances...\n")
print(paste("Rscript spatial_states.R", discovery, recovery, fractions, coo))
PushToJobQueue((paste("Rscript spatial_states.R", discovery, recovery, fractions, coo))) 
RunJobQueue()
cat("\nCalculating ecotype abundances...\n")
PushToJobQueue((paste("Rscript spatial_ecotypes.R", discovery, recovery, fractions, coo))) 
RunJobQueue()

cat("\nPlotting cell state heatmaps...\n")
key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{
    n_states = key[key[,1] == cell_type, 2]
    PushToJobQueue((paste("Rscript spatial_plot_states.R", discovery, recovery, fractions, coo, cell_type)))
}   
RunJobQueue()

PushToJobQueue((paste("Rscript spatial_plot_ecotypes.R", discovery, recovery, fractions, coo))) 
RunJobQueue()

cat("\nCopying EcoTyper results to the output folder!\n")

if(file.exists(final_output) && length(list.files(final_output)) > 0)
{
	old_results_folder = paste0(final_output, format(Sys.time(), " %a %b %d %X %Y"))
	dir.create(old_results_folder, recursive = T, showWarnings = F)
	warning(paste0("The output folder contains files from a previous run. Moving those files to: '", old_results_folder, "'"))	
	system(paste0("mv -f ", final_output, "/* '", old_results_folder, "'"))
}

dir.create(final_output, recursive = T, showWarnings = F)

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
for(cell_type in key[,1])
{			
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, paste0(cell_type, "_spatial_heatmaps.pdf")), final_output))
	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, paste0(cell_type, "_spatial_heatmaps.png")), final_output))
}	

system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, "state_abundances.txt"), final_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, "ecotype_abundances.txt"), final_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, "Ecotype_spatial_heatmaps.pdf"), final_output))
system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, "Ecotype_spatial_heatmaps.png"), final_output))

end = Sys.time()
cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 0), "\n"))

