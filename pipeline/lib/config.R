library(data.table)
check_discovery_configuration <- function(config){
	input_mat = config$Input$"Expression matrix"
	discovery = config$Input$"Discovery dataset name"
	annotation = config$Input$"Annotation file"
	output_dir = file.path("datasets/discovery", discovery)

	if(!file.exists(input_mat))
	{
		stop(paste0("Error: Input file '", input_mat, "' is missing!"))
	}
	
	mat = fread(input_mat, sep = "\t")
	if(ncol(mat) < 2)
	{
		stop(paste0(ncol(mat), " columns detected in file '", input_mat, "'. Please make sure that the file is tab-delimited!"))
	}	
	
	dir.create(output_dir, recursive = T, showWarning = F)
	system(paste0("ln -sf ", normalizePath(input_mat), " ", file.path(output_dir, "data.txt")))
	if(!is.null(annotation))
	{
		if(!file.exists(annotation))
		{
			stop(paste0("Error: Annotation file '", annotation, "' is missing!"))
		}
		system(paste0("ln -sf ", normalizePath(annotation), " ", file.path(output_dir, "annotation.txt")))
	}
} 

check_discovery_configuration_presorted <- function(config){
	input_path = config$Input$"Expression matrices"
	discovery = config$Input$"Discovery dataset name"
	annotation = config$Input$"Annotation file"

	if(config$"Pipeline settings"$"Filter non cell type specific genes")
	{
		fractions = "scRNA_specific_genes"
	}else{
		fractions = "scRNA_all_genes"
	}

	output_dir = file.path("datasets/discovery", discovery)
	csx_dir = file.path("CIBERSORTx/hires", discovery, fractions)
	dir.create(output_dir, recursive = T,  showWarning = F)
	dir.create(csx_dir, recursive = T,  showWarning = F)

	if(!is.null(annotation))
	{
		if(!file.exists(annotation))
		{
			stop(paste0("Error: Annotation file '", annotation, "' is missing!"))
		}
		system(paste0("ln -sf ", normalizePath(annotation), " ", file.path(output_dir, "annotation.txt")))
	}

	classes = NULL
	for(file in list.files(input_path))
	{
		cell_type = gsub(".txt$", "",  file)
		input_mat = file.path(input_path, file)
		if(!file.exists(input_mat))
		{
			stop(paste0("Error: Input file '", input_mat, "' is missing!"))
		}
	
		mat = fread(input_mat, sep = "\t")
		if(ncol(mat) < 2)
		{
			stop(paste0(ncol(mat), " columns detected in file '", input_mat, "'. Please make sure that the file is tab-delimited!"))
		}	
		classes = rbind(classes, data.frame(x = cell_type))
		system(paste0("ln -sf ", normalizePath(input_mat), " ", file.path(csx_dir, paste0(cell_type, ".txt"))))
	}	
	write.table(t(classes), file.path(csx_dir, "classes.txt"), sep = "\t", row.names = F, col.names = F)
	
} 