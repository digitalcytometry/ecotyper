library(data.table)
check_discovery_configuration <- function(config){
	input_mat = config$Input$"Expression matrix"
	discovery = config$Input$"Discovery dataset name"
	annotation = config$Input$"Annotation file"	
	output_dir = file.path("datasets/discovery", discovery)

	if(!file.exists(input_mat))
	{
		stop(paste0("Input format error: Input file '", input_mat, "' is missing!"))
	}
	
	mat = fread(input_mat, sep = "\t", nrows = 5)
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
			stop(paste0("Input format error: Annotation file '", annotation, "' is missing!"))
		}
		system(paste0("ln -sf ", normalizePath(annotation), " ", file.path(output_dir, "annotation.txt")))
	}
} 

check_discovery_configuration_scRNA <- function(config){
	input_mat = config$Input$"Expression matrix"
	discovery = config$Input$"Discovery dataset name"
	annotation = config$Input$"Annotation file"
	p_value_cutoff = as.numeric(as.character(config$"Pipeline settings"$"Jaccard matrix p-value cutoff"))
	output_dir = file.path("datasets/discovery", discovery)

	if(!file.exists(input_mat))
	{
		stop(paste0("Input format error: Input file '", input_mat, "' is missing!"))
	}
	
	if(length(p_value_cutoff) == 0 || is.na(p_value_cutoff) || p_value_cutoff <= 0 || p_value_cutoff >1)
	{
		
		stop(paste0("The p-value cutoff in field 'Jaccard matrix p-value cutoff' needs to be a number in the iterval (0,1]. Value provided:", p_value_cutoff, "."))
	}

	mat = fread(input_mat, sep = "\t", nrows = 5)
	if(ncol(mat) < 2)
	{
		stop(paste0(ncol(mat), " columns detected in file '", input_mat, "'. Please make sure that the file is tab-delimited!"))
	}	
	
	dir.create(output_dir, recursive = T, showWarning = F)
	system(paste0("ln -sf ", normalizePath(input_mat), " ", file.path(output_dir, "data.txt")))
	
	if(!file.exists(annotation))
	{
		stop(paste0("Input format error: Annotation file '", annotation, "' is missing! This file needs to be provided for discovery in scRNA-seq data!"))
	}

	ann = read.delim(annotation)
	if(!all(c("ID", "CellType", "Sample") %in% colnames(ann)))
	{
		stop(paste0("Input format error: Annotation file '", annotation, "' does not contain columns: 'ID', 'CellType' or 'Sample'! All three columns are required for discovery in scRNA-seq data!"))	
	}
	
	if(!all(as.character(ann$ID) == make.names(as.character(ann$ID))))
	{
		stop(paste0("Input format error: The values in column 'ID' of the annotation file '", annotation, "' contain special characters (e.g. ' ', '-'), or start with digits. Please make sure that this column and the column names of the expression matrix do not contain values modified by the R function 'make.names'."))	
	}

	if(!all(colnames(mat)[-1] %in% ann$ID))
	{
		stop(paste0("Input format error: The following ids present in the column names of the expression matrix are missing from the annotation file (column 'ID'): '", paste(colnames(mat)[-1][!colnames(mat)[-1] %in% ann$ID], collapse = "', '"), "'."))	
	}

	system(paste0("ln -sf ", normalizePath(annotation), " ", file.path(output_dir, "annotation.txt")))

} 

check_discovery_configuration_presorted <- function(config){
	input_path = config$Input$"Expression matrices"
	discovery = config$Input$"Discovery dataset name"
	annotation = config$Input$"Annotation file"

	if(config$"Pipeline settings"$"Filter non cell type specific genes")
	{
		fractions = "Cell_type_specific_genes"
	}else{
		fractions = "All_genes"
	}

	output_dir = file.path("datasets/discovery", discovery)
	csx_dir = file.path("CIBERSORTx/hires", discovery, fractions)
	dir.create(output_dir, recursive = T,  showWarning = F)
	dir.create(csx_dir, recursive = T,  showWarning = F)

	if(!is.null(annotation))
	{
		if(!file.exists(annotation))
		{
			stop(paste0("Input format error: Annotation file '", annotation, "' is missing!"))
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
			stop(paste0("Input format error: Input file '", input_mat, "' is missing!"))
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