suppressPackageStartupMessages({
library(data.table)
source("lib/misc.R")
})

args = c("discovery_scRNA_CRC", "Cell_type_specific_genes", "B.cells") 
args = c("discovery_scRNA_CRC", "Top_1000", "B.cells") 
args = commandArgs(T)  

dataset = args[1]
fractions = args[2]
cell_type = args[3]
scaling_column = args[4] 

dataset_type = "discovery"
 
input_dir = file.path("../datasets", dataset_type, dataset) 
output_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "Cell_type_specific_genes") 
dir.create(output_dir, recursive = T, showWarning = F)

annotation = read.delim(file.path(input_dir, paste0("annotation.txt")))

raw_input = fread(file.path(input_dir, paste0("data.txt")), data.table  = F)
rownames(raw_input) =  raw_input[,1]
raw_input = raw_input[,-1]
annotation = annotation[annotation$ID %in% colnames(raw_input),]

tb = as.data.frame(table(as.character(annotation$CellType)))

if(!cell_type %in% tb[,1])
{
	stop(paste0("Cell type ", cell_type, " not found in annotation file!"))
}
if(tb[tb[,1] == cell_type,2] < 50)
{
	warning(paste0("Only ", nrow(annotation), " single cells are available for cell type: ", cell_type, ". At least 50 are required. Skipping this cell type!\n"))			
	quit(status=0, save='no')
}

annotation = annotation[annotation$CellType %in% tb[,1],]
splits = split(annotation, annotation$CellType)

set.seed(1234)
annotation = do.call(rbind, lapply(splits, function(spl){
	spl[sample(1:nrow(spl), min(500, nrow(spl))),]
	}))

raw_input = raw_input[,colnames(raw_input) %in% annotation$ID]
#raw_input = raw_input[apply(raw_input,  1, var) > 0 && apply(raw_input, 1, function(x) sum(x > 0) >= 3),] 

log_data = log2(raw_input + 1)

clinical = read_clinical(colnames(log_data), dataset = dataset, dataset_type = "discovery")
if(is.na(scaling_column)) 
{
	by = NULL
}else{
	by = as.character(clinical[,scaling_column])
}

scaled_data = scale_data(log_data, by = by)
scaled_data[is.na(scaled_data)] = 0
if(fractions == "Cell_type_specific_genes")
{
	gene_info = doDE(scaled_data, clinical$CellType, cell_type)
	colnames(gene_info)[2] = "CellType"
	write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes_raw.txt")), sep = "\t", row.names = F)
	gene_info = gene_info[gene_info$Q <= 0.05 & gene_info$FC > 0,]	
	write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes.txt")), sep = "\t", row.names = F)
}else{
	if(grepl("Top", fractions))
	{
		n_genes = as.integer(gsub("Top_", "", fractions))
		if(is.na(n_genes))
		{
			stop("The number of genes to be selected is not an integer!")
		}
		gene_info = doCV(log_data, clinical$CellType, cell_type)
		colnames(gene_info)[2] = "CellType"
		if(n_genes > nrow(gene_info))
		{
			warning(paste0("Could not select top ", n_genes, ". Only ", nrow(gene_info), " availale!"))
		}
		write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes_raw.txt")), sep = "\t", row.names = F)
		gene_info = gene_info[1:(min(nrow(gene_info), n_genes)),]	
		write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes.txt")), sep = "\t", row.names = F)
	}	
}

