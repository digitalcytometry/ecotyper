suppressPackageStartupMessages({
library(data.table)
source("lib/misc.R")
})

args = c("PresortedDiscovery", "scRNA_all_genes", "B.cells") 
args = commandArgs(T)  
discovery = args[1]
fractions = args[2]
to_test_cell_type = args[3]

dataset_type = "discovery"
 
input_dir = file.path("../CIBERSORTx/hires", discovery, fractions) 
output_dir = file.path("../EcoTyper", discovery, fractions, "Analysis", "Cell_type_specific_genes") 
dir.create(output_dir, recursive = T, showWarning = F)

classes_path = file.path("../CIBERSORTx/hires", discovery, fractions, "classes.txt")
classes = read.delim(classes_path)

comb_mat = NULL
annotation = NULL
for(cell_type in colnames(classes))
{
	raw_input = fread(file.path(input_dir, paste0(cell_type, ".txt")), data.table  = F)
	rownames(raw_input) = raw_input[,1]	
	colnames(raw_input)[-1] = paste0(cell_type, "_", colnames(raw_input)[-1])
	colnames(raw_input)[1]= "Gene"
	if(is.null(comb_mat)){
		comb_mat = raw_input
	}else{
		comb_mat= merge(comb_mat, raw_input, by = "Gene", all = T)
	}
	ann = data.frame(ID = colnames(raw_input)[-1], CellType = cell_type)
	annotation  = rbind(annotation,  ann)
}

raw_input = comb_mat
rownames(raw_input)= raw_input[,1]
raw_input = raw_input[,-1]
annotation = annotation[annotation$ID %in% colnames(raw_input),]

splits = split(annotation, annotation$CellType)

set.seed(1234)
annotation = do.call(rbind, lapply(splits, function(spl){
	spl[sample(1:nrow(spl), min(500, nrow(spl))),]
	}))

raw_input = raw_input[,colnames(raw_input) %in% annotation$ID]
#raw_input = raw_input[apply(raw_input,  1, var) > 0 && apply(raw_input, 1, function(x) sum(x > 0) >= 3),] 

log_data = log2(raw_input + 1)

scaled_data = scale_data(log_data, by = NULL)
#scaled_data[is.na(scaled_data)] = 0

gene_info = doDE(scaled_data, annotation$CellType, to_test_cell_type)
tmp = log_data[,annotation$CellType == to_test_cell_type]
df  = data.frame(Gene = rownames(log_data), Mean = apply(tmp, 1, mean,  na.rm = T))
df  = df[match(gene_info$Gene, df$Gene),]

colnames(gene_info)[2] = "CellType"
write.table(gene_info, file.path(output_dir, paste0(to_test_cell_type, "_cell_type_specific_genes_raw.txt")), sep = "\t", row.names = F)
gene_info = gene_info[(gene_info$Q <= 0.05 & gene_info$FC > 0) | (is.na(gene_info$FC) & !is.na(df$Mean)),]
write.table(gene_info, file.path(output_dir, paste0(to_test_cell_type, "_cell_type_specific_genes.txt")), sep = "\t", row.names = F)
