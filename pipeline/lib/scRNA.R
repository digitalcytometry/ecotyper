suppressWarnings({
library(data.table)
})

read_scRNA <- function(dataset, cell_type, log2 = F, scale_rows = F)
{
	#print(dataset)
	annotation = read.delim(file.path("../datasets/scRNA", dataset, "annotation.txt"))
	annotation = annotation[annotation$CellType == cell_type, ]
	if(nrow(annotation)==0)
	{
		#warning(paste0("There are no cells annotated to cell type '", cell_type, "' in dataset '", dataset, 
		#	"'. Please make sure that column 'CellType' in file '", file.path("../datasets/scRNA", dataset, "annotation.txt"), "' contains value '", cell_type, "'"))
		return (NULL)
	}

	expr = fread(file.path("../datasets/scRNA", dataset, "data.txt"), data.table = F)
	expr = expr[!is.na(expr[,1]),]
	expr = expr[!duplicated(expr[,1]),]

	rownames(expr) = expr[,1]
	expr = expr[,-1]

	annotation = annotation[annotation$ID %in% colnames(expr),]
	if(nrow(annotation)==0)
	{
		stop(paste0("The cell IDs do not match between the column names of expression matrix and column 'CellType' of the annotation file, for cell type '", cell_type, "'."))
	}
	rownames(annotation) = annotation$ID
	expr = expr[,match( annotation$ID, colnames(expr))]

	if(log2)
	{
		if(!all(expr >= 0))
		{
			stop(paste0("Not all values in the expression matrix of dataset '", dataset, "' are positive. Could not take log2!"))
		}
		expr = log2(expr + 1)
	}

	if(scale_rows)
	{
		tmp_cp = colnames(expr)
		expr = t(apply(expr, 1, scale))
		colnames(expr) = tmp_cp
	}

	return (list(expr = expr, annotation = annotation))
}