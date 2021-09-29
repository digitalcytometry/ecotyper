suppressPackageStartupMessages({
library(reshape2)
})

args = c("Carcinoma", "VisiumBreast", "Carcinoma_Fractions", "Epithelial.cells")
args = commandArgs(T) 

discovery = args[1]
recovery = args[2]
fractions = args[3]
malignant_cell = args[4]

ecotypes_dir = file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery") 
output_dir = file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery) 
dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))

casted = read.delim(file.path(output_dir, "state_abundances.txt"))

mapping = read.delim(file.path(ecotypes_dir, "ecotypes.txt"))
if(all(paste0("CE", 1:10) %in% mapping$Ecotype))
{
	mapping$Ecotype = factor(as.character(mapping$Ecotype), levels = paste0("CE", 1:10))
}else{	
	mapping$Ecotype = as.factor(as.character(mapping$Ecotype))
}
mapping = mapping[order(mapping$Ecotype),]

splits = split(casted, casted$Sample)

ecotype_data = do.call(rbind, lapply(splits, function(spl)
{
	tmp = spl[,match(mapping$ID, colnames(spl))]	
	tmp[is.na(tmp)] = 0
	tmp[tmp > 1] = 1

	ecotypes = sapply(levels(mapping$Ecotype), function(x)
	{
		apply(tmp[,mapping$Ecotype == x], 1, mean)
	})
	
	ecotypes[is.na(ecotypes)] = 0
	ecotypes = apply(ecotypes, 2, function(x) x / quantile(ecotypes, .99, na.rm = T))
	ecotypes[ecotypes>1] = 1

	df = data.frame(spl[,1:5], (ecotypes))	
	df$Malignant[is.na(df$Malignant)] = 0
	#df$Malignant = df$Malignant * df$Malignant
	df$Malignant = df$Malignant - quantile(df$Malignant, .01, na.rm = T)
	df$Malignant = ifelse(df$Malignant < 0, 0, df$Malignant)
	df$Malignant = df$Malignant / quantile(df$Malignant, .99, na.rm = T)
	df$Malignant = ifelse(df$Malignant > 1, 1, df$Malignant)
	
	df
}))

melted = melt(as.data.frame(ecotype_data), id.vars = c("ID",  "X", "Y", "Sample", "Malignant"))
colnames(melted) = c("ID",  "X", "Y", "Sample", "Malignant", "Ecotype", "Abundance")

write.table(melted, file.path(output_dir, "ecotype_abundances_long.txt"), sep = "\t") 
write.table(ecotype_data, file.path(output_dir, "ecotype_abundances.txt"), sep = "\t") 

