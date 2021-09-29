suppressPackageStartupMessages({
library(data.table)
library(ggplot2)
source("lib/misc.R")
source("lib/heatmaps.R")
})

args = c("MyDiscovery", "Carcinoma")
args = commandArgs(T) 
dataset = args[1]
fractions = args[2]

key_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
output_dir = file.path("../EcoTyper", dataset, fractions,  "Analysis", "advanced_QC_filter")
dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path(key_dir, "rank_data.txt"))

all_df = NULL
tmp = apply(key, 1, function(x){
	x<<-x 
	#print(x)
	cell_type = x[[1]]
	cluster = gsub(" ","",x[[2]])

	input_dir = file.path("../CIBERSORTx/hires", dataset, fractions) 
	csx_data = fread(file.path(input_dir, paste0(cell_type, ".txt")), data.table = F)
	rownames(csx_data) = csx_data[,1]
	csx_data = csx_data[,-1]
	
	all_coef = NULL
	#print(cluster) 
	input_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type, cluster) 
	
	gene_info = read.delim(file.path(input_dir, "gene_info.txt"))
	classes = read.delim(file.path(input_dir, "state_assignment.txt"))
	
	data = csx_data[match(gene_info$Gene, rownames(csx_data)),]
	data = data[,match(classes$ID, colnames(data))]

	var1 = apply(log2(data), 1, mean)
	var2 = apply(data, 1, function(x) rev(sort(table(unlist(x))))[1]) / ncol(data)

	small_genes_table = gene_info
	small_genes_table$Expression = var1
	small_genes_table$Dropout = var2
	
	splits = split(small_genes_table, as.character(small_genes_table$State))
	for(spl in splits)
	{
		
		df = data.frame(State = spl$State[1], CellType = cell_type,
			Dropout = mean(spl$Dropout), Expression = mean(spl$Expression))
		all_df <<- rbind(all_df, df)
	}
})

all_df$ID = paste0(all_df$CellType, "_", all_df$State)
write.table(all_df, file.path(output_dir, "dropout_scores.txt"), sep = "\t")

all_df = read.delim(file.path(output_dir, "dropout_scores.txt"))

all_df$Dropout_index = scale((as.vector(scale(all_df$Dropout)) + as.vector(scale(all_df$Expression))))
all_df$Dropout_index_discrete = as.factor(ifelse(all_df$Dropout_index < 1.96, "Non-outlier", "Outlier"))
pdf(file.path(output_dir, "dropout_scores.pdf"), width = 5, height = 5, family = "Helvetica", useDingbats = F)
g <- ggplot(all_df, aes(x = Dropout_index)) +
	geom_density()+
	geom_vline(aes(xintercept = 1.96), lty = 2, color = "red") + 
	geom_jitter(aes(y = 0.05, color = Dropout_index_discrete), height = 0.025, width = 0, size = 2) + 
	theme_bw() + 
	theme(panel.grid = element_blank()) + 
	theme(aspect.ratio = 1) + 
	labs(x = "Dropout score", y = "Density", color = "Status") + 
	scale_color_manual(values = c("gray", "red")) 
plot(g)
tmp = dev.off()

png(file.path(output_dir, "dropout_scores.png"), width = 5, height = 5, units = "in", res = 400, family = "Helvetica")
g <- ggplot(all_df, aes(x = Dropout_index)) +
	geom_density()+
	geom_vline(aes(xintercept = 1.96), lty = 2, color = "red") + 
	geom_jitter(aes(y = 0.05, color = Dropout_index_discrete), height = 0.025, width = 0, size = 2) + 
	theme_bw() + 
	theme(panel.grid = element_blank()) + 
	theme(aspect.ratio = 1) + 
	labs(x = "Dropout score", y = "Density", color = "Status") + 
	scale_color_manual(values = c("gray", "red")) 
plot(g)
tmp = dev.off()

write.table(all_df, file.path(output_dir, "dropout_scores_with_outliers.txt"), sep = "\t")