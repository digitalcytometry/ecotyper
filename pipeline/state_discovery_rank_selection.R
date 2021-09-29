suppressPackageStartupMessages({
library(NMF)
library(reshape2)
library(ggplot2)
library(cowplot)
source("lib/misc.R")
})
 
args = c("discovery_cross_cor", "discovery_scRNA_CRC", "scRNA_specific_genes", "20", "0.95")
args = commandArgs(T)   
dataset_type=args[1]
dataset = args[2] 
fractions = args[3]
max_n_clusters = as.integer(as.character(args[4]))
min_cophenetic = as.numeric(as.character(args[5]))

base_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", dataset_type)
output_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
dir.create(output_dir, recursive = T, showWarning = F)

cell_types = list.files(base_dir)

all_data = NULL
for(cell_type in cell_types)
{
	for(n_clusters in 2:max_n_clusters)
	{
		#print(n_clusters)
		input_dir = file.path(base_dir, cell_type, n_clusters) 
		
		if(!file.exists(file.path(input_dir, "rank_data.txt")))
		{
			next
		}
		all_coef = read.delim(file.path(input_dir, "rank_data.txt"))
		all_coef$CellType = cell_type
		all_data = rbind(all_data, all_coef)
	}
}

find_threshold <- function(spl, min_cophenetic)
{	
	cross = c(F, F, sapply(3:nrow(spl), function(i) spl[i,"Cophenetic"] < min_cophenetic & spl[i - 1, "Cophenetic"] >= min_cophenetic & spl[i - 2, "Cophenetic"] >= min_cophenetic))
	if(length(which(cross)) == 0)
	{
		idx = max_n_clusters
	}else{
		last_cross <<- which(cross)[length(which(cross))]
		idx = last_cross - 2 + which.min(c(abs(min_cophenetic - spl[last_cross - 1, "Cophenetic"]), abs(min_cophenetic - spl[last_cross,"Cophenetic"])))		
	}

	data.frame(CellType = spl$CellType[1], Chosen_Rank = spl[idx, "n_clusters"])
}

biggest_drop <- function(spl)
{
	casted = dcast(all_data, n_clusters~CellType, value.var = "Cophenetic")
	cross = c(NA, sapply(2:nrow(spl), function(i) spl[i - 1, "Cophenetic"] - spl[i,"Cophenetic"]))	
	spl[which.max(cross)-1, "n_clusters"]
}

splits = split(all_data, as.character(all_data$CellType))
nclst = do.call(rbind, lapply(splits, function(spl) find_threshold(spl, min_cophenetic)))

if(!all(!is.na(nclst$Chosen_Rank)))
{
	#nclst = sapply(splits, biggest_drop)	
	cat(paste0("************************************************\n",
	"WARNING: ",	
	"The Cophenetic coefficient cutoff provided in the config file (", min_cophenetic,
	 ") did not work for the following cell types: ", paste0(as.character(nclst[is.na(nclst$Chosen_Rank), "CellType"]), collapse = ", "), " ",
	"For these cell types, the Cophenetic coefficient corresponding to the biggest drop has been selected. ",
	"Please inspect the final output file 'rank_plot.pdf' to ensure that the selected number of states is appropriate! ",
	"If the values are not appropriate, please adjust the cutoff in the config file and re-run EcoTyper from this step onwards!\n",
	"************************************************\n"))	
	drops = sapply(splits, biggest_drop)
	nclst$Redefined = ifelse(is.na(nclst$Chosen_Rank), T, F)
	nclst$Chosen_Rank = ifelse(is.na(nclst$Chosen_Rank), drops, nclst$Chosen_Rank)
}else{
	cat(paste0(
		"Number of states has been successfully selected for each cell type!",
		" Please inspect the final output file 'rank_plot.pdf' to ensure that the selected number of states is appropriate!\n"
		))
	nclst$Redefined = F
}


#nclst = do.call(rbind, lapply(splits, find_threshold))

ord = sort(unique(as.character(all_data$CellType)))
ord = ord[ord %in% nclst$CellType]
nclst = nclst[match(ord, nclst$CellType),]

write.table(nclst, file.path(output_dir, "rank_data.txt"), row.names = F, sep = "\t")
write.table(nclst, file.path(output_dir, paste0("rank_data_cutoff_", min_cophenetic,".txt")), row.names = F, sep = "\t")

nclst$min_cophenetic = min_cophenetic	
nclst$Alpha = factor(ifelse(nclst$Redefined, "", "Cophenetic cutoff supplied"), levels = c("", "Cophenetic cutoff supplied"))
scale_alpha = c(0, 1)
if(length(unique(nclst$Alpha)) == 1)
{
	scale_alpha = 1
}
pdf(file.path(output_dir, "rank_plot.pdf"), family = "Helvetica", useDingbats = F,
	width = 3 * sqrt(length(table(all_data$CellType))), height = 2.5 * sqrt(length(table(all_data$CellType))))
g <- ggplot(all_data, aes(x = n_clusters, y = Cophenetic)) + 
	geom_point() +
	geom_line() + 
	geom_vline(data = nclst, aes(xintercept = Chosen_Rank, col = "Number of states selected by EcoTyper"), lty = 2) + 
	geom_hline(data = nclst, aes(yintercept = min_cophenetic, alpha = Alpha), col = "black", lty = 2) + 
	theme_publication() + 
	theme(aspect.ratio = 1) + 
	scale_color_manual(values = c("red")) + 
	scale_alpha_manual(values = scale_alpha) + 
	labs(y = "Cophenetic coefficient", x = "Number of states", color = NULL, alpha = NULL) + 
	theme(legend.position = "bottom") + 
	facet_wrap(.~CellType, scales = "free") 
	
plot(g)
tmp = dev.off()

pdf(file.path(output_dir, paste0("rank_plot_cutoff_", min_cophenetic,".pdf")), family = "Helvetica", useDingbats = F,
	width = 3 * sqrt(length(table(all_data$CellType))), height = 2.5 * sqrt(length(table(all_data$CellType))))
g <- ggplot(all_data, aes(x = n_clusters, y = Cophenetic)) + 
	geom_point() +
	geom_line() + 
	geom_vline(data = nclst, aes(xintercept = Chosen_Rank, col = "Number of states selected by EcoTyper"), lty = 2) + 
	geom_hline(data = nclst, aes(yintercept = min_cophenetic, alpha = Alpha), col = "black", lty = 2) + 
	theme_publication() + 
	theme(aspect.ratio = 1) + 
	scale_color_manual(values = c("red")) + 
	scale_alpha_manual(values = scale_alpha) + 
	labs(y = "Cophenetic coefficient", x = "Number of states", color = NULL, alpha = NULL) + 
	theme(legend.position = "bottom") + 
	facet_wrap(.~CellType, scales = "free") 
	
plot(g)
tmp = dev.off()

png(file.path(output_dir, "rank_plot.png"), family = "Helvetica",
	width = 3 * sqrt(length(table(all_data$CellType))), height = 2.5 * sqrt(length(table(all_data$CellType))), units = "in", res = 150)
plot(g)
tmp = dev.off()

