suppressPackageStartupMessages({
library(ggplot2)
library(data.table)
library(reshape2)
library(ComplexHeatmap)
source("lib/misc.R")
})

args = c("Carcinoma", "VisiumBreast", "Carcinoma_Fractions", "Epithelial.cells")
args = commandArgs(T)

discovery = args[1]
recovery = args[2]
fractions = args[3]
malignant_cell = args[4]

new_dataset_type = "visium"

output_dir = file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery) 
dir.create(output_dir, recursive = T, showWarning = F) 

annotation = read.delim(file.path("../datasets/visium", recovery, "annotation.txt"))

cell_fractions = read_fractions(file.path("../CIBERSORTx/fractions/visium", recovery, fractions))
colnames(cell_fractions) = make.names(colnames(cell_fractions))

write.table(cell_fractions, file.path(output_dir, paste0(recovery, "_CSx_fractions.txt")), sep = "\t", row.names = F)

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
all_melted = NULL
for(cell_type in key[,1])
{
	cat(paste0("Calculating cell state abundance for ", cell_type, "...\n"))
	
	clusters = key[key[,1] == cell_type, 2]
	
	mapping = read.delim(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, clusters, "mapping_to_initial_states.txt"))
	
	states_dir = file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, clusters)
	H = read.delim(file.path(states_dir, "initial_state_abundances.txt"))
	
	H = apply(H, 2, function(x) (x == max(x)) + 0) 
	H = H[match(mapping[,2], rownames(H)),]
	rownames(H) = mapping[,1]
	
	H = t(H)
	
	tmp = cell_fractions[,c(colnames(cell_fractions)[1], cell_type, malignant_cell)]
	colnames(tmp)[2:3] = c(cell_type, "Malignant")

	merged = merge(H, tmp, by.x = "row.names", by.y = colnames(cell_fractions)[1], all.x = T)
	colnames(merged)[1] = "ID"
	
	melted = melt(merged, id.vars = c("ID", cell_type, "Malignant"))
		
	colnames(melted) = c("ID", "CSx_fraction", "Malignant", "State", "NMF")
	melted$CellType = cell_type
	melted$Sample = recovery
	melted$NMF_weighted_raw = melted$NMF * melted$CSx_fraction
	melted$Abundance = melted$NMF_weighted_raw / quantile(melted$NMF_weighted_raw, .99, na.rm = T)
	melted$Abundance  = ifelse(melted$Abundance > 1, 1, melted$Abundance)
	
	all_melted = rbind(all_melted, melted)
}

all_melted$X = annotation[match(all_melted$ID, annotation$ID),"X"]
all_melted$Y = annotation[match(all_melted$ID, annotation$ID),"Y"]

write.table(all_melted, file.path(output_dir, "state_abundances_long.txt"), sep = "\t", row.names = F)

all_melted$State_ID = paste0(all_melted$CellType, "_", all_melted$State)
casted_data = dcast(all_melted, ID + X + Y + Sample + Malignant  ~ State_ID, value.var = "Abundance")
write.table(casted_data, file.path(output_dir, "state_abundances.txt"), sep = "\t", row.names = F)

