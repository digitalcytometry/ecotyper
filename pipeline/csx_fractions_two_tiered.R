suppressPackageStartupMessages({
source("lib/misc.R")
})

args = commandArgs(T)
dataset_type = args[1]
dataset = args[2]
sig1_name = args[3]
sig1_mode = args[4]
sig2_name = args[5]
sig2_mode = args[6]
output_signature = args[7]

sig1_folder = file.path("../CIBERSORTx/fractions", dataset_type, dataset, sig1_name, sig1_mode)
sig1 = read_fractions(sig1_folder)

sig2_folder = file.path("../CIBERSORTx/fractions", dataset_type, dataset, sig2_name, sig2_mode)
sig2 = read_fractions(sig2_folder)

output_dir = file.path("../CIBERSORTx/fractions", dataset_type, dataset, output_signature) 
dir.create(output_dir, recursive = T, showWarning = F)

tmp_sig = merge(sig1, sig2, by = "Mixture")

if(output_signature == "Carcinoma_Fractions")
{
	mappings = read.delim("../utils/Carcinoma_mapping.txt", header = F)
}else{
	mappings = read.delim("../utils/Lymphoma_mapping.txt", header = F)
}

mappings[,1] = make.names(mappings[,1])
mappings[,2] = make.names(mappings[,2])
chunks = split(mappings, mappings[,2])

sig = do.call(cbind, lapply(chunks, function(x){
	apply(tmp_sig[,as.character(x[,1]),drop = F], 1, sum) * tmp_sig[, "CD45"]
}))
sig = cbind(tmp_sig[,!(colnames(tmp_sig) %in% as.character(mappings[,1]) | colnames(tmp_sig) %in% as.character(mappings[,2]))], sig)
sig = sig[,colnames(sig) != "CD45"]
sig = sig[,colnames(sig) != "Eos"]

if(output_signature == "Carcinoma_Fractions")
{
	colnames(sig)[which(colnames(sig) == "EPCAM")] = "Epithelial.cells"	
}else{
	sig = sig[,colnames(sig) != "EPCAM"]
}
colnames(sig)[which(colnames(sig) == "CD10")] = "Fibroblasts"
colnames(sig)[which(colnames(sig) == "CD31")] = "Endothelial.cells"

sig[,-1] = t(apply(sig[, -1], 1, function(x) x / sum(x)))

write.table(sig, file.path(output_dir, "CIBERSORTx_Results.txt"), sep = "\t", quote = F, row.names = F)

