suppressPackageStartupMessages({
library(data.table)
source("lib/multithreading.R")
source("lib/misc.R")
})

args = c("discovery", "Lung", "TR12", "", "", 10)
args = commandArgs(T)  
dataset_type = args[1]
dataset = args[2]
fractions = args[3]
username = args[4]
token = args[5]
n_threads = as.integer(args[6])

output_dir = file.path("../CIBERSORTx/hires", dataset, fractions) 
dir.create(output_dir, recursive = T, showWarning = F)

input_mixture = fread(file.path("../datasets/", dataset_type, dataset, "data.txt"), data.table = F)
fractions_dir = file.path("../CIBERSORTx/fractions", dataset_type, dataset, fractions)
csx_fractions = read_fractions(fractions_dir)

#csx_fractions = csx_fractions[,1:(ncol(csx_fractions) - 3)]
csx_fractions = csx_fractions[match(colnames(input_mixture)[-1], csx_fractions$Mixture),]

write.table(csx_fractions, file.path(output_dir, "cibresults.txt"), sep = "\t", row.names = F)
classes = colnames(csx_fractions)[2:ncol(csx_fractions)]
write.table(t(classes), file.path(output_dir, "classes.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

n_workers = ceiling(nrow(input_mixture) * ncol(input_mixture) / 500000)
spl = rep_len(1:n_workers, length.out = nrow(input_mixture))

splits = split(input_mixture, spl)
tmp = lapply(1:length(splits), function(x){
	sub_dir = file.path(output_dir, paste0("worker_", x))
	dir.create(sub_dir, recursive = T, showWarning = F)
	fwrite(splits[[x]], file.path(sub_dir, paste0("input.txt")), sep = "\t", row.names = F, quote = F)
	PushToJobQueue(paste("Rscript csx_hires_worker.R", dataset, fractions, x, 1, username, token))
})
RunJobQueue()
