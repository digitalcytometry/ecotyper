library(data.table)
source("lib/misc.R")

args = c("discovery", "Lung", "LM22", "B_mode", "", "", F)
args = commandArgs(T) 

dataset_type = args[1]
dataset = args[2]
sigmatrix = args[3]
id = args[4]
username = args[5]
token = args[6]
filt_0 = as.logical(args[7])

if(is.na(filt_0))
{
	filt_0 = F
}

if(!id %in% c("no_batch", "B_mode", "S_mode"))
{
	stop(paste('Argument 4 of script "csx_fractions.R" is restricted to: "no_batch", "B_mode", "S_mode". The value supplied is:', id))
}

input_mixture = file.path("../datasets", dataset_type, dataset, "data.txt")
sigmatrix_path = file.path("../utils/signature_matrices", sigmatrix, sigmatrix)
refsample_path = file.path("../utils/signature_matrices", sigmatrix, "data.txt")
sourceGEPs_path = file.path("../utils/signature_matrices", sigmatrix, paste0(sigmatrix, "_sourceGEP.txt"))

output_dir = file.path("../CIBERSORTx/fractions/", dataset_type, dataset, sigmatrix, id) 
dir.create(output_dir, recursive = T, showWarning = F)

if(!filt_0)
{
	err = system(paste0("cp -fL '", normalizePath(input_mixture), "' '", file.path(normalizePath(output_dir), "mixture.txt'")))

	if(err != 0)
	{
		stop()
	}
}else{
	cat("Warning: Running CIBERSORTx on the subset of ST spots that contain expression for at least one gene in the signature matrix. This is necessary to prevent matrix singularity issues observed when running CIBERSORTx on very sparse data, as is the case with ST arrays!\n")
	
	data = fread(input_mixture, data.table = F)
	sig = read.delim(file.path(sigmatrix_path))
	small_data = data[data[,1] %in% sig[,1],]
	data = data[,c(T, apply(small_data[,-1], 2, function(x) sum(x > 0) > 0))]

	write.table(data, file.path(normalizePath(output_dir), "mixture.txt"), row.names = F, sep = "\t")
}
err = system(paste0("cp -fL '", normalizePath(sigmatrix_path), "' '", file.path(normalizePath(output_dir), "sigmatrix.txt'")))
if(err != 0)
{
	stop()
}
err = system(paste0("cp -fL '", normalizePath(sourceGEPs_path), "' '", file.path(normalizePath(output_dir), "sourceGEPs.txt'")))
if(err != 0)
{
	stop()
}

if(id == "no_batch")
{
	cat(paste0("Running CIBERSORTx fractions with no batch correction, using the '", sigmatrix , "' signature matrix and '", dataset, "' dataset.\n"))
	err = system(paste("Rscript csx_fractions_no_batch_mode.R ", output_dir, username, token))
	if(err != 0)
	{
		stop()
	}
}
if(id == "B_mode")
{
	cat(paste0("Running CIBERSORTx fractions with B-mode batch correction, using the '", sigmatrix , "' signature matrix and '", dataset, "' dataset.\n"))
	err = system(paste("Rscript csx_fractions_B_mode.R ", output_dir, username, token))
	if(err != 0)
	{
		stop()
	}
}
if(id == "S_mode")
{
	cat(paste0("Running CIBERSORTx fractions with S-mode batch correction, using the '", sigmatrix , "' signature matrix and '", dataset, "' dataset.\n"))
	err = system(paste("Rscript csx_fractions_S_mode.R ", output_dir, username, token))
	if(err != 0)
	{
		stop()
	}
}
