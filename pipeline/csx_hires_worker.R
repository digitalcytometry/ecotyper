suppressPackageStartupMessages({
source("lib/misc.R")
})

#args = c("discovery", "Lung", "TR12", "two_tiered", "1", "1")
args = commandArgs(T)
dataset = args[1]
fractions = args[2]
worker = args[3]
threads = args[4]
username = args[5]
token = args[6]

if(is.null(threads) | is.na(threads))
{
	threads = 1
}

input_dir = file.path("../CIBERSORTx/hires", dataset, fractions) 
output_dir = file.path("../CIBERSORTx/hires", dataset, fractions, paste0("worker_", worker))
dir.create(output_dir, recursive = T, showWarning = F)

system(paste0("cp -fL ", normalizePath(file.path(input_dir, "classes.txt")), " ", file.path(output_dir, "classes.txt")))
system(paste0("cp -fL ", normalizePath(file.path(input_dir, "cibresults.txt")), " ", file.path(output_dir, "cibresults.txt")))

cmd_line = "docker run \\
				-v <output_dir>:/src/data \\
				-v <output_dir>:/src/outdir \\
				cibersortx/hires \\
				--username <username> --token <token>  \\
				--mixture /src/data/input.txt \\
				--sigmatrix /src/data/classes.txt \\
				--classes /src/data/classes.txt \\
				--cibresults /src/data/cibresults.txt \\
				--heatmap FALSE \\
				--threads <threads>
				"
cmd_line = gsub("<output_dir>", normalizePath(output_dir), cmd_line)
cmd_line = gsub("<threads>", threads, cmd_line)
cmd_line = gsub("<username>", username, cmd_line)
cmd_line = gsub("<token>", token, cmd_line)

#print(cmd_line)
err = system(cmd_line)
if(err != 0)
{
	stop()
}

