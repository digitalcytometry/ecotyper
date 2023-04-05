args = commandArgs(T)
output_dir = normalizePath(args[1])
username = args[2]
token = args[3]
image_path = args[4]
id = "B_mode"

if(is.na(image_path) || is.null(image_path) || image_path == "NULL" || image_path == "'NULL'" || image_path == "NA")
{
	cat("Running on docker...\n")
	cmd_line = "docker run \\
				-v '<output_dir>':/src/data \\
				-v '<output_dir>':/src/outdir \\
				cibersortx/fractions \\
				--username <username> --token <token>  \\
				--mixture /src/data/mixture.txt \\
				--sigmatrix /src/data/sigmatrix.txt \\
				--sourceGEPs /src/data/sourceGEPs.txt \\
				--rmbatchBmode TRUE \\
				--verbose TRUE
				"
	cmd_line = gsub("<output_dir>", output_dir, cmd_line)
	cmd_line = gsub("<username>", username, cmd_line)
	cmd_line = gsub("<token>", token, cmd_line)	

}else{
	cat("Running on singularity...\n")
	cmd_line = "singularity exec -c \\
				-B '<output_dir>/':/src/data \\
				-B '<output_dir>/':/src/outdir \\
				'<image>' \\
				/src/CIBERSORTxFractions --username <username> --token <token>  \\
				--mixture /src/data/mixture.txt \\
				--sigmatrix /src/data/sigmatrix.txt \\
				--sourceGEPs /src/data/sourceGEPs.txt \\
				--rmbatchBmode TRUE \\
				--verbose TRUE
				"
	cmd_line = gsub("<output_dir>", output_dir, cmd_line)
	cmd_line = gsub("<image>", image_path, cmd_line)
	cmd_line = gsub("<username>", username, cmd_line)
	cmd_line = gsub("<token>", token, cmd_line)
}

#print(cmd_line)
err = system(cmd_line)
if(err != 0)
{
	stop()
}


