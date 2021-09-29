args = commandArgs(T)
output_dir = normalizePath(args[1])
username = args[2]
token = args[3]
id = "S_mode"

cmd_line = "docker run \\
				-v <output_dir>:/src/data \\
				-v <output_dir>:/src/outdir \\
				cibersortx/fractions \\
				--username <username> --token <token>  \\
				--mixture /src/data/mixture.txt \\
				--sigmatrix /src/data/sigmatrix.txt \\
				--refsample /src/data/refsample.txt \\				
				--rmbatchSmode TRUE \\
				--verbose TRUE
				"
cmd_line = gsub("<output_dir>", output_dir, cmd_line)
cmd_line = gsub("<username>", username, cmd_line)
cmd_line = gsub("<token>", token, cmd_line)
#print(cmd_line)
err = system(cmd_line)
if(err != 0)
{
	stop()
}


