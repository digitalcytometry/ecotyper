library(parallel)

job_queue = c()
PushToJobQueue <- function(cmd){
	job_queue <<- c(job_queue, cmd)
}

RunJobQueue <- function()
{
	if(length(job_queue) == 0)
	{
		return(NULL)
	}
	res = mclapply(job_queue, FUN = system, mc.cores = n_threads)	
	job_queue <<- c()
	errors = sum(unlist(res))
	if(errors > 0)
	{
		stop("EcoTyper failed. Please check the error message above!")
	}
}