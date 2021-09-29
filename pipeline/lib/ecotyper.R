library(NMF)

NMFpredict <- function(W, new_data)
{
	trainig_gene_set = gsub("__pos", "", rownames(W)[grepl("__pos", rownames(W))])
	new_data = new_data[match(trainig_gene_set, rownames(new_data)),]
	rownames(new_data) = trainig_gene_set
	new_data[is.na(new_data)] = 0
	to_predict = as.matrix(posneg(new_data))

	my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...)
	{
		w <- .basis(x)
		h <- .coef(x)
		nb <- nbterms(x)
		nc <- ncterms(x)
		h <- NMF:::std.divergence.update.h(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
		#w <- NMF:::std.divergence.update.w(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
		if (i%%10 == 0) {
			h <- pmax.inplace(h, eps, icterms(x))
			#w <- pmax.inplace(w, eps, ibterms(x))
		}
		if (copy) {
			#.basis(x) <- w
			.coef(x) <- h
		}
		return(x)
	}
	
	ws = W
	ws <- ws[apply(to_predict, 1, function(x) var(x) > 0),]
	to_predict = to_predict[apply(to_predict, 1, function(x) var(x) > 0),]
	ws = as.matrix(ws)

	dummy = rnmf(ncol(W), to_predict)
	
	my.seeding.method <- function(model, target){
		basis(model) <- ws #estim.r@fit@W
		# initialize H randomly
		coef(model) <- dummy@H 
		# return updated object
		return(model)
	}
	
	nmf_method <- NMFStrategy('my-method', 'brunet', Update = my_method, objective = 'KL', Stop='connectivity')

	new_nmf = nmf(to_predict, ncol(W), nrun = 1, method = nmf_method, seed = my.seeding.method, .opt='P1')
	new_nmf
}


averge_FC_marker_genes <- function(data, gene_info, assignment) 
{
	data = data[match(gene_info$Gene, rownames(data)),match(assignment$ID, colnames(data))]
	FC = sapply(levels(as.factor(gene_info$State)), function(x){
		if(sum(assignment$State == x) == 0)
		{
			return(rep(NA, nrow(data)))
		}
		apply(data[,assignment$State == x,drop = F], 1, function(y) mean(y, na.rm = T)) - apply(data[,assignment$State != x,drop = F], 1, function(y) mean(y, na.rm = T))
		})

	splits = split(as.data.frame(FC), gene_info$State)

	sapply(levels(as.factor(gene_info$State)), function(x){
		x <<- x
		if((!x %in% names(splits)) || (nrow(splits[[x]]) < 2))
		{
			return(NA)
		}

		mean(splits[[x]][,x], na.rm = T)
	})
}

calculate_recovery_significance <- function(data, initial_gene_info, initial_assignment, output_dir, n_bootstraps = 30)
{
	ref_average = averge_FC_marker_genes(data, initial_gene_info, initial_assignment)

	cat("Running recovery significance calculation. This step might take a long time...")
	coef_null = NULL
	boot_averages = NULL
	bootstrap_output = file.path(output_dir, "bootstrap_data")
	dir.create(bootstrap_output, recursive = T, showWarning = F)
	lapply(1:n_bootstraps, function(x){		
		set.seed(100 + x)

		rand_data = t(apply(data,1,sample))
		colnames(rand_data) = colnames(data)
		write.table(rand_data, file.path(bootstrap_output, paste0("bootstrap_", x, "_expression_matrix.txt")), sep = "\t")	
		
		nmf <- NMFpredict(W, rand_data)
		save(nmf, file = file.path(bootstrap_output, paste0("bootstrap_", x, ".RData")))
		coef <- apply(nmf@fit@H, 2, function(x) (x/sum(x)))
		assignment = data.frame(ID = colnames(rand_data), State = apply(coef, 2, function(x) colnames(W)[which.max(x)]))
		rownames(coef) = sprintf("IS%02d.B%02d", 1:n_states, x)
		coef_null <<- rbind(coef_null, coef)

		boot_average = averge_FC_marker_genes(rand_data, initial_gene_info, assignment)
		boot_averages <<- rbind(boot_averages, boot_average)
		})
	write.table(coef_null, file.path(output_dir, "coef_null.txt"), sep = "\t")

	means = apply(boot_averages, 2, function(x) mean(x, na.rm = T))
	sds = apply(boot_averages, 2, function(x) sd(x, na.rm = T))

	Z = ( ref_average - means) / sds
	df = as.data.frame(ref_average)
	zs_corr = data.frame(State = rownames(df), Z = Z, RefMean = df[,1],  BootstrapMean = means, BootstrapSD = sds) 

	write.table(zs_corr, file.path(output_dir, "recovery_z_scores_initial_states.txt"), sep = "\t", row.names = F)

	filt_zs_corr = zs_corr[match(mapping[,2], zs_corr$State),]
	filt_zs_corr$State = as.character(mapping[,1])
	write.table(filt_zs_corr, file.path(output_dir, "recovery_z_scores.txt"), sep = "\t", row.names = F)
}


