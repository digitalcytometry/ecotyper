library(RColorBrewer)
library(plyr)

abspath <- function(files)file.path(normalizePath(dirname(files)), basename(files))

ecotype_to_factor <- function(vector){
	unq = sort(unique(as.character(vector)))
	
	suppressWarnings({
	numbers = gsub("CE", "", unq)
	numbers = as.integer(numbers)

	if(all(!is.na(numbers)))
	{
		vector = factor(as.character(vector), levels = paste0("CE", sort(numbers)))		
		return(vector)
	}else{
		numbers = gsub("E", "", unq)
		numbers = as.integer(numbers)		
	}
		
	if(all(!is.na(numbers)))
	{
		vector = factor(as.character(vector), levels = paste0("E", sort(numbers)))
	}else{
		vector = as.factor(as.character(vector))
	}
	}) 
	
	vector
}

read_fractions <- function(input_dir){
	if(file.exists(file.path(input_dir, "CIBERSORTx_Adjusted.txt")))
	{
		cat(paste0("Reading adjusted fractions from folder: ", input_dir, "\n"))
		data = read.delim(file.path(input_dir, "CIBERSORTx_Adjusted.txt"))
	}else{
		if(file.exists(file.path(input_dir, "CIBERSORTx_Results.txt")))
		{
			cat(paste0("Reading fractions from folder: ", input_dir, "\n"))
			data = read.delim(file.path(input_dir, "CIBERSORTx_Results.txt"))
		}
		else{
			stop(paste0("File ", file.path(input_dir, "CIBERSORTx_Adjusted.txt"), " does not exist! Please make sure you ran CIBERSORTx fractions!"))
		}
	}
	if("P.value" %in% make.names(colnames(data)))
	{
		data = data[1:(which(make.names(colnames(data)) == "P.value") - 1)]
	}
	data
}

# --------------------------------------------------
# get top 1000 variable genes
# --------------------------------------------------
# m: normalized matrix
get_variable_genes <-function(m)
{    
	if(nrow(m) < 1000)
	{
		return(m)
	}
	m = t(m)
    df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
    df$dispersion<-with(df,var/mean)
    df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,unique(quantile(mean,seq(0.1,1,0.05)),Inf))))
    var_by_bin<-ddply(df,"mean_bin",function(x) {
        data.frame(bin_median=median(x$dispersion),
        bin_mad=mad(x$dispersion))
    })
    df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
    df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
    df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  
    n_genes = min(sum(!is.nan(df$dispersion_norm)), 1000)
	if(n_genes == 0)
	{
		m_n_1000 <- data.frame()
		return (NULL)
	}
    disp_cut_off<-rev(sort(df$dispersion_norm))[n_genes]
    df$used<-df$dispersion_norm >= disp_cut_off
    set.seed(0) 
    m_n_1000 <- m[,head(order(-df$dispersion_norm),n_genes)]
    t(m_n_1000)
}

scale_data <- function(data, by = NULL)
{
	if(is.null(by))
	{
		variable = rep("all samples", ncol(data))
	}else{
		variable = by
	}
	
	df = data.frame(ID = colnames(data), Variable = variable)
	splits = split(df, as.character(df$Variable))
	reconstituted_data = NULL
	for(spl in splits)
	{
		if(nrow(spl) < 2)
		{
			next
		}
	 
		#print(paste0("Scaling across: ", spl$Variable[1]))
		tmp_data = data[,colnames(data) %in% spl$ID, drop = F]
		scaled = t(apply(tmp_data, 1, scale))
		colnames(scaled) = colnames(tmp_data)

		if(is.null(reconstituted_data))
		{
			reconstituted_data = scaled
		}else{
			reconstituted_data = cbind(reconstituted_data, scaled) 
		}
	}
	rownames(reconstituted_data) = rownames(data)
	data = reconstituted_data[,match(colnames(data), colnames(reconstituted_data))]
	data
}


read_clinical <- function(sample_list, dataset = "Lung", dataset_type = "discovery", row_names = T, first_10 = F)
{
	
	if(!file.exists(file.path("../datasets", dataset_type, dataset, "annotation.txt")))
	{
		cat(paste0("No annotation file found for dataset: ", dataset, ". Creating empty annotation.\n"))
		clinical = data.frame(ID = sample_list)
		return(clinical)
	}
	clinical = read.delim(file.path("../datasets",dataset_type, dataset, "annotation.txt"))
	
	barcodes = sample_list
	clinical = clinical[match(barcodes, clinical$ID),]
	
	
	clinical$ID = sample_list
	
	if(row_names)
	{
		rownames(clinical) = make.names(clinical$ID, T)
	}
	return(clinical)
}

get_colors <- function(count, column = NULL, palette = NULL, dataset = "Carcinoma")
{
	if(!is.null(column))
	{
		if(column %in% c("Cluster", "State"))
		{
			palette = 1
		}
		if(column %in% c("Tissue", "SourceTissue"))
		{
			palette = 2
		}
		if(column %in% c("Dataset"))
		{
			palette = 3
		}
		if(column %in% c("acronym", "Histology"))
		{
			palette = 4
		}
		if(column %in% c("CellType", "Cell Type",  "Cell type"))
		{
			palette = 5
		}
		if(column %in% c("Filtered"))
		{
			palette = 6
		}
		if(column %in% c("CAF"))
		{
			palette = 7
		}
		if(tolower(column) %in% c("metacluster", "ecotype"))
		{
			palette = 8 
		}
		if(tolower(column) %in% c("coo"))
		{
			palette = 9
		}
	}
	
	if(is.null(palette))
	{
		palette = 8
	}

	if(palette == 1)
	{
		cols = c(brewer.pal(9, "Set1"), brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)] , brewer.pal(8, "Set3"))
		cols[6] = "gold"
		a1 <- col2rgb(cols)
		a2 <- rgb2hsv(a1)
		cols=hsv(a2[1,] * 0.93, a2[2,], a2[3,])
		if(dataset == "Lymphoma")
		{
			cols = c("#EB7D5B", "#FED23F", "#B5D33D","#6CA2EA", "#442288", cols) 
		}
	}else{
		if(palette == 2)
		{
			cols = rep(c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), 50)	
		}else{
			if(palette == 3)
			{
				cols = rep(c(colorRampPalette(brewer.pal(8, "Accent"))(8)[c(3,5,6,4,7,8,1,2)], colorRampPalette(brewer.pal(8, "Set2"))(8), colorRampPalette(brewer.pal(8, "Pastel1"))(8)), 10000)
				cols[6] = "lightblue"
			}else{ 
				if(palette == 4) 
				{
					cols = rep(c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")[-c(5,6,7)] , rev(brewer.pal(8, "Set2")), rev(brewer.pal(8, "Accent"))), 50)
				}else{
					if(palette == 5) 
					{
						cols = rep(rev(c(brewer.pal(8, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)])), 50)
					}else{
						if(palette == 6) 
						{
							cols = c("#016450", "#fff7fb")
						}else{
							if(palette == 7) 
							{
								cols = rev(c("blue", "hotpink4","salmon", "darkorchid", "deeppink1", "lightpink","orange2", "aquamarine", "brown1", "aquamarine4"))
							}else{
								if(palette == 8)  
								{
									cols = rep(c(brewer.pal(12, "Paired")[-c(2, 4, 6, 8, 10, 11, 12)], brewer.pal(8, "Dark2")[-c(5,7,8)] , brewer.pal(8, "Set3")), 50)
									cols[6] = "gold"
									aux = cols[10]
									cols[10] = cols[2]
									cols[2] = aux
									aux = cols[5]
									cols[5] = cols[9]
									cols[9] = aux
									aux = cols[3]
									cols[3] = cols[8]
									cols[8] = aux
									aux = cols[3]
									cols[3] = cols[7]
									cols[7] = aux
									aux = cols[4]
									cols[4] = cols[7]
									cols[7] = aux
									cols_copy = cols 

									cols = c(rgb(214, 55, 46, maxColorValue = 255), rgb(250, 221, 75, maxColorValue = 255), rgb(112, 180, 96, maxColorValue = 255), rgb(230, 144, 193, maxColorValue = 255), rgb(152, 94, 168, maxColorValue = 255), rgb(163, 163, 163, maxColorValue = 255), rgb(183, 211, 229, maxColorValue = 255), rgb(230, 216, 194, maxColorValue = 255), rgb(240, 143, 53, maxColorValue = 255), rgb(81, 137, 187, maxColorValue = 255))
									cols = c(cols, cols_copy)
									cols = cols[-c(11:17)]
									if(dataset == "Lymphoma")
									{
										cols = c("#3E2883", "#76A2E4", "#8E549F", "#458833", "#BBD058", "#FFFD61", "#F8D25D", "#F3A83B", "#EC5428", cols)
									}
								}else{
									if(palette == 9)  
									{
										cols = rep(c("#73BBCE", "#DD7A1E", "#7F7F7F"), 100)
									}else{
										cols = rep(c(colorRampPalette(brewer.pal(8, "Accent"))(limit), colorRampPalette(brewer.pal(8, "Set2"))(limit), colorRampPalette(brewer.pal(8, "Pastel1"))(limit)), 10000)
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cols[1:count] 
}


collapse_samples <- function(data, variable)
{
	agg_data = do.call(cbind, lapply(levels(variable), function(x){
		if(ncol(data[, variable == x, drop = F]) == 0)
		{
			return (rep(NA, nrow(data)))
		}
		apply(data[, variable == x, drop = F], 1, function(x) mean(x, na.rm = T))
	}))
	colnames(agg_data) = levels(variable)
	agg_data
}


get_FC_for_clusters <- function(data, variable)
{
	logFC_data = do.call(cbind, lapply(levels(as.factor(as.character(variable))), function(cluster){
			clust_samples = variable == cluster
			clust_mean_exp = apply(data[, clust_samples,drop = F], 1, function(x){
				mean(x, na.rm = T)
			})
			non_clust_mean_exp = apply(data[,!clust_samples,drop = F], 1, function(x){
				mean(x, na.rm = T)
			})
			clust_mean_exp - non_clust_mean_exp
		}))
	colnames(logFC_data) = levels(as.factor(as.character(variable)))
	logFC_data = as.data.frame(logFC_data)

	tmp = do.call(rbind, (apply(logFC_data, 1, function(x) {
			x <<- x
			maxes = which.max(x)
			if(length(maxes) == 0)
			{
				return (list(NA, NA))
			}
			if(length(maxes) > 1)
			{
				ret = maxes[which.max(x[maxes])][1]
			}else
			{
				ret = which.max(x)
			}
			list(colnames(logFC_data)[ret], x[[ret]])
		})))
	tmp[,1] = unlist(tmp[,1])
	tmp[,2] = unlist(tmp[,2])

	logFC_data = data.frame(Gene = rownames(logFC_data), logFC_data, State = unlist(tmp[,1]), MaxFC = unlist(tmp[,2]))
	logFC_data = logFC_data[order(logFC_data$State, -logFC_data$MaxFC),]
	logFC_data
}

doDE <- function(data, variable, value = NULL)
{
	library(matrixTests)
	library(parallel)

	vals = levels(as.factor(as.character(variable)))
	if(!is.null(value))
	{
		vals = value
	}
	logFC_data = do.call(rbind, lapply(vals, function(cluster){
		clust_samples = variable == cluster
		clust_mean_exp = apply(data[,clust_samples,drop = F], 1, function(x){
			mean(x, na.rm = T)
		})
		non_clust_mean_exp = apply(data[,!(clust_samples),drop = F], 1, function(x){
			mean(x, na.rm = T)
		})

		p = row_wilcoxon_twosample(data[,clust_samples,drop = F], data[,!clust_samples,drop = F])$pvalue

		data.frame(Gene = rownames(data), State= cluster, FC = clust_mean_exp - non_clust_mean_exp, P = p, Q = p.adjust(p, method = "BH"))
	}))
	logFC_data$P = ifelse(is.na(logFC_data$P), 1, logFC_data$P)
	logFC_data$Q = ifelse(is.na(logFC_data$Q), 1, logFC_data$Q)
	logFC_data = logFC_data[order(logFC_data$State, -logFC_data$FC, logFC_data$P),]
	logFC_data 
}

lookup_celltype <- function(values)
{
	unlist(lapply(as.character(values), function(x) gsub(".", " ", gsub(".and.", "/", x, fixed = T), fixed = T)))
}

lookup_celltype_short <- function(values)
{
	prefix_length = 3:max(unlist(lapply(as.character(values), function(x) nchar(x))))

	for(pl in prefix_length)
	{
		vals = unlist(lapply(as.character(values), function(x){
			l = nchar(x)
			m = min(l, pl)
			substr(x, 1, m)
			}))
		if(sum(duplicated(vals)) == 0)
		{
			break
		}
	}
	vals
}

theme_publication <- function (base_size = 11, base_family = "", base_line_size = base_size/33, 
    base_rect_size = base_size/22) 
{
    half_line <- base_size/2 
    t <- theme(line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"), 
    	rect = element_rect(fill = NA, colour = "black", size = base_rect_size, linetype = 1), 
        text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, lineheight = 0.9, 
            hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE), 
        axis.line = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.x = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.x.top = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.x.bottom = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.y = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.y.left = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.line.y.right = element_line(colour = "black", size = base_line_size, lineend = "round"),  
        axis.text = element_text(size = rel(0.8), colour = "black"), 
        axis.text.x = element_text(margin = margin(t = 0.8 * half_line/2), vjust = 1), 
        axis.text.x.top = element_text(margin = margin(b = 0.8 * half_line/2), vjust = 0), 
        axis.text.y = element_text(margin = margin(r = 0.8 * half_line/2), hjust = 1), 
        axis.text.y.right = element_text(margin = margin(l = 0.8 * half_line/2), hjust = 0), 
        axis.ticks = element_line(colour = "black", lineend = "round"), 
        axis.ticks.length = unit(half_line * 3/4, "pt"), 
        axis.ticks.length.x = NULL,  
        axis.ticks.length.x.top = NULL,
        axis.ticks.length.x.bottom = NULL, 
        axis.ticks.length.y = NULL,
        axis.ticks.length.y.left = NULL, 
        axis.ticks.length.y.right = NULL, 
        axis.title.x = element_text(margin = margin(t = half_line/2), vjust = 1), 
        axis.title.x.top = element_text(margin = margin(b = half_line/2), vjust = 0), 
        axis.title.y = element_text(angle = 90, margin = margin(r = half_line/2), vjust = 1), 
        axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2), vjust = 0),
        legend.background = element_rect(colour = NA), 
        legend.spacing = unit(2 * half_line, "pt"), 
        legend.spacing.x = NULL, 
        legend.spacing.y = NULL, 
        legend.margin = margin(half_line, half_line, half_line, half_line),
        legend.key = element_blank(), 
        legend.key.size = unit(1.2, "lines"), 
        legend.key.height = NULL, 
        legend.key.width = NULL, 
        legend.text = element_text(size = rel(0.8)), 
        legend.text.align = NULL, 
        legend.title = element_text(hjust = 0), 
        legend.title.align = NULL, 
        legend.position = "right", 
        legend.direction = NULL, 
        legend.justification = "center", 
        legend.box = NULL, 
        legend.box.margin = margin(0, 0, 0, 0, "cm"),
        legend.box.background = element_blank(), 
        legend.box.spacing = unit(2 * half_line, "pt"), 
        panel.background = element_rect(fill = NA, colour = NA), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(half_line, "pt"), 
        panel.spacing.x = NULL, 
        panel.spacing.y = NULL, 
        panel.ontop = FALSE, 
        strip.background = element_rect(fill = NA, colour = NA), 
        strip.text = element_text(colour = "black", size = rel(0.8), margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)), 
        strip.text.x = NULL,
        strip.text.y = element_text(angle = -90), 
        strip.text.y.left = element_text(angle = 90),
        strip.placement = "inside", 
        strip.placement.x = NULL, strip.placement.y = NULL, strip.switch.pad.grid = unit(half_line/2, "pt"),
        strip.switch.pad.wrap = unit(half_line/2, "pt"), 
        plot.background = element_rect(colour = NA), 
        plot.title = element_text(size = rel(1.2), hjust = 0, vjust = 1, margin = margin(b = half_line)),
        plot.title.position = "panel", 
        plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = half_line)), 
        plot.caption = element_text(size = rel(0.8), hjust = 1, vjust = 1, margin = margin(t = half_line)), 
        plot.caption.position = "panel", 
        plot.tag = element_text(size = rel(1.2), hjust = 0.5, vjust = 0.5), 
        plot.tag.position = "topleft",
        plot.margin = margin(half_line, half_line, half_line, half_line),
        complete = TRUE)
    #ggplot_global$theme_all_null %+replace% t
    t
}
