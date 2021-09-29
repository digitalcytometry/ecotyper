library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

heatmap_color_annotation <- function(data, columns, palettes = NULL, dataset = "Carcinoma", name = "ann", which = "column", annotation_legend_param = list(), show_annotation_name = T)
{
	#print(columns)
	df = data[,match(columns, colnames(data)),drop = F]
	rownames(df) = rownames(data)
	colnames(df) = columns
	
	top_colors = list()
	for(col in columns)
	{
		if(!is.factor(df[,col]))
		{
			df[,col] = as.factor(as.character(df[,col]))
		}
		if(all(is.na(df[,col])))
		{
			df[,col] = as.factor(as.character(paste('random_', rep(sample(1:2, 1), nrow(df)))))
		}
		if(is.null(palettes))
		{
			type_cols = get_colors(length(levels(df[,col])), column = col, dataset = dataset)
		}else{
			type_cols = get_colors(length(levels(df[,col])), palette = palettes[which(columns == col)], dataset = dataset)
		}
		names(type_cols) = levels(df[,col])
		top_colors[[col]] = type_cols
	}
	if(which == "row")
	{
		annotation_name_side = "top"
	}else{
		annotation_name_side = "right"
	}
	HeatmapAnnotation(df = df, name = name, col = top_colors, annotation_name_side = annotation_name_side, annotation_legend_param = annotation_legend_param,
		show_annotation_name = show_annotation_name, which = which, border = T)
}


heatmap_simple <- function(mat, 
	left_annotation = NULL, top_annotation = NULL, bottom_annotation = NULL, right_annotation = NULL,
	left_columns = NULL, top_columns = NULL, bottom_columns = NULL, right_columns = NULL, dataset = "Carcinoma",	
	cluster_rows = F, cluster_columns = F,
	show_row_dend = F, show_column_dend = F,
	show_row_names = F, show_column_names = F,
	scale_rows = F, scale_columns = F, show_annotation_name = T,
	sort_rows_by = NULL, sort_columns_by = NULL,
	legend_name = "",
	col = NULL,
	color_range = c(-2, 0, 2), color_palette = c("gray", "#762a83", "black", "yellow"), name = "hmap", ...)
{
	if(nrow(mat) == 0 || ncol(mat) == 0)
	{
		return (NULL)
	}
	
	if(scale_rows && ncol(mat) > 1)
	{
		tmp <- colnames(mat)
		mat <- t(apply(mat, 1, function(x) scale(x)))
		colnames(mat) <- tmp
	}
	
	if(scale_columns && nrow(mat) > 1)
	{
		tmp <- rownames(mat)
		mat <- (apply(mat, 2, function(x) scale(x)))
		rownames(mat) <- tmp
	}
	
	if(!is.null(left_annotation) && class(left_annotation) != "HeatmapAnnotation")
	{
		left_annotation = left_annotation[match(rownames(mat), rownames(left_annotation)),]
		nrows <- as.integer(toupper(max(1, max(unlist(lapply(left_columns, function(x) sqrt(length(unique(left_annotation[,x])))))))))
		nrows <- max(1, min(nrows, min(unlist(lapply(left_columns, function(x) length(unique(left_annotation[,x])))))))
		left_annotation = heatmap_color_annotation(left_annotation, left_columns, name = "left_ann", which = "row", dataset = dataset, 
			show_annotation_name = show_annotation_name,
			annotation_legend_param = list(legend_direction = "horizontal", nrow = nrows, by_row = F,
				title_position = "topcenter", title_gp = gpar(fontsize = 11)))

	}
	if(!is.null(right_annotation) && class(right_annotation) != "HeatmapAnnotation")
	{
		right_annotation = right_annotation[match(rownames(mat), rownames(right_annotation)),]
		nrows <- as.integer(toupper(max(1, max(unlist(lapply(right_columns, function(x) sqrt(length(unique(right_annotation[,x])))))))))
		nrows <- max(1, min(nrows, min(unlist(lapply(right_columns, function(x) length(unique(right_annotation[,x])))))))
		right_annotation = heatmap_color_annotation(right_annotation, right_columns, name = "right_ann", which = "row", dataset = dataset, 
			show_annotation_name = show_annotation_name,
			annotation_legend_param = list(legend_direction = "horizontal", nrow = nrows, by_row = F,
				title_position = "topcenter", title_gp = gpar(fontsize = 11)))
	}
	if(!is.null(top_annotation) && class(top_annotation) != "HeatmapAnnotation")
	{
		top_annotation = top_annotation[match(colnames(mat), rownames(top_annotation)),]
		nrows <- as.integer(toupper(max(1, max(unlist(lapply(top_columns, function(x) sqrt(length(unique(top_annotation[,x])))))))))
		nrows <- max(1, min(nrows, min(unlist(lapply(top_columns, function(x) length(unique(top_annotation[,x])))))))
		top_annotation = heatmap_color_annotation(top_annotation, top_columns, name = "top_ann", dataset = dataset, 
			show_annotation_name = show_annotation_name,
			annotation_legend_param = list(legend_direction = "horizontal", nrow = nrows, by_row = F,
				title_position = "topcenter", title_gp = gpar(fontsize = 11)))
	}
	if(!is.null(bottom_annotation) && class(bottom_annotation) != "HeatmapAnnotation")
	{ 
		bottom_annotation = bottom_annotation[match(colnames(mat), rownames(bottom_annotation)),]
		nrows <- as.integer(toupper(max(1, max(unlist(lapply(bottom_columns, function(x) sqrt(length(unique(bottom_annotation[,x])))))))))
		nrows <- max(1, min(nrows, min(unlist(lapply(bottom_columns, function(x) length(unique(bottom_annotation[,x])))))))
		bottom_annotation = heatmap_color_annotation(bottom_annotation, bottom_columns, name = "bottom_ann", dataset = dataset, 
			show_annotation_name = show_annotation_name,
			annotation_legend_param = list(legend_direction = "horizontal", nrow = nrows, by_row = F,
			title_position = "topcenter", title_gp = gpar(fontsize = 11)))
	}
	if(is.null(col))
	{
		#print(color_palette)
		col = colorRamp2(breaks = color_range, colors = color_palette[-1])
		na_col = color_palette[1]
	}
	#print(dim(mat))
	
	tryCatch({
		ht_opt$message = F
	}, error = function(x){})
	ht_opt$verbose = F

	Heatmap(as.matrix(mat), top_annotation = top_annotation, left_annotation = left_annotation, 
		bottom_annotation = bottom_annotation, right_annotation = right_annotation,
		col = col, na_col = na_col, name = name, border = T,
	cluster_rows = cluster_rows, cluster_columns = cluster_columns,
	show_row_dend = show_row_dend, show_column_dend = show_column_dend, 
	show_row_names = show_row_names, show_column_names = show_column_names, #clustering_distance_columns = distfunc, clustering_distance_rows = distfunc, 
	heatmap_legend_param = list(legend_direction = "horizontal", color_bar = "continuous", #border = "black", 
		title = legend_name, title_position = "topcenter", title_gp = gpar(fontsize = 11), legend_width = unit(1.5, "in")),
	use_raster = T, raster_device = c("png"), ...
	)
}


top_markers_left <- function(data, top_n = 5, ord = NULL, ord_col = NULL, cex = NULL, exhaustion_markers = T)
{
	if(exhaustion_markers)
	{
		markers = read.delim("../utils/naive_cytotoxicity_exhaustion_simp.txt", header = F)
		markers[,2] = ""
	}else{
		markers = NULL
	}
	
	if(!is.null(ord))
	{
		spl = split(ord, ord[,ord_col])
		markers2 = do.call(rbind, lapply(spl, function(x){
			df = data.frame(V1 = rownames(x)[1:(min(top_n, nrow(x) - 1))], V2 = x[1, ord_col])
		}))
		markers = rbind(markers, markers2)	
	}

	markers = markers[!duplicated(markers[,1]),]
	labels = markers[match(rownames(data), markers[,1]),]

	labels = apply(labels, 1, function(x) if(is.na(x[[1]])) "" else ifelse(x[[2]] == "", paste0(x[[1]]), paste0(x[[1]])))	
	dup = !duplicated(ord[,ord_col])	
	dup = c(dup[-1], F)
	labels[dup] <- " "	

	tmp = labels[labels != ""]
	cols = ifelse(tmp %in% c("", " "), "white", "black")
	n_labels = length(which(labels != ""))
	
	if(is.null(cex))
	{
		cex = 30 * 0.65 / max(n_labels, 25)
	}

	rowAnnotation(link = anno_mark(at = which(labels != ""), 
			labels_gp = gpar(cex = cex, fontface = "italic"), 
			link_gp = gpar(col = cols),
			side = "left",
			extend = c(0, 0.3),
			link_width = unit(0.5, "in"),
			labels = labels[labels != ""]), width = unit(.1, "in") + max_text_width(labels))
}

rectangle_annotation_coordinates <- function(rows, columns)
{
	levels = union(levels(rows), levels(columns))
	rows = factor(as.character(rows), levels = levels)
	columns = factor(as.character(columns), levels = levels)

	dup = sapply(levels, function(lvl)which(rows == lvl)[1]-1)	
	fract = dup / length(rows)
	
	dup = sapply(levels, function(lvl)which(columns == lvl)[1]-1)	
	fract2 = dup / length(columns)
	
	for(i in length(levels):1)
	{
		
		if(i == 1)
		{
			if(is.na(fract[i]))
			{
				fract[i] = 0
			}
			if(is.na(fract2[i]))
			{
				fract2[i] = 0
			}
		}else{
			if(i == length(levels))
			{
				if(is.na(fract[i]))
				{
					fract[i] = 1
				}
				if(is.na(fract2[i]))
				{
					fract2[i] = 1
				}
			}else{
				if(is.na(fract[i]))
				{
					fract[i] = fract[i+1]
				}
				if(is.na(fract2[i]))
				{
					fract2[i] = fract2[i+1]
				}
			}
		}		
	}	
	height =  c(fract[-1], 1) - fract
	width =  c(fract2[-1], 1) - fract2 
	list(y=1-fract, h=height, x=fract2, w=width)  
}

