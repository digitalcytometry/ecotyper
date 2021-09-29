suppressPackageStartupMessages({
library(ComplexHeatmap)
library(colorspace) 
library(reshape2)
library(ggplot2)
library(circlize)
source("lib/misc.R")
})

args = c("Carcinoma", "VisiumBreast", "Carcinoma_Fractions", "Epithelial.cells")
args = commandArgs(T) 

discovery = args[1]
recovery = args[2]
fractions = args[3]
malignant_cell = args[4]
 
ecotypes_dir = file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery") 
output_dir = file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery) 
dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))

ecotypes = read.delim(file.path(ecotypes_dir, "ecotypes.txt"))
ecotype_data = read.delim(file.path(output_dir, "ecotype_abundances_long.txt")) 
ecotype_data$Ecotype = ecotype_to_factor(ecotype_data$Ecotype)

tb = as.data.frame(table(ecotypes$Ecotype))
ecotype_data$Abundance = ecotype_data$Abundance * tb[match(as.character(ecotype_data$Ecotype), tb[,1]), 2]
ecotype_data$Abundance = ecotype_data$Abundance / quantile(ecotype_data$Abundance, .99, na.rm = T)
	
cat(paste0("Plotting ecotype heatmaps...\n"))

spatial_heatmap_hexagon_confocal <- function(data, columns, colors, name = name, spot_col = "black", newpage = F, ...){
	
	mat = dcast(data, X ~ Y, value.var = columns[1]) 
	rn = as.character(mat[,1])
	rownames(mat) = paste0("R", as.character(mat[,1]))
	mat = mat[,-1]
	mat = rbind(rep(NA, ncol(mat)), mat)
	mat = rbind(rep(NA, ncol(mat)), mat)
	mat = rbind(mat, rep(NA, ncol(mat)))
	mat = rbind(mat, rep(NA, ncol(mat)))
	rownames(mat) = paste0("R", rownames(mat))
	rownames(mat)[3:(nrow(mat) - 2)] = rn

	mat = cbind(rep(NA, nrow(mat)), mat)
	mat = cbind(rep(NA, nrow(mat)), mat)
	mat = cbind(mat, rep(NA, nrow(mat)))
	mat = cbind(mat, rep(NA, nrow(mat)))

	n = nrow(mat) 
	m = ncol(mat) 
	width = unit(2, "in") 
	height = unit(4 * (3/4) * (2*sqrt(3.0)/3) * (n + 1/4) / ((m + 0.5 / (2*sqrt(3.0)/3))) , "in") 

	hmp = Heatmap(as.matrix(mat), name = paste0(paste(columns, collapse = "_"), name), cluster_rows = F, cluster_columns = F,
		na_col = "white", col = colorRamp2(breaks =  c(0, 1), colors = c("white", "white")),
		show_row_names = F, show_column_names = F, show_heatmap_legend = F,
		width = width,
		height = height,
		use_raster = T, raster_quality= 1,
		cell_fun = function(jj, ii, x, y, w, h, col) {
			ii <<- ii
			jj <<- jj	
			if(is.na(mat[ii,jj]))
			{ 
				return()
			}	
			
	    	for(column in columns)
	    	{
	    		col = colors[which(columns == column)]
	    		alpha = data[as.character(data$Y) == (colnames(mat)[jj]) & as.character(data$X) == (rownames(mat)[ii]),column]
	    			    			    		
	    		if(length(alpha) == 0)
	    		{
	    			next 
	    		}
	    		
	    		cX = eval(parse(text = x))
	    		cY = eval(parse(text = y))
	    		radius = eval(parse(text = w)) * 0.99 
	    			    		
	    		angle = 2 * pi / 6
	    		lx = c()
	    		ly = c()
	    		for(pos in 0:5)
				{
				    px = cX + radius * (2*sqrt(3.0)/3) * sin(pos * angle)
				    py = cY + radius * cos(pos * angle)
				    lx = c(lx, px)
				    ly = c(ly, py)
				}
	    		grid.polygon(x = lx, y = ly, gp = gpar(col = col, fill = col, alpha = alpha, lwd = 0))	    		
	    	}
	    
		}, ...)

	draw(hmp, newpage = newpage)

	decorate_heatmap_body(paste0(paste(columns, collapse = "_"), name), {
	grid.rect(x = unit(1/2, "npc"), y = unit(1/2, "npc"),
          width = width, height = height,
          gp = gpar(col = "black", lwd = 1, fill = NA))
	})
	return(1)
}

background = "#CCCCCC"
charcoal = "#36454f"

splits = split(ecotype_data, ecotype_data$Ecotype)

maxc = 5
nr = ceiling(length(splits)/maxc)
if(nr == 1)
{
	nc = length(splits)
}else{
	nc = maxc
}

pdf(file.path(output_dir, paste0("Ecotype", "_spatial_heatmaps.pdf")), height = 2.6 * nr + 1, width = nc * 2.3)		
pushViewport(viewport(layout = grid.layout(nc = 1, nr = 2, heights = c(unit(2.6 * nr, "in"), unit(0.25, "in")))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(layout = grid.layout(nc = nc, nr = nr)))

r = 1
c = 1

for(spl in splits)
{ 
	spl$Malignant = spl$Malignant - 0.1
	spl$Malignant = ifelse(spl$Malignant < 0.1, 0.1, spl$Malignant)
	spl$Malignant = spl$Malignant / max(spl$Malignant, na.rm = T)
		
	column = "Abundance"
	columns = c("Malignant", column)
	
	colors = c(background, charcoal)

	title = as.character(spl$Ecotype)[1]
	
	spl[,columns][is.na(spl[,columns])] = 0
	
	spl[,column] = spl[,column] / 0.9	
	spl[,column] = ifelse(spl[,column] > 1, 1, spl[,column])
				
	pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))	
	spatial_heatmap_hexagon_confocal(spl, columns, colors, name = title, column_title = title, spot_col = "black")
	upViewport() 
		
	c = c + 1
	if(c > maxc)
	{
		c = 1
		r = r + 1
	}
}	
upViewport() 
upViewport() 
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
g <- Heatmap(as.matrix(c(paste0("Cell of origin (", malignant_cell, ")"), "Ecotype")), col = c(background, charcoal), name = "hmap", width = 0, height = 0, 
	heatmap_legend_param = list(title = NULL, legend_direction = "horizontal", nrow = 1, title_position = "leftcenter"))
draw(g, newpage = F, heatmap_legend_side = "bottom")
upViewport() 
tmp = dev.off()	

png(file.path(output_dir, paste0("Ecotype", "_spatial_heatmaps.png")), height = 2.8 * nr + 1, width = nc * 2.3, units = "in", res = 300)		
pushViewport(viewport(layout = grid.layout(nc = 1, nr = 2, heights = c(unit(2.6 * nr, "in"), unit(0.25, "in")))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(layout = grid.layout(nc = nc, nr = nr)))

r = 1
c = 1

for(spl in splits)
{ 
	spl$Malignant = spl$Malignant - 0.1
	spl$Malignant = ifelse(spl$Malignant < 0.1, 0.1, spl$Malignant)
	spl$Malignant = spl$Malignant / max(spl$Malignant, na.rm = T)
		
	column = "Abundance"
	columns = c("Malignant", column)
	
	colors = c(background, charcoal)

	title = as.character(spl$Ecotype)[1]
	
	spl[,columns][is.na(spl[,columns])] = 0
	
	spl[,column] = spl[,column] / 0.9	
	spl[,column] = ifelse(spl[,column] > 1, 1, spl[,column])
				
	pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))	
	spatial_heatmap_hexagon_confocal(spl, columns, colors, name = title, column_title = title, spot_col = "black")
	upViewport() 
		
	c = c + 1
	if(c > maxc)
	{
		c = 1
		r = r + 1
	}
}	
upViewport() 
upViewport() 
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
g <- Heatmap(as.matrix(c(paste0("Cell of origin (", malignant_cell, ")"), "Ecotype")), col = c(background, charcoal), name = "hmap", width = 0, height = 0,
	heatmap_legend_param = list(title = NULL, legend_direction = "horizontal", nrow = 1, title_position = "leftcenter"))
draw(g, newpage = F, heatmap_legend_side = "bottom")
upViewport() 

tmp = dev.off()	



