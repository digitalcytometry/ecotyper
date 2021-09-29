library(Matrix)
args = commandArgs(T)
#args = c("VisiumBreast", "../example_data/VisiumBreast") 

id = args[1]
matrix_dir = args[2]

barcode.path = file.path(matrix_dir, "barcodes.tsv.gz")
features.path = file.path(matrix_dir, "features.tsv.gz")
matrix.path = file.path(matrix_dir, "matrix.mtx.gz")
positions.path = file.path(matrix_dir, "tissue_positions_list.csv")

if(!file.exists(barcode.path))
{
	stop(paste0("File 'barcodes.tsv.gz' is missing from the input dir '", matrix_dir, "'"))
}
if(!file.exists(features.path))
{
	stop(paste0("File 'features.tsv.gz' is missing from the input dir '", matrix_dir, "'"))
}
if(!file.exists(matrix.path))
{
	stop(paste0("File 'matrix.mtx.gz' is missing from the input dir '", matrix_dir, "'"))
}
if(!file.exists(positions.path))
{
	stop(paste0("File 'tissue_positions_list.csv' is missing from the input dir '", matrix_dir, "'"))
}

output_dir = file.path("../datasets/visium/", id)
dir.create(output_dir, recursive = T, showWarning = F)

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

data = as.matrix(mat)
data = data.frame(Gene = feature.names[match(rownames(data), feature.names$V1),]$V2, data)
data = data[!duplicated(data$Gene),]

spatial_coords = read.delim(positions.path, sep = ",", header = F)
spatial_coords$ID = paste0(make.names(spatial_coords[,1]))

spatial_coords = spatial_coords[match(colnames(data)[-1], make.names(spatial_coords[,1])),]

colnames(data)[-1] = spatial_coords$ID

clinical = data.frame(ID = colnames(data)[-1])
clinical$X = as.numeric(as.character(spatial_coords[,3]))
clinical$Y = as.numeric(as.character(spatial_coords[,4]))

write.table(data, file.path(output_dir, "data.txt"), row.names = F, sep = "\t")
write.table(clinical, file.path(output_dir, "annotation.txt"), row.names = F, sep = "\t")

