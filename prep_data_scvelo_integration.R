library(Seurat)
library(tidyverse)
library(data.table)
library(magrittr)
source('~/R/RNAseq/finetuned/read_tools.R')

setwd("C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/R_references/Seurat_integration")

layer_files = dir()
layer_files = grep("matrix.csv", layer_files, value = T)
layer_names = gsub("_matrix.csv","",layer_files)

metadata = read.csv("all_metadata.csv", header = T, row.names = "X")
varnames = read.csv("varnames.csv")


layers = list()
for (i in 1:length(layer_files)) {
  matrix <- fread(layer_files[i], header = T)
  matrix = transpose(matrix[,-1])
  matrix <- as.data.frame(matrix)
  rownames(matrix) = varnames[,1]
  colnames(matrix) = rownames(metadata)
  seurat_obj <- CreateSeuratObject(matrix, meta.data = metadata)
  layers[[i]]  <- seurat_obj
}

integrated <- list()
for (i in 1:length(layer_names)) {
  split <- SplitObject(layers[[i]], split.by = "dataset")
  integ <- integration_basic(split, names(split), number = 6000)
  integrated[[i]] <- integ
}
remove(list = c("integ", "split", "matrix","seurat_obj"))
names(integrated) <- layer_names
saveRDS(integrated, "integrated_list.rds")



var_list <- list()
for (i in names(integrated)) {
  var_list[[i]] <- integrated[[i]]@assays$integrated@var.features
}
vars_to_transfer <- Reduce(intersect, var_list[c("unspliced", "spliced")] )

twice_integrated <- integrated[[1]] 
twice_integrated <- RenameAssays(twice_integrated, integrated = "ambiguous")
twice_integrated[["RNA"]] <- NULL
twice_integrated[["SCT"]] <- NULL
twice_integrated[["unspliced"]] <- integrated[["unspliced"]]@assays$integrated
twice_integrated[["spliced"]] <- integrated[["spliced"]]@assays$integrated

remove(integrated)

for (i in layer_names) {
  twice_integrated <- twice_integrated %>% RunPCA(reduction.name = paste0(i, "pca"), assay = i, verbose = F)
}


twice_integrated <- FindMultiModalNeighbors(twice_integrated, reduction.list = paste0(c("spliced", "unspliced"), "pca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "Splicing.weight", weighted.graph = T, k.nn = 50
)

twice_integrated <- RunUMAP(twice_integrated, nn.name = "weighted.nn",
                            reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                            umap.method = "uwot-learn", repulsion.strength = 0.8, local.connectivity = 3,
                            negative.sample.rate = 8, learning.rate = 0.7)

DimPlot(twice_integrated, reduction = "wnn.umap", cols = as.factor(twice_integrated@meta.data$seurat_clusters))

saveRDS(twice_integrated, "twice_integrated.rds")

fwrite(as.data.frame(
  twice_integrated@assays$spliced@scale.data[vars_to_transfer,]
  ),row.names = T, file = "spliced_seur.csv")

fwrite(as.data.frame(
  twice_integrated@assays$unspliced@scale.data[vars_to_transfer,]
  ), row.names = T, file = "unspliced_seur.csv")

meta <- twice_integrated@meta.data
umap <- twice_integrated@reductions$wnn.umap@cell.embeddings
meta <- meta %>% bind_cols(as_tibble(umap))

fwrite(meta, row.names = T, file = "meta_seur.csv")
fwrite(twice_integrated@neighbors$weighted.nn@nn.idx, row.names = T, file = "neighbors_seur.csv")
fwrite(twice_integrated@neighbors$weighted.nn@nn.dist, row.names = T, file = "dist_seur.csv")




setwd("~/R/RNAseq")
