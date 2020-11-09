library(Seurat)
library(stringr)
library(phateR)



location <- "~/R/RNAseq/finetuned/RNAseq_scadden/whole_BM/"
setwd(location)
files <- grep("sub.+.rds",dir(), value = T)

outut_loc = "~/R/RNAseq/scientificProject/data/Scadden/R_references/subsetting_EC/"

for (i in 1:length(files)) {
  substd <- readRDS(paste0(location, files[i]))
  name <- stringr::str_extract(files[i],"_(\\w+)_") %>% str_remove_all("_")
  write.csv(paste0(outut_loc, name, ".csv"), x =  substd$seurat_clusters)
}
remove(substd)

scadden_preprocessed_for_EC <- readRDS("~/R/RNAseq/finetuned/RNAseq_scadden/whole_BM/scadden_preprocessed_for_EC.rds")

origin <- scadden_preprocessed_for_EC$origin_f
origin <- origin %>% str_extract("_(\\w+)_") %>% str_remove_all("_")

cells <-  names(scadden_preprocessed_for_EC$orig.ident)
cells <- cells %>% str_extract("[ACTG]{5,}")
clusters <- scadden_preprocessed_for_EC$seurat_clusters

umap <- scadden_preprocessed_for_EC@reductions$umap@cell.embeddings

phate <- phate(scadden_preprocessed_for_EC@reductions$pca@cell.embeddings[,1:20],
      gamma = -1,mds.solver = "smacof", knn = 50, decay = 200, t = 50)

meta <- tibble(origin, cells, clusters)
meta <- meta %>% bind_cols(as_tibble(umap),as_tibble(phate$embedding))


outut_loc = "~/R/RNAseq/scientificProject/data/Scadden/R_references/"

write.csv(paste0(outut_loc, "seurat_meta.csv"), x =  meta)



