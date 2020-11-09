import scvelo as scv
import pandas as pd
import anndata
import os
import re
import sys
from pathlib import Path
import numpy as np
import cellrank as cr
import scipy


def subset_anndata(loom, csvfile):
    adata = scv.read(loom)
    adata.var_names_make_unique()
    adata.obs.index = [re.search("[ACTG]{6,}", x).group(0) for x in adata.obs.index]
    my_pd = pd.read_csv(csvfile + ".csv", index_col=0, encoding="utf8")
    intersected = adata.obs.index.intersection(my_pd.index)
    adata = adata[intersected, :]
    return adata


os.chdir(r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/VeloData")
csv_loc = r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/R_references/subsetting_EC/"
scv.settings.set_figure_params('scvelo')
file_list = os.listdir()

con_dir = {}
for file in file_list:
    name = re.sub("_.+", "", file)
    con_dir[name] = subset_anndata(file, csv_loc + name)

concat = anndata.concat(con_dir, axis=0, label="dataset")
path = Path(r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/" + r"Concat_raw.h5ad")
concat.write_h5ad(filename=path,)
del concat
del con_dir
del file_list


loc = r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/Raw_based/"
os.chdir(loc)
adata = scv.read(path)
adata.obs.dataset = [x for x in adata.obs.dataset]
new_index = []
for ob in range(len(adata.obs.index)):
    cell = adata.obs.index[ob]
    dataset = adata.obs["dataset"][ob]
    n_ind = cell + "_" + dataset
    new_index.append(n_ind)
adata.obs.index = new_index

# adding seurat data
csv_loc = r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/R_references/"
seurat_meta = pd.read_csv(csv_loc + "seurat_meta.csv", index_col=0)

seurat_meta[["origin", "cells", "clusters"]] = seurat_meta[["origin", "cells", "clusters"]].astype("string")
seurat_meta["rec_cells"] = seurat_meta["cells"].str.cat(seurat_meta["origin"], sep="_")
rec_cels = seurat_meta["rec_cells"].to_list()
seurat_meta.index = rec_cels
seurat_meta = seurat_meta.loc[adata.obs.index.values]

adata.obs['seurat_clusters'] = seurat_meta.clusters.to_list()
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

adata.obsm['X_umap_seurat'] = np.asanyarray(seurat_meta[["UMAP_1", "UMAP_2"]])
adata.obsm['X_PHATE_seurat'] = np.asanyarray(seurat_meta[["PHATE1", "PHATE2"]])


# saveing data for R
csv_loc = r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/R_references/Seurat_integration/"
layers = ["spliced", "unspliced", "ambiguous"]
for layer in layers:
    layer_matrix = scipy.sparse.csc_matrix.todense(adata.layers[layer])
    layer_matrix = pd.DataFrame(data=layer_matrix)
    layer_matrix.index = adata.obs.index
    layer_matrix.to_csv(csv_loc + layer + "_matrix.csv")

adata.obs.to_csv(csv_loc + "all_metadata.csv")
adata.var_names.to_series().to_csv(csv_loc + "varnames.csv")

# basic scvelo computation
scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=25, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)
scv.tl.umap(adata)
scv.tl.velocity_pseudotime(adata)
scv.tl.velocity_confidence(adata)
scv.tl.terminal_states(adata)
scv.tl.velocity_clusters(adata)
adata.write_h5ad(filename=Path(loc + r"Raw_computed.h5ad"))


adata = scv.read(Path(loc + r"Raw_computed.h5ad"))


scv.pl.velocity_embedding_stream(adata, basis="PHATE_seurat", color='latent_time', density=1,
                                 dpi=200, n_neighbors=3, min_mass=1)

scv.pl.velocity_embedding_stream(adata, basis="umap_seurat", color='seurat_clusters', density=1,
                                 dpi=200, min_mass=0)

scv.pl.velocity_embedding_stream(adata, basis="umap", color='seurat_clusters', dpi=300, min_mass=0, density=1)
scv.pl.velocity_embedding_stream(adata, basis="umap", color='latent_time', dpi=300, min_mass=0,
                                 density=1, color_map='gnuplot')

scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[30, 95], dpi=300)

scv.pl.velocity_embedding_stream(adata, basis="PHATE_seurat", color='root_cells',
                                 dpi=300, min_mass=0, density=1, n_neighbors=10)

scv.pl.velocity_embedding_stream(adata, basis="umap", color='root_cells', dpi=300, min_mass=0,
                                 density=1, color_map='gnuplot')


scv.pl.velocity_embedding(adata, basis="umap", color='velocity_pseudotime', dpi=800, color_map='gnuplot')


scv.pl.velocity_embedding_stream(adata, basis="umap_seurat", color='dataset', dpi=800, min_mass=0,
                                 density=1, color_map='gnuplot', size=10, alpha=1)


x, y = scv.utils.get_cell_transitions(adata, basis='PHATE_seurat', starting_cell=1000)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False, basis="PHATE_seurat")
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, basis='PHATE_seurat')



adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='seurat_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


scv.pl.paga(adata, basis='umap_seurat', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, dpi=800)