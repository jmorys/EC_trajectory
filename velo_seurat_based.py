import scvelo as scv
import pandas as pd
import anndata
import os
import re
import sys
from pathlib import Path
import numpy as np
import scipy

os.chdir(r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/Seurat_based")
seur_loc = r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/R_references/Seurat_integration/"

namez = pd.read_csv(seur_loc + 'spliced_seur.csv', index_col=0).index
namez = namez.to_list()
s = scv.read(seur_loc + 'spliced_seur.csv', cache=True, first_column_names=True, index_col=0)
s = s.transpose()
u = scv.read(seur_loc + 'unspliced_seur.csv', cache=True, first_column_names=True)
u = u.transpose()

adata = s
adata.layers['spliced'] = s.X
adata.layers['unspliced'] = u.X

n = pd.read_csv(seur_loc + 'neighbors_seur.csv', index_col=0)
n = n-1
adata.uns.neighbors = n
adata.var_names = namez
scv.pp.moments(adata, n_neighbors=15)
scv.tl.recover_dynamics(adata, fit_steady_states=False, fit_connected_states=True)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)
scv.tl.velocity_pseudotime(adata)
scv.tl.velocity_confidence(adata)
scv.tl.terminal_states(adata)
scv.tl.velocity_clusters(adata)
scv.tl.umap(adata)
adata.write_h5ad(filename=r"Seur_computed.h5ad")


scv.pl.velocity_embedding(adata, basis="umap", color='velocity_pseudotime', dpi=800, color_map='gnuplot')
scv.pl.velocity_embedding_stream(adata, basis="umap", color='velocity_pseudotime', dpi=800, color_map='gnuplot')
scv.pl.velocity_embedding_stream(adata, basis="umap", color='latent_time', dpi=800, color_map='gnuplot')

scv.pl.velocity_embedding_stream(adata, basis="umap", color='root_cells', dpi=300, min_mass=0,
                                 density=1, color_map='gnuplot')


scv.tl.paga(adata, groups='velocity_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')



scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, dpi=800)


scv.pl.velocity(adata, "Ly6c1", dpi=800)

scv.tl.rank_velocity_genes(adata, min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()