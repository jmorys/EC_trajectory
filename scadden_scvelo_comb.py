import scvelo as scv
import pandas as pd
import numpy as np
from scipy import sparse
import os

os.chdir(r"C:\Users\USER\Documents\R\RNAseq\scientificProject\data\combined")

scv.settings.set_figure_params('scvelo')

s = scv.read('sub_spliced.csv', cache=True, first_column_names=True)
s = s.transpose()
u = scv.read('sub_unspliced.csv', cache=True)
u = u.transpose()

adata = s
adata.layers['spliced'] = s.X
adata.layers['unspliced'] = u.X

m = pd.read_csv("combined_meta.csv")
adata.obs['cell_cluster'] = list(m['seurat_clusters'])
adata.obs['cell_cluster'] = adata.obs['cell_cluster'].astype('category')

UMAP_D = m[['UMAP_1', 'UMAP_2']]
adata.obsm['X_umap_ori'] = np.asanyarray(UMAP_D)


del u
del s
del m

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

scv.pp.moments(adata, n_pcs=20, n_neighbors=50)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)
scv.tl.umap(adata)

scv.tl.velocity_clusters(adata)
adata.write_h5ad(filename='S_scvelo_EC.h5ad')

adata = scv.read("S_scvelo_EC.h5ad")
scv.tl.velocity_pseudotime(adata)
scv.tl.velocity_confidence(adata)
scv.tl.terminal_states(adata)
adata.write_h5ad(filename='S_scvelo_EC.h5ad')

scv.pl.velocity_embedding_stream(adata, basis="umap", color='cell_cluster',
                                 dpi=1000, save="velocity_umap_cell_cluster.png")
scv.pl.velocity_embedding_stream(adata, basis="umap_ori", color='cell_cluster',
                                 dpi=1000, save="velocity_umap_ori_cell_cluster.png")

scv.pl.velocity_embedding_stream(adata, basis="umap", color='velocity_clusters', smooth=0.4, min_mass=2, dpi=1000,
                                 save="velocity_umap_velocity_clusters.png")
scv.pl.velocity_embedding_stream(adata, basis="umap_ori", color='velocity_clusters', smooth=0.4, min_mass=2, dpi=1000,
                                 save="velocity_umap_ori_velocity_clusters.png")

scv.pl.velocity_embedding(adata, basis="umap", color='velocity_clusters', arrow_length=2, arrow_size=1.5, dpi=1000,
                          save="velocity_umap_velocity_clusters3.png")


scv.pl.velocity_embedding_stream(adata, basis="umap", color='root_cells',
                                 dpi=1000, save="velocity_umap_root_cells.png")
scv.pl.velocity_embedding_stream(adata, basis="umap", color='end_points',
                                 dpi=1000, save="velocity_umap_end_points.png")

scv.pl.velocity_embedding_stream(adata, basis="umap_ori", color='latent_time',
                                 dpi=1000, save="velocity_umap_ori_latent_time.png")

scv.pl.velocity_embedding_stream(adata, basis="umap_ori", color='velocity_pseudotime',
                                 dpi=1000, save="velocity_umap_ori_velocity_pseudotime.png")

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]

scv.pl.heatmap(adata, var_names=top_genes, tkey='velocity_pseudotime',
               n_convolve=100, col_color='cell_cluster', save="heatmap_velocity_pseudotime_top300.png")
scv.pl.heatmap(adata, var_names=top_genes[:9], tkey='velocity_pseudotime',
               n_convolve=100, col_color='cell_cluster', save="heatmap_velocity_pseudotime_top9.png")

scv.pl.scatter(adata, basis=top_genes[:9], size=10, ncols=3, fontsize=10, color='latent_time', alpha=0.8,
               legend_loc='none', frameon=False,
               dpi=1000, save="scatter_top9.png")

scv.pl.scatter(adata, y=top_genes[:9], size=10, ncols=3, fontsize=10, x='velocity_pseudotime', alpha=0.8,
               legend_loc='none', frameon=False, color='latent_time',
               dpi=1000, save="scatter_velocity_pseudotime_top9.png")

scv.pl.scatter(adata, y=top_genes[:9], size=10, ncols=3, fontsize=10, x='latent_time', alpha=0.8, legend_loc='none',
               frameon=False, color='latent_time',
               dpi=1000, save="scatter_latent_time_top9.png")

scv.pl.scatter(adata, color='velocity_confidence', perc=[2, 98],
               dpi=1000, save="velocity_umap_velocity_confidence.png")

scv.tl.tsne(adata)
scv.pl.velocity_embedding_stream(adata, basis="tsne", color='cell_cluster')

c_data = pd.DataFrame(adata.obs[['velocity_pseudotime', 'velocity_clusters']])
c_data.to_csv("cells_data.csv")


trans_graph = adata.uns['velocity_graph']
sparse.save_npz("trans_graph", trans_graph)
