import scvelo as scv
import pandas as pd
import os
import numpy as np

os.chdir(r"C:/Users/USER/Documents/R/RNAseq/scientificProject/data/Scadden/Magic_based")
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
adata.var_names = namez
scv.pp.moments(adata, n_neighbors=15)


s = scv.read(seur_loc + 'spliced_magic.csv', cache=True, first_column_names=True)
u = scv.read(seur_loc + 'unspliced_magic.csv', cache=True, first_column_names=True)

adata.layers['Ms'] = s.X
adata.layers['Mu'] = u.X

Ms = adata.layers["Ms"]
for col in range(Ms.shape[1]):
    gene = Ms[:, col]
    ss = gene[gene < 0].mean()
    Ms[:, col] = gene - ss
adata.layers["Ms"] = Ms

Mu = adata.layers["Mu"]
for col in range(Mu.shape[1]):
    gene = Mu[:, col]
    ss = gene[gene < 0].mean()
    Mu[:, col] = gene - ss
adata.layers["Mu"] = Mu


scv.tl.recover_dynamics(adata, fit_steady_states=False)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)
scv.tl.velocity_pseudotime(adata)
scv.tl.velocity_confidence(adata)
scv.tl.terminal_states(adata)
scv.tl.velocity_clusters(adata)
scv.tl.umap(adata)
adata.write_h5ad(filename=r"Magic_computed.h5ad")

adata = scv.read(filename=r"Magic_computed.h5ad")

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


scv.pl.velocity(adata, "Msn", dpi=800)

scv.tl.rank_velocity_genes(adata, min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()


meta_seur = pd.read_csv(seur_loc + 'meta_seur.csv', index_col=0)

adata.obs['seurat_clusters'] = [str(x) + "cluster" for x in meta_seur.seurat_clusters.to_list()]
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

adata.obsm['X_umap_wnn'] = np.asanyarray(meta_seur[["wnnUMAP_1", "wnnUMAP_2"]])

kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
top_genes = top_genes.to_list()
scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby='seurat_clusters')
scv.pl.scatter(adata, basis=top_genes[0:20], ncols=5, add_outline='fit_diff_kinetics', dpi=1000, show=False,
               save="fitdifkin.png")


scv.get_df(adata[:, top_genes], ['fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
scv.pl.scatter(adata, y=top_genes[:9], size=10, ncols=3, fontsize=10, x='velocity_pseudotime', alpha=0.8,
               legend_loc='none', frameon=False, color='latent_time',
               dpi=1000, save="scatter_velocity_pseudotime_top9.png")

scv.pl.velocity_embedding_stream(adata, basis="umap_wnn", color='latent_time', dpi=800, color_map='gnuplot')


x, y = scv.utils.get_cell_transitions(adata, basis='umap_wnn', starting_cell=1103, n_steps=200)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False, basis="umap_wnn", threshold=0.1)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
