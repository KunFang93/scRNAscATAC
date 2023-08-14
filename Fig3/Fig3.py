import scvelo as scv
import scanpy as sc

datdir = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/scvelo'
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
adata = sc.read_h5ad(f'{datdir}/my_data.h5ad')
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')
# load loom files for spliced/unspliced matrices for each sample:
ldata_pt1 = scv.read(f'{datdir}/PT1.loom', cache=True)
ldata_pt2 = scv.read(f'{datdir}/PT2.loom', cache=True)
ldata_pt3 = scv.read(f'{datdir}/PT5.loom', cache=True)
ldata_rt1 = scv.read(f'{datdir}/RT3.loom', cache=True)
ldata_rt2 = scv.read(f'{datdir}/RT4.loom', cache=True)
ldata_rt3 = scv.read(f'{datdir}/RT6.loom', cache=True)

# rename barcodes in order to merge:
pt1_barcodes = [bc.split(':')[1] for bc in ldata_pt1.obs.index.tolist()]
pt1_barcodes = ['PT1_' + bc[0:len(bc)-1] + '-1' for bc in pt1_barcodes]
pt2_barcodes = [bc.split(':')[1] for bc in ldata_pt2.obs.index.tolist()]
pt2_barcodes = ['PT2_' + bc[0:len(bc)-1] + '-1' for bc in pt2_barcodes]
pt3_barcodes = [bc.split(':')[1] for bc in ldata_pt3.obs.index.tolist()]
pt3_barcodes = ['PT3_' + bc[0:len(bc)-1] + '-1' for bc in pt3_barcodes]
rt1_barcodes = [bc.split(':')[1] for bc in ldata_rt1.obs.index.tolist()]
rt1_barcodes = ['RT1_' + bc[0:len(bc)-1] + '-1' for bc in rt1_barcodes]
rt2_barcodes = [bc.split(':')[1] for bc in ldata_rt2.obs.index.tolist()]
rt2_barcodes = ['RT2_' + bc[0:len(bc)-1] + '-1' for bc in rt2_barcodes]
rt3_barcodes = [bc.split(':')[1] for bc in ldata_rt3.obs.index.tolist()]
rt3_barcodes = ['RT3_' + bc[0:len(bc)-1] + '-1' for bc in rt3_barcodes]
ldata_pt1.obs.index = pt1_barcodes
ldata_pt2.obs.index = pt2_barcodes
ldata_pt3.obs.index = pt3_barcodes
ldata_rt1.obs.index = rt1_barcodes
ldata_rt2.obs.index = rt2_barcodes
ldata_rt3.obs.index = rt3_barcodes

# make variable names unique
ldata_pt1.var_names_make_unique()
ldata_pt2.var_names_make_unique()
ldata_pt3.var_names_make_unique()
ldata_rt1.var_names_make_unique()
ldata_rt2.var_names_make_unique()
ldata_rt3.var_names_make_unique()

# extract cancer cells
cancer_meta = adata.obs
cancer_meta['sample'] = cancer_meta['barcode'].str.split('_',expand=True)[0]
ldata_pt1_cancer = ldata_pt1[cancer_meta.loc[cancer_meta['sample']=='PT1','barcode'].values,]
ldata_pt2_cancer = ldata_pt2[cancer_meta.loc[cancer_meta['sample']=='PT2','barcode'].values,]
ldata_pt3_cancer = ldata_pt3[cancer_meta.loc[cancer_meta['sample']=='PT3','barcode'].values,]
ldata_rt1_cancer = ldata_rt1[cancer_meta.loc[cancer_meta['sample']=='RT1','barcode'].values,]
ldata_rt2_cancer = ldata_rt2[cancer_meta.loc[cancer_meta['sample']=='RT2','barcode'].values,]
ldata_rt3_cancer = ldata_rt3[cancer_meta.loc[cancer_meta['sample']=='RT3','barcode'].values,]

# concatenate the three loom
ldata = ldata_pt1_cancer.concatenate([ldata_pt2_cancer, ldata_pt3_cancer, ldata_rt1_cancer,
                                      ldata_rt2_cancer, ldata_rt3_cancer])
# merge matrices into the original adata object
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', show=False, save='_celltypes.pdf')
# compute RNA velocity using the steady-state model (stochastic option)
scv.pl.proportions(adata, groupby='seurat_clusters', show=False, save='pie.pdf')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata,n_pcs=30, n_neighbors=30)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata, n_jobs=10)

# Visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters',
                               show=False,save='embedding_grid.pdf',scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['seurat_clusters'], save='embedding_stream.pdf', show=False, title='')
# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['ESR1','ERBB2','LDHB','EDN1'],legend_loc='best',
                color='seurat_clusters', save="genes_profiling.pdf")

# Downstream
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
rank_velo_genes = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])

kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline='0,1,2,3,4,5,6,7,8,9,10,11,12')
for i in [str(i) for i in range(13)]:
    scv.pl.scatter(adata, rank_velo_genes[i][:5], ylabel=i, frameon=False, color='seurat_clusters',
                   legend_loc='right margin',size=10, linewidth=1.5, save=f"{i}_top5_veloGene.pdf")


# velocity confidence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="confidence.pdf")

# graph
scv.pl.velocity_graph(adata, threshold=.1, color='seurat_clusters',save="graph.pdf")

# draw descendents/anscestors coming from a specified cell
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell='TATTGGGCATAC')
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

# pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',save="pseudotime.pdf")

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='seurat_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save="PAGA.pdf")

adata.write_loom(f'{datdir}/adata_process.loom')

# dynamic model
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

scv.tl.recover_dynamics(adata,n_jobs=10)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata,n_jobs=10)

scv.pl.velocity_embedding_stream(adata, basis='umap', save='dynamic_embedding_stream.png',color='seurat_clusters')

adata_var_df = adata.var
adata_var_df = adata_var_df[(adata_var_df['fit_likelihood'] > .1) & adata_var_df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(adata_var_df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(adata_var_df['fit_beta'] * adata_var_df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(adata_var_df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="latent_time.pdf")

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, save='dyn_heatmap.pdf')

scv.pl.velocity(adata, ['PROCR', 'CD44'],  color='seurat_clusters',
                legend_loc='best', ncols=1, mode='dynamical')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)
#Top-likelihood genes
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False, color='seurat_clusters',save='toplk_genes.pdf')
# Cluster-specific top-likelihood genes
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
lk_df = scv.get_df(adata, 'rank_dynamical_genes/names')
lk_df.head(5)

for cluster in [str(i) for i in range(13)]:
    scv.pl.scatter(adata, lk_df[cluster][:5], ylabel=cluster, color='seurat_clusters',frameon=False, save=f'{cluster}_lk.pdf')

scv.tl.paga(adata, groups='seurat_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save="Dyn_PAGA.pdf")

adata.write(f'{datdir}/adata_dynamic.h5ad', compression='gzip')
# adata = scv.read(f'{datdir}/adata_dynamic.h5ad')
