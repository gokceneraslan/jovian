# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # jovian - 03 Differential Expression

# %%
import scanpy as sc

# Plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# numpy et al.
import numpy as np
import scipy.sparse as sp
import pandas as pd

# R integration
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, numpy2ri, r
from rpy2.robjects.vectors import StrVector, FloatVector, ListVector
import rpy2.robjects as ro
import anndata2ri

from pathlib import Path
import math
from tqdm.auto import tqdm
import warnings
import shelve

# %%
# %load_ext rpy2.ipython

# %%
sc.set_figure_params(dpi=100, fontsize=12)
matplotlib.rcParams['font.sans-serif'] = matplotlib.rcParamsDefault['font.sans-serif']

sc.settings.verbosity = 'hint'

# %% [markdown]
# ## Load the session

# %%
with shelve.open('session.pkl', protocol=4, writeback=False) as db:
    for k in db.keys():
        globals()[k] = db[k]

# %%
del conf_samples, conf_samples_processed

# %% [markdown]
# ***

# %% [markdown]
# ## Parameters

# %%
par_de_group = 'leiden'
par_de_n_genes = 2000
par_de_method = 't-test_overestim_var'

par_per_group_de = True
par_group_key = 'tissue'

par_save_filename_de = 'de-genes.xlsx'
par_save_filename_de_group = 'de-genes-%s.xlsx'

# %% [markdown]
# ***

# %% [markdown]
# ## Read dataset

# %%
adata = sc.read(par_save_filename)
adata

# %% [markdown]
# ## Differential expression

# %%
# %%time

sc.tl.rank_genes_groups(adata, par_de_group, n_genes=par_de_n_genes, method=par_de_method)

# %% [markdown]
# ## Visualization

# %%
umap_point_size = np.maximum(120000/adata.n_obs, 2)

# %%
f, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color='leiden', legend_loc='on data', legend_fontoutline=3, legend_fontsize=14, legend_fontweight='normal', title='Clusters', ax=ax, show=False, size=umap_point_size);

# %%
sc.pl.rank_genes_groups_dotplot(adata, standard_scale='var', mean_only_expressed=True)

# %%
sc.pl.rank_genes_groups(adata, n_genes=30, ncols=6, fontsize=7, sharey=False)


# %% [markdown]
# ## Cell typing

# %%
def predict_cell_types(adata, use_raw=True, species='mouse', cluster_key='leiden', numCores=1, **kwds):

    if species == 'mouse':
        ref_names = ['MouseRNAseqData','ImmGenData']
    else:
        ref_names = ['HumanPrimaryCellAtlasData', 'MonacoImmuneData']

    s = importr('SingleR')
    refs = [s.__dict__[ref_name]() for ref_name in ref_names]
    ref_genes = StrVector(set.intersection(*[set(r['rownames'](d)) for d in refs]))
    ref = r('cbind')(*[r('`[`')(d, ref_genes) for d in refs]) # merge references

    ad = adata.raw if use_raw else adata
    obs = adata.obs

    common_genes = sorted(list(set(ad.var_names) & set(ref_genes)))
    assert len(common_genes) > 10000, 'Not enough genes overlapping with ref SingleR datasets...'

    ref = r('`[`')(ref, ro.vectors.StrVector(common_genes))

    mat = ad[:, common_genes].X.T.copy()

    if sp.issparse(mat):
        mat = anndata2ri.scipy2ri.py2rpy(sp.csr_matrix(mat))
    else:
        mat = numpy2ri.py2rpy(mat)

    mat = r("`rownames<-`")(mat, ro.vectors.StrVector(ad[:, common_genes].var_names))
    mat = r("`colnames<-`")(mat, ro.vectors.StrVector(obs.index)) # TODO: really needed?
    clusters = ro.vectors.StrVector(obs[cluster_key].values.tolist())

    if numCores > 1:
        par = r('BiocParallel::MulticoreParam')(workers = numCores)
        kwds['BPPARAM'] = par

    labels = s.SingleR(test=mat,
                       ref=ref,
                       labels = r('`$`')(ref, 'label.fine'), # use label.main too
                       method='cluster',
                       clusters=clusters, **kwds)

    labels = pandas2ri.rpy2py(r('as.data.frame')(labels))

    adata.obs['predicted_cell_types'] = ''
    for cluster in adata.obs[cluster_key].cat.categories:
        adata.obs.loc[adata.obs[cluster_key] == cluster, 'predicted_cell_types'] = labels.loc[cluster]['pruned.labels']

    adata.uns['cell_type_prediction'] = labels['pruned.labels'].to_dict()


# %% [markdown]
# ## Save markers

# %%
def save_markers(adata, filename, group_key='leiden'):

    with pd.ExcelWriter(filename) as writer:

        for cluster in tqdm(adata.obs[group_key].cat.categories):
            marker_df = sc.get.rank_genes_groups_df(adata, cluster)
            expr = adata.raw[adata.obs[group_key] == cluster][:, marker_df.names.values].X.A
            marker_df['percent_cell_expressed'] = (expr != 0).sum(0) / expr.shape[0]
            marker_df['mean_normalized_expression'] = np.ma.array(expr, mask=(expr==0)).mean(0).filled()

            if 'cell_type_prediction' in adata.uns_keys():
                pred = adata.uns["cell_type_prediction"][cluster]
                pred = pred.translate(str.maketrans({'\\': '_',
                                                     '/': '_',
                                                     '*': '_',
                                                     '[': '_',
                                                     ']': '_',
                                                     ':': '_',
                                                     '?': '_'}))
                sheet_name = f'{cluster} ({pred})'
            else:
                sheet_name = cluster
            marker_df.to_excel(writer, sheet_name=sheet_name, index=False)


# %%
save_markers(adata, par_save_filename_de)

# %% [markdown]
# ## Differential expression per group

# %%
# %%time

from IPython.core.display import display, HTML
sc.settings.verbosity = 0

if par_per_group_de:
    ascat = pd.Categorical(adata.obs[par_group_key])
    for group in tqdm(ascat.categories):
        display(HTML(f'<h1>{par_group_key.capitalize()}: {group}</h1>'))

        ad = sc.read(par_save_filename_group % group)
        sc.tl.rank_genes_groups(ad, par_de_group, n_genes=par_de_n_genes, method=par_de_method)

        f, ax = plt.subplots(figsize=(6, 6))
        sc.pl.umap(ad, color='leiden', legend_loc='on data', legend_fontoutline=3, legend_fontsize=14, legend_fontweight='normal', title='Clusters', ax=ax, show=False, size=umap_point_size);

        sc.pl.rank_genes_groups_dotplot(ad, standard_scale='var', mean_only_expressed=True)
        sc.pl.rank_genes_groups(ad, n_genes=30, ncols=6, fontsize=7, sharey=False)

        predict_cell_types(ad, species=par_species)

        f, ax = plt.subplots(figsize=(6, 6))
        sc.pl.umap(ad, color='predicted_cell_types', title='Predicted cell types', ax=ax, show=False, size=umap_point_size);
        display(f)

        save_markers(ad, par_save_filename_de_group % group)
        ad.write(par_save_filename_group % group)

# %% [markdown]
# ## Transfer predicted annotations

# %%
adata.obs['predicted_cell_types'] = ''
adata.obs['sample_leiden'] = ''

# %%
for group in adata.obs[par_group_key].astype('category').cat.categories:
    ad = sc.read(par_save_filename_group % group)
    adata.obs['predicted_cell_types'].update(ad.obs['predicted_cell_types'].astype(str))
    adata.obs['sample_leiden'].update(ad.obs['leiden'].astype(str))

# %%
adata.obs

# %%
adata.obs['predicted_cell_types_coarse'] = [x.split(':')[0] for x in adata.obs['predicted_cell_types']]

# %%
plt.rcParams['figure.figsize'] = (7, 7)
sc.pl.umap(adata, color=['predicted_cell_types_coarse',
                         'predicted_cell_types'], ncols=1, legend_loc='on data', legend_fontsize=7, legend_fontweight='normal', legend_fontoutline=2);

# %%
sc.tl.draw_graph(adata)

# %%
plt.rcParams['figure.figsize'] = (7, 7)
sc.pl.draw_graph(adata, color=['predicted_cell_types_coarse',
                               'predicted_cell_types'], ncols=1, legend_loc='on data', legend_fontsize=7, legend_fontweight='normal', legend_fontoutline=2);

# %%
adata.write(par_save_filename)

# %% [markdown]
# ## Serialize session

# %%
k = None
var = None

with shelve.open('session.pkl', protocol=4) as db:
    for k, var in globals().items():
        if k.startswith('par_') or k.startswith('conf_'):
            print(f'Storing {k}...')
            db[k] = var
db.close()
