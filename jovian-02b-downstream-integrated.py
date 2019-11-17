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
# # jovian - 02b Downstream

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
from rpy2.robjects import pandas2ri
import rpy2.robjects.numpy2ri as numpy2ri
import anndata2ri

numpy2ri.activate()
pandas2ri.activate()
anndata2ri.activate()

from pathlib import Path
import math
from tqdm.auto import tqdm
import warnings
import pickle
from urllib.request import urlopen
import shelve

# %%
# %load_ext rpy2.ipython

# %%
sc.set_figure_params(dpi=100, fontsize=12)
matplotlib.rcParams['font.sans-serif'] = matplotlib.rcParamsDefault['font.sans-serif']

sc.settings.verbosity = 'hint'

# %% [markdown]
# ***

# %% [markdown]
# ## Parameters

# %% tags=["parameters"]
par_save_filename = 'adata.h5ad'
par_save_filename_group = 'adata-group-%s.h5ad'

par_remove_doublets = True

par_generate_plots_per_group = True
par_group_key = 'tissue'

par_merge_type = 'inner'
par_batch_key = 'sample_name'

# %% [markdown]
# ***

# %% [markdown]
# ## Load session

# %%
with shelve.open('session.pkl', protocol=4, writeback=False) as db:
    for k in db.keys():
        globals()[k] = db[k]

# %%
del conf_samples

for sample, ad in conf_samples_processed.items():
    ad.X = ad.layers['counts'].copy()
    del ad.layers['counts']

    display(ad)

# %% [markdown]
# ## Merge samples

# %%
# %%time

batch_categories, ads = zip(*conf_samples_processed.items())

adata = sc.AnnData.concatenate(*ads, join=par_merge_type, batch_key=par_batch_key, batch_categories=batch_categories)
del conf_samples_processed

adata

# %%
sc.pp.filter_genes(adata, min_counts=1)
sc.pp.filter_cells(adata, min_counts=1)

adata

# %%
adata.obs.head()

# %% [markdown]
# ## Highest expression

# %%
sc.pl.highest_expr_genes(adata);

# %% [markdown]
# ## Normalization and log transform

# %%
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=par_preprocessing_target_sum)
sc.pp.log1p(adata)
adata.raw = adata

# %%
adata

# %% [markdown]
# # Downstream

# %%
# %%time

sc.pp.highly_variable_genes(adata, n_top_genes=par_downstream_n_top_genes)
sc.pp.pca(adata, n_comps=par_downstream_n_pcs, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=par_downstream_n_neighbors, metric=par_downstream_neighbor_metric)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=par_downstream_louvain_resolution)
sc.tl.diffmap(adata)

adata

# %% [markdown]
# ## Store the AnnData file

# %%
adata.write(par_save_filename)

# %%
adata.obs.head()

# %% [markdown]
# ## Percent mito UMI

# %%
f, ax = plt.subplots(1, 4, figsize=(20, 4))
sc.pl.scatter(adata,
              x='n_umis',
              y='n_genes',
              color='mt_frac',
              ax=ax[0],
              show=False,
              right_margin=2.85,
              title='Percent mitochondrial UMIs')
ax[0].set_xscale('log')
ax[0].set_yscale('log')

sc.pl.scatter(adata, x='n_umis', y='mt_frac', ax=ax[1], show=False)
ax[1].set_xscale('log')
plt.subplots_adjust(wspace=0.5)

sc.pl.violin(adata, keys='log10_n_umis', groupby='sample_name', rotation=90, ax=ax[2], show=False)
sc.pl.violin(adata, keys='log10_n_umis', groupby='sample_name', rotation=90, ax=ax[3], show=False)

# %% [markdown]
# # Visualization

# %% [markdown]
# ## PCA

# %%
sc.pl.pca_variance_ratio(adata)

# %%
sc.pl.pca_loadings(adata, components=range(1, 6))

# %%
sc.pl.pca(adata, color='leiden', components=['1,2', '3,2', '3,4', '5,4'], legend_loc=None)

# %%
feat = sorted(list(set(adata.obs.columns) & set(conf_visualize_features + conf_sample_features)))

sc.pl.pca(adata, color=feat, wspace=0.55, ncols=4)

# %% [markdown]
# ### Embedding & clustering

# %%
umap_point_size = np.maximum(120000/adata.n_obs, 2)

# %%
f, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color='leiden', legend_loc='on data', legend_fontoutline=3, legend_fontsize=14, legend_fontweight='normal', title='Clusters', ax=ax, show=False);

# %% [markdown]
# ### QC on embeddings

# %%
sc.pl.umap(adata, color=feat, wspace=0.7, ncols=4)

# %%
sc.pl.umap(adata, color='scrublet_score', vmax=0.5, cmap='Reds')

# %% [markdown]
# ### Samples

# %%
f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=True, sharex=True)
axs = axs.flatten()
for sample, ax in zip(adata.obs.sample_name.cat.categories, axs):
    sc.pl.umap(adata, color='sample_name', groups=sample, ax=ax, show=False, legend_loc=None, title=sample)

plt.subplots_adjust()

# %% [markdown]
# ## Remove doublets

# %%
if par_remove_doublets:

    adata = adata[~adata.obs.scrublet].copy()

    sc.pp.pca(adata, n_comps=par_downstream_n_pcs, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=par_downstream_n_neighbors, metric=par_downstream_neighbor_metric)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=par_downstream_louvain_resolution)
    sc.tl.diffmap(adata)

    f, ax = plt.subplots(figsize=(10, 10))
    sc.pl.umap(adata, color='leiden', legend_loc='on data', legend_fontoutline=3, legend_fontsize=14, legend_fontweight='normal', title='Clusters', ax=ax, show=False);

    ### QC on embeddings
    sc.pl.umap(adata, color=feat, wspace=0.7, ncols=4)
    sc.pl.umap(adata, color='scrublet_score', vmax=0.5, cmap='Reds')

    ### Samples
    f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                          figsize=(conf_plotting_width,
                                   conf_plotting_n_rows*conf_plotting_height_per_row),
                          sharey=True, sharex=True)
    axs = axs.flatten()
    for sample, ax in zip(adata.obs.sample_name.cat.categories, axs):
        sc.pl.umap(adata, color='sample_name', groups=sample, ax=ax, show=False, legend_loc=None, title=sample)
    plt.subplots_adjust()

    adata.write(par_save_filename)

# %% [markdown]
# ## Plots per group

# %%
# %%time

from IPython.core.display import display, HTML
sc.settings.verbosity = 0

if par_generate_plots_per_group:

    ascat = pd.Categorical(adata.obs[par_group_key])
    for group in ascat.categories:

        display(HTML(f'<h1>{par_group_key.capitalize()}: {group}</h1>'))

        ad = adata[ascat == group].copy()

        assert np.min(ad.X) >= 0.0
        ad._inplace_subset_var(ad.X.sum(0).A1 > 0)
        ad._inplace_subset_obs(ad.X.sum(1).A1 > 0)

        sc.pp.highly_variable_genes(ad, n_top_genes=par_downstream_n_top_genes)

        sc.pp.pca(ad, n_comps=par_downstream_n_pcs, svd_solver='arpack')
        sc.pp.neighbors(ad, n_neighbors=par_downstream_n_neighbors, metric=par_downstream_neighbor_metric)
        sc.tl.umap(ad)
        sc.tl.leiden(ad, resolution=par_downstream_louvain_resolution)
        sc.tl.diffmap(ad)

        ##### Visualize
        f, ax = plt.subplots(1, 4, figsize=(20, 4))
        sc.pl.scatter(ad,
                      x='n_umis',
                      y='n_genes',
                      color='mt_frac',
                      ax=ax[0],
                      show=False,
                      right_margin=2.85,
                      title='Percent mitochondrial UMIs')
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')

        sc.pl.scatter(ad, x='n_umis', y='mt_frac', ax=ax[1], show=False)
        ax[1].set_xscale('log')
        plt.subplots_adjust(wspace=0.5)

        sc.pl.violin(ad, keys='log10_n_umis', groupby='sample_name', rotation=90, ax=ax[2], show=False)
        sc.pl.violin(ad, keys='log10_n_umis', groupby='sample_name', rotation=90, ax=ax[3], show=False)

        f, ax = plt.subplots(figsize=(10, 10))
        sc.pl.umap(ad, color='leiden', legend_loc='on data', legend_fontoutline=3, legend_fontsize=14, legend_fontweight='normal', title='Clusters', ax=ax, show=False)

        sc.pl.umap(ad, color=feat, wspace=0.7, ncols=4)
        sc.pl.umap(ad, color='scrublet_score', vmax=0.5, cmap='Reds')

        ##### Save
        ad.write(par_save_filename_group % group)

# %% [markdown]
# ## Serialize the session

# %%
k = None
var = None

with shelve.open('session.pkl', protocol=4) as db:
    for k, var in globals().items():
        if k.startswith('par_') or k.startswith('conf_'):
            print(f'Storing {k}...')
            db[k] = var
db.close()
