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
# # jovian - 02 Downstream QC

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
# ## Load session

# %%
with shelve.open('session.pkl', protocol=4, writeback=False) as db:
    for k in db.keys():
        globals()[k] = db[k]

# %% [markdown]
# ## Parameters

# %% {"tags": ["parameters"]}
par_cutoff_min_counts = 200
par_cutoff_min_genes  = 200
par_cutoff_max_genes  = None
par_final_empty_drops_fdr_cutoff = 0.01

par_remove_mito_genes = True
par_mito_cutoff = 'mito_cutoff_'
par_remove_sex_genes = False

par_preprocessing_target_sum = 10000
par_regress_out_variables = []
par_regress_out_n_jobs = 6

# downstream parameters
par_downstream_n_top_genes = 3000
par_downstream_n_pcs = 50
par_downstream_n_neighbors = 15
par_downstream_louvain_resolution = 1.5
par_downstream_neighbor_metric = 'euclidean'

par_save_filename_sample = 'adata-sample-%s.h5ad'

# %% [markdown]
# ***

# %%
conf_samples

# %% [markdown]
# ## Filtering and cell calling

# %%
for k, ad in tqdm(list(conf_samples.items())):
    if par_cutoff_min_genes: sc.pp.filter_cells(ad, min_genes=par_cutoff_min_genes)
    if par_cutoff_max_genes: sc.pp.filter_cells(ad, max_genes=par_cutoff_max_genes)
    if par_cutoff_min_counts: sc.pp.filter_cells(ad, min_counts=par_cutoff_min_counts)
    if par_final_empty_drops_fdr_cutoff: ad._inplace_subset_obs(ad.obs.empty_drops_FDR < par_final_empty_drops_fdr_cutoff)
    if par_mito_cutoff:
        if isinstance(par_mito_cutoff, str) and par_mito_cutoff.endswith('_'):
            ad._inplace_subset_obs(ad.obs.mt_frac < ad.obs[par_mito_cutoff])
        else:
            ad._inplace_subset_obs(ad.obs.mt_frac < par_mito_cutoff)

    # remove all-zero genes
    if sp.issparse(ad.X):
        ad._inplace_subset_var(ad.X.sum(0).A1 > 0)
    else:
        ad._inplace_subset_var(ad.X.sum(0) > 0)

    display(ad)

# %% [markdown]
# ## Doublet detection

# %%
import scrublet as scr

def run_scrublet(ad):
    scrub = scr.Scrublet(ad.X)
    n_prin_comps = min(min(ad.X.shape[0], ad.X.shape[1]), 30)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=n_prin_comps)

    ad.obs['scrublet'] = predicted_doublets
    ad.obs['scrublet_score'] = doublet_scores


# %%
for adata in tqdm(list(conf_samples.values())):
    run_scrublet(adata)

# %% [markdown]
# ## Predict sex

# %%
conf_xist_gene_name = 'Xist' if par_species == 'mouse' else 'XIST'

for adata in tqdm(list(conf_samples.values())):
    if conf_xist_gene_name is in adata.var_names:
        adata.obs['predicted_sex'] = ['female' if x else 'male' for x in adata.obs_vector(conf_xist_gene_name) > 0]
    else:
        adata.obs['predicted_sex'] = 'male'

# %% [markdown]
# ## Highest expression

# %%
f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=False, sharex=False)
axs = axs.flatten()
for (sample, adata), ax in zip(list(conf_samples.items()), axs):
    sc.pl.highest_expr_genes(adata, ax=ax, show=False, n_top=15);
    ax.set_title(sample)

plt.subplots_adjust(wspace=0.6, hspace=0.4)

# %% [markdown]
# ### Exclude mito genes and & sex genes

# %%
if par_remove_mito_genes:
    for adata in tqdm(list(conf_samples.values())):
        adata._inplace_subset_var(~adata.var_names.str.startswith(conf_mt_prefix))
        display(adata)

# %%
if par_remove_sex_genes:
    from pyannotables import tables

    for sample in tqdm(list(conf_samples.keys())):
        genes = tables['mus_musculus-ensembl95-GRCm38'] if par_species == 'mouse' else tables['homo_sapiens-ensembl95-GRCh38']
        sex_genes = genes.gene_name[genes.contig.isin(['X', 'Y'])].values

        conf_samples[sample] = conf_samples[sample][:, ~conf_samples[sample].var_names.isin(sex_genes)].copy()
        display(conf_samples[sample])

    del tables

# %%
f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=False, sharex=False)
axs = axs.flatten()
for (sample, adata), ax in zip(list(conf_samples.items()), axs):
    sc.pl.highest_expr_genes(adata, ax=ax, show=False, n_top=15);
    ax.set_title(sample)

plt.subplots_adjust(wspace=0.6, hspace=0.4)

# %% [markdown]
# ## Normalization and log transform

# %%
for adata in tqdm(list(conf_samples.values())):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=par_preprocessing_target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    display(adata)

# %% [markdown]
# ## Cell cycle scores

# %%
gene_list_url = 'https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt'

if par_species == 'mouse':
    cell_cycle_genes = [str(x.strip(), 'utf-8').capitalize() for x in urlopen(gene_list_url)] # capitalize = shame
else:
    cell_cycle_genes = [str(x.strip(), 'utf-8') for x in urlopen(gene_list_url)]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

for adata in tqdm(list(conf_samples.values())):
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# %% [markdown]
# ## Regress out

# %%
# %%time

for adata in tqdm(list(conf_samples.values())):
    for v in par_regress_out_variables:
        sc.pp.regress_out(adata, v)

# %% [markdown]
# ## Downstream

# %%
sc.settings.verbosity = 0

for sample, ad in tqdm(list(conf_samples.items())):

    assert np.min(ad.X) >= 0.0
    ad._inplace_subset_var(ad.X.sum(0).A1 > 0)
    ad._inplace_subset_obs(ad.X.sum(1).A1 > 0)

    sc.pp.highly_variable_genes(ad, n_top_genes=par_downstream_n_top_genes)

    sc.pp.pca(ad, n_comps=par_downstream_n_pcs, svd_solver='arpack')
    sc.pp.neighbors(ad, n_neighbors=par_downstream_n_neighbors, metric=par_downstream_neighbor_metric)
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=par_downstream_louvain_resolution)
    sc.tl.diffmap(ad)

    ad._sanitize()

# %% [markdown]
# # Visualization

# %%
conf_visualize_features = ['log10_n_umis', 'log10_n_genes', 'neg_log10_empty_drops_FDR',
                           'mt_frac', '10x_cell_calling', 'predicted_sex',
                           'phase', 'scrublet']

# %%
# %%time

from IPython.core.display import display, HTML

for sample, ad in tqdm(list(conf_samples.items())):

    display(HTML(f'<h1>Sample: {sample}</h1>'))

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

    sc.pl.umap(ad, color=conf_visualize_features, wspace=0.7, ncols=4)
    sc.pl.umap(ad, color='scrublet_score', vmax=0.5, cmap='Reds')

    ##### Save
    ad.write(par_save_filename_sample % sample)

# %% [markdown]
# ## Serialize the session

# %%
conf_samples_processed = conf_samples
del conf_samples

# %%
k = None
var = None

with shelve.open('session.pkl', protocol=4) as db:
    for k, var in globals().items():
        if k.startswith('par_') or k.startswith('conf_'):
            print(f'Storing {k}...')
            db[k] = var
db.close()
