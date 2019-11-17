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
# # jovian - 01 Upstream QC

# %%
import scanpy as sc

# Plotting
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# numpy et al.
import numpy as np
import scipy.sparse as sp
import scipy
import pandas as pd

# R integration
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.vectors import StrVector, FloatVector, ListVector
from rpy2.robjects import r
import rpy2.robjects as ro
import anndata2ri # scipy.sparse + AnnData support

numpy2ri.activate()
pandas2ri.activate()
anndata2ri.activate()

from pathlib import Path
import math
from tqdm.auto import tqdm
import warnings
import shelve

# %%
# %load_ext rpy2.ipython

# %%
sc.settings.verbosity = 'hint'
plt.rcParams['figure.dpi'] = 150


# %% [markdown]
# ## Parameters

# %% {"tags": ["parameters"]}
par_species = 'human' # mouse or human
par_data_dir = 'data'

# Unfiltered QC
par_initial_umi_cutoff = 10
par_initial_gene_cutoff = 10

# EmptyDrops
par_empty_drops_lower_umi_cutoff = 200
par_empty_drops_ignore_cutoff = 10
par_empty_drops_niters = 10000
par_empty_drops_fdr_cutoff = 0.01
par_empty_drops_retain = 800

# %%
conf_mt_prefix = 'MT-' if par_species == 'human' else 'mt-'

# %% [markdown]
# ***

# %% [markdown]
# ## Read the sample sheet and data sets

# %%
conf_sample_sheet = pd.read_csv(Path(par_data_dir) / 'samples.csv').set_index('h5ad_or_h5_path')
conf_sample_sheet

# %%
conf_n_samples = len(conf_sample_sheet)

conf_plotting_n_cols = 5
conf_plotting_width = 20
conf_plotting_height_per_row = 3.5

conf_plotting_n_rows = math.ceil(conf_n_samples/conf_plotting_n_cols)
conf_sample_features = conf_sample_sheet.columns.tolist()

# %%
assert 'sample_name' in conf_sample_sheet.columns, 'Sample sheet must have a unique sample_name column'
assert 'raw' in conf_sample_sheet.columns, 'Sample sheet must have a raw (True or False) column'
assert conf_sample_sheet.sample_name.nunique() == conf_n_samples, 'Sample sheet must have a unique sample_name column'


# %% [markdown]
# ## Tools for prefiltered data

# %%
def emptydrops(adata, lower, niters, ignore, retain=None, **kwargs):
    du = importr('DropletUtils')
    ed = du.emptyDrops(adata.X.T,
                       lower=lower,
                       niters=niters,
                       ignore=ignore,
                       retain=retain,
                       **kwargs)
    fdr = ed.slots['listData'].rx2('FDR')
    adata.obs['empty_drops_FDR'] = fdr

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fdr[fdr == 0.0] = 1.0/niters

        adata.obs['neg_log10_empty_drops_FDR'] = -np.log10(fdr)


# %%
# TODO: loading an R package that uses reticulate causes a segfault, see https://github.com/rstudio/reticulate/issues/208
# Seurat 3.1.0 loads leiden package that causes this. Use Seurat 3.0.2 for now.
def soupX(path_list, sample_names, soup_range=(0, 10), keep_droplets=True, n_top_cont_genes=20):

    sx = importr('SoupX')
    scl = sx.load10X(StrVector(path_list),
                     channelNames=StrVector(sample_names),
                     soupRange=FloatVector(soup_range),
                     keepDroplets=keep_droplets)

    scl = sx.inferNonExpressedGenes(scl)
    for sample in sample_names:
        df = scl.rx2('channels').rx2(sample).rx2('nonExpressedGenes')
        cont_genes = df[df.isUseful == 1][:n_top_cont_genes].index.values
        sx_cont_genes = ListVector(dict(c_genes=cont_genes.tolist()))
        scl = sx.calculateContaminationFraction(scl, sample, cont_genes)
        scl = sx.interpolateCellContamination(scl, sample, useGlobal=True)

    return scl


# %% [markdown]
# ## Store QC info and run tools on prefiltered data

# %%
# %%time

conf_samples = {}

for sample in tqdm(list(conf_sample_sheet.itertuples())):
    file = Path(par_data_dir) / sample.Index

    if str(sample.raw).lower() == 'true':
        ad = sc.read_10x_h5(file).copy()
    else:
        ad = sc.read(file)

    ad.var_names_make_unique()
    display(ad)

    # Read 10X filtered file
    file_path = Path(file)
    filtered_h5_path = file_path.parent / file_path.name.replace('raw', 'filtered')
    if filtered_h5_path.exists():
        cell_barcodes = sc.read_10x_h5(filtered_h5_path).obs_names.values
        ad.obs['10x_cell_calling'] = False
        ad.obs.loc[cell_barcodes] = True
    else:
        ad.obs['10x_cell_calling'] = np.nan

    for sample_feature in conf_sample_features:
        ad.obs[sample_feature] = sample._asdict()[sample_feature]

    if str(sample.raw).lower() == 'true':
        # EmptyDrops
        emptydrops(ad,
                   lower=par_empty_drops_lower_umi_cutoff,
                   niters=par_empty_drops_niters,
                   ignore=par_empty_drops_ignore_cutoff,
                   retain=par_empty_drops_retain)
    else:
        ad.obs['empty_drops_FDR'] = np.nan
        ad.obs['neg_log10_empty_drops_FDR'] = np.nan

    ad.obs['n_umis']  = ad.X.sum(1)
    ad.obs['n_genes'] = (ad.X != 0).sum(1).A1
    ad.obs['log10_n_umis'] = np.log10(ad.X.sum(1))
    ad.obs['log10_n_genes'] = np.log10((ad.X != 0).sum(1).A1)

    # Total UMI/gene cutoffs
    ad = ad[ad.obs.n_umis  > par_initial_umi_cutoff]
    ad = ad[ad.obs.n_genes > par_initial_gene_cutoff]
    display(ad)

    # Save barcode rank
    ad.obs['barcode_rank'] = scipy.stats.rankdata(-ad.obs['n_umis'])

    # MT
    mt_gene_mask = ad.var_names.str.startswith(conf_mt_prefix)
    assert mt_gene_mask.sum() > 0, 'Wrong mt prefix'
    ad.obs['mt_frac'] = ad.X[:, mt_gene_mask].sum(1).A1 / ad.obs['n_umis']

    conf_samples[sample.sample_name] = ad.copy()

# %%
conf_samples

# %% [markdown]
# ## UMI/Gene Rank Plots

# %%
f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=True, sharex=True)
axs = axs.flatten()

for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
    ax.scatter(ad.obs.barcode_rank.values,
               ad.obs.n_umis.values,
               s=1, alpha=0.5)

    ax.axhline(par_empty_drops_ignore_cutoff, color='red')
    ax.axhline(par_empty_drops_lower_umi_cutoff, color='red')

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlabel('Barcode rank')
    ax.set_ylabel('nUMIs')

    ax.set_title(sample)

plt.subplots_adjust()

# %% [markdown]
# ## nGenes vs nUMIs and distributions

# %%
f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=True, sharex=True)
axs = axs.flatten()

for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
    ax.scatter(ad.obs.n_umis.values+1,
               ad.obs.n_genes.values+1,
               alpha=0.5,
               s=0.1)
    ax.axvline(par_empty_drops_ignore_cutoff, color='red')
    ax.axvline(par_empty_drops_lower_umi_cutoff, color='red')

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title(sample)

    ax.set_xlabel('nUMIs')
    ax.set_ylabel('nGenes')


f.suptitle('UMIs vs genes', fontsize=16)
plt.subplots_adjust(top=0.9)

## UMI distr.

f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=True, sharex=True)
axs = axs.flatten()

for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
    ax.hist(ad.obs.log10_n_umis.values, bins=100)
    ax.axvline(np.log10(par_empty_drops_ignore_cutoff), color='red')
    ax.axvline(np.log10(par_empty_drops_lower_umi_cutoff), color='red')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_title(sample)
    ax.set_xlabel('nUMIs (log10)')

f.suptitle('UMI distributions', fontsize=16)
plt.subplots_adjust(top=0.9)

## Gene distr.

f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                      figsize=(conf_plotting_width,
                               conf_plotting_n_rows*conf_plotting_height_per_row),
                      sharey=True, sharex=True)
axs = axs.flatten()

for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
    ax.hist(ad.obs.log10_n_genes.values, bins=100)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.set_title(sample)
    ax.set_xlabel('nGenes (log10)')

f.suptitle('Gene distributions', fontsize=16)
plt.subplots_adjust(top=0.9)

# %% [markdown]
# ## nGenes vs nUMIs and distributions (with EmptyDrops)

# %%
from mpl_toolkits.axes_grid1 import make_axes_locatable

if (conf_sample_sheet.raw.astype(str).str.lower() == 'true').any():

    f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                          figsize=(conf_plotting_width,
                                   conf_plotting_n_rows*conf_plotting_height_per_row),
                          sharey=True, sharex=True)
    axs = axs.flatten()

    for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
        if np.any(np.isfinite(ad.obs.neg_log10_empty_drops_FDR)):
            idx = np.argsort(ad.obs.neg_log10_empty_drops_FDR)[::-1]
            ad = ad[idx]

            pc = ax.scatter(ad.obs.n_umis.values+1,
                            ad.obs.n_genes.values+1,
                            s=0.1,
                            c=ad.obs.neg_log10_empty_drops_FDR)
            ax.axvline(par_empty_drops_ignore_cutoff, color='red')
            ax.axvline(par_empty_drops_lower_umi_cutoff, color='red')

            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_title(sample)

            ax.set_xlabel('nUMIs')
            ax.set_ylabel('nGenes')

            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            f.colorbar(pc, cax=cax, orientation='vertical')

    f.suptitle('UMIs vs genes (colored by EmptyDrops FDR)', fontsize=16)
    plt.subplots_adjust(wspace=0.3, hspace=0.3, top=0.9)

    ## UMI distr.

    f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                          figsize=(conf_plotting_width,
                                   conf_plotting_n_rows*conf_plotting_height_per_row),
                          sharey=True, sharex=True)
    axs = axs.flatten()

    for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
        if np.any(np.isfinite(ad.obs.neg_log10_empty_drops_FDR)):
            ad = ad[ad.obs.neg_log10_empty_drops_FDR > -np.log10(par_empty_drops_fdr_cutoff)]
            ax.hist(ad.obs.log10_n_umis.values, bins=100)
            ax.axvline(np.log10(par_empty_drops_ignore_cutoff), color='red')
            ax.axvline(np.log10(par_empty_drops_lower_umi_cutoff), color='red')

            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.set_title(sample)
            ax.set_xlabel('nUMIs (log10)')

    f.suptitle('EmptyDrops-selected barcodes (UMI distributions)', fontsize=16)
    plt.subplots_adjust(top=0.9)

    # Gene distr.

    f, axs = plt.subplots(conf_plotting_n_rows, conf_plotting_n_cols,
                          figsize=(conf_plotting_width,
                                   conf_plotting_n_rows*conf_plotting_height_per_row),
                          sharey=True, sharex=True)
    axs = axs.flatten()

    for sample, ad, ax in tqdm(list(zip(conf_samples.keys(), conf_samples.values(), axs))):
        if np.any(np.isfinite(ad.obs.neg_log10_empty_drops_FDR)):
            ad = ad[ad.obs.neg_log10_empty_drops_FDR > -np.log10(par_empty_drops_fdr_cutoff)]
            x = ad.obs.log10_n_genes.values
            ax.hist(x, bins=100)

            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.set_title(sample)
            ax.set_xlabel('nGenes (log10)')

    f.suptitle('EmptyDrops-selected barcodes (gene distributions)', fontsize=16)
    plt.subplots_adjust(top=0.9)

# %% [markdown]
# ## Serialize the session

# %%
k = None
var = None

with shelve.open('session.pkl', protocol=4, flag='n') as db:
    for k, var in globals().items():
        if k.startswith('par_') or k.startswith('conf_'):
            print(f'Storing {k}...')
            db[k] = var
db.close()
