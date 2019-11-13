ARG BASE_CONTAINER=jupyter/scipy-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Gokcen Eraslan <geraslan@broadinstitute.org>"

USER root

# Ubuntu packages needed for R and Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends ffmpeg cmake libigraph-dev && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID

# Install R and python packages through conda-forge
RUN conda install --quiet --yes -c bioconda \
    plotly \
    plotnine \
    seaborn \
    'pytables=3.6*' \
    'numba=0.46*' \
    'python-igraph=0.7*' \
    leidenalg \
    louvain \
    'r-base=3.6.1' \
    'r-devtools=2.2*' \
    'r-irkernel=1.0*' \
    'r-plyr=1.8*' \
    'r-ggplot2=3.2*' \
    'r-rcurl=1.95*' \
    'r-biocmanager' \
    'r-seurat=3.0.2' \
    'r-huge=1.3*' \
    'r-psych=1.8*' \
    'rpy2=3.1*' && \
    conda clean --all -f -y

# Install R packages
RUN R -e 'devtools::install_github("constantAmateur/SoupX", upgrade=F)' \
 && R -e 'BiocManager::install(c("SingleR", "DropletUtils", "scater", "scran", "scRNAseq"))' \
 && R -e 'install.packages(c("bootnet"), repos = "http://cran.us.r-project.org")' \
 && R -e 'IRkernel::installspec()' \
 && fix-permissions /home/$NB_USER

# Install python3 packages
RUN pip install git+https://github.com/theislab/scanpy.git && \
    pip install scvelo scrublet fa2 mnnpy MulticoreTSNE scplot \
                openpyxl scvi cellxgene skggm pyannotables \
                papermill && \
    pip install git+https://github.com/flying-sheep/anndata2ri.git && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

RUN ipython profile create && \
    echo "c.InlineBackend.figure_format = 'retina'" >> ~/.ipython/profile_default/ipython_kernel_config.py && \
    echo "c.InteractiveShell.cache_size = 0" >> ~/.ipython/profile_default/ipython_kernel_config.py

USER root

# install arbitrary Ubuntu packages here to speed things up
RUN apt-get update && \
    apt-get install -y --no-install-recommends htop vim && \
    rm -rf /var/lib/apt/lists/*

# Configure container startup
ENV JUPYTER_ENABLE_LAB=1
ENV GRANT_SUDO=1
ENV OPENBLAS_NUM_THREADS=1
ENTRYPOINT ["start-notebook.sh"]
