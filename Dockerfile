ARG BASE_CONTAINER=jupyter/scipy-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Gokcen Eraslan <geraslan@broadinstitute.org>"

USER root

# Ubuntu packages needed for R and Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends ffmpeg cmake less libigraph-dev && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get install -y wget build-essential curl libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev git unzip && \
    apt-get install --no-install-recommends -y lsb-release gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=290.0.0-0 && \
    apt-get -qq -y autoremove && \
    apt-get clean

USER $NB_UID

# to avoid https://github.com/conda/conda/issues/9681
RUN conda install --quiet --yes conda=4.8.3

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install R and python packages through conda-forge
RUN conda install --quiet --yes \
    plotly \
    plotnine \
    seaborn \
    'pytables=3.6*' \
    'numba=0.48*' \
    python-igraph \
    leidenalg \
    louvain \
    jupytext \
    graph-tool \
    hvplot

RUN conda install --quiet --yes \
    'r-base=3.6*' \
    'r-devtools=2.3*' \
    'r-plyr=1.8*' \
    'r-ggplot2=3.3*' \
    'r-rcurl=1.95*' \
    'r-biocmanager' \
    'r-seurat=3.1.5' \
    'r-rgl=0.100*'
    
RUN conda install --quiet --yes -c pytorch \
    'pytorch=1.5*' cpuonly \
    torchvision && \
    conda clean --all -f -y

ENV TAR=/bin/tar
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"

# Install R packages
RUN R -e 'BiocManager::install(c("SingleR", "DropletUtils", "scater", "scran", "scRNAseq", "MAST", "multtest"))' \
 && R -e 'devtools::install_github("constantAmateur/SoupX", upgrade=F)' \
 && R -e 'install.packages(c("lme4", "ggpubr", "IRkernel"), repos = "http://cran.us.r-project.org")' \
 && R -e 'IRkernel::installspec()' \
 && fix-permissions /home/$NB_USER
 
RUN R -e 'install.packages(c("Formula", "maxLik"), repos = "http://cran.us.r-project.org")' \
 && R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/DirichletReg/DirichletReg_0.6-3.1.tar.gz", repos = NULL, type = "source")'

# Install python3 packages
RUN pip install scanpy anndata -U && \
    pip install scvelo scrublet fa2 mnnpy MulticoreTSNE scplot \
                openpyxl scvi cellxgene skggm pyannotables  \
                papermill rpy2 harmony-pytorch && \
    pip install git+https://github.com/flying-sheep/anndata2ri.git && \
    pip install git+https://github.com/broadinstitute/CellBender.git && \
    pip install black -U && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

RUN ipython profile create && \
    echo "c.InlineBackend.figure_format = 'retina'" >> ~/.ipython/profile_default/ipython_kernel_config.py && \
    echo "c.InteractiveShell.cache_size = 0" >> ~/.ipython/profile_default/ipython_kernel_config.py

RUN jupyter labextension install @jupyterlab/toc \
    && fix-permissions /home/$NB_USER

USER root

# install arbitrary Ubuntu packages here to speed things up
RUN apt-get update && \
    apt-get install -y --no-install-recommends htop vim aria2 && \
    rm -rf /var/lib/apt/lists/*

# Configure container startup
ENV JUPYTER_ENABLE_LAB=1
ENV GRANT_SUDO=1
ENV OPENBLAS_NUM_THREADS=1
ENV JUPYTER_TOKEN=aacb5622c49616c257f8975e6edc6c2be307330b0ac6c6931a098c7ed3f7afa5
ENTRYPOINT ["start-notebook.sh"]
