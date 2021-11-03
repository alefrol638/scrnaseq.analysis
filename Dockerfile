FROM bioconductor/bioconductor_docker:RELEASE_3_13
LABEL MAINTAINER Aleksej Frolov
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
RUN LANG=C.UTF-8
RUN LC_ALL=C.UTF-8
RUN export DEBIAN_FRONTEND="noninteractive"
RUN export DEBCONF_NONINTERACTIVE_SEEN=true
RUN apt-get update -q && \
apt-get install -q -y --no-install-recommends apt-utils \
libudunits2-dev \
libmysqlclient-dev \
libhdf5-dev \
libgdal-dev \
gdal-bin \
proj-bin \
libproj-dev \
libgsl-dev \
libigraph0-dev \
zlib1g-dev \
libtool \
bison \
flex \
automake \
autoconf \
libpng*-dev \
libglpk-dev \
xorg \
libx11-dev \
libglu1-mesa-dev \
libfreetype*-dev \
p7zip \
build-essential \
libssl-dev \
libffi-dev \
libxslt1-dev \
python3.8 \
python-dev \
python3-dev \
python3.8-dev \
python3-pip \
python3-venv \
libxt-dev \
libgtk2.0-dev \
libcairo2-dev \
xvfb \
xauth \
xfonts-base \
python-opengl \
python-pyrex \
libgle3 \
python3-venv \
libxslt1-dev \
libldap2-dev \
libsasl2-dev \
libpq-dev \
pigz \
lbzip2 \
texlive \
texlive-latex-extra \
texlive-fonts-extra \
texlive-bibtex-extra \
texlive-science \
texi2html \
texinfo \
bzip2 \
ca-certificates \
git \
libglib2.0-0 \
libsm6 \
libxext6 \
libxrender1 \
mercurial \
openssh-client \
procps \
subversion \
wget \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*
RUN PATH=/opt/conda/bin:$PATH
ARG CONDA_VERSION=py39_4.10.3
RUN set -x && \
UNAME_M="$(uname -m)" && \
if [ "${UNAME_M}" = "x86_64" ]; then \
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"; \
SHA256SUM="1ea2f885b4dbc3098662845560bc64271eb17085387a70c2ba3f29fff6f8d52f"; \
elif [ "${UNAME_M}" = "s390x" ]; then \
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-s390x.sh"; \
SHA256SUM="1faed9abecf4a4ddd4e0d8891fc2cdaa3394c51e877af14ad6b9d4aadb4e90d8"; \
elif [ "${UNAME_M}" = "aarch64" ]; then \
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-aarch64.sh"; \
SHA256SUM="4879820a10718743f945d88ef142c3a4b30dfc8e448d1ca08e019586374b773f"; \
elif [ "${UNAME_M}" = "ppc64le" ]; then \
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-ppc64le.sh"; \
SHA256SUM="fa92ee4773611f58ed9333f977d32bbb64769292f605d518732183be1f3321fa"; \
fi && \
wget "${MINICONDA_URL}" -O miniconda.sh -q && \
echo "${SHA256SUM} miniconda.sh" > shasum && \
if [ "${CONDA_VERSION}" != "latest" ]; then sha256sum --check --status shasum; fi && \
mkdir -p /opt && \
sh miniconda.sh -b -p /opt/conda && \
rm miniconda.sh shasum && \
ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
echo "conda activate base" >> ~/.bashrc && \
find /opt/conda/ -follow -type f -name '*.a' -delete && \
find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
/opt/conda/bin/conda clean -afy
RUN chown 5532:1001 /opt/conda
RUN R -e 'devtools::install_version("future", version = "1.12.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE)'
RUN R -e 'BiocManager::install("BiocStyle",Ncpus=future::availableCores())'
RUN R -e 'install.packages("devtools",Ncpus=future::availableCores())'
RUN R -e 'install.packages("pkgbuild",Ncpus=future::availableCores())'
RUN R -e 'install.packages("cowplot",Ncpus=future::availableCores())'
RUN R -e 'install.packages("readxl",Ncpus=future::availableCores())'
RUN R -e 'install.packages("doMC",Ncpus=future::availableCores())'
RUN R -e 'install.packages("git2r",Ncpus=future::availableCores())'
RUN R -e 'install.packages("renv",Ncpus=future::availableCores())'
RUN R -e 'install.packages("qpcR",Ncpus=future::availableCores())'
RUN R -e 'install.packages("DESeq2",Ncpus=future::availableCores())'
RUN R -e 'install.packages("here",Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("tidyverse", version = "1.3.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("future.apply", version = "1.2.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("reticulate", version = "1.18", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("patchwork", version = "1.1.1", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("RColorBrewer", version = "1.1-2", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("pheatmap", version = "1.0.12", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("ggrepel", version = "0.9.1", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("ggthemes", version = "4.2.4", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("gplots", version = "3.1.1", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("ggpubr", version = "0.4.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("VennDiagram", version = "1.6.20", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("R.utils", version = "2.10.1", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("Seurat", version = "4.0.1", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("SoupX", version = "1.5.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("stringdist", version = "0.9.6.3", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("clustree", version = "0.4.3", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("NbClust", version = "3.0", repos = "http://cran.us.r-project.org",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("remotes", version = "2.3.0", repos = "http://cran.us.r-project.org",Ncpus=future::availableCores())'
RUN R -e 'install.packages("NMF")'
RUN R -e 'devtools::install_github("jokergoo/circlize")'
RUN R -e 'install.packages(c("clue","expm"),Ncpus = future::availableCores())'
RUN R -e 'install.packages("rjson",Ncpus = future::availableCores())'
RUN R -e 'devtools::install_github("jokergoo/GetoptLong")'
RUN R -e 'devtools::install_github("jokergoo/ComplexHeatmap")'
RUN R -e 'install.packages("rngtools",Ncpus = future::availableCores())'
RUN R -e 'devtools::install_version("SingleCellExperiment", version="1.14.1", repos=c("https://bioconductor.org/packages/3.13/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_github("chris-mcginnis-ucsf/DoubletFinder",auth_token = "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")'
RUN R -e 'devtools::install_github("satijalab/seurat-wrappers",auth_token = "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")'
RUN R -e 'devtools::install_github("sqjin/CellChat",auth_token = "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")'
RUN R -e 'devtools::install_github("immunogenomics/harmony",auth_token = "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")'
RUN R -e 'devtools::install_github("jlmelville/uwot",auth_token = "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")'
RUN R -e 'devtools::install_github("mojaveazure/seurat-disk",auth_token = "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")'
RUN R -e 'devtools::install_github("barkasn/fastSave",auth_token = "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")'
RUN R -e 'devtools::install_version("clusterProfiler", version="4.0.5", repos=c("https://bioconductor.org/packages/3.13/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("DOSE", version="3.16.0", repos=c("https://bioconductor.org/packages/3.12/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("tximport", version="1.18.0", repos=c("https://bioconductor.org/packages/3.12/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("biomaRt", version="2.46.3", repos=c("https://bioconductor.org/packages/3.12/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("LoomExperiment", version="1.8.0", repos=c("https://bioconductor.org/packages/3.12/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("EnhancedVolcano", version="1.10.0", repos=c("https://bioconductor.org/packages/3.13/bioc","http://cran.us.r-project.org"),upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("org.Hs.eg.db", version="3.12.0", repos="https://bioconductor.org/packages/3.12/data/annotation",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
RUN R -e 'devtools::install_version("org.Mm.eg.db", version="3.13.0", repos="https://bioconductor.org/packages/3.13/data/annotation",upgrade = "never", force = TRUE,Ncpus=future::availableCores())'
ADD scrnaseq.analysis /tmp/scrnaseq.analysis
RUN R -e 'devtools::install("/tmp/scrnaseq.analysis/")'
RUN R -e 'reticulate::conda_create("scvelo")'
RUN R -e 'reticulate::conda_install("scvelo","scvelo",pip=T,python_version = 3.7)'
RUN R -e 'reticulate::conda_create("totalVI")'
RUN R -e 'reticulate::conda_install("totalVI","scvi-tools",channel = c("bioconda","conda-forge"),python_version = 3.7)'
RUN R -e 'reticulate::conda_create("scirpy")'
RUN R -e 'reticulate::conda_install("scirpy","scvi-tools",channel = c("bioconda","conda-forge"),python_version = 3.7)'

CMD exec /bin/bash "$@"

