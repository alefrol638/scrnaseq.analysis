# Singularity 

All the required packages are installed in singularity container: alefrol/default/bioconductor_3.13:fullR.
To start a R studio server session execute this on your remote server (f.e 10.0.161.2): 


```{bash}
###create config folders and files, to be able to run locally (first fill in individual information)

mkdir -p run var-lib-rstudio-server

printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf

## run container 
start container
PASSWORD=test123 singularity exec 
--nv --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf 
--bind server_location:container_location library://alefrol/default/bioconductor_3.13fullR rserver --www-address=10.0.161.2 --auth-none=0 
--auth-pam-helper-path=pam-helper --secure-cookie-key-file ~/tmp/r-server --server-data-dir ~/var/run/rstudio-server
--www-port=<portofchoice>&

```

# scRNAseq Analysis

Wrapper for multiple R packages for performing scRNAseq analysis on the data format of the PRECISE output. Includes generation of Volcano
Plots, GSEA, cluster abundances and other useful plots.

## Installation

You can install this package (including unresolved dependecies) using the remotes package and a deploy token registered for this repository:

```{r}

BiocManager::install("EnhancedVolcano",Ncpus=100)


###requires remotes version 2.3.0
remotes::install_git("https://gitlab.dzne.de/frolova/scrnaseq.analysis.git",
                     credentials=git2r::cred_user_pass("gitlab+deploy-token-8", "Li18nfcNd5bTZBgJkNzm",build_vignettes = T))

```

if that does not work you probably need to install remotes version 2.3.0:

```{r}

require(devtools)
install_version("remotes", version = "2.3.0", repos = "http://cran.us.r-project.org")

```

## Usage

you can browse the vignettes using:

```{r}
### check all installed vignettes
browseVignettes("scRNAseq.analysis")
###have a look at the tutorial html file 
vignette("Tutorial_scRNAseq_analysis")

###or get the source code for the tutorial vignette

edit(vignette("Tutorial_scRNAseq_analysis"))


```