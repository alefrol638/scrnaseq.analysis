# scRNAseq Analysis

Wrapper for multiple R packages for performing scRNAseq analysis on the data format of the PRECISE output. Includes generation of Volcano
Plots, GSEA, cluster abundances and other useful plots.

## Singularity 

All the required packages are installed in docker image: alefrol94/scrnaseq.analysis.
To start a R studio server session execute this on your remote server (f.e 10.0.161.2): 


```{bash}
###create config folders and files, to be able to run locally (first fill in individual information)

mkdir -p run var-lib-rstudio-server

printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf

## run container 

PASSWORD=test123 singularity exec 
--nv --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf 
--bind server_location:container_location docker:alefrol94/scrnaseq.analysis rserver --www-address=10.0.161.2 --auth-none=0 
--auth-pam-helper-path=pam-helper --secure-cookie-key-file ~/tmp/r-server --server-data-dir ~/var/run/rstudio-server
--www-port=<portofchoice>&

```


## Installation

You can install this package outside of the docker container (including unresolved dependecies) using the remotes package and a deploy token registered for this repository. However, 
you might need to install more packages which are not listed here:

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

Please download the vignette from this repository under ./vignettes/Tutorial_scRNAseq_analysis.html and open it 
with your internet browser of choice.