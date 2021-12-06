# scRNAseq Analysis

Wrapper for multiple R packages for performing scRNAseq analysis on the data format of the DropSeq output. Includes generation of Volcano
Plots, GOEA, cluster abundances and other useful plots.

## Singularity 

All the required packages are installed in docker image: alefrol94/scrnaseq.analysis.
To start a R studio server session execute this on your remote server:

```{bash}
###create config folders and files, to be able to run locally (first fill in individual information)

mkdir -p run var-lib-rstudio-server

printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf

## run container 

PASSWORD=<yourPassword> singularity exec 
--nv --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf 
--bind server_location:container_location docker:alefrol94/scrnaseq.analysis rserver --www-address=<yourServer> --auth-none=0 
--auth-pam-helper-path=pam-helper --secure-cookie-key-file ~/tmp/r-server --server-data-dir ~/var/run/rstudio-server
--www-port=<portofchoice>&

```
Alternatively in Docker: 
```{bash}

docker run -d -p <localport>:<remote> -e PASSWORD=<yourPassword> -v = /home/user:/home/user alefrol94/scrnaseq.analysis:reticulate rserver --www-port=8787 --secure-cookie-key-file /home/user/tmp/r-server --server-daemonize=0

```

This image also contains a full latex installation and miniconda, if you would like to document your work or use python packages via 
reticulate.


start using the environment, f.e scvelo:

```{r}
reticulate::use_condaenv("scvelo")
```

The following environments are preinstalled: 

- scvelo (Velocity analysis)
- scirpy (TCR clonotypes)
- totalVI (RNA + Protein seq)



## Installation

You can install this package outside of the docker container (including unresolved dependecies) using the remotes package and a deploy token registered for this repository. However, 
you might need to install more packages which are not listed here:

```{r}

BiocManager::install("EnhancedVolcano",Ncpus=100)


###requires remotes version 2.3.0
remotes::install_git("https://gitlab.dzne.de/frolova/scrnaseq.analysis.git",
                     credentials=git2r::cred_user_pass("gitlab+deploy-token-8", "Li18nfcNd5bTZBgJkNzm"))

```

if that does not work you probably need to install remotes version 2.3.0:

```{r}

require(devtools)
install_version("remotes", version = "2.3.0", repos = "http://cran.us.r-project.org")

```

## Usage

Please download the vignette from this repository under ./vignettes/Tutorial_scRNAseq_analysis.html and open it 
with your internet browser of choice.
