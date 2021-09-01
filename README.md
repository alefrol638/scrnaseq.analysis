# scRNAseq Analysis

Wrapper for multiple R packages for performing scRNAseq analysis on the data format of the PRECISE output. Includes generation of Volcano
Plots, GSEA, cluster abundances and other useful plots.

## Installation

You can install this package (including unresolved dependecies) using the remotes package and a deploy token registered for this repository:

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

you can browse the vignettes using:

```{r}
### check all installed vignettes
browseVignettes("scRNAseq.analysis")
###have a look at the tutorial html file 
vignette("Tutorial_scRNAseq_analysis")

###or get the source code for the tutorial vignette

edit(vignette("Tutorial_scRNAseq_analysis"))


```