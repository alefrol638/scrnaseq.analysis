# scRNAseq Analysis

Wrapper for multiple R packages for performing scRNAseq analysis on the data format of the PRECISE output. Includes generation of Volcano
Plots, GSEA, cluster abundances and other useful plots.

## Installation

You can install this package using the remotes package and a deploy token registered for this repository:

```{r}

###requires remotes version 2.3.0
remotes::install_git("https://gitlab.dzne.de/frolova/scrnaseq.analysis.git",
                     credentials=git2r::cred_user_pass("gitlab+deploy-token-8", "Li18nfcNd5bTZBgJkNzm"))

```