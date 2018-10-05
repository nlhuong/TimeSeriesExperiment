[![Build Status](https://travis-ci.com/nlhuong/vistimeseq.svg?branch=master)](https://travis-ci.com/nlhuong/vistimeseq)

# vistimeseq
Analysis and visualization of short time course data

## Overview

vistimeseq is a package for visualization and analysis of short time-course datasets.
The package is a comprehensive toolbox built on S4 `vistimeseq` data structure designed
specifically to handle time-course data. Functions included allow user to perform the
analysis efficiently. Apart from native functions, `vistimeseq` also includes wrappers 
around functions from other package to let the user easily work on `vistimeseq`
object without having to switch between different frameworks. The package
is mostly indended for gene expression data, but can be applied to any
datasets with observations taken over time.

`vistimeseq` performs:

- data normalization and transformation
- heatmap plotting with colorbars
- PCA projection of samples 
- PCA projection of features with corresponding time-series overlayed
- Time-series features clustering
- Finding genes differentially expressed between two conditions/groups at specific timepoints
- Finding genes with differential trajectories between two conditions/groups
- Pathway enrichment analysis. Looking for over-representation of differential 
genes in GO or KEGG pathways.

To install the package do the following:

```{r}
install.packages("devtools")
devtools::install_github("nlhuong/vistimeseq")
```
