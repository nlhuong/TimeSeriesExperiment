[![Build Status](https://travis-ci.com/nlhuong/TimeSeriesExperiment.svg?branch=master)](https://travis-ci.com/nlhuong/TimeSeriesExperiment)

# TimeSeriesExperiment
Analysis and visualization of short time series data

## Overview

`TimeSeriesExperiment` is a package for visualization and analysis of short, 
regular time-series datasets. The package is a comprehensive toolbox built on 
`TimeSeriesExperiment` class (extension of `SummarizedExperiment`) designed
specifically to handle time-course data. Functions included allow user to 
perform the analysis efficiently. Apart from native functions, 
`TimeSeriesExperiment` also includes wrappers around functions from other 
package to let the user easily work on `TimeSeriesExperiment` object without 
having to switch between different frameworks. The package is mostly indended 
for gene expression data, but can be applied to any datasets with observations 
taken over time.

`TimeSeriesExperiment` performs:

- data normalization and transformation
- heatmap plotting with colorbars
- PCA projection of samples
- PCA projection of features with corresponding time-series overlayed
- Time-series features clustering
- Finding genes differentially expressed between two conditions/groups at 
specific timepoints
- Finding genes with differential trajectories between two conditions/groups
- Pathway enrichment analysis. Looking for over-representation of differential
genes in GO or KEGG pathways.

To install the package do the following:

```{r}
install.packages("devtools")
devtools::install_github("nlhuong/TimeSeriesExperiment")
```

The package is currently submitted to Bioconductor and will be available
in the future from:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("TimeSeriesExperiment")
```