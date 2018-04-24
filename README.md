# MetICA: Independent component analysis for high-resolution mass-spectrometry based metabolomics

### Context
ICA is an important alternative to classical statistical approaches for non-targeted metabolomics data. It extends the concept of regular correlation (e.g. in PCA, ASCA and PLS-DA) to statistical dependance by capturing higher order dependencies. However, its algorithm instability (output variations between different algorithm runs) and the biological validity of components have been overlooked when applied to complex metabolomics data. MetICA adresses these problems by gathering ICs estimated from multiple algorithm runs and from bootstrapped datasets, clustering them so as to find the most representative components. While evaluating the algorithmic stability, MetICA also suggests multiple criteria to select the correct number of components and to rank the extracted components.

## Installation from Github using R (with devtools)

```R
library("devtools")
install_github("MetICA2","daniellyz",dependencies=TRUE)
library("MetICA2")
```


