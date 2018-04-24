# MetICA: Independent component analysis for high-resolution mass-spectrometry based metabolomics

## Context
ICA is an important alternative to classical statistical approaches for non-targeted metabolomics data. It extends the concept of regular correlation (e.g. in PCA, ASCA and PLS-DA) to statistical dependance by capturing higher order dependencies. However, its algorithm instability (output variations between different algorithm runs) and the biological validity of components have been overlooked when applied to complex metabolomics data. MetICA adresses these problems by gathering ICs estimated from multiple algorithm runs and from bootstrapped datasets, clustering them so as to find the most representative components. While evaluating the algorithmic stability, MetICA also suggests multiple criteria to select the correct number of components and to rank the extracted components.

## Installation from Github using R (with devtools)

```R
library(devtools)
install_github("daniellyz/MetICA2")
library(MetICA)
```

## An example of data analysis using MetICA

### Load yeast metabolomics data:

```R
data(yeast_metabolome) 
# Check what is inside the example data:
yeast_metabolome$feature[1:10,]  # Metabolic features (m/z values and ids)
yeast_metabolome$X[1:10,1:10] # samples x metabolic features data matrix
```
### Alternative way to load example data from .csv file:

```R

yeast_metabolom = read.csv()

```


### Load yeast metabolomics data:


