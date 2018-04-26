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

```{r}
data(yeast_metabolome) 
# Check what is inside the example data:
yeast_metabolome$features[1:10,]  # Metabolic features (m/z values and ids)
yeast_metabolome$X[1:10,1:10] # samples x metabolic features data matrix
X = yeast_metabolome$X
```
### Also possible to load example data from .csv file:

```{r}
yeast_metabolome = read.csv("https://raw.githubusercontent.com/daniellyz/MetICA2/master/inst/Yeast-metabolome.csv")
features = yeast_metabolome[,c("ID","Mass")]
X = yeast_metabolome[,3:ncol(yeast_metabolome)] # Only keep intensity data for MetICA
rownames(X) = features$ID
X = t(X) # Transpose the data since MetICA accepts samples x variables matrices
```

### Run MetICA simulations:

```{r}
# Begin a MetICA simulation with 2000 estimated components in total. The samples are not time-dependent, so trend = FALSE. Numbers of clusters are evaluated between 2 and 10:
M1=MetICA(X,pcs = 10,max_iter = 200,boot.prop = 0.3,trends = F, max.cluster = 10)
# Users can confirm the number of pcs used for denoising if they think enough variance is explained, they can modify the number of pcs at this moment as well:
```
![choose](inst/Launch_MetICA.JPG)

```{r}
# Validation plot to decide the number of clusters
results=validationPlot(M1)


