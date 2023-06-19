# scVar
A wrapper method to measure cell-to-cell variability in scRNA-seq data.

Each method/metric follows the vignette or paper description. Based on our evaluation [paper](https://www.biorxiv.org/content/10.1101/2022.11.24.517880v1.full.pdf+html), we suggest data has at least 50 cells per cell type to obtain reliable estimation.

### Installation 
```
devtools::install_github("huiwenzh/scVar")
library(scVar)
```

### Usage 
```
dat <- matrix(runif(n=1000), nrow=5000)
# use CV
dat_CVvar <- scVar(dat,norm=T,metric='CV')
# use scran
dat_scranvar <- scVar(dat,norm=F,metric='scran')
```


