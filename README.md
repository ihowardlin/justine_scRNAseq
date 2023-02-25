Test

```
# Load libraries
library(DropletUtils)
library(readr)
library(scran)
library(scater)
library(celldex)
library(ensembldb)
library(pheatmap)
library(forcats)

# Plotting functions
library(ggplot2) 

# The main class we will use for Single Cell data
library(SingleCellExperiment) 

# Setting the seed for reproducibility
set.seed(12345)
```

```
#read in data and create a SCE object for downstream analysis
p0_sce <- read10xCounts("C:/Users/ihowa/Dropbox/justine/P0", type="sparse")
p5_sce <- read10xCounts("C:/Users/ihowa/Dropbox/justine/P5", type="sparse")
```
