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
p56_sce <- read10xCounts("C:/Users/ihowa/Dropbox/justine/P56", type="sparse")
```

<details><summary> P0 QC </summary>
  
```
# create a table of statistics using emptyDropsCellRanger
droplet_df <- DropletUtils::emptyDropsCellRanger(counts(p0_sce))

# view rows where FDR is not `NA`
droplet_df[!is.na(droplet_df$FDR), ]

# filter droplets using `which` to prevent NA trouble
cells_to_retain <- which(droplet_df$FDR < 0.01)
filtered_sce <- p0_sce[, cells_to_retain]

filtered_sce

```
</details>
  
```
# read in a table of mitochondrial genes and extract ids
mito_genes <- readr::read_tsv("hs_mitochondrial_genes.tsv") |>
  # filter to only the genes that are found in our dataset
  dplyr::filter(gene_id %in% rownames(filtered_sce)) |>
  # create a vector from the gene_id column
  dplyr::pull(gene_id)
```
