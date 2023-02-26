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
library(miQC)

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

# read in a table of mitochondrial genes and extract ids
mito_genes <- readr::read_tsv("hs_mitochondrial_genes.tsv") |>
     # filter to only the genes that are found in our dataset
     dplyr::filter(symbol %in% rownames(filtered_sce)) |>
     # create a vector from the gene_id column
     dplyr::pull(symbol)

filtered_sce <- scuttle::addPerCellQC(filtered_sce,
                                      subsets = list(mito = mito_genes))
head(colData(filtered_sce))

# use miQC::plotMetrics()
miQC::plotMetrics(filtered_sce) + theme_bw()

# fit the miQC model
miqc_model <- miQC::mixtureModel(filtered_sce) 

# plot the miQC model
miQC::plotModel(filtered_sce, miqc_model) + 
  theme_bw()
  
# look at miQC filtering
miQC::plotFiltering(filtered_sce, miqc_model, 
                    posterior_cutoff = 0.75) + 
  theme_bw()

# perform miQC filtering
qcfiltered_sce <- miQC::filterCells(filtered_sce, 
                                    model = miqc_model)
                  
# filter cells by unique gene count (`detected`)
qcfiltered_sce <- qcfiltered_sce[, which(qcfiltered_sce$detected >= 200)]                  
```
</details>
  
![mitochondrial](https://user-images.githubusercontent.com/56315895/221388029-14fc8b01-ed9a-4cc9-8f41-f8db53bb4b35.jpeg)
![compromised](https://user-images.githubusercontent.com/56315895/221388032-fd3a81f7-92c7-4d94-9987-34564051b2f3.jpeg)
![kept](https://user-images.githubusercontent.com/56315895/221388035-7084d5c7-f1b8-413e-97df-c2f6898c653f.jpeg)
