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
  
# perform normalization using scaling factors
# and save as a new SCE object
normalized_sce <- scuttle::logNormCounts(qcfiltered_sce)
# identify 2000 genes 
num_genes <- 2000
 
# model variance, partitioning into biological and technical variation
gene_variance <- scran::modelGeneVar(normalized_sce)
 
# get the most variable genes
hv_genes <- scran::getTopHVGs(gene_variance,
                               n = num_genes)
# calculate and save PCA results
normalized_sce <- scater::runPCA(
     normalized_sce,
     ncomponents = 50, # how many components to keep
     subset_row = hv_genes # use only the variable genes we chose
                  )

normalized_sce

# extract the PCA matrix
pca_matrix <- reducedDim(normalized_sce, "PCA")
 
# look at the shape of the matrix
dim(pca_matrix)

p0_normalized_sce <- normalized_sce
colnames(p0_normalized_sce) <- colData(p0_normalized_sce)$Barcode
readr::write_rds(p0_normalized_sce, file = "p0_output_sce_file.rds", compress = "gz")  
  
```
  
![mitochondrial](https://user-images.githubusercontent.com/56315895/221388029-14fc8b01-ed9a-4cc9-8f41-f8db53bb4b35.jpeg)
![compromised](https://user-images.githubusercontent.com/56315895/221388032-fd3a81f7-92c7-4d94-9987-34564051b2f3.jpeg)
![kept](https://user-images.githubusercontent.com/56315895/221388035-7084d5c7-f1b8-413e-97df-c2f6898c653f.jpeg)

</details>  

<details><summary> P5 QC </summary>
  
```
# create a table of statistics using emptyDropsCellRanger
droplet_df <- DropletUtils::emptyDropsCellRanger(counts(p5_sce))

# view rows where FDR is not `NA`
droplet_df[!is.na(droplet_df$FDR), ]

# filter droplets using `which` to prevent NA trouble
cells_to_retain <- which(droplet_df$FDR < 0.01)
filtered_sce <- p5_sce[, cells_to_retain]

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
  
# perform normalization using scaling factors
# and save as a new SCE object
normalized_sce <- scuttle::logNormCounts(qcfiltered_sce)
# identify 2000 genes 
num_genes <- 2000
 
# model variance, partitioning into biological and technical variation
gene_variance <- scran::modelGeneVar(normalized_sce)
 
# get the most variable genes
hv_genes <- scran::getTopHVGs(gene_variance,
                               n = num_genes)
# calculate and save PCA results
normalized_sce <- scater::runPCA(
     normalized_sce,
     ncomponents = 50, # how many components to keep
     subset_row = hv_genes # use only the variable genes we chose
  )
                  
normalized_sce

# extract the PCA matrix
pca_matrix <- reducedDim(normalized_sce, "PCA")
 
# look at the shape of the matrix
dim(pca_matrix)

p5_normalized_sce <- normalized_sce
colnames(p5_normalized_sce) <- colData(p5_normalized_sce)$Barcode
readr::write_rds(p5_normalized_sce, file = "p5_output_sce_file.rds", compress = "gz")  
  
```

![mitochondrial](https://user-images.githubusercontent.com/56315895/221389459-1353ae72-a7c6-4454-9bd8-c93eeb95c996.jpeg)
![compromised](https://user-images.githubusercontent.com/56315895/221389840-b2a651a1-6e70-4948-8339-3e74ac98ef19.jpeg)
![kept](https://user-images.githubusercontent.com/56315895/221389841-d0c09463-3c3a-4eeb-9ca5-83110af3e441.jpeg)

</details>  
  

<details><summary> P56 QC </summary>
  
```
# create a table of statistics using emptyDropsCellRanger
droplet_df <- DropletUtils::emptyDropsCellRanger(counts(p56_sce))

# view rows where FDR is not `NA`
droplet_df[!is.na(droplet_df$FDR), ]

# filter droplets using `which` to prevent NA trouble
cells_to_retain <- which(droplet_df$FDR < 0.01)
filtered_sce <- p56_sce[, cells_to_retain]

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
  
# perform normalization using scaling factors
# and save as a new SCE object
normalized_sce <- scuttle::logNormCounts(qcfiltered_sce)
# identify 2000 genes 
num_genes <- 2000
 
# model variance, partitioning into biological and technical variation
gene_variance <- scran::modelGeneVar(normalized_sce)
 
# get the most variable genes
hv_genes <- scran::getTopHVGs(gene_variance,
                               n = num_genes)
# calculate and save PCA results
normalized_sce <- scater::runPCA(
     normalized_sce,
     ncomponents = 50, # how many components to keep
     subset_row = hv_genes # use only the variable genes we chose
  )
                  
normalized_sce

# extract the PCA matrix
pca_matrix <- reducedDim(normalized_sce, "PCA")
 
# look at the shape of the matrix
dim(pca_matrix)

p56_normalized_sce <- normalized_sce
colnames(p56_normalized_sce) <- colData(p56_normalized_sce)$Barcode
readr::write_rds(p56_normalized_sce, file = "p56_output_sce_file.rds", compress = "gz")  
  
```
![mitochondrial](https://user-images.githubusercontent.com/56315895/221391134-e8b120bf-f573-4f20-a34e-7eb7d5863a92.jpeg)
![compromised](https://user-images.githubusercontent.com/56315895/221391135-0eca2aaf-37b5-40cd-b4df-90cb20d133eb.jpeg)
![kept](https://user-images.githubusercontent.com/56315895/221391136-4a5c6e5d-d017-4f4a-90c2-f90c600b737f.jpeg)

</details>    

```  
# Define directory to save integrated SCE object to
output_dir <- file.path("C:/Users/ihowa/Dropbox/justine/")  
              
# Define output file name for the integrated object
integrated_sce_file <- file.path(output_dir, "rms_integrated_subset.rds")
                       
# Vector of all the samples to read in:
sample_names <- c("p0",
                  "p5",
                  "p6")
                       
sce_paths <- c("p0_output_sce_file.rds",
                  "p5_output_sce_file.rds",
                  "p56_output_sce_file.rds")

# Use purrr::map() to read all files into a list at once
sce_list <- purrr::map(
  sce_paths, 
  readr::read_rds)
            
# Assign the sample names as the names for sce_list
names(sce_list) <- sample_names

# Use `purrr::map()` to quickly extract rowData column names for all SCEs
purrr::map(sce_list,
           \(x) colnames(rowData(x))) 
  
# Define vector of shared genes
shared_genes <- sce_list |>
# get rownames (genes) for each SCE in sce_list
purrr::map(rownames) |>
# reduce to the _intersection_ among lists
purrr::reduce(intersect)    

purrr::map(sce_list,
           \(x) colnames(colData(x)) )
 
format_sce <- function(sce, sample_name) {
  # Input arguments:
  ## sce: An SCE object to format
  ## sample_name: The SCE objects name
  # This function returns a formatted SCE object.
  
  ###### Ensure that we can identify the originating sample information ######
  # Add a column called `sample` that stores this information
  # This will be stored in `colData`
  sce$Sample <- sample_name
  
  
  ###### Ensure cell ids will be unique ######
  # Update the SCE object column names (cell ids) by prepending `sample_name`
  colnames(sce) <- glue::glue("{sample_name}-{colnames(sce)}")
        
  # Return the formatted SCE object
  return(sce)
}
   
# We can use `purrr::map2()` to loop over two list/vector arguments
sce_list_formatted <- purrr::map2(
  # Each "iteration" will march down the first two 
  #  arguments `sce_list` and `names(sce_list)` in order
  sce_list,
  names(sce_list),
  # Name of the function to run
  format_sce
)

sce_list_formatted <- sce_list_formatted |> purrr::map(\(sce) sce[shared_genes, ])	
	
# Print resulting list
sce_list_formatted
  
# Merge SCE objects 
merged_sce <- do.call(cbind, sce_list_formatted)

# Print the merged_sce object
merged_sce  
 
# Specify the number of genes to identify
num_genes <- 2000

# Calculate variation for each gene
gene_variance <- scran::modelGeneVar(merged_sce,
                                     # specify the grouping column:
                                     block = merged_sce$Sample)

# Get the top `num_genes` high-variance genes to use for dimension reduction
hv_genes <- scran::getTopHVGs(gene_variance,
                              n = num_genes)

# Use batchelor to calculate PCA for merged_sce, considering only
#  the high-variance genes
# We will need to include the argument `preserve.single = TRUE` to get
#  a single matrix with all samples and not separate matries for each sample
merged_pca <- batchelor::multiBatchPCA(merged_sce,
                                       subset.row = hv_genes,
                                       batch = merged_sce$Sample,
                                       preserve.single = TRUE)
              
# add PCA results to merged SCE object 
reducedDim(merged_sce, "merged_PCA") <- merged_pca[[1]]

# add merged_UMAP from merged_PCA
merged_sce <- scater::runUMAP(merged_sce,
                              dimred = "merged_PCA",
                              name = "merged_UMAP")
              
# UMAPs scaled together when calculated from the merged SCE
scater::plotReducedDim(merged_sce,
                       dimred = "merged_UMAP",
                       colour_by = "Sample",
                       # Some styling to help us see the points:
                       point_size = 0.5,
                       point_alpha = 0.2) +
  # Modify the legend key so its points are larger and easier to see
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  # Add a plot title
  ggtitle("UMAP calculated on merged_sce")
 
# integrate with fastMNN, again specifying only our high-variance genes
integrated_sce <- batchelor::fastMNN(
     merged_sce,
     batch = merged_sce$Sample, 
     subset.row = hv_genes
 ) 
              
# Make a new reducedDim named fastmnn_PCA from the corrected reducedDim in integrated_sce
reducedDim(merged_sce, "fastmnn_PCA") <- reducedDim(integrated_sce, "corrected")

# Calculate UMAP
merged_sce <- scater::runUMAP(
  merged_sce, 
  dimred = "fastmnn_PCA", 
  name = "fastmnn_UMAP"
)

scater::plotReducedDim(merged_sce,
                        # plot the fastMNN coordinates
                        dimred = "fastmnn_UMAP",
                        # color by Sample
                        colour_by = "Sample",
                        # Some styling to help us see the points:
                        point_size = 0.5,
                        point_alpha = 0.2) +
     # Modify legend so they key is larger and easier to see
     guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
     # add plot title
     ggtitle("UMAP after integration with fastMNN")                
 
# integrate with harmony, setting the argument `do_pca = FALSE`
# since we are providing a PCA matrix directly
harmony_pca <- harmony::HarmonyMatrix(
    data_mat = reducedDim(merged_sce, "merged_PCA")[,1:10],
    do_pca = FALSE,
    meta_data = colData(merged_sce),
    vars_use = 'Sample'
)
	       
# As before, calculate UMAP on this PCA matrix with appropriate names
merged_sce <- scater::runUMAP(merged_sce, 
                               dimred = "harmony_PCA", 
                               name   = "harmony_UMAP")
scater::plotReducedDim(merged_sce,
                        dimred = "harmony_UMAP",
                        colour_by = "Sample",
                        point_size = 0.5,
                        point_alpha = 0.2) +
     ggtitle("UMAP after integration with harmony") +
     guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))  
# Export to RDS file with "gz" compression
readr::write_rds(merged_sce, integrated_sce_file, compress = "gz")	       
```  
![umap_merged](https://user-images.githubusercontent.com/56315895/221396715-43b1d1a6-6531-49aa-801a-e3c5c06dee99.jpeg)
![umap_after_integration](https://user-images.githubusercontent.com/56315895/221396717-f96e0d20-be94-4406-ad93-81a441dd165a.jpeg)
![umap with harmony](https://user-images.githubusercontent.com/56315895/221400172-a3d4eba8-d4b9-43af-bb26-658741c03fe8.jpeg)
