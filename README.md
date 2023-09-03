
# scDecouple
**`scDecouple`** is an R package that aims to decouple true cellular response from infected proportion bias in scCRISPR-seq. It models the distribution of gene expression profiles in perturbed cells and then iteratively finds the maximum likelihood of cell cluster proportions as well as the cellular response for each gRNA. **`scDecouple`** comes equipped with an integrated one-click data preprocessing pipeline, infected proportion bias quantification, cellular response estimation, and downstream analytical tools.

For more information, please refer to the [manuscript](https://doi.org/10.1101/2023.01.31.526445)

## Installation
Users can install the developmental version from [GitHub](https://github.com/MengQiuchen/scDecouple)

```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("MengQiuchen/scDecouple", build_vignettes = TRUE)
```

To load the installed **`scDecouple`** in R:

```
library(scDecouple)
```

## Input
**`scDecouple`** takes the following inputs:

1. `Y`: The perturbation group, cell-by-gene matrix.
2. `Z`: The control group, cell-by-gene matrix.
3. `using_rank` (optional): A logical value indicating whether to use rank-based selection of principal components (PCs) or value-based selection. The default is TRUE.
4. `top.sdev` (optional): When using_rank=TRUE, this parameter sets the threshold of standard deviations of PCs. The default is 30.
5. `top.multimodal` (optional): When using_rank=TRUE, this parameter sets the threshold of dip statistics. The default is 3.
6. `th.sdev` (optional): When using_rank=FALSE, this parameter sets the threshold of standard deviations of PCs. The default is 1.
7. `th.multimodal` (optional): When using_rank=FALSE, this parameter sets the threshold of dip statistics. The default is 0.0125.
8. `degenes.num` (optional): The number of differentially expressed genes (DE genes) for GO enrichment calculation. The default is 1500.
9. `seed_` (optional): The random seed for the EM algorithm. The default is 0.

## Usage
To use **`scDecouple`**, you can call the `run_scDecouple` function with your input data. Here's an example of how to use it:

```R
# Load the scDecouple package
library(scDecouple)

# Provide your input data Y and Z
# Example: run_scDecouple(Y, Z)
result <- run_scDecouple(Y, Z)

# Access the results in the 'result' object
# For example: result$processed_data, result$selected_pcs, etc.
```

For more details on the function and its parameters, please refer to the package documentation and vignettes.

## Documentation
You can access the detailed documentation and vignettes by running the following command in R:

```R
# Load the scDecouple package
library(scDecouple)

# View the documentation
?run_scDecouple

```


## Demo Example
Here's a demo example to demonstrate how to use **`scDecouple`**:

```R
# load basic R library
suppressPackageStartupMessages({
    library(tidyverse)
    library(enrichr)
    library(mixtools)
    library(diptest)
    library(Seurat)
})

# demo data
data(test_data)

# load the control group Z and the perturbation group Y
# Z represents the control group: Cells by Genes
# Y represents the perturbation group: Cells by Genes
ls()

# one-stop

res <- run_scDecouple(Y, Z)
```
Or you can use **`scDecouple`** step-by-step and customize the inputs and parameters for each step:

```R
# run step by step

## step1: preprocessing, to get log-normalized data and highly variable genes
res.step1 <- data_preprocessing(Y, Z)

## step2: get PCs which may have high infection bias
res.step2 <- pc_selection(Z.pca.mat = res.step1$Z.pca.results$x, Z.pca.sdev = res.step1$Z.pca.results$sdev, using_rank = TRUE)

## step3: decouple cellular response and infection bias
res.step3 <- scDecouple(Z.pca.mat = res.step1$Z.pca.results$x, Z.pca.rotation = res.step1$Z.pca.results$rotation, Y.norm.log.sub = res.step1$Y.norm.log.sub, select.pcs = res.step2)

## step4: downstream analysis: gene ranking and pathway enrichment
res.step4 <- downstream_analysis(Z.pca.rotation = res.step1$Z.pca.results$rotation, beta.PC.scDecouple = res.step3$beta.PC.scDecouple, Z.norm.log.sub = res.step1$Z.norm.log.sub, Y.norm.log.sub = res.step1$Y.norm.log.sub, genes.variable = res.step1$genes.variable)
```

This demo example illustrates the typical workflow of using **`scDecouple`** for scCRISPR-seq data analysis, including data preprocessing, PC selection, decoupling cellular response, and downstream analysis.

## Citation
If you use **`scDecouple`** in your research, please cite the corresponding manuscript:
[Link to manuscript](https://doi.org/10.1101/2023.01.31.526445)

## License
This package is distributed under the GNU General Public License.

For any questions or issues, please feel free to [report them on GitHub](https://github.com/MengQiuchen/scDecouple/issues).

Enjoy using **`scDecouple`** for your scCRISPR-seq data analysis!
