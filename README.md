
<!-- README.md is generated from README.Rmd. Please edit that file -->

# omnideconv

<!-- badges: start -->

[![R-CMD-check](https://github.com/icbi-lab/omnideconv/workflows/R-CMD-check/badge.svg)](https://github.com/icbi-lab/omnideconv/actions)
<!-- badges: end -->

The goal of omnideconv is to unify second generation immune
deconvolution methods in an R package.

## Installation

Install the CRAN package devtools and use it to install omnideconv from
[GitHub](https://github.com/):

``` r
install.packages("devtools")
devtools::install_github("icbi-lab/omnideconv")
```

## Usage

The main functions of this package are to build a signature matrix and
deconvolute bulk RNA-seq data based on this signature matrix.

The basic workflow is the following:

### 1. Build a Signature Matrix

``` r
omnideconv::build_model(single_cell_data, cell_type_annotations, 
    method = "dwls")
```

In this method, single\_cell\_data is a matrix of single cell RNA-seq
data with genes in rows and cells in columns, cell\_type\_annotations a
vector that contains an annotation for each cell in your
single\_cell\_data, and the possible methods are:

-   Bisque (“bisque”)
-   DWLS (“dwls”)
-   MOMF (“momf”)
-   Scaden (“scaden”)
-   CibersortX (“cibersortx”)

### 2. Deconvolute

``` r
omnideconv::deconvolute(bulk, signature_matrix, method = "dwls")
```

Here, bulk is your bulk RNA-seq data as a matrix with genes in rows and
samples in column, signature\_matrix is the signature matrix you created
in the previous step and the method can again be one of the four methods
listed above.

This is, what the cell type properties in your bulk RNA-seq data
computed in the deconvolution step could look like:

|               |  B cell | Macrophage | Monocyte conventional | Monocyte non-conventional | NK cell | T cell CD4 | T cell CD8 | T cell dividing | T cell regulatory |
|:--------------|--------:|-----------:|----------------------:|--------------------------:|--------:|-----------:|-----------:|----------------:|------------------:|
| HD30\_PBMC\_0 | 0.12816 |          0 |               0.51934 |                   0.13378 | 0.18987 |    0.00000 |    0.00561 |               0 |           0.02325 |
| HD30\_PBMC\_1 | 0.09559 |          0 |               0.67510 |                   0.00696 | 0.19640 |    0.00000 |    0.00000 |               0 |           0.02595 |
| HD30\_PBMC\_3 | 0.17181 |          0 |               0.51368 |                   0.05744 | 0.16574 |    0.07258 |    0.00000 |               0 |           0.01876 |
| HD30\_PBMC\_7 | 0.28474 |          0 |               0.44150 |                   0.03053 | 0.17612 |    0.05122 |    0.00000 |               0 |           0.01590 |
| HD31\_PBMC\_0 | 0.08652 |          0 |               0.71617 |                   0.00000 | 0.05495 |    0.10743 |    0.02555 |               0 |           0.00937 |
| HD31\_PBMC\_1 | 0.09561 |          0 |               0.68077 |                   0.00000 | 0.04447 |    0.14894 |    0.02204 |               0 |           0.00818 |
| HD31\_PBMC\_3 | 0.05961 |          0 |               0.71674 |                   0.00000 | 0.05728 |    0.12130 |    0.02990 |               0 |           0.01516 |
| HD31\_PBMC\_7 | 0.06928 |          0 |               0.69114 |                   0.00000 | 0.07708 |    0.08902 |    0.04345 |               0 |           0.03003 |

### Learn More

For more information and an example workflow see the vignette of this
package.

## Available Methods

More information about the second generation deconvolution metods
unified in this package as well as their citations can be obtained from
the list below.

-   Bisque: Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A.,
    Garske, K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., &
    Halperin, E. (2020). Publisher Correction: Accurate estimation of
    cell composition in bulk expression through robust integration of
    single-cell information. Nature Communications, 11(1), 2891.
    <https://doi.org/10.1038/s41467-020-16607-9>
-   CibersortX: Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J.,
    Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015).
    Robust enumeration of cell subsets from tissue expression profiles.
    Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>
-   DWLS: Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan,
    G.-C. (2019). Accurate estimation of cell-type composition from gene
    expression data. Nature Communications, 10(1), 2975.
    <https://doi.org/10.1038/s41467-019-10802-z>
-   MOMF: Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and
    flexible method for deconvoluting bulk RNAseq data with single-cell
    RNAseq data, 2019, DOI: 10.5281/zenodo.3373980
-   Scaden: Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K.,
    Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition
    analysis from tissue expression profiles.
    <https://doi.org/10.1101/659227>
