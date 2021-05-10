
<!-- README.md is generated from README.Rmd. Please edit that file -->

# immunedeconv2

<!-- badges: start -->

[![R-CMD-check](https://github.com/icbi-lab/immunedeconv2/workflows/R-CMD-check/badge.svg)](https://github.com/icbi-lab/immunedeconv2/actions)
<!-- badges: end -->

The goal of immunedeconv2 is to unify second generation immune
deconvolution methods in an R package.

## Installation

Install the CRAN package devtools and use it to install immunedeconv2
from [GitHub](https://github.com/):

``` r
install.packages("devtools")
devtools::install_github("icbi-lab/immunedeconv2")
```

## Usage

The main functions of this package are to build a signature matrix and
deconvolute bulk RNA-seq data based on this signature matrix.

The basic workflow is the following:

### 1\. Build a Signature Matrix

``` r
immunedeconv2::build_model(single_cell_data, cell_type_annotations, 
    method = "dwls")
```

In this method, single\_cell\_data is a matrix of single cell RNA-seq
data with genes in rows and cells in columns, cell\_type\_annotations a
vector that contains an annotation for each cell in your
single\_cell\_data, and the possible methods are:

  - Bisque (“bisque”)
  - DWLS (“dwls”)
  - MOMF (“momf”)
  - Scaden (“scaden”)

### 2\. Deconvolute

``` r
immunedeconv2::deconvolute(bulk, signature_matrix, method = "dwls")
```

Here, bulk is your bulk RNA-seq data as a matrix with genes in rows and
samples in column, signature\_matrix is the signature matrix you created
in the previous step and the method can again be one of the four methods
listed above.

This is, what the cell type properties in your bulk RNA-seq data
computed in the deconvolution step could look like:

|                            | HD30\_PBMC\_0 | HD30\_PBMC\_1 | HD30\_PBMC\_3 | HD30\_PBMC\_7 | HD31\_PBMC\_0 | HD31\_PBMC\_1 | HD31\_PBMC\_3 | HD31\_PBMC\_7 |
| :------------------------- | ------------: | ------------: | ------------: | ------------: | ------------: | ------------: | ------------: | ------------: |
| T\_cell\_CD4               |       0.00000 |       0.00000 |       0.07258 |       0.05122 |       0.10743 |       0.14894 |       0.12130 |       0.08902 |
| Macrophage                 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |
| T\_cell\_dividing          |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |
| T\_cell\_CD8               |       0.00561 |       0.00000 |       0.00000 |       0.00000 |       0.02555 |       0.02204 |       0.02990 |       0.04345 |
| Monocyte\_conventional     |       0.51934 |       0.67510 |       0.51368 |       0.44150 |       0.71617 |       0.68077 |       0.71674 |       0.69114 |
| Monocyte\_non-conventional |       0.13378 |       0.00696 |       0.05744 |       0.03053 |       0.00000 |       0.00000 |       0.00000 |       0.00000 |
| T\_cell\_regulatory        |       0.02325 |       0.02595 |       0.01876 |       0.01590 |       0.00937 |       0.00818 |       0.01516 |       0.03003 |
| B\_cell                    |       0.12816 |       0.09559 |       0.17181 |       0.28474 |       0.08652 |       0.09561 |       0.05961 |       0.06928 |
| NK\_cell                   |       0.18987 |       0.19640 |       0.16574 |       0.17612 |       0.05495 |       0.04447 |       0.05728 |       0.07708 |

### Learn More

For more information and an example workflow see the vignette of this
package.

## Available Methods

More information about the second generation deconvolution metods
unified in this package as well as their citations can be obtained from
the list below.

  - Bisque: Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A., Garske,
    K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., & Halperin, E.
    (2020). Publisher Correction: Accurate estimation of cell
    composition in bulk expression through robust integration of
    single-cell information. Nature Communications, 11(1), 2891.
    <https://doi.org/10.1038/s41467-020-16607-9>
  - DWLS: Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan,
    G.-C. (2019). Accurate estimation of cell-type composition from gene
    expression data. Nature Communications, 10(1), 2975.
    <https://doi.org/10.1038/s41467-019-10802-z>
  - MOMF: Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and
    flexible method for deconvoluting bulk RNAseq data with single-cell
    RNAseq data, 2019, DIO: 10.5281/zenodo.3373980
  - Scaden: Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K.,
    Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition
    analysis from tissue expression profiles.
    <https://doi.org/10.1101/659227>
