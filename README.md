
<!-- README.md is generated from README.Rmd. Please edit that file -->

# omnideconv

<!-- badges: start -->

[![R-CMD-check](https://github.com/icbi-lab/omnideconv/workflows/R-CMD-check/badge.svg)](https://github.com/icbi-lab/omnideconv/actions)
<!-- badges: end -->

The goal of omnideconv is to unify second generation immune
deconvolution methods in an R package.

## Installation

There are two ways to install it, the long one, installing all packages
needed for all methods and the short one, installing only the ones for
basic functionality.

If the short install is chosen, the required packages for methods are
installed when the methods are used. The same goes for required python
packages.

Install the CRAN package devtools and use it to install omnideconv from
[GitHub](https://github.com/):

``` r
install.packages("devtools")

# long/complete install
devtools::install_github("icbi-lab/omnideconv", dependencies = c("Imports",
    "Suggests"))
omnideconv::install_all_python()

# short/quick install
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

-   AutoGeneS (“autogenes”)
-   Bisque (“bisque”)
-   BSeq-sc (“bseqsc”)
-   CDSeq (“cdseq”)
-   CibersortX (“cibersortx”)
-   CPM (“cpm”)
-   DWLS (“dwls”)
-   MOMF (“momf”)
-   MuSiC (“music”)
-   Scaden (“scaden”)
-   SCDC (“scdc”)

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

|               |    B cell | Macrophage | Monocyte conventional | Monocyte non-conventional |   NK cell | T cell CD4 | T cell CD8 | T cell dividing | T cell regulatory |
|:--------------|----------:|-----------:|----------------------:|--------------------------:|----------:|-----------:|-----------:|----------------:|------------------:|
| HD30\_PBMC\_0 | 0.0800780 |      0e+00 |             0.2764050 |                 0.1427114 | 0.1385384 |  0.3316996 |  0.0277395 |       0.0026353 |         0.0001929 |
| HD30\_PBMC\_1 | 0.0735586 |      1e-06 |             0.3554604 |                 0.1236851 | 0.1253645 |  0.2968113 |  0.0243933 |       0.0007228 |         0.0000029 |
| HD30\_PBMC\_3 | 0.0862939 |      0e+00 |             0.2591685 |                 0.1246382 | 0.1305280 |  0.3595342 |  0.0360420 |       0.0037421 |         0.0000530 |
| HD30\_PBMC\_7 | 0.1156300 |      0e+00 |             0.3316376 |                 0.0226567 | 0.1140384 |  0.3596168 |  0.0545500 |       0.0017350 |         0.0001355 |
| HD31\_PBMC\_0 | 0.1141260 |      3e-07 |             0.4713667 |                 0.0024376 | 0.0511600 |  0.3026050 |  0.0571770 |       0.0010128 |         0.0001147 |
| HD31\_PBMC\_1 | 0.1143393 |      0e+00 |             0.4559677 |                 0.0000000 | 0.0408233 |  0.3235003 |  0.0651854 |       0.0001089 |         0.0000752 |
| HD31\_PBMC\_3 | 0.0865631 |      0e+00 |             0.4423767 |                 0.0000000 | 0.0486223 |  0.3520335 |  0.0686691 |       0.0015950 |         0.0001403 |
| HD31\_PBMC\_7 | 0.0930244 |      0e+00 |             0.4262222 |                 0.0000000 | 0.0583644 |  0.3587579 |  0.0632944 |       0.0000000 |         0.0003367 |

### Learn More

For more information and an example workflow see the vignette of this
package.

## Requirements

Most methods do not require additional software/tokens, but there are a
few exceptions:

-   A working version of Docker is required for CibersortX
-   A token for CibersortX is required from this website:
    <https://cibersortx.stanford.edu/>
-   The Cibersort source code is required for BSeq-sc (see tutorial in
    ?omnideconv::bseqsc\_config)

## Available methods, Licenses, Citations

Note that, while *omnideconv* itself is free ([GPL
3.0](https://github.com/icbi-lab/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

> CITATION

| method                                                 | license                                                                             | citation                                                                                                                                                                                                                                                                                                                                                                                      |
|--------------------------------------------------------|-------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [AutoGeneS](https://github.com/theislab/AutoGeneS/)    | free ([MIT](https://github.com/theislab/AutoGeneS/blob/master/LICENSE))             | Aliee, H., & Theis, F. (2021). AutoGeneS: Automatic gene selection using multi-objective optimization for RNA-seq deconvolution. <https://doi.org/10.1101/2020.02.21.940650>                                                                                                                                                                                                                  |
| [Bisque](https://github.com/cozygene/bisque)           | free ([GPL 3.0](https://github.com/cozygene/bisque/blob/master/DESCRIPTION))        | Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A., Garske, K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., & Halperin, E. (2020). Publisher Correction: Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nature Communications, 11(1), 2891. <https://doi.org/10.1038/s41467-020-16607-9>                            |
| [BSeq-sc](https://github.com/shenorrLab/bseqsc)        | free ([GPL 2.0](https://github.com/shenorrLab/bseqsc/blob/master/DESCRIPTION))      | Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. In Cell Systems (Vol. 3, Issue 4, pp. 346–360.e4). <https://doi.org/10.1016/j.cels.2016.08.011> |
| [CDSeq](https://github.com/kkang7/CDSeq_R_Package)     | free ([GPL 3.0](https://github.com/kkang7/CDSeq_R_Package/blob/master/DESCRIPTION)) | Kang, K., Huang, C., Li, Y. et al. CDSeqR: fast complete deconvolution for gene expression data from bulk tissues. BMC Bioinformatics 22, 262 (2021). <https://doi.org/10.1186/s12859-021-04186-5>                                                                                                                                                                                            |
| [CibersortX](https://cibersortx.stanford.edu/)         | free for non-commerical use only                                                    | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>                                                                                                                        |
| [CPM](https://github.com/amitfrish/scBio)              | free ([GPL 2.0](https://github.com/amitfrish/scBio/blob/master/DESCRIPTION))        | Frishberg, A., Peshes-Yaloz, N., Cohn, O., Rosentul, D., Steuerman, Y., Valadarsky, L., Yankovitz, G., Mandelboim, M., Iraqi, F. A., Amit, I., Mayo, L., Bacharach, E., & Gat-Viks, I. (2019). Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16(4), 327–332. <https://doi.org/10.1038/s41592-019-0355-5>                                                 |
| [DWLS](https://bitbucket.org/yuanlab/dwls/src/master/) | free ([GPL](https://bitbucket.org/yuanlab/dwls/src/master/DESCRIPTION))             | Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan, G.-C. (2019). Accurate estimation of cell-type composition from gene expression data. Nature Communications, 10(1), 2975. <https://doi.org/10.1038/s41467-019-10802-z>                                                                                                                                                             |
| [MOMF](https://github.com/sqsun/MOMF)                  | free ([GPL 3.0](https://github.com/sqsun/MOMF/blob/master/LICENSE.md))              | Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and flexible method for deconvoluting bulk RNAseq data with single-cell RNAseq data, 2019, DIO: 10.5281/zenodo.3373980                                                                                                                                                                                                                  |
| [MuSiC](https://github.com/xuranw/MuSiC/)              | free ([GPL 3.0](https://github.com/xuranw/MuSiC/blob/master/LICENSE))               | Wang, X., Park, J., Susztak, K., Zhang, N. R., & Li, M. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1), 380. <https://doi.org/10.1038/s41467-018-08023-x>                                                                                                                                                      |
| [Scaden](https://github.com/KevinMenden/scaden)        | free ([MIT](https://github.com/KevinMenden/scaden/blob/master/LICENSE))             | Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K., Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition analysis from tissue expression profiles. <https://doi.org/10.1101/659227>                                                                                                                                                                                   |
| [SCDC](https://github.com/meichendong/SCDC)            | NA                                                                                  | Dong, M., Thennavan, A., Urrutia, E., Li, Y., Perou, C. M., Zou, F., & Jiang, Y. (2020). SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references. Briefings in Bioinformatics. <https://doi.org/10.1093/bib/bbz166>                                                                                                                                        |
