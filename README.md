# omnideconv

[![R-CMD-check](https://github.com/omnideconv/omnideconv/actions/workflows/test.yml/badge.svg)](https://github.com/omnideconv/omnideconv/actions/workflows/test.yml)
[![license](https://img.shields.io/badge/license-GPL3-blue.svg)](https://github.com/omnideconv/omnideconv/blob/master/LICENSE.md)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://omnideconv.github.io/omnideconv)
[![Codecov test coverage](https://codecov.io/gh/omnideconv/omnideconv/branch/master/graph/badge.svg)](https://app.codecov.io/gh/omnideconv/omnideconv?branch=master)

The goal of omnideconv is to unify second generation cell-type deconvolution methods in an R package.

## Installation

There are two ways to install `omnideconv`:

- The _minimal_ installation installs only the dependencies required for the basic functionalities. All deconvolution
  methods need to be installed on-demand.
- The _complete_ installation installs all dependencies including all deconvolution methods. This may take
  a considerable time.

Since not all dependencies are on CRAN or Bioconductor, `omnideconv` is available from GitHub only.
We recommend installing it through the [pak](https://github.com/r-lib/pak) package manager:

```r
# install the `pak` package manager
install.packages("pak")

# minimal installation
pak::pkg_install("omnideconv/omnideconv")

# complete installation, including Python dependencies
pak::pkg_install("omnideconv/omnideconv", dependencies = TRUE)
omnideconv::install_all_python()
```

## Usage

The main functions of this package are to build a signature matrix and
deconvolute bulk RNA-seq data based on this signature matrix.

The basic workflow (in this example using the method “dwls”) is the
following:

### 1. Build a Signature Matrix

```r
omnideconv::build_model(single_cell_data_1, cell_type_annotations_1,
    method = "dwls")
```

In this method, single_cell_data is a matrix of single cell RNA-seq
data with genes in rows and cells in columns, cell_type_annotations a
vector that contains an annotation for each cell in your
single_cell_data, and the possible methods are:

- AutoGeneS (“autogenes”)
- Bisque (“bisque”)
- BSeq-sc (“bseqsc”)
- CDSeq (“cdseq”)
- CIBERSORTx (“cibersortx”)
- CPM (“cpm”)
- DWLS (“dwls”)
- MOMF (“momf”)
- MuSiC (“music”)
- Scaden (“scaden”)
- SCDC (“scdc”)

### 2. Deconvolute

```r
omnideconv::deconvolute(bulk, signature_matrix, method = "dwls")
```

Here, bulk is your bulk RNA-seq data as a matrix with genes in rows and
samples in column, signature_matrix is the signature matrix you created
in the previous step and the method can again be one of the four methods
listed above.

This is, what the cell type properties in your bulk RNA-seq data
computed in the deconvolution step could look like:

|             |     B | CD4 T | CD8 T |    DC |  Mono |    NK |
| :---------- | ----: | ----: | ----: | ----: | ----: | ----: |
| HD30_PBMC_0 | 0.087 | 0.444 | 0.318 | 0.047 | 0.037 | 0.066 |
| HD30_PBMC_1 | 0.086 | 0.439 | 0.311 | 0.052 | 0.039 | 0.073 |
| HD30_PBMC_3 | 0.085 | 0.540 | 0.244 | 0.044 | 0.032 | 0.055 |
| HD30_PBMC_7 | 0.091 | 0.472 | 0.295 | 0.048 | 0.028 | 0.067 |
| HD31_PBMC_0 | 0.080 | 0.617 | 0.167 | 0.041 | 0.045 | 0.049 |
| HD31_PBMC_1 | 0.081 | 0.566 | 0.214 | 0.041 | 0.042 | 0.056 |
| HD31_PBMC_3 | 0.080 | 0.525 | 0.243 | 0.043 | 0.044 | 0.065 |
| HD31_PBMC_7 | 0.053 | 0.851 | 0.000 | 0.000 | 0.059 | 0.037 |

### Learn More

For more information and an example workflow see the vignette of this
package.

## Requirements

Most methods do not require additional software/tokens, but there are a
few exceptions:

- A working version of Docker is required for CIBERSORTx
- A token for CIBERSORTx is required from this website:
  <https://cibersortx.stanford.edu/>
- The CIBERSORT source code is required for BSeq-sc (see tutorial in
  ?omnideconv::bseqsc_config)

## Available methods, Licenses, Citations

Note that, while _omnideconv_ itself is free ([GPL
3.0](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

> CITATION

| method                                                 | license                                                                             | citation                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------------------------------------------------ | ----------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [AutoGeneS](https://github.com/theislab/AutoGeneS/)    | free ([MIT](https://github.com/theislab/AutoGeneS/blob/master/LICENSE))             | Aliee, H., & Theis, F. (2021). AutoGeneS: Automatic gene selection using multi-objective optimization for RNA-seq deconvolution. <https://doi.org/10.1101/2020.02.21.940650>                                                                                                                                                                                                                  |
| [Bisque](https://github.com/cozygene/bisque)           | free ([GPL 3.0](https://github.com/cozygene/bisque/blob/master/DESCRIPTION))        | Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A., Garske, K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., & Halperin, E. (2020). Publisher Correction: Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nature Communications, 11(1), 2891. <https://doi.org/10.1038/s41467-020-16607-9>                            |
| [BSeq-sc](https://github.com/shenorrLab/bseqsc)        | free ([GPL 2.0](https://github.com/shenorrLab/bseqsc/blob/master/DESCRIPTION))      | Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. In Cell Systems (Vol. 3, Issue 4, pp. 346–360.e4). <https://doi.org/10.1016/j.cels.2016.08.011> |
| [CDSeq](https://github.com/kkang7/CDSeq_R_Package)     | free ([GPL 3.0](https://github.com/kkang7/CDSeq_R_Package/blob/master/DESCRIPTION)) | Kang, K., Huang, C., Li, Y. et al. CDSeqR: fast complete deconvolution for gene expression data from bulk tissues. BMC Bioinformatics 22, 262 (2021). <https://doi.org/10.1186/s12859-021-04186-5>                                                                                                                                                                                            |
| [CIBERSORTx](https://cibersortx.stanford.edu/)         | free for non-commerical use only                                                    | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>                                                                                                                        |
| [CPM](https://github.com/amitfrish/scBio)              | free ([GPL 2.0](https://github.com/amitfrish/scBio/blob/master/DESCRIPTION))        | Frishberg, A., Peshes-Yaloz, N., Cohn, O., Rosentul, D., Steuerman, Y., Valadarsky, L., Yankovitz, G., Mandelboim, M., Iraqi, F. A., Amit, I., Mayo, L., Bacharach, E., & Gat-Viks, I. (2019). Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16(4), 327–332. <https://doi.org/10.1038/s41592-019-0355-5>                                                 |
| [DWLS](https://bitbucket.org/yuanlab/dwls/src/master/) | free ([GPL](https://bitbucket.org/yuanlab/dwls/src/master/DESCRIPTION))             | Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan, G.-C. (2019). Accurate estimation of cell-type composition from gene expression data. Nature Communications, 10(1), 2975. <https://doi.org/10.1038/s41467-019-10802-z>                                                                                                                                                             |
| [MOMF](https://github.com/sqsun/MOMF)                  | free ([GPL 3.0](https://github.com/sqsun/MOMF/blob/master/LICENSE.md))              | Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and flexible method for deconvoluting bulk RNAseq data with single-cell RNAseq data, 2019, DIO: 10.5281/zenodo.3373980                                                                                                                                                                                                                  |
| [MuSiC](https://github.com/xuranw/MuSiC/)              | free ([GPL 3.0](https://github.com/xuranw/MuSiC/blob/master/LICENSE))               | Wang, X., Park, J., Susztak, K., Zhang, N. R., & Li, M. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1), 380. <https://doi.org/10.1038/s41467-018-08023-x>                                                                                                                                                      |
| [Scaden](https://github.com/KevinMenden/scaden)        | free ([MIT](https://github.com/KevinMenden/scaden/blob/master/LICENSE))             | Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K., Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition analysis from tissue expression profiles. <https://doi.org/10.1101/659227>                                                                                                                                                                                   |
| [SCDC](https://github.com/meichendong/SCDC)            | ([MIT](https://github.com/meichendong/SCDC/blob/master/README.md))                  | Dong, M., Thennavan, A., Urrutia, E., Li, Y., Perou, C. M., Zou, F., & Jiang, Y. (2020). SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references. Briefings in Bioinformatics. <https://doi.org/10.1093/bib/bbz166>                                                                                                                                        |
