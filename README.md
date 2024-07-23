# omnideconv

[![R-CMD-check](https://github.com/omnideconv/omnideconv/actions/workflows/test.yml/badge.svg)](https://github.com/omnideconv/omnideconv/actions/workflows/test.yml)
[![license](https://img.shields.io/badge/license-GPL3-blue.svg)](https://github.com/omnideconv/omnideconv/blob/master/LICENSE.md)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://omnideconv.github.io/omnideconv)
[![Codecov test coverage](https://codecov.io/gh/omnideconv/omnideconv/branch/main/graph/badge.svg)](https://app.codecov.io/gh/omnideconv/omnideconv?branch=main)

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

Upon the first loading, miniconda will be installed if not already present. A dedicated conda environment will be created to host the python-based methods.

## Available methods

The methods currently implemented in omnideconv are:

- AutoGeneS (“autogenes”)
- Bisque (“bisque”)
- BayesPrism ("bayesprism")
- BSeq-sc (“bseqsc”)
- CDSeq (“cdseq”)
- CIBERSORTx (“cibersortx”)
- CPM (“cpm”)
- DWLS (“dwls”)
- MOMF (“momf”)
- MuSiC (“music”)
- Scaden (“scaden”)
- SCDC (“scdc”)

## General usage

All the deconvolution methods included in omnideconv can be run in one step,
trough the function `deconvolute`, which takes in input the matrix of bulk RNAseq to be deconvolved (bulk_gene_expression), along with the training single cell expression matrix (single_cell_object) with the cell type annotations and sample information.

```r
deconvolution <- omnideconv::deconvolute(bulk_gene_expression, method,
                                         single_cell_object, cell_type_annotations, batch_ids)
```

## Signature matrix/model building

The methods AutoGeneS, BSeq-Sc, DWLS, CIBERSORTx, MOMF and Scaden first optimize their internal model, for example building a signature matrix, and then use this model to perform deconvolution. For these methods, the `build_model` function can be used. The obtained model can then be given in input to the `deconvolute` function, omitting the single cell data.

```r
signature <- omnideconv::build_model(single_cell_object, cell_type_annotations,
                                     batch_ids, method, bulk_gene_expression)

deconvolution <- omnideconv::deconvolute(bulk_gene_expression, signature)
```

The `deconvolute` function returns a sample x cell type matrix with the estimated cell fractions

## Input data

Different methods have different requirements in terms of input data.
This list has been compiled considering the methods documentation, described data procssing or
authors recommendation

| Method     | Single cell normalization | Bulk normalization |
| ---------- | ------------------------- | ------------------ |
| AutogeneS  | CPM                       | TPM                |
| BayesPrism | Counts                    | Counts             |
| Bisque     | Counts                    | Counts             |
| Bseq-Sc    | Counts                    | TPM                |
| CDseqR     | Counts                    | Counts             |
| CIBERSORTx | CPM                       | TPM                |
| CPM        | Counts                    | Counts             |
| DWLS       | Counts                    | TPM                |
| MOMF       | Counts                    | Counts             |
| MuSiC      | Counts                    | TPM                |
| Scaden     | Counts                    | TPM                |
| SCDC       | Counts                    | TPM                |

### Learn More

For more information and an example workflow see the vignette of this
package.

## Requirements

Most methods do not require additional software/tokens, but there are a
few exceptions:

- A working version of Docker or Singularity is required for CIBERSORTx
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

> Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Dietrich, Alexander and Merotto, Lorenzo and Pelz, Konstantin and Eder, Bernhard and Zackl, Constantin and Reinisch, Katharina and Edenhofer, Frank and Marini, Federico and Sturm, Gregor and List, Markus and Finotello, Francesca. (2024) https://doi.org/10.1101/2024.06.10.598226 

| method                                                 | license                                                                                    | citation                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [AutoGeneS](https://github.com/theislab/AutoGeneS/)    | free ([MIT](https://github.com/theislab/AutoGeneS/blob/master/LICENSE))                    | Aliee, H., & Theis, F. (2021). AutoGeneS: Automatic gene selection using multi-objective optimization for RNA-seq deconvolution. <https://doi.org/10.1101/2020.02.21.940650>                                                                                                                                                                                                                  |
| [BayesPrism](https://github.com/Danko-Lab/BayesPrism)  | free ([GPL 3.0](https://github.com/Danko-Lab/BayesPrism/blob/main/BayesPrism/DESCRIPTION)) | Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022). <https://doi.org/10.1038/s43018-022-00356-3>                                                                                                                      |
| [Bisque](https://github.com/cozygene/bisque)           | free ([GPL 3.0](https://github.com/cozygene/bisque/blob/master/DESCRIPTION))               | Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A., Garske, K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., & Halperin, E. (2020). Publisher Correction: Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nature Communications, 11(1), 2891. <https://doi.org/10.1038/s41467-020-16607-9>                            |
| [BSeq-sc](https://github.com/shenorrLab/bseqsc)        | free ([GPL 2.0](https://github.com/shenorrLab/bseqsc/blob/master/DESCRIPTION))             | Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. In Cell Systems (Vol. 3, Issue 4, pp. 346–360.e4). <https://doi.org/10.1016/j.cels.2016.08.011> |
| [CDSeq](https://github.com/kkang7/CDSeq_R_Package)     | free ([GPL 3.0](https://github.com/kkang7/CDSeq_R_Package/blob/master/DESCRIPTION))        | Kang, K., Huang, C., Li, Y. et al. CDSeqR: fast complete deconvolution for gene expression data from bulk tissues. BMC Bioinformatics 22, 262 (2021). <https://doi.org/10.1186/s12859-021-04186-5>                                                                                                                                                                                            |
| [CIBERSORTx](https://cibersortx.stanford.edu/)         | free for non-commerical use only                                                           | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>                                                                                                                        |
| [CPM](https://github.com/amitfrish/scBio)              | free ([GPL 2.0](https://github.com/amitfrish/scBio/blob/master/DESCRIPTION))               | Frishberg, A., Peshes-Yaloz, N., Cohn, O., Rosentul, D., Steuerman, Y., Valadarsky, L., Yankovitz, G., Mandelboim, M., Iraqi, F. A., Amit, I., Mayo, L., Bacharach, E., & Gat-Viks, I. (2019). Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16(4), 327–332. <https://doi.org/10.1038/s41592-019-0355-5>                                                 |
| [DWLS](https://bitbucket.org/yuanlab/dwls/src/master/) | free ([GPL](https://bitbucket.org/yuanlab/dwls/src/master/DESCRIPTION))                    | Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan, G.-C. (2019). Accurate estimation of cell-type composition from gene expression data. Nature Communications, 10(1), 2975. <https://doi.org/10.1038/s41467-019-10802-z>                                                                                                                                                             |
| [MOMF](https://github.com/sqsun/MOMF)                  | free ([GPL 3.0](https://github.com/sqsun/MOMF/blob/master/LICENSE.md))                     | Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and flexible method for deconvoluting bulk RNAseq data with single-cell RNAseq data, 2019, DIO: 10.5281/zenodo.3373980                                                                                                                                                                                                                  |
| [MuSiC](https://github.com/xuranw/MuSiC/)              | free ([GPL 3.0](https://github.com/xuranw/MuSiC/blob/master/LICENSE))                      | Wang, X., Park, J., Susztak, K., Zhang, N. R., & Li, M. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1), 380. <https://doi.org/10.1038/s41467-018-08023-x>                                                                                                                                                      |
| [Scaden](https://github.com/KevinMenden/scaden)        | free ([MIT](https://github.com/KevinMenden/scaden/blob/master/LICENSE))                    | Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K., Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition analysis from tissue expression profiles. <https://doi.org/10.1101/659227>                                                                                                                                                                                   |
| [SCDC](https://github.com/meichendong/SCDC)            | ([MIT](https://github.com/meichendong/SCDC/blob/master/README.md))                         | Dong, M., Thennavan, A., Urrutia, E., Li, Y., Perou, C. M., Zou, F., & Jiang, Y. (2020). SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references. Briefings in Bioinformatics. <https://doi.org/10.1093/bib/bbz166>                                                                                                                                        |
