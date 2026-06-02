# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Build and check

```r
# Full R CMD check (equivalent to CI)
R CMD check .

# Or within R:
devtools::check()

# Document (regenerate man/ from roxygen2 comments)
devtools::document()
```

### Tests

```r
# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-b-buildmodel.R")

# Run tests matching a pattern
devtools::test(filter = "deconvolution")
```

Test files are prefixed to indicate order: `test-a-python`, `test-b-buildmodel`, `test-c-deconvolution`, `test-d-dataconversion`. Small test data lives in `inst/small_test_data/`.

### Linting / style

Pre-commit hooks enforce `styler` (tidyverse style) and other checks. To run manually:

```r
styler::style_pkg()          # format all R files
lintr::lint_package()        # lint
```

Pre-commit must be installed separately (`pre-commit install`) for hooks to run on commit.

## Architecture

### Package structure

`omnideconv` is an R package that wraps 12 second-generation cell-type deconvolution methods under a unified interface. The two public entry points are:

- **`build_model(single_cell_object, cell_type_annotations, method, ...)`** — builds a signature matrix or model from single-cell data. Only relevant for methods that have a separate training step: AutoGeneS, BSeq-sc, DWLS, CIBERSORTx, MOMF, Scaden.
- **`deconvolute(bulk_gene_expression, model, method, single_cell_object, ...)`** — runs deconvolution on bulk RNA-seq. For single-step methods (BayesPrism, Bisque, CDSeq, CPM, MuSiC, SCDC) `model` is ignored.

Both functions dispatch via `switch(method, ...)` in `R/deconvolution_algorithms.R`. Each method has its own file (`R/AutoGeneS.R`, `R/BayesPrism.R`, etc.) containing `build_model_<method>` and `deconvolute_<method>` functions.

### Input / output conventions

- All inputs are normalized to a plain R `matrix` internally via `convert_to_matrix()` in `R/data_processing.R`. Accepted input types: `matrix`, `data.frame`, `SingleCellExperiment`, `AnnData` (R6).
- Matrices are always **rows = genes, columns = cells/samples**.
- `deconvolute()` returns a matrix with **rows = samples, columns = cell types**, sorted alphabetically.
- Special characters in gene/cell-type names are escaped before passing to methods and restored after — see `escape_special_chars()` / `deescape_special_chars()` in `R/helper.R`.

### Python-based methods

AutoGeneS and Scaden use Python via `reticulate`. On package load (`R/zzz.R`), a conda environment named `r-omnideconv` (Python 3.8) is created/activated automatically. To install all Python dependencies at once:

```r
omnideconv::install_all_python()
```

### CIBERSORTx

Requires a running Docker daemon or Apptainer/Singularity. The token and email must be obtained from https://cibersortx.stanford.edu/ and passed as arguments. Check availability with `omnideconv::check_container()`.

### On-demand dependency installation

Most method-specific R packages (BayesPrism, DWLS, MuSiC, etc.) are listed in `Suggests`, not `Imports`. `check_and_install()` in `R/deconvolution_algorithms.R` prompts the user to install them the first time a method is called. For a complete installation:

```r
pak::pkg_install("omnideconv/omnideconv", dependencies = TRUE)
```

### CI

CI runs only on Ubuntu 22.04 (Mac and Windows are explicitly unsupported due to dependency constraints). Workflows are in `.github/workflows/`. The `test.yml` workflow installs method packages from their respective GitHub forks under the `omnideconv` organization.
