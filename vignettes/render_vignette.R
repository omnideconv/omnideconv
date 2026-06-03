# Run this script from RStudio to pre-compile the vignette.
# Working directory must be the package root (omnideconv/).
#
# In RStudio: Session > Set Working Directory > To Project Directory, then source this file.

rmarkdown::render(
  input       = "vignettes/omnideconv_example.Rmd",
  output_file = "omnideconv_example.html",
  output_dir  = "vignettes",
  envir       = new.env()
)
