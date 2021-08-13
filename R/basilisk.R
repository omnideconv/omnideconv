# scaden_env <- basilisk::BasiliskEnvironment(envname="scaden_env",
#                                pkgname="omnideconv",
#                                packages=c("scaden==1.1.2")
# )


autogenes_env <- basilisk::BasiliskEnvironment(
  envname = "autogenes_env",
  pkgname = "omnideconv",
  packages = c("anndata==0.7.6"),
  channels = c("bioconda", "conda-forge"),
  pip = c("autogenes==1.0.4")
)
