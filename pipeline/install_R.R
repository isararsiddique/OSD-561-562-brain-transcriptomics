# Run once in the same R environment used by the pipeline.
if (getRversion() < "4.5.0" || getRversion() >= "4.6.0") {
  stop("The recorded study environment uses R 4.5.x / Bioconductor 3.22.")
}
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}
BiocManager::install(version="3.22", ask=FALSE)
BiocManager::install(c("statmod", "locfit"), ask=FALSE, update=FALSE)
for (record in list(c("limma", "3.66.0"), c("edgeR", "4.8.2"))) {
  package = record[1]; version = record[2]
  if (requireNamespace(package, quietly=TRUE) && as.character(packageVersion(package)) == version) next
  url = paste0("https://bioconductor.org/packages/3.22/bioc/src/contrib/", package, "_", version, ".tar.gz")
  install.packages(url, repos=NULL, type="source")
}
stopifnot(as.character(packageVersion("limma")) == "3.66.0",
          as.character(packageVersion("edgeR")) == "4.8.2")
