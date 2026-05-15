#!/usr/bin/env Rscript
# Gene selection module (scrapper-backed) for omnibenchmark.
#
# Selects highly variable genes using scrapper::modelGeneVariances() on a
# TENx-format normalized expression matrix.

suppressPackageStartupMessages({
  library(HDF5Array)
  library(scrapper)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))

main <- function() {
  args <- parse_select_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5", "number_selected")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  mat <- TENxMatrix(args$input_h5, group = "matrix")
  mat <- as(mat, "dgCMatrix")
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(mat), ncol(mat)))

  gene_var <- modelGeneVariances(mat)
  hvgs <- chooseHighlyVariableGenes(
    gene_var$statistics$residuals,
    top = args$number_selected
    )
  sel_feats <- rownames(mat)[hvgs]
  cat(sprintf("  selected %d features\n", length(sel_feats)))

  out <- file.path(args$output_dir, paste0(args$name, "_normalized_selected.h5"))
  cat("output_file:", out, "\n")
  writeTENxMatrix(mat[sel_feats, ], out, group = "matrix")
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
