#!/usr/bin/env Rscript
# Gene selection module (scrapper-backed) for omnibenchmark.
#
# Selects highly variable genes using scrapper::modelGeneVariances() on a
# TENx-format normalized expression matrix.

suppressPackageStartupMessages({
  library(HDF5Array)
  library(scrapper)
})


# arg parsing
source("src/common/cli.R")
p <- arg_parser("FEAT module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "FEAT")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--number_selected", type = "integer", help = "number of PCs")
args <- parse_args(p)                    # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  mat <- TENxMatrix(args$normalized_h5, group = "matrix")
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
