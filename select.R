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
p <- add_argument(p, "--backend", type = "character", default = "memory",
                  help = "counts backend: memory (anndataR, in-RAM) or delayed (H5SparseMatrix, out-of-core)")
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
  # TENxMatrix is already a DelayedArray; "memory" is the coercion, not the read.
  if (args$backend == "memory") mat <- as(mat, "dgCMatrix")
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
  # level=0: no gzip. Left at the default NULL, writeTENxMatrix uses
  # getHDF5DumpCompressionLevel() == 6, paid on every write here and again on
  # every read downstream. These are pipeline intermediates, not artifacts to
  # ship, so disk is the cheap axis to spend.
  writeTENxMatrix(mat[sel_feats, ], out, group = "matrix", level = 0)
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
