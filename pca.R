#!/usr/bin/env Rscript
# PCA module (scrapper-backed) for omnibenchmark.
#
# Output format: see docs/pca_output.md (neutral HDF5, format_version "1").
#
# Implementation notes
# --------------------
# - scrapper::runPca operates on the gene-by-cell matrix directly (it
#   internally centers/scales rows). No explicit scale step here.
# - Subsetting to selected genes happens here for now; this responsibility
#   should move to a dedicated upstream cleanup stage. See load_subset_matrix.

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(scrapper)
  library(BiocSingular)
  library(data.table)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))

run_pca <- function(X, args) {
  # X: gene-by-cell sparse matrix (rows = genes).
  set.seed(args$random_seed)

  if (args$solver == "irlba") {
    pca <- runPca(X, number = args$n_components, num.threads = 1L)
    # scrapper::runPca returns components (n_components x n_cells), rotation (n_genes x n_components)
    embedding <- t(pca$components)
    loadings  <- pca$rotation
    variance  <- as.numeric(pca$variance.explained)
    total_var <- if (!is.null(pca$total.variance)) as.numeric(pca$total.variance) else sum(variance)
  } else {
    # random / exact via BiocSingular::runSVD on the transposed (cells x genes) matrix
    bsparam <- switch(args$solver,
      random = RandomParam(),
      exact  = ExactParam(),
      stop("unknown solver: ", args$solver)
    )
    # runSVD expects cells-as-rows; center across genes (i.e. center=TRUE centers columns)
    svd <- runSVD(t(X), k = args$n_components, center = TRUE, BSPARAM = bsparam)
    embedding <- svd$u %*% diag(svd$d)        # (n_cells, n_components)
    loadings  <- svd$v                         # (n_genes, n_components)
    n         <- ncol(X)
    variance  <- svd$d^2 / (n - 1)
    # total variance: sum of per-gene variances, computed sparse-safe
    rs2 <- Matrix::rowSums(X^2)
    rs1 <- Matrix::rowSums(X)
    total_var <- sum((rs2 - rs1^2 / n) / (n - 1))
  }

  variance_ratio <- variance / total_var
  # decorate embeddings w/ row/colnames
  rownames(embedding) <- colnames(X) 
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))

  # loadings, etc are here in case needed as output later
  list(
    embedding      = embedding, 
    loadings       = matrix(as.double(loadings),  nrow = nrow(loadings)),
    variance       = as.double(variance),
    variance_ratio = as.double(variance_ratio)
  )
}

main <- function() {
  args <- parse_pca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5",
              "solver", "n_components", "random_seed")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix") # read into memory
  cat(sprintf("  matrix (genes x cells): %d x %d\n",
              nrow(m), ncol(m)))

  res <- run_pca(m, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(res$embedding), ncol(res$embedding),
              nrow(res$loadings),  ncol(res$loadings)))

  out <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
  cat("output_file:", out, "\n")
  fwrite(data.frame(cell_id = rownames(x), res$embedding), out, 
         sep = "\t", quote = FALSE, row.names = TRUE)
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
