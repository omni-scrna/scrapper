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

# arg parsing
source("src/common/cli.R")
p <- arg_parser("PCA module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "PCA")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--solver", type = "character", help = "name of solver")
p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")
args <- parse_args(p)                    # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


run_pca <- function(X, args) {
  # X: gene-by-cell sparse matrix (rows = genes).
  set.seed(args$random_seed)

  if (args$solver == "irlba") {
    pca <- runPca(X, number = args$n_components, num.threads = 1L, seed = args$random_seed)
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
  # decorate embeddings/loadings w/ row/colnames
  rownames(embedding) <- colnames(X) 
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))
  rownames(loadings)  <- rownames(X)
  colnames(loadings)  <- paste0("PC", seq_len(ncol(loadings)))

  # loadings, etc are here in case needed as output later
  list(
    embedding      = embedding, 
    loadings       = loadings,
    variance       = as.double(variance),
    variance_ratio = as.double(variance_ratio)
  )
}


main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  m <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
  m <- as(m, "dgCMatrix")
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))

  res <- run_pca(m, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
    nrow(res$embedding), ncol(res$embedding),
    nrow(res$loadings),  ncol(res$loadings)))

  out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_embedding.tsv", args$name))
  fwrite(data.frame(cell_id = rownames(res$embedding), res$embedding), out_embeddings_tsv,
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_embeddings_tsv))

  out_loadings_tsv <- file.path(args$output_dir, sprintf("%s_loadings.tsv", args$name))
  fwrite(data.frame(gene = rownames(res$loadings), res$loadings), out_loadings_tsv,
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_loadings_tsv))
}

if (sys.nframe() == 0L) {
  main()
}
