#!/usr/bin/env Rscript
# PCA module for Omnibenchmark using scrapper and BiocSingular.
#
# Input:
#   Normalized, HVG-selected HDF5 matrix (genes x cells).
#
# Implementations:
#   - irlba:  scrapper::runPca() using implicitly restarted Lanczos bidiagonalization.
#   - random: BiocSingular::runSVD() using randomized SVD.
#
#
# Preprocessing:
#   PCA is performed on centered, unscaled expression data.
#
# Parameter mapping for the two PCA solvers:
# - IRLBA (scrapper::runPca):
#     n_iter -> iterations (maximum Lanczos restart iterations)
#     n_oversamples -> extra.work (additional working-subspace dimensions)
# - Randomized SVD (BiocSingular::runSVD):
#     n_iter -> q (number of power iterations)
#     n_oversamples -> p (number of oversampling vectors)
#
# Outputs:
#   - <name>_embedding.tsv:     PC scores (cells x components).
#   - <name>_loadings.tsv: PC loadings (genes x components).
#   - <name>_runtime.csv:  PCA runtime and execution parameters.
#
# The module preserves cell and gene identifiers for downstream comparisons.
suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(scrapper)
  library(BiocSingular)
  library(BiocParallel)
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
p <- add_argument(p, "--n_iter", type = "integer", default = 2L, help = "number of power iterations for randomized SVD")
p <- add_argument(p, "--n_oversamples", type = "integer", default = 10L, help = "number of oversamples for randomized SVD")
p <- add_argument(p, "--random_seed", type = "integer", default = 42L, help = "seed")
# Thread count. scrapper's PCA is a C++ irlba over Eigen (assorthead), built without
# EIGEN_USE_BLAS. Neither the sparse nor the dense path ever enters R's BLAS.
# Snakemake exports OMP_NUM_THREADS = <rule threads> (default 1), which we inherit.
p <- add_argument(p, "--num_threads", type = "integer",
                  default = as.integer(Sys.getenv("OMP_NUM_THREADS", "1")),
                  help = "threads for runPca (default: OMP_NUM_THREADS)")
# Dense vs sparse control: switches the tatami representation handed to irlba.
p <- add_argument(p, "--dense", type = "character", default = "false",
                  help = "materialise the matrix dense before PCA (true/false)")
args <- parse_args(p)                    # argparser's own parser

#thread control
make_bpparam <- function(n_threads) {
  if (n_threads < 1L) {stop("--num_threads must be >= 1")}
  if (n_threads == 1L) {SerialParam()} else if (.Platform$OS.type == "unix") {MulticoreParam(workers = n_threads, progressbar = FALSE)}
  else {SnowParam(workers = n_threads, progressbar = FALSE)}}

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
    pca <- runPca(X, number = args$n_components, num.threads = args$num_threads)
    # scrapper::runPca returns components (n_components x n_cells), rotation (n_genes x n_components)
    embedding <- t(pca$components)
    loadings  <- pca$rotation
    variance  <- as.numeric(pca$variance.explained)
    total_var <- if (!is.null(pca$total.variance)) as.numeric(pca$total.variance) else sum(variance)
  } else {
    # random / exact via BiocSingular::runSVD on the transposed (cells x genes) matrix
    bsparam <- switch(args$solver,
      random = RandomParam(p = args$n_oversamples, q = args$n_iter),
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
  # Coerce straight from the DelayedArray to the target representation, so peak
  # memory reflects the format under test rather than a sparse copy plus a
  # dense one.
  if (identical(args$dense, "true")) {
    m <- as(m, "matrix")
    cat(sprintf("LOG: dense matrix: %.0f MB\n", as.numeric(object.size(m)) / 1e6))
  } else {
    m <- as(m, "dgCMatrix")
  }
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))
# runtime recording 
  pca_start <- Sys.time()
  res <- run_pca(m, args)
  pca_end <- Sys.time()
  pca_runtime_s <- as.numeric(difftime(pca_end, pca_start, units = "secs"))
  cat(sprintf("  PCA runtime: %.6f seconds\n", pca_runtime_s))
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
#save runtime 
  out_runtime_csv <- file.path(args$output_dir, sprintf("%s_runtime.csv", args$name))
  fwrite(data.frame(name = args$name, solver = args$solver, random_seed = args$random_seed, n_components = args$n_components, n_iter = args$n_iter,
      n_oversamples = args$n_oversamples, num_threads = args$num_threads, pca_start = format(pca_start, "%Y-%m-%dT%H:%M:%OS3%z"),
      pca_end = format(pca_end, "%Y-%m-%dT%H:%M:%OS3%z"), runtime_s = pca_runtime_s), out_runtime_csv)
  cat(sprintf("  wrote: %s\n", out_runtime_csv))}

if (sys.nframe() == 0L) {
  main()
}
