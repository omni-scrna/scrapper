#!/usr/bin/env Rscript
# PCA module (scrapper-backed) for omnibenchmark.
#
# scrapper wraps libscran, which ships one PCA: a C++ irlba over a sparse matrix.


suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(scrapper)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("PCA module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "PCA")            # the stage I/O contract
p <- add_argument(p, "--solver", type = "character", default = "irlba",
                  help = "name of solver (irlba)")
p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")

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

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


run_pca <- function(X, args) {
  # X: gene-by-cell sparse matrix (rows = genes).
  # runPca ignores R's RNG (its Rcpp export is rng=false); the seed it takes is for
  # irlba's initial random vector, so it goes in as an argument, not via set.seed().
  seed <- if (is.na(args$random_seed)) 5489L else args$random_seed  # 5489 = runPca's default

  if (!identical(args$solver, "irlba")) {
    stop("unknown solver: ", args$solver,
         " (this module only ships libscran's irlba)")
  }

  cat(sprintf("LOG: runPca seed = %d, num.threads = %d\n", seed, args$num_threads))
  pca <- runPca(X, number = args$n_components, seed = seed,
                num.threads = args$num_threads)
  # scrapper::runPca returns components (n_components x n_cells), rotation (n_genes x n_components)
  embedding <- t(pca$components)
  loadings  <- pca$rotation
  variance  <- as.numeric(pca$variance.explained)
  total_var <- if (!is.null(pca$total.variance)) as.numeric(pca$total.variance) else sum(variance)

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

  res <- run_pca(m, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
    nrow(res$embedding), ncol(res$embedding),
    nrow(res$loadings),  ncol(res$loadings)))

  out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
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
