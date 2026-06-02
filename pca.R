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
  library(rhdf5)
  library(yaml)
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

run_pca_per_batch <- function(m, args) {
  # load cell ids and batch assignments from raw data
  cell_ids_all <- as.character(h5read(args$rawdata_h5ad, "obs/_index"))
  batch_raw <- h5read(args$rawdata_h5ad, paste0("obs/", args$batch_variable))|>
    as.character()
  
  # keep only cells present in normalized data
  cell_ids_keep <- colnames(m)
  idx_keep <- match(cell_ids_keep, cell_ids_all)
  batch_keep <- batch_raw[idx_keep]

  # get all unique batch labels
  batches <- unique(batch_keep)

  cat(sprintf("  per-batch PCA: %d batches (%s)\n",
    length(batches), paste(batches, collapse = ", ")))

  # run PCA separately for each batch (we should consider parallelizing this)
  pca_res_ls <- batches |> lapply(function(b) {
    cells_batch <- which(batch_keep == b)
    m_batch <- m[, cells_batch]
    cat(sprintf("    batch '%s': %d cells\n", b, ncol(m_batch)))
    tmp_res <- run_pca(m_batch, args)
    list(
      embedding = data.table(
        cell_id = colnames(m_batch),
        tmp_res$embedding,
        batch_id = b
      ),
      loadings = data.table(
        gene   = rownames(tmp_res$loadings),
        tmp_res$loadings,
        batch_id = b
      )
    )
  })
  
  # concatenate all embeddings and loadings into a single data.table
  embeddings_dt <- rbindlist(lapply(pca_res_ls, `[[`, "embedding"))
  loadings_dt   <- rbindlist(lapply(pca_res_ls, `[[`, "loadings"))
  list( # could also consider saving variances
    embeddings = embeddings_dt,
    loadings   = loadings_dt
  )
}

main <- function() {
  args <- parse_pca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5",
              "solver", "n_components", "random_seed", "per_batch")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix")
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))

  if (args$per_batch) {
    res <- run_pca_per_batch(m, args)

    out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_pcas_per_batch.tsv", args$name))
    fwrite(res$embeddings, out_embeddings_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_embeddings_tsv))

    out_loadings_tsv  <- file.path(args$output_dir, sprintf("%s_loadings_per_batch.tsv", args$name))
    fwrite(res$loadings, out_loadings_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_loadings_tsv))

  } else {
    res <- run_pca(m, args)
    cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
      nrow(res$embedding), ncol(res$embedding),
      nrow(res$loadings),  ncol(res$loadings)))

    out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
    fwrite(data.frame(cell_id = rownames(res$embedding), res$embedding), out_embeddings_tsv,
      sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_embeddings_tsv))

    out_loadings_tsv <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
    fwrite(data.frame(gene = rownames(res$loadings), res$loadings), out_loadings_tsv,
      sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  wrote: %s\n", out_loadings_tsv))
  }
}

if (sys.nframe() == 0L) {
  main()
}
