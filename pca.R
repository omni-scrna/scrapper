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
  library(rhdf5)
  library(Matrix)
  library(scrapper)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))

OUTPUT_FORMAT_VERSION <- "1"
TOOL <- "scrapper"


load_selected_genes <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con)
  trimws(lines[nzchar(trimws(lines))])
}


load_subset_matrix <- function(h5_path, selected_genes) {
  # Load TENx HDF5 (genes x cells) and subset to selected_genes.
  # Returns list(X = dgCMatrix genes-as-rows, gene_ids, cell_ids).
  data    <- as.numeric(h5read(h5_path, "matrix/data"))
  indices <- as.integer(h5read(h5_path, "matrix/indices"))
  indptr  <- as.integer(h5read(h5_path, "matrix/indptr"))
  shape   <- as.integer(h5read(h5_path, "matrix/shape"))  # c(n_genes, n_cells)

  gene_ids <- tryCatch(
    as.character(h5read(h5_path, "matrix/features/id")),
    error = function(e) tryCatch(
      as.character(h5read(h5_path, "matrix/genes")),
      error = function(e) sprintf("gene_%d", seq_len(shape[1]) - 1L)
    )
  )
  cell_ids <- tryCatch(
    as.character(h5read(h5_path, "matrix/barcodes")),
    error = function(e) sprintf("cell_%d", seq_len(shape[2]) - 1L)
  )

  # TENx HDF5 stores data in CSC order (indptr = column pointers),
  # i.e. columns are cells. dgCMatrix is also CSC; the matrix below is
  # therefore genes-as-rows, cells-as-columns.
  m <- sparseMatrix(
    i = indices + 1L, p = indptr, x = data,
    dims = shape, index1 = TRUE, repr = "C"
  )

  missing <- setdiff(selected_genes, gene_ids)
  if (length(missing) > 0) {
    sample <- head(sort(missing), 10)
    stop(sprintf(
      "%d selected gene(s) not present in normalized.h5; first %d: %s",
      length(missing), length(sample), paste(sample, collapse = ", ")
    ))
  }
  mask <- gene_ids %in% selected_genes
  list(X = m[mask, , drop = FALSE],
       gene_ids = gene_ids[mask],
       cell_ids = cell_ids)
}


run_pca <- function(X, gene_ids, cell_ids, args) {
  # X: gene-by-cell sparse matrix (rows = genes).
  set.seed(args$random_seed)
  pca <- runPca(X, number = args$n_components, num.threads = 1L)

  # scrapper::runPca returns components as (n_components × n_cells)
  # and rotation as (n_genes × n_components).
  embedding <- t(pca$components)              # (n_cells, n_components)
  loadings  <- pca$rotation                   # (n_genes, n_components)
  variance  <- as.numeric(pca$variance.explained)
  total_var <- if (!is.null(pca$total.variance)) {
    as.numeric(pca$total.variance)
  } else {
    sum(variance)
  }
  variance_ratio <- variance / total_var

  list(
    embedding = matrix(as.double(embedding), nrow = nrow(embedding)),
    loadings  = matrix(as.double(loadings),  nrow = nrow(loadings)),
    variance  = as.double(variance),
    variance_ratio = as.double(variance_ratio)
  )
}


write_output <- function(path, res, cell_ids, gene_ids, args) {
  if (file.exists(path)) file.remove(path)
  h5createFile(path)

  h5write(res$embedding,       path, "embedding")
  h5write(res$loadings,        path, "loadings")
  h5write(res$variance,        path, "variance")
  h5write(res$variance_ratio,  path, "variance_ratio")
  h5write(as.character(cell_ids), path, "cell_ids")
  h5write(as.character(gene_ids), path, "gene_ids")

  fid <- H5Fopen(path)
  on.exit(H5Fclose(fid), add = TRUE)

  scrapper_version <- as.character(packageVersion("scrapper"))
  h5writeAttribute(OUTPUT_FORMAT_VERSION, fid, "format_version")
  h5writeAttribute(TOOL,                  fid, "tool")
  h5writeAttribute(scrapper_version,      fid, "tool_version")
  h5writeAttribute(args$solver,           fid, "solver")
  h5writeAttribute(as.integer(args$n_components), fid, "n_components")
  h5writeAttribute(as.integer(args$random_seed),  fid, "random_seed")
}


main <- function() {
  args <- parse_pca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "normalized_h5", "selected_genes",
              "solver", "n_components", "random_seed")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  selected <- load_selected_genes(args$selected_genes)
  cat(sprintf("  selected genes: %d\n", length(selected)))

  loaded <- load_subset_matrix(args$normalized_h5, selected)
  cat(sprintf("  matrix (genes x cells): %d x %d\n",
              nrow(loaded$X), ncol(loaded$X)))

  res <- run_pca(loaded$X, loaded$gene_ids, loaded$cell_ids, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(res$embedding), ncol(res$embedding),
              nrow(res$loadings),  ncol(res$loadings)))

  out <- file.path(
    args$output_dir,
    sprintf("%s_%s_n_%d.h5", args$name, args$solver, args$n_components)
  )
  write_output(out, res, loaded$cell_ids, loaded$gene_ids, args)
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
