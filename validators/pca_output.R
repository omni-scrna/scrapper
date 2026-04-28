#!/usr/bin/env Rscript
# Validate a PCA output HDF5 file against the format-version 1 spec.
#
# Usage
# -----
#     Rscript validators/pca_output.R path/to/name_solver_n_K.h5
#
# Exit codes: 0 = valid, 1 = validation failure, 2 = file not found / IO error.

suppressPackageStartupMessages(library(rhdf5))


REQUIRED_DATASETS <- list(
  embedding      = list(dtype = "double", ndim = 2L),
  loadings       = list(dtype = "double", ndim = 2L),
  variance       = list(dtype = "double", ndim = 1L),
  variance_ratio = list(dtype = "double", ndim = 1L),
  cell_ids       = list(dtype = NULL,     ndim = 1L),
  gene_ids       = list(dtype = NULL,     ndim = 1L)
)

REQUIRED_ATTRS <- c("format_version", "tool", "tool_version", "solver",
                    "n_components", "random_seed")

VALID_SOLVERS_BY_TOOL <- list(
  scanpy   = c("arpack", "randomized"),
  scrapper = c("irlba")
)


read_attrs <- function(path) {
  fid <- H5Fopen(path)
  on.exit(H5Fclose(fid))
  attrs <- list()
  n <- H5Aget_num_attrs(fid)
  if (n > 0) {
    for (i in seq_len(n) - 1L) {
      aid <- H5Aopen_by_idx(fid, ".", n = i)
      name <- H5Aget_name(aid)
      attrs[[name]] <- H5Aread(aid)
      H5Aclose(aid)
    }
  }
  attrs
}


list_dataset_names <- function(path) {
  ls <- h5ls(path, recursive = FALSE)
  ls$name[ls$otype == "H5I_DATASET"]
}


dataset_info <- function(path, name) {
  ls <- h5ls(path, recursive = FALSE)
  row <- ls[ls$name == name & ls$otype == "H5I_DATASET", ]
  if (nrow(row) == 0) return(NULL)
  dims <- as.integer(strsplit(row$dim[1], " x ")[[1]])
  list(dtype = row$dclass[1], ndim = length(dims), dims = dims)
}


validate <- function(path) {
  errors <- character(0)

  if (!file.exists(path)) {
    return(sprintf("file not found: %s", path))
  }

  attrs <- tryCatch(read_attrs(path),
                    error = function(e) {
                      return(NULL)
                    })
  if (is.null(attrs)) {
    return(sprintf("cannot open HDF5 file: %s", path))
  }

  # root attributes
  missing_attrs <- setdiff(REQUIRED_ATTRS, names(attrs))
  for (a in sort(missing_attrs)) {
    errors <- c(errors, sprintf("missing root attribute: %s", a))
  }

  if (!is.null(attrs$format_version) && as.character(attrs$format_version) != "1") {
    errors <- c(errors, sprintf(
      "unexpected format_version: %s (expected '1')",
      attrs$format_version
    ))
  }

  tool <- attrs$tool
  solver <- attrs$solver
  if (!is.null(tool) && !(tool %in% names(VALID_SOLVERS_BY_TOOL))) {
    errors <- c(errors, sprintf(
      "unknown tool: '%s' (valid: %s)",
      tool, paste(sort(names(VALID_SOLVERS_BY_TOOL)), collapse = ", ")
    ))
  } else if (!is.null(tool) && !is.null(solver) &&
             !(solver %in% VALID_SOLVERS_BY_TOOL[[tool]])) {
    errors <- c(errors, sprintf(
      "solver '%s' not valid for tool '%s' (valid: %s)",
      solver, tool, paste(sort(VALID_SOLVERS_BY_TOOL[[tool]]), collapse = ", ")
    ))
  }

  # datasets
  for (name in names(REQUIRED_DATASETS)) {
    spec <- REQUIRED_DATASETS[[name]]
    info <- dataset_info(path, name)
    if (is.null(info)) {
      errors <- c(errors, sprintf("missing dataset: /%s", name))
      next
    }
    # rhdf5's h5ls reports class names like "FLOAT", "INTEGER", "STRING"
    if (!is.null(spec$dtype)) {
      if (spec$dtype == "double" && !(info$dtype %in% c("FLOAT"))) {
        errors <- c(errors, sprintf(
          "/%s: expected float dtype, got %s", name, info$dtype
        ))
      }
    }
    if (info$ndim != spec$ndim) {
      errors <- c(errors, sprintf(
        "/%s: expected %dD array, got %dD", name, spec$ndim, info$ndim
      ))
    }
  }

  # cross-dataset shape consistency
  shapes <- lapply(names(REQUIRED_DATASETS), function(n) dataset_info(path, n))
  names(shapes) <- names(REQUIRED_DATASETS)
  if (!any(vapply(shapes, is.null, logical(1)))) {
    # h5ls reports dims in HDF5 order; rhdf5 transposes on read so file
    # shape "n_cells x n_components" appears as "n_components x n_cells".
    # We check via expected file-order consistency.
    emb_dims <- shapes$embedding$dims     # (n_components, n_cells) in file order
    lod_dims <- shapes$loadings$dims      # (n_components, n_genes)
    n_comps <- emb_dims[1]
    n_cells <- emb_dims[2]
    n_genes <- lod_dims[2]

    if (lod_dims[1] != n_comps) {
      errors <- c(errors, sprintf(
        "/loadings first dim %d != n_components %d from /embedding",
        lod_dims[1], n_comps
      ))
    }
    if (shapes$variance$dims[1] != n_comps) {
      errors <- c(errors, sprintf(
        "/variance length %d != n_components %d",
        shapes$variance$dims[1], n_comps
      ))
    }
    if (shapes$variance_ratio$dims[1] != n_comps) {
      errors <- c(errors, sprintf(
        "/variance_ratio length %d != n_components %d",
        shapes$variance_ratio$dims[1], n_comps
      ))
    }
    if (shapes$cell_ids$dims[1] != n_cells) {
      errors <- c(errors, sprintf(
        "/cell_ids length %d != n_cells %d",
        shapes$cell_ids$dims[1], n_cells
      ))
    }
    if (shapes$gene_ids$dims[1] != n_genes) {
      errors <- c(errors, sprintf(
        "/gene_ids length %d != n_genes %d",
        shapes$gene_ids$dims[1], n_genes
      ))
    }

    if (!is.null(attrs$n_components)) {
      stored <- as.integer(attrs$n_components)
      if (stored != n_comps) {
        errors <- c(errors, sprintf(
          "n_components attribute (%d) != actual embedding components (%d)",
          stored, n_comps
        ))
      }
    }
  }

  # value sanity checks
  vr <- tryCatch(h5read(path, "variance_ratio"), error = function(e) NULL)
  if (!is.null(vr)) {
    if (any(vr < 0) || any(vr > 1)) {
      errors <- c(errors, "/variance_ratio values must be in [0, 1]")
    }
    if (sum(vr) > 1.0 + 1e-6) {
      errors <- c(errors, sprintf(
        "/variance_ratio sums to %.6f, expected <= 1.0", sum(vr)
      ))
    }
  }
  v <- tryCatch(h5read(path, "variance"), error = function(e) NULL)
  if (!is.null(v) && any(v < 0)) {
    errors <- c(errors, "/variance contains negative values")
  }

  errors
}


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    cat("usage: validators/pca_output.R <pca_output.h5>\n", file = stderr())
    quit(status = 2)
  }
  path <- args[1]
  result <- validate(path)
  if (is.character(result) && length(result) == 1 &&
      grepl("^(file not found|cannot open)", result)) {
    cat(sprintf("INVALID: %s\n  - %s\n", path, result))
    quit(status = 2)
  }
  if (length(result) > 0) {
    cat(sprintf("INVALID: %s\n", path))
    for (e in result) cat(sprintf("  - %s\n", e))
    quit(status = 1)
  }
  cat(sprintf("OK: %s\n", path))
  attrs <- read_attrs(path)
  emb <- dataset_info(path, "embedding")$dims
  lod <- dataset_info(path, "loadings")$dims
  cat(sprintf("  cells=%d  genes=%d  components=%d\n",
              emb[2], lod[2], emb[1]))
  cat(sprintf("  tool=%s %s  solver=%s\n",
              attrs$tool, attrs$tool_version, attrs$solver))
  quit(status = 0)
}

if (sys.nframe() == 0L) {
  main()
}
