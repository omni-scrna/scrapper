#!/usr/bin/env Rscript
# Argument parser for omnibenchmark scrapper modules.
#
# Reusable across stage entrypoints (pca.R, and future umap.R, cluster.R, ...).
#
# Conventions:
# - All arguments are required. No defaults — callers (omnibenchmark configs)
#   must pass everything explicitly so runs are fully reproducible from the
#   invocation line.

suppressPackageStartupMessages(library(optparse))

build_pca_parser <- function() {
  option_list <- list(
    make_option("--output_dir", type = "character",
                help = "Output directory for results"),
    make_option("--name", type = "character",
                help = "Module name/identifier"),
    make_option("--normalized_selected.h5", type = "character",
                help = "TENx-format HDF5 of normalized expression (genes x cells)"),
    #make_option("--selected.genes", type = "character",
    #            help = "Gzipped text file of selected gene ids (one per line)"),
    make_option("--solver", type = "character",
                help = "PCA solver (scrapper: irlba, random, exact)"),
    make_option("--n_components", type = "integer",
                help = "Number of principal components to compute"),
    make_option("--random_seed", type = "integer",
                help = "Seed for randomized solvers (and for reproducibility)")
  )

  OptionParser(
    option_list = option_list,
    description = "OmniBenchmark PCA module (scrapper)"
  )
}

parse_pca_args <- function() {
  parser <- build_pca_parser()
  raw <- parse_args(parser)

  args <- list(
    output_dir      = raw$output_dir,
    name            = raw$name,
    input_h5        = raw[["normalized_selected.h5"]],
    #selected_genes  = raw[["selected.genes"]],
    solver          = raw$solver,
    n_components    = raw$n_components,
    random_seed     = raw$random_seed
  )

  required <- names(args)
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  valid_solvers <- c("irlba", "random", "exact")
  if (!(args$solver %in% valid_solvers)) {
    stop("Invalid --solver: ", args$solver,
         " (valid: ", paste(valid_solvers, collapse = ", "), ")")
  }

  args

