#!/usr/bin/env Rscript
# Normalization (scrapper-backed) for omnibenchmark.
#
# Normalizes scRNA-seq data using scrapper::normalizeRnaCounts.se() on a
# cell-subsetted H5AD object

suppressPackageStartupMessages({
  library(HDF5Array)
  library(scrapper)
  library(SingleCellExperiment)
  library(anndataR)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("NORM module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "NORM")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
#p <- add_argument(p, "--param", type = "integer", help = "number of PCs")
args <- parse_args(p)                      # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args))
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
cat(sprintf("----------------------------------\n"))

# TODO: throw error when args are not right

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# read cellids to subset on
cellids <- readLines(gzfile(args$filtered_cellids))
cat("length(cellids):", length(cellids), "\n")

# read featureids to subset on
featureids <- readLines(gzfile(args$filtered_featureids))
cat("length(featureids):", length(featureids), "\n")

# read H5AD into SCE
sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")
sce <- sce[featureids,cellids]
#sce <- normalizeRnaCounts.se(sce)
#d <- Matrix(logcounts(sce), sparse = TRUE)

# only really valid for 4.5.x
cnts <- counts(sce)
d <- normalizeCounts(cnts, centerSizeFactors(colSums(cnts)))

output_file <- file.path(args$output_dir, paste0(args$name, "_normalized.h5"))
writeTENxMatrix(d, output_file, group="matrix")
cat(sprintf("wrote: %s\n", output_file))
file.info(output_file)[,c("size", "ctime")]

