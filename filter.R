#!/usr/bin/env Rscript
# Filtering cells (scrapper-backed) for omnibenchmark.
#
# Does cell-level filtering of raw single-cell RNA-seq data,
# hopefully adapting to the 'Sample' column of the metadata

suppressPackageStartupMessages({
  library(HDF5Array)
  library(scrapper)
  library(SingleCellExperiment)
  library(DelayedArray)
  library(argparser)
  library(jsonlite)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("FILT module")
p <- add_base_args(p)                   # --output_dir, --name
p <- add_stage_args(p, "filter")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
#p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
args <- parse_args(p)                        # argparser's own parser

# logging
cat(sprintf("LOG: command line args\n----------------------------------\n"))
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
for (k in args[-1]) {
  cat(sprintf("  %s: %s\n", k, args[[k]]))
}
cat(sprintf("LOG: command line args\n----------------------------------\n"))

# TODO: throw error when args are not right

# create output dir
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# read input H5AD
sce <- read_h5ad(args$input_h5, as = "SingleCellExperiment")

# do sample-wise filtering scrapper-style
is.mito <- grepl("^[Mm][Tt]-", rownames(sce))
rna.qc.metrics <- computeRnaQcMetrics(assay(sce), 
                                      subsets = list(mt = is.mito))

rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics, block = args$batch_variable)
keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics, block = args$batch_variable)

# write selected cellids
output_file <- file.path(args$output_dir, paste0(args$name, "_cellids.txt.gz"))
writeLines(colnames(sce)[keep], gzfile(output_file))
file.info(output_file)[,c("size", "ctime")]

