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
  library(anndataR)
  library(yaml)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("FILT module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "FILT")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
#p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
args <- parse_args(p)                    # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

# TODO: throw error when args are not right

# create output dir
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# read input H5AD
sce <- read_h5ad(args$rawdata.h5ad, as = "SingleCellExperiment")

# read input H5AD
props <- read_yaml(args$properties.info)
batch <- props$batch_variable

# do sample-wise filtering scrapper-style
is.mito <- grepl("^[Mm][Tt]-", rownames(sce))
rna.qc.metrics <- computeRnaQcMetrics(assay(sce), 
                                      subsets = list(mt = is.mito))

rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics, block = batch)
keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics, block = batch)

# write selected cellids
output_file <- file.path(args$output_dir, paste0(args$name, "_cellids.txt.gz"))
writeLines(colnames(sce)[keep], gzfile(output_file))
file.info(output_file)[,c("size", "ctime")]

