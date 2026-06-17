#!/usr/bin/env Rscript
# Normalization (scrapper-backed) for omnibenchmark.
#
# Normalizes scRNA-seq data using scrapper::normalizeRnaCounts.se() on a
# cell-subsetted H5AD object

suppressPackageStartupMessages({
  library(HDF5Array)
  library(scrapper)
  library(SingleCellExperiment)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("FILT module")
p <- add_base_args(p)                   # --output_dir, --name
p <- add_stage_args(p, "three-normalize")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
#p <- add_argument(p, "--param", type = "integer", help = "number of PCs")
args <- parse_args(p)                        # argparser's own parser

print(args)

# logging
cat(sprintf("LOG: command line args\n----------------------------------\n"))
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
for (k in args[-1]) {
  cat(sprintf("  %s: %s\n", k, args[[k]]))
}
cat(sprintf("LOG: command line args\n----------------------------------\n"))

# TODO: throw error when args are not right

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# read cellids to subset on
cellids <- readLines(gzfile(args$filtered.cellids))
cat("length(cellids):", length(cellids), "\n")

# read H5AD into SCE
sce <- read_h5ad(args$input_h5, as = "SingleCellExperiment")
sce <- sce[,cellids]
sce <- normalizeRnaCounts.se(sce)
d <- Matrix(logcounts(sce), sparse = TRUE)

output_file <- file.path(args$output_dir, paste0(args$name, "_normalized.h5"))
writeTENxMatrix(d, output_file, group="matrix")
cat(sprintf("wrote: %s\n", output_file))
file.info(output_file)[,c("size", "ctime")]

