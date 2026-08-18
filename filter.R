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
p <- add_argument(p, "--min_cells", type = "integer", default = 5,
                   help = "minimum number of cells a feature must be detected in to be kept")
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
sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")

# read input H5AD
props <- read_yaml(args$properties_info)
batch <- props$batch_variable

# do sample-wise filtering scrapper-style
is.mito <- grepl("^[Mm][Tt]-", rownames(sce))
rna.qc.metrics <- computeRnaQcMetrics(assay(sce), 
                                      subsets = list(mt = is.mito))

rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics, block = batch)
keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics, block = batch)

# do feature-wise filtering
keep_features <- rowSums(assay(sce[, keep], "counts") > 0) >= args$min_cells
cat(sprintf("LOG: keeping %d / %d features (min_cells = %d)\n",
            sum(keep_features), length(keep_features), args$min_cells))

# calculate size factors from the filtered raw counts
X_filtered <- assay(sce[keep_features, keep], "counts")
library_size <- colSums(X_filtered)
mean_library_size <- mean(library_size)

if (!is.finite(mean_library_size) || mean_library_size <= 0) {
  stop("mean library size must be finite and greater than zero")
}

size_factor <- library_size / mean_library_size
cat(sprintf("LOG: size factor range: [%g, %g]\n", min(size_factor), max(size_factor)))

# write selected cellids
output_file <- file.path(args$output_dir, paste0(args$name, "_cellids.txt.gz"))
writeLines(colnames(sce)[keep], gzfile(output_file))

# write selected feature ids
features_output_file <- file.path(args$output_dir, paste0(args$name, "_featureids.txt.gz"))
writeLines(rownames(sce)[keep_features], gzfile(features_output_file))

# write size factors
size_factor_output_file <- file.path(args$output_dir, paste0(args$name, "_size_factors.tsv"))
write.table(data.frame(cell_id = colnames(sce)[keep],
                       size_factor = as.numeric(size_factor)),
            file = size_factor_output_file, sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("LOG: cellids output file info\n")
print(file.info(output_file)[,c("size", "ctime")])

cat("LOG: featureids output file info\n")
print(file.info(features_output_file)[,c("size", "ctime")])
