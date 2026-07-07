#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(scrapper)
})

# read_neighbors(): load a CSR neighbors graph from a *_neighbors.h5 file.
source("src/read_neighbors.R")

# arg parsing
source("src/common/cli.R")
p <- arg_parser("CLUST module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "CLUST")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--flavor", type = "character", help = "Clustering algorithm (multilevel|leiden_modularity|leiden_cpm)")
p <- add_argument(p, "--resolution", type = "numeric", help = "Clustering resolution")
p <- add_argument(p, "--random_seed", type = "integer", help = "Random seed")
args <- parse_args(p)                      # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


# Reproducibility
set.seed(args$random_seed)


# Load neighbors graph
neighbors <- read_neighbors(args$neighbors_h5, group = "connectivities")

neighbors_graph <- neighbors$graph
cell_id <- neighbors$cell_ids


# Run clustering
if (args$flavor == "multilevel") {
  
  clustering_result <- scrapper::clusterGraph(
    neighbors_graph,
    method = "multilevel",
    multilevel.resolution = args$resolution,
    seed = args$random_seed
  )
  
} else if (args$flavor == "leiden_modularity") {
  
  clustering_result <- scrapper::clusterGraph(
    neighbors_graph,
    method = "leiden",
    leiden.resolution = args$resolution,
    leiden.objective = "modularity",
    seed = args$random_seed
  )
  
}  else if (args$flavor == "leiden_cpm") {
    
    clustering_result <- scrapper::clusterGraph(
      neighbors_graph,
      method = "leiden",
      leiden.resolution = args$resolution,
      leiden.objective = "cpm",
      seed = args$random_seed
    )
  
}  else {
  stop("Unknown --flavor: ", args$flavor)
}


# Extract clusters matrix
m_clusters <- data.frame(
  cell_id = cell_id,
  cluster = as.character(clustering_result$membership)
)


cat("Cluster matrix dimensions:\n")
print(dim(m_clusters))
cat("\n")


# Save cluster matrix as .tsv
output_file <- file.path(
  args$output_dir,
  paste0(args$name, "_clusters.tsv")
)

cat("Writing output to:\n")
cat(output_file, "\n\n")

write.table(
  m_clusters,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(file.info(output_file)[, c("size", "ctime")])
