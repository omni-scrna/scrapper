# Read a CSR-stored graph from a *_neighbors.h5 file into the graph-list format
# expected by scrapper::clusterGraph(). The graph lives under `group` as scipy
# CSR triples (data/indices/indptr, 0-based) with cell barcodes at /cell_ids.
# The column indices stored in the H5 file are 0-based, as in Python/scipy CSR,
# so we add 1 to convert them to R 1-based vertex indices. The row
# indices are reconstructed from `indptr`.
#
# This function converts the CSR representation into an edge list:
# each non-zero matrix entry becomes a candidate edge from row vertex `from` to
# column vertex `to`, with its corresponding weight from `data`.
# We keep only entries with from < to. This removes duplicated edges and self-loops.
#
# The returned object is a list with:
#   graph$vertices = number of cells/vertices
#   graph$edges    = integer vector of vertex pairs
#   graph$weights  = numeric vector with one weight per edge
#   cell_ids       = original cell barcodes, kept separately because the scrapper
#                    graph stores only numeric vertex indices, not cell names.

read_neighbors <- function(path, group = "connectivities") {
  ids <- as.character(h5read(path, "cell_ids"))
  n <- length(ids)
  
  indptr <- h5read(path, paste0(group, "/indptr"))
  indices <- h5read(path, paste0(group, "/indices"))  # 0-based column indices
  data <- as.numeric(h5read(path, paste0(group, "/data")))

  from <- rep.int(seq_len(n), diff(indptr)) # Convert rows to 1-based R indices

  to <- as.integer(indices) + 1L # Convert columns from 0-based Python indices to 1-based R indices
  
  keep <- from < to # Keep only one copy of each undirected edge
  
  graph <- list(
    vertices = as.integer(n),
    edges = as.integer(c(rbind(from[keep], to[keep]))),
    weights = as.numeric(data[keep])
  )
  
  list(
    graph = graph,
    cell_ids = ids
  )
}