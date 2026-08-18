# Counts loading for FILT/NORM, switchable between an in-memory read and an
# out-of-core view of the same h5ad.
#
# Why "delayed" exists: anndataR::read_h5ad() materialises the CSR indptr as an
# R integer vector. Past 2^31 nonzeros passes at 2.62e9 and that overflows to NA,
# and anndataR drops layers/counts with a *warning* rather than an error,
# handing back an object with no counts.
# H5SparseMatrix() reads the same group as a DelayedArray and never builds that
# vector, so it is the only path that reaches the top of a cell-count ladder.
#
# Cost: block reads instead of RAM. Expect out-of-core to be markedly slower;
# it is a separate operating point, not a free switch.

suppressPackageStartupMessages({
  library(HDF5Array)
  library(SingleCellExperiment)
  library(rhdf5)
})

read_counts <- function(path, backend = c("memory", "delayed")) {
  backend <- match.arg(backend)
  if (backend == "memory") {
    return(anndataR::read_h5ad(path, as = "SingleCellExperiment"))
  }
  # h5ad holds CSR(cells x genes), bit-identical to CSC(genes x cells), so
  # H5SparseMatrix reports genes x cells -- the same orientation the memory
  # path yields as an SCE. Callers need no special casing.
  m <- H5SparseMatrix(path, "layers/counts")
  dimnames(m) <- list(as.character(h5read(path, "var/_index")),
                      as.character(h5read(path, "obs/_index")))
  SingleCellExperiment(list(counts = m))
}
