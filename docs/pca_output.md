# PCA output format (version 1)

File: `{output_dir}/{name}_{solver}_n_{n_components}.h5`

Flat HDF5 — intentionally not h5ad — so R (`rhdf5` / `HDF5Array`) and Julia (`HDF5.jl`) can read it without an anndata dependency. The schema is shared with the scanpy-backed module; files produced by either tool are interchangeable.

## Datasets

| Path | dtype | Shape | Description |
|---|---|---|---|
| `/embedding` | float64 | `(n_cells, n_components)` | PCA scores, cells-as-rows |
| `/loadings` | float64 | `(n_genes, n_components)` | Gene loadings (rotation / components), genes-as-rows |
| `/variance` | float64 | `(n_components,)` | Variance explained by each PC |
| `/variance_ratio` | float64 | `(n_components,)` | Fraction of total variance explained per PC |
| `/cell_ids` | UTF-8 string | `(n_cells,)` | Cell barcodes, aligned with rows of `/embedding` |
| `/gene_ids` | UTF-8 string | `(n_genes,)` | Gene identifiers (selected/HVG subset), aligned with rows of `/loadings` |

**Note for R consumers:** `/embedding` is cells-as-rows, matching `SingleCellExperiment::reducedDim` orientation. No transpose needed in most cases; use `t()` only if your downstream expects cells-as-columns.

## Root attributes

| Attribute | Type | Description |
|---|---|---|
| `format_version` | string | Always `"1"` for this spec |
| `tool` | string | Producing tool: `scanpy` or `scrapper` |
| `tool_version` | string | Version of the producing tool |
| `solver` | string | PCA solver: `arpack`/`randomized` (scanpy) or `irlba` (scrapper) |
| `n_components` | int | Number of PCs computed |
| `random_seed` | int | Random seed passed to the solver |

## Preprocessing

`scrapper::runPca` operates on the gene-by-cell matrix directly and centers/scales internally. The input matrix is subsetted to the selected genes before being passed to the solver.

## Validation

```bash
Rscript validators/pca_output.R path/to/name_irlba_n_K.h5
```

Exit codes: `0` = valid, `1` = validation failure, `2` = IO / usage error.

## Reading in R

```r
library(rhdf5)

emb <- h5read("name_irlba_n_K.h5", "embedding")     # n_cells × n_components
lod <- h5read("name_irlba_n_K.h5", "loadings")      # n_genes × n_components
var <- h5read("name_irlba_n_K.h5", "variance")
vr  <- h5read("name_irlba_n_K.h5", "variance_ratio")
cells <- h5read("name_irlba_n_K.h5", "cell_ids")
genes <- h5read("name_irlba_n_K.h5", "gene_ids")
```

## Reading in Python

```python
import h5py

with h5py.File("name_irlba_n_K.h5", "r") as h5:
    embedding      = h5["embedding"][:]       # (n_cells, n_components)
    loadings       = h5["loadings"][:]        # (n_genes, n_components)
    variance       = h5["variance"][:]
    variance_ratio = h5["variance_ratio"][:]
    cell_ids       = h5["cell_ids"].asstr()[:]
    gene_ids       = h5["gene_ids"].asstr()[:]
    tool           = h5.attrs["tool"]
    solver         = h5.attrs["solver"]
    n_components   = int(h5.attrs["n_components"])
```
