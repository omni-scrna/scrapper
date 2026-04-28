# Scrapper

Scrapper-backed (libscran) module for omnibenchmark scRNA pipelines.

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` loads all runtime libraries and prints `OK`. Run it after install to confirm the environment is healthy.

## Usage

### PCA

```sh
pixi run Rscript pca.R \
  --output_dir <dir> \
  --name <name> \
  --normalized.h5 <normalized.h5> \
  --selected.genes <selected.genes.gz> \
  --solver <irlba> \
  --n_components <int> \
  --random_seed <int>
```

Output: `<output_dir>/<name>_<solver>_n_<n_components>.h5` — see [`docs/pca_output.md`](docs/pca_output.md) for the full format spec. The schema is shared with the scanpy-backed module.

#### Validation

```sh
pixi run validate <output_dir>/<name>_<solver>_n_<n_components>.h5
```

Exit codes: `0` = valid, `1` = validation failure, `2` = IO / usage error.

## Conda environment export

```sh
pixi run export-env
```

Exports the resolved environment to `envs/scrapper.yml`. The environment is named after the repo root folder.

## Citation

If you use this module in your research, please cite it using the information in `CITATION.cff`.
