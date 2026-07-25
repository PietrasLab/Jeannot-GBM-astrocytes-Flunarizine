# External data

Data in this folder is **not** produced by this project and is **not** version-controlled
here. It is downloaded from a public repository at analysis time.

## Murine glioblastoma 10x Flex scRNA-seq (Lindgren / Rosberg et al.)

Used by [`scripts/external-data-figures.R`](../../scripts/external-data-figures.R) to compute
in vivo TGF-beta hallmark activity in reactive astrocytes, as independent support for the
in vitro TGF-beta results in this manuscript.

- **Record:** https://zenodo.org/records/17922854
- **DOI:** [10.5281/zenodo.17922854](https://doi.org/10.5281/zenodo.17922854)
- **Title:** Processed single-cell and spatial transcriptomics Seurat objects for murine glioblastoma (v1.0, 2025-12-13)
- **Licence:** CC-BY-4.0
- **Required file:** `seurat_flex_filtered_v1.0.rds` (520 MB)

Only that one file is needed; the rest of the record (raw Flex object, Visium and Visium HD
objects) is not used here.

### How to get it

Download `seurat_flex_filtered_v1.0.rds` from the record and place it in this directory:

```sh
curl -L -o data/external/seurat_flex_filtered_v1.0.rds \
  "https://zenodo.org/records/17922854/files/seurat_flex_filtered_v1.0.rds?download=1"
```

Alternatively, if you already hold the file elsewhere, point the script at it instead of
copying:

```sh
export FLEX_SEURAT_RDS=/path/to/seurat_flex_filtered_v1.0.rds
Rscript scripts/external-data-figures.R
```
