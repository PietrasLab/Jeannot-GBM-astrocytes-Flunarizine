# JEANNOT-GBM-ASTROCYTES-FLUNARIZINE

**Repository:** https://github.com/PietrasLab/JEANNOT-GBM-ASTROCYTES-FLUNARIZINE  
**Preprint:** Jeannot P. *et al.* bioRxiv (2025)  
https://doi.org/10.1101/2025.07.12.664538

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
<!-- Add Zenodo badge once you mint a DOI:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)
-->

---

## Overview
This repository contains analysis code and figures presented in the the study:

> Jeannot P., Rosberg R., Lindgren D., Dawson J.C., Pracucci E., Börjesson-Freitag C., Martinez J., Pantazopoulou V., Malmberg M., Smolag K.I., Manou D., Elliott R.J.R., Ceberg C., Berg T.J., Ahlenius H., Carragher N.O., Pietras A. (2025).  
> *Phenotypic Screening Identifies Flunarizine as an Inhibitor of Radiotherapy-Induced Astrocyte Reactivity with Therapeutic Potential in Glioblastoma.*  
> **bioRxiv**. https://doi.org/10.1101/2025.07.12.664538

---

## Data availability
- **Raw data**: ArrayExpress accession [E-MTAB-15215](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15215)  
- **Processed data objects**: DESeq2 `.dds` and `.rds` objects in `data/` generated from paired-end RNA-seq reads processed using the **nf-core/rnaseq** pipeline (v3.14.0). 
- **Results**: contrasts, PCA, heatmaps, and figure-ready CSV/PDF outputs in `results/`

---


## Repository contents
- `data/` — DESeq2 objects and metadata  
- `results/` — figures and figure data of figures in manuscript
- `scripts/` — analysis scripts to reproduice manuscript figures (e.g., `manuscript-figures.R`)  
- `LICENSE`, `CITATION.cff`, `README.md`  

---

Objects included in `data/`:
- `deseq2.dds.RData` — DESeq2 dataset object  
- `deseq2.filtered.rds` — filtered DESeq2 object  
- `fdata_full.rds`, `pdata.rds`, `SampleSheet-2024_169.csv` — metadata inputs  

---
