# PPGLomics

An interactive web application for exploring pheochromocytoma and paraganglioma (PPGL) transcriptomic data.

## Overview

PPGLomics provides a user-friendly interface for querying and visualizing gene expression data from two publicly available PPGL cohorts:

- **TCGA-PCPG** (n=160): Multiple driver genotypes including SDHx, VHL, RET, EPAS1, and others
- **A5 Consortium** (n=91): SDHB germline mutation carriers with detailed clinical annotation

## Features

- Single-gene expression analysis by genotype, cluster, or clinical variables
- Hierarchical clustering heatmaps
- Pairwise gene correlation analysis
- Genome-wide correlation search
- Differential expression with interactive volcano plots
- Kaplan-Meier survival analysis (A5 dataset)
- OncoPrint visualization of driver alterations

## Run Locally

```r
# Install dependencies
install.packages(c("shiny", "ggplot2", "plotly", "survival", "survminer"))
BiocManager::install("ComplexHeatmap")

# Clone and run
shiny::runApp()
```

## Data Sources

- Fishbein L, et al. (2017). Comprehensive molecular characterization of pheochromocytoma and paraganglioma. *Cancer Cell*. [PMID: 28162975](https://pubmed.ncbi.nlm.nih.gov/28162975/)
- Flynn A, et al. (2021). The genomic landscape of phaeochromocytoma. *Nature Communications*. [PMID: 33479239](https://pubmed.ncbi.nlm.nih.gov/33479239/)

## License

MIT

## Contact

Hussam Alkaissi, MD MS
National Institutes of Health
email: hussam.alkaissi@nih.gov
