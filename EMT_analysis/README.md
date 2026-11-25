# snRNA-seq and Visium analysis

Analysis of the 12 ovarian carcinosarcoma (OCS) PDX tumours profiled by 10x Flex single-nucleus RNA-seq (snRNA-seq) and Visium CytAssist spatial transcriptomics.

## Workflow

Scripts, run in order:

| Script | Description | Paper |
|---|---|---|
| `00_Install_Packages.R` | R package installation | |
| `01_snRNAseq_Per_Sample.R` | Per-sample QC, normalisation and clustering of the snRNA-seq samples | Supp Table 3 |
| `02_Visium_Per_Sample.R` | Per-sample QC, normalisation and clustering of the Visium samples | Supp Table 3 |
| `03_OCS_Signatures.R` | Genes associated with the pan-cancer EMT score in the four mixed tumours (edgeR quasi-likelihood), OCS-specific epithelial and mesenchymal signatures | Fig 2G; Supp Tables 10-13 |
| `04_Merge_Samples.R` | Samples of each platform merged without integration, UMAP, pan-cancer signature scores, OCS signature scores (AUCell) | Supp Figs 1B-C, 3 |
| `05_EMT_Bins_Pseudobulk.R` | Five EMT bins from the OCS EMT score, pseudo-bulk DE of each bin versus the other four bins, KEGG analysis of the bin markers | Fig 3A-C, 4; Supp Fig 4; Supp Tables 14-33 |

- `setup.R` loads the packages and defines the annotation directory and colour palettes;
- `config.R` defines paths, the sample sheet, the pan-cancer signatures and helper functions.

The figure panels are made from the output of these scripts by the scripts in [`../Figures`](../Figures).

The sample sheet [`../annotation/sample_info.tsv`](../annotation/sample_info.tsv) gives, for each sample, the name used in the paper, the location of its Cell Ranger and Space Ranger output, and the QC thresholds.

## Running the pipeline

Set `CELLRANGER_DIR` and `SPACERANGER_DIR` in `config.R` to the Cell Ranger and Space Ranger output directories, then

```r
setwd("path/to/this/repository/EMT_analysis")

source("00_Install_Packages.R")
source("01_snRNAseq_Per_Sample.R")
source("02_Visium_Per_Sample.R")
source("03_OCS_Signatures.R")
source("04_Merge_Samples.R")
source("05_EMT_Bins_Pseudobulk.R")
```

The full R analysis pipeline and figure generation are estimated to take approximately 2–6 hours on a modern workstation with 64 GB RAM and SSD storage, excluding package installation and Cell Ranger/Space Ranger processing.

Notes:
- AUCell breaks ties in the gene rankings at random, so the AUCell scores, and hence the EMT bins, depend on the random seed set in `04_Merge_Samples.R`.
