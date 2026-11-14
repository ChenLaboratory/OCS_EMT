# Figures

Scripts that regenerate the snRNA-seq and spatial transcriptomics (ST) figure panels of the paper from the output of [`../EMT_analysis`](../EMT_analysis), except Fig 4 and Supp Fig 4, which are drawn from the KEGG results of the published analysis in [`published_results/`](published_results).

## Scripts

| Script | Panels |
|---|---|
| `figS1BC_figS3_merged_samples.R` | Supp Fig 1B-C (UMAP of the merged samples), Supp Fig 3A-B (pan-cancer signature scores) |
| `fig3_EMT_bins.R` | Fig 3A (OCS EMT scores and bin cut-offs), Fig 3B (percentage of nuclei/spots in each EMT bin), Fig 3C (EMT bins on the H&E images) |
| `fig4_figS4_KEGG.R` | Fig 4A-E (snRNA-seq), Supp Fig 4A-E (ST): top 10 KEGG pathways of the markers of each EMT bin |

- `setup.R` loads the packages and defines the annotation directory and colour palettes;
- `config.R` defines paths, the sample sheet and helper functions;
- `published_results/` contains the KEGG results of the published analysis used for Fig 4 and Supp Fig 4.

## Running the scripts

Run the scripts in [`../EMT_analysis`](../EMT_analysis) first (not needed for `fig4_figS4_KEGG.R`), then

```r
setwd("path/to/this/repository/Figures")
source("figS1BC_figS3_merged_samples.R")
source("fig3_EMT_bins.R")
source("fig4_figS4_KEGG.R")
```
