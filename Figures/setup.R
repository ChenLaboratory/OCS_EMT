## setup.R

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

annotation_dir <- "../annotation"  # path to annotation files

## ---- Colour palettes -------------------------------------------------
## Samples, in the order of the sample sheet
col.sample <- RColorBrewer::brewer.pal(12, "Set3")

## EMT bins: Fig 3A-B, and the Loupe Browser colours of Fig 3C
col.bin <- c(EMT1 = "#F8A27E", EMT2 = "#F8E07F", EMT3 = "#8DCA77", EMT4 = "#B8C4E4", EMT5 = "#CE9AC6")
col.bin.spatial <- c(EMT1 = "#E03410", EMT2 = "#E8D107", EMT3 = "#80DD0E", EMT4 = "#4A90E2", EMT5 = "#BD10E0")

## Platform names used in the paper
platform_name <- c(snRNAseq = "snRNA-seq", Visium = "ST")
