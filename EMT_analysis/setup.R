## setup.R

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(edgeR)
})

annotation_dir <- "../annotation"  # path to annotation files

## EMT bin colours (as in Fig 3C)
col.bin <- c(EMT1 = "#E03410", EMT2 = "#E8D107", EMT3 = "#80DD0E", EMT4 = "#4A90E2", EMT5 = "#BD10E0")
