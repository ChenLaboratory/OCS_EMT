## 00_Install_Packages.R

# CRAN
install.packages(c(
    "Seurat",
    "hdf5r",
    "ggplot2",
    "patchwork",
    "RColorBrewer"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "edgeR",
    "limma",
    "AUCell",
    "org.Hs.eg.db"
))
