## 02_Visium_Per_Sample.R

source("setup.R")
source("config.R")

pdf(file.path(OUTPUT_DIR, "plots", "02_Visium_per_sample.pdf"), height = 5, width = 15)
for (id in samples) {
    p <- sample_info[id, ]

## ---- Read in and QC (Supp Table 3) -----------------------------------
so <- read_visium(id)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

par(mfrow = c(1, 3))
hist(so$nFeature_Spatial, breaks = 100, col = "dodgerblue", xlab = "", main = paste(p$paper_name, "- # genes"))
abline(v = p$visium_min_genes, lwd = 2, lty = 2, col = "red")
hist(so$nCount_Spatial, breaks = 100, col = "dodgerblue", xlab = "", main = "# UMIs")
abline(v = p$visium_min_umi, lwd = 2, lty = 2, col = "red")
hist(so$percent.mt, breaks = 100, col = "dodgerblue", xlab = "", main = "% mitochondrial UMIs")
abline(v = p$visium_max_mt, lwd = 2, lty = 2, col = "red")

## Spots with few UMIs are removed first, then genes detected in at most 3
## spots or with at most 10 UMIs, then spots with few genes or high mito
so <- so[, so$nCount_Spatial > p$visium_min_umi]
counts <- LayerData(so, layer = "counts")
so <- so[Matrix::rowSums(counts > 0) > 3 & Matrix::rowSums(counts) > 10, ]
keep.spot <- so$nFeature_Spatial > p$visium_min_genes & so$percent.mt < p$visium_max_mt
print(table(keep.spot))
so <- so[, keep.spot]

## ---- Normalisation and clustering ------------------------------------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 1000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
print(SpatialDimPlot(so, pt.size.factor = 1.4) + ggtitle(p$paper_name) + DimPlot(so, label = TRUE))

saveRDS(so, file.path(OUTPUT_DIR, "Visium", paste0(id, ".rds")))
}
dev.off()

print(sessionInfo())
