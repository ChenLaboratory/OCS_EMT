## 01_snRNAseq_Per_Sample.R

source("setup.R")
source("config.R")

pdf(file.path(OUTPUT_DIR, "plots", "01_snRNAseq_per_sample.pdf"), height = 5, width = 15)
for (id in samples) {
    p <- sample_info[id, ]

## ---- Read in and QC (Supp Table 3) -----------------------------------
## Genes detected in fewer than 3 nuclei and nuclei with fewer than 200 genes
## are removed when the object is created.
so <- read_flex(id)
so[["percent_mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

par(mfrow = c(1, 3))
hist(so$nFeature_RNA, breaks = 100, col = "dodgerblue", xlab = "", main = paste(p$paper_name, "- # genes"))
abline(v = p$flex_min_genes, lwd = 2, lty = 2, col = "red")
hist(so$nCount_RNA, breaks = 100, col = "dodgerblue", xlab = "", main = "# UMIs")
abline(v = p$flex_max_umi, lwd = 2, lty = 2, col = "red")
hist(so$percent_mt, breaks = 100, col = "dodgerblue", xlab = "", main = "% mitochondrial UMIs")
abline(v = p$flex_max_mt, lwd = 2, lty = 2, col = "red")

keep.cell <- so$nFeature_RNA > p$flex_min_genes & so$nCount_RNA < p$flex_max_umi &
    so$percent_mt < p$flex_max_mt
print(table(keep.cell))
so <- so[, keep.cell]

## ---- Normalisation and clustering ------------------------------------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 1500, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, seed.use = 8888, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.2, verbose = FALSE)
print(DimPlot(so, label = TRUE) + ggtitle(p$paper_name))

saveRDS(so, file.path(OUTPUT_DIR, "snRNAseq", paste0(id, ".rds")))
}
dev.off()

print(sessionInfo())
