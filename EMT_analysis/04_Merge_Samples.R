## 04_Merge_Samples.R

source("setup.R")
source("config.R")
library(AUCell)

ocs <- list(
    Epithelial = readLines(file.path(OUTPUT_DIR, "signatures", "OCS_epithelial_signature.txt")),
    Mesenchymal = readLines(file.path(OUTPUT_DIR, "signatures", "OCS_mesenchymal_signature.txt")))

for (platform in c("snRNAseq", "Visium")) {

## ---- Merge without integration ---------------------------------------
## Raw counts of the nuclei/spots retained by the per-sample QC. For Visium,
## only genes present in all samples are kept.
so_list <- lapply(samples, function(id) {
    so <- if (platform == "snRNAseq") read_flex(id) else read_visium(id)
    so[, colnames(read_sample(platform, id))]
})
common <- Reduce(intersect, lapply(so_list, rownames))
if (platform == "Visium") so_list <- lapply(so_list, function(so) so[common, ])
so <- merge(so_list[[1]], so_list[-1], add.cell.ids = samples)
rm(so_list)
so$orig.ident <- factor(so$orig.ident, levels = samples)
so$sample <- factor(sample_info[as.character(so$orig.ident), "paper_name"], levels = sample_info$paper_name)
print(table(so$sample))

## ---- UMAP (Supp Fig 1B-C) --------------------------------------------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = if (platform == "snRNAseq") 1500 else 1000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
so <- JoinLayers(so)

## ---- Pan-cancer signature scores (Supp Fig 3) ------------------------
## Mean log-normalised expression over the signature genes present in all samples
score <- signature_means(LayerData(so, layer = "data")[common, ], pan_cancer)
so$pan_Epithelial <- score[, "Epithelial"]
so$pan_Mesenchymal <- score[, "Mesenchymal"]

## ---- OCS signature scores (AUCell) -----------------------------------
## AUCell breaks ties in the gene rankings at random, hence the seed.
## OCS EMT score = mesenchymal AUC - epithelial AUC.
set.seed(2025)
auc <- getAUC(AUCell_run(LayerData(so, layer = "counts"), ocs))
so$OCS_Epithelial <- auc["Epithelial", colnames(so)]
so$OCS_Mesenchymal <- auc["Mesenchymal", colnames(so)]
so$EMT <- so$OCS_Mesenchymal - so$OCS_Epithelial

saveRDS(so, file.path(OUTPUT_DIR, paste0(platform, "_merged.rds")))
}

print(sessionInfo())
