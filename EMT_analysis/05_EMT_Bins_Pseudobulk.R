## 05_EMT_Bins_Pseudobulk.R

source("setup.R")
source("config.R")
library(org.Hs.eg.db)

pdf(file.path(OUTPUT_DIR, "plots", "05_pseudobulk.pdf"), height = 8, width = 10)
for (platform in c("snRNAseq", "Visium")) {
    so <- readRDS(file.path(OUTPUT_DIR, paste0(platform, "_merged.rds")))

## ---- EMT bins (Fig 3A-C) ---------------------------------------------
## snRNA-seq nuclei with < 3000 UMIs are removed. The bin cut-offs are the
## quintiles of the OCS EMT score in PH142, the sample with the widest spread.
if (platform == "snRNAseq") so <- so[, so$nCount_RNA > 3000]
cutoffs <- c(min(so$EMT), quantile(so$EMT[so$sample == "PH142"], probs = seq(0.2, 0.8, by = 0.2)), max(so$EMT))
write.csv(cutoffs, file.path(OUTPUT_DIR, "tables", paste0(platform, "_EMT_bin_cutoffs.csv")))
so$bin <- cut(so$EMT, breaks = cutoffs, include.lowest = TRUE, labels = names(col.bin))

tab <- table(so$sample, so$bin)
print(tab)
write.csv(as.data.frame.matrix(tab), file.path(OUTPUT_DIR, "tables", paste0(platform, "_EMT_bin_counts.csv")))

## Bin of each nucleus/spot, used by ../Figures (the published Fig 3C was
## made by importing this table into Loupe Browser)
write.csv(data.frame(sample = so$orig.ident, barcode = sub("^.*_", "", colnames(so)), bin = so$bin),
    file.path(OUTPUT_DIR, "tables", paste0(platform, "_barcode_EMT_bins.csv")), row.names = FALSE)

## ---- Pseudo-bulk samples ---------------------------------------------
## UMI counts summed by sample and bin. Pseudo-bulk samples with small
## library sizes are removed.
y <- Seurat2PB(so, sample = "orig.ident", cluster = "bin", assay = DefaultAssay(so))
rm(so)
min.lib.size <- if (platform == "snRNAseq") 5e5 else 3e6
y <- y[, y$samples$lib.size > min.lib.size]
y$samples$group <- factor(y$samples$cluster)
if (platform == "snRNAseq") {
    keep <- filterByExpr(y)
} else {
    keep <- filterByExpr(y, min.count = 100, min.total.count = 1000)
}
print(table(keep))
y <- y[keep, ]
y <- normLibSizes(y)
y$genes <- data.frame(SYMBOL = rownames(y),
    ENTREZID = AnnotationDbi::mapIds(org.Hs.eg.db, rownames(y), column = "ENTREZID", keytype = "SYMBOL"))
plotMDS(y, col = col.bin[as.character(y$samples$group)], main = platform)

## ---- Bin markers and KEGG (Fig 4, Supp Fig 4, Supp Tables 14-33) -----
## voom with sample as a blocking factor. Markers of a bin are up-regulated
## relative to each of the other four bins, with FDR < 0.05 from the F-test
## of the four contrasts.
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
v <- voomLmFit(y, design, block = y$samples$sample, plot = TRUE)
for (k in 1:5) {
    contr <- makeContrasts(contrasts = paste0("EMT", k, "-EMT", setdiff(1:5, k)), levels = design)
    fit <- eBayes(contrasts.fit(v, contr))
    tt <- topTable(fit, coef = 1:4, number = Inf)
    up <- rowSums(fit$coefficients[rownames(tt), ] > 0) == 4
    markers <- tt[up & tt$adj.P.Val < 0.05, ]
    cat(platform, "EMT", k, "markers:", nrow(markers), "\n")
    write.csv(markers, file.path(OUTPUT_DIR, "tables", paste0(platform, "_EMT", k, "_markers.csv")), row.names = FALSE)
    kegg <- topKEGG(kegga(markers$ENTREZID, species = "Hs"), number = Inf)
    write.csv(kegg, file.path(OUTPUT_DIR, "tables", paste0(platform, "_EMT", k, "_KEGG.csv")))
}
}
dev.off()

print(sessionInfo())
