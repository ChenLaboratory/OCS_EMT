## 03_OCS_Signatures.R

source("setup.R")
source("config.R")

## ---- Genes associated with the pan-cancer EMT score ------------------
## Pan-cancer EMT score = mean log-normalised expression of the mesenchymal
## signature minus that of the epithelial signature. Each gene is regressed on
## the score with edgeR quasi-likelihood, separately for each mixed tumour and
## platform.
tabs <- list()
for (platform in c("snRNAseq", "Visium")) {
    for (id in mixed_samples) {
        so <- read_sample(platform, id)
        score <- signature_means(LayerData(so, layer = "data"), pan_cancer)
        EMT <- score[, "Mesenchymal"] - score[, "Epithelial"]

        y <- DGEList(counts = LayerData(so, layer = "counts"))
        y <- normLibSizes(y)
        design <- model.matrix(~ EMT)
        fit <- glmQLFit(y, design)
        res <- glmQLFTest(fit, coef = 2)
        print(summary(decideTests(res)))

        k <- paste(platform, id, sep = "_")
        tabs[[k]] <- topTags(res, n = Inf)$table
        write.csv(tabs[[k]], file.path(OUTPUT_DIR, "signatures", paste0("EMT_genes_", k, ".csv")))
    }
}

## ---- OCS signatures (Supp Tables 10-11) ------------------------------
## Genes with FDR < 0.01 and the same sign in all eight tables, restricted to
## genes detected in all 12 snRNA-seq samples.
detected <- Reduce(intersect, lapply(samples, function(id) rownames(read_sample("snRNAseq", id))))
down <- lapply(tabs, function(x) rownames(x)[x$FDR < 0.01 & x$logFC < 0])
up <- lapply(tabs, function(x) rownames(x)[x$FDR < 0.01 & x$logFC > 0])
ocs <- list(Epithelial = intersect(Reduce(intersect, down), detected),
    Mesenchymal = intersect(Reduce(intersect, up), detected))
print(lengths(ocs))
writeLines(ocs$Epithelial, file.path(OUTPUT_DIR, "signatures", "OCS_epithelial_signature.txt"))
writeLines(ocs$Mesenchymal, file.path(OUTPUT_DIR, "signatures", "OCS_mesenchymal_signature.txt"))

## Genes shared with the pan-cancer signatures (Fig 2G, Supp Tables 12-13)
shared <- Map(intersect, ocs, pan_cancer)
print(lengths(shared))

print(sessionInfo())
