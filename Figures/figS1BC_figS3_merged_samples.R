source("setup.R")
source("config.R")

figS1 <- c(snRNAseq = "FigS1B", Visium = "FigS1C")
figS3 <- c(snRNAseq = "FigS3A", Visium = "FigS3B")
for (platform in c("snRNAseq", "Visium")) {
    so <- readRDS(file.path(data_dir, paste0(platform, "_merged.rds")))
    df <- data.frame(sample = so$sample, barcode = colnames(so), Embeddings(so, "umap"),
        pan_Epithelial = so$pan_Epithelial, pan_Mesenchymal = so$pan_Mesenchymal)
    rm(so)

    ### Supplementary Figure 1B (snRNA-seq) and 1C (ST): UMAP of the merged samples
    p <- ggplot(df, aes(umap_1, umap_2, colour = sample)) + geom_point(size = 0.1) +
        scale_colour_manual(values = col.sample) + coord_fixed() +
        labs(x = "UMAP_1", y = "UMAP_2", colour = NULL, title = platform_name[platform]) + theme_classic() +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    ggsave(file.path(out_dir, paste0(figS1[platform], "_UMAP_", platform, ".pdf")), p, height = 7, width = 9)
    save_source(df[, 1:4], paste0(figS1[platform], "_UMAP_", platform))

    ### Supplementary Figure 3A (snRNA-seq) and 3B (ST): pan-cancer signature scores
    p <- ggplot(df, aes(pan_Epithelial, pan_Mesenchymal, colour = sample)) + geom_point(size = 0.1) +
        scale_colour_manual(values = col.sample) + facet_wrap(~ sample, ncol = 4) +
        labs(x = "Epithelial score", y = "Mesenchymal score") + theme_bw() + theme(legend.position = "none")
    ggsave(file.path(out_dir, paste0(figS3[platform], "_pan_cancer_scores_", platform, ".pdf")), p, height = 7, width = 7)
    save_source(df[, c(1:2, 5:6)], paste0(figS3[platform], "_pan_cancer_scores_", platform))
}

sessionInfo()
