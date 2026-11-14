source("setup.R")
source("config.R")

p3A <- p3B <- list()
for (platform in c("snRNAseq", "Visium")) {
    so <- read_merged(platform)
    so <- so[, !is.na(so$bin)]  # snRNA-seq nuclei with > 3000 UMIs
    cutoffs <- read.csv(file.path(data_dir, "tables", paste0(platform, "_EMT_bin_cutoffs.csv")))$x

    ### Figure 3A: OCS EMT scores, with the EMT bin cut-offs as dashed lines
    p3A[[platform]] <- RidgePlot(so, features = "EMT", group.by = "sample", cols = col.sample) +
        geom_vline(xintercept = cutoffs[2:5], linetype = "dashed") +
        labs(x = "EMT score", y = NULL, title = platform_name[platform]) + NoLegend()
    save_source(data.frame(sample = so$sample, barcode = colnames(so), EMT = so$EMT, bin = so$bin),
        paste0("Fig3A_EMT_scores_", platform))

    ### Figure 3B: percentage of nuclei/spots in each EMT bin. Samples in the
    ### reverse order of Fig 3A. PH003mixed was removed from the downstream
    ### analysis of the snRNA-seq data (373 nuclei), so its bar is left empty.
    tab <- table(so$sample, so$bin)
    if (platform == "snRNAseq") tab["PH003mixed", ] <- NA
    prop <- as.data.frame(prop.table(tab, 1) * 100)
    colnames(prop) <- c("sample", "bin", "percent")
    prop$sample <- factor(prop$sample, levels = rev(levels(prop$sample)))
    save_source(prop, paste0("Fig3B_EMT_bin_percentages_", platform))
    p3B[[platform]] <- ggplot(prop, aes(sample, percent, fill = bin)) +
        geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.3, na.rm = TRUE) +
        scale_fill_manual(values = col.bin) +
        labs(x = NULL, y = paste("Percentage (%) of", c(snRNAseq = "nuclei", Visium = "barcoded spots")[platform]),
            fill = NULL, title = platform_name[platform]) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    rm(so)
}
ggsave(file.path(out_dir, "Fig3A_EMT_score_ridges.pdf"), patchwork::wrap_plots(p3A, ncol = 2), height = 6, width = 11)
ggsave(file.path(out_dir, "Fig3B_EMT_bin_percentages.pdf"),
    patchwork::wrap_plots(p3B, ncol = 2, guides = "collect"), height = 5, width = 11)

### Figure 3C: EMT bins of the barcoded spots on the H&E image of each ST section
bins <- read.csv(file.path(data_dir, "tables", "Visium_barcode_EMT_bins.csv"))
p3C <- lapply(samples, function(id) {
    so <- readRDS(file.path(data_dir, "Visium", paste0(id, ".rds")))
    b <- bins[bins$sample == id, ]
    so$bin <- factor(b$bin[match(colnames(so), b$barcode)], levels = names(col.bin.spatial))
    suppressMessages(SpatialDimPlot(so, group.by = "bin", stroke = NA, combine = FALSE)[[1]] +
        scale_fill_manual(values = col.bin.spatial, drop = FALSE) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        labs(fill = NULL, title = sample_info[id, "paper_name"]))
})
## One legend, from PH142, which has all five bins
keep <- which(sample_info$paper_name == "PH142")
p3C[-keep] <- lapply(p3C[-keep], "+", NoLegend())
ggsave(file.path(out_dir, "Fig3C_spatial_EMT_bins.pdf"),
    patchwork::wrap_plots(p3C, ncol = 4, guides = "collect"), height = 13, width = 17)

sessionInfo()
