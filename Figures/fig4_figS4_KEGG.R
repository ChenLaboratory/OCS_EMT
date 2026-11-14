source("setup.R")
source("config.R")

### Figure 4A-E (snRNA-seq) and Supplementary Figure 4A-E (ST): top 10 KEGG
### pathways of the markers of each EMT bin, ordered by the number of markers
### in the pathway. P-values are the unadjusted kegga P-values.
### The KEGG results of the published analysis are used, because kegga uses the
### current KEGG release. For the results of ../EMT_analysis, set
### kegg_dir <- file.path(data_dir, "tables").
kegg_dir <- published_dir
fig <- c(snRNAseq = "Fig4", Visium = "FigS4")
for (platform in names(fig)) {
    pdf(file.path(out_dir, paste0(fig[platform], "_KEGG_", platform, ".pdf")), height = 4, width = 8)
    for (k in 1:5) {
        kegg <- read.csv(file.path(kegg_dir, paste0(platform, "_EMT", k, "_KEGG.csv")), row.names = 1)
        top <- head(kegg, 10)  # topKEGG output is sorted by P.DE
        save_source(top, paste0(fig[platform], LETTERS[k], "_KEGG_EMT", k, "_", platform), row.names = TRUE)
        print(ggplot(top, aes(DE, reorder(Pathway, DE), size = DE, colour = -log10(P.DE))) +
            geom_point(alpha = 0.7) + scale_colour_gradient(low = "blue", high = "red") +
            guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
            labs(x = "Gene Count", y = "KEGG Pathway", size = "Gene Count", colour = "-log10(P-value)") +
            theme_minimal() + theme(axis.title = element_text(face = "bold")))
    }
    dev.off()
}

sessionInfo()
