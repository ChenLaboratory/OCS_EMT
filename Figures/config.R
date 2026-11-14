## config.R

if (!exists("annotation_dir")) source("setup.R")

## ---- Directories -----------------------------------------------------
data_dir      <- "../EMT_analysis/output"  # output of the analysis scripts
published_dir <- "published_results"       # KEGG results of the published analysis
out_dir       <- "output"
dir.create(file.path(out_dir, "source_data"), showWarnings = FALSE, recursive = TRUE)

## ---- Sample sheet ----------------------------------------------------
sample_info <- read.delim(file.path(annotation_dir, "sample_info.tsv"))
rownames(sample_info) <- sample_info$sample_id
samples <- sample_info$sample_id
names(col.sample) <- sample_info$paper_name

## ---- Helpers ---------------------------------------------------------
## Merged object of a platform (04) with the EMT bin of each nucleus/spot (05).
## snRNA-seq nuclei with < 3000 UMIs were not binned (bin = NA).
read_merged <- function(platform) {
    so <- readRDS(file.path(data_dir, paste0(platform, "_merged.rds")))
    bins <- read.csv(file.path(data_dir, "tables", paste0(platform, "_barcode_EMT_bins.csv")))
    so$bin <- factor(bins$bin[match(colnames(so), paste(bins$sample, bins$barcode, sep = "_"))],
        levels = names(col.bin))
    so
}

## Write the numbers plotted in a panel to out_dir/source_data/<name>.csv
save_source <- function(x, name, row.names = FALSE)
    write.csv(x, file.path(out_dir, "source_data", paste0(name, ".csv")), row.names = row.names)
