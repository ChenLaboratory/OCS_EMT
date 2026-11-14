## config.R

if (!exists("annotation_dir")) source("setup.R")

## ---- Directories -----------------------------------------------------
CELLRANGER_DIR  <- "Data/CellRanger"   # Cell Ranger output of the snRNA-seq (10x Flex) samples
SPACERANGER_DIR <- "Data/SpaceRanger"  # Space Ranger output of the Visium samples
OUTPUT_DIR      <- "output"

for (d in c("snRNAseq", "Visium", "signatures", "tables", "plots")) {
    dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
}

## ---- Sample sheet ----------------------------------------------------
sample_info <- read.delim(file.path(annotation_dir, "sample_info.tsv"))
rownames(sample_info) <- sample_info$sample_id
samples <- sample_info$sample_id

## The four tumours with both carcinoma and sarcoma (#233, #105, PH142, PH419),
## used to derive the OCS signatures
mixed_samples <- samples[sample_info$mixed]

## Pan-cancer signatures (Supp Tables 4-5)
pan_cancer <- list(
    Epithelial = readLines(file.path(annotation_dir, "pan_cancer_epithelial_signature.txt")),
    Mesenchymal = readLines(file.path(annotation_dir, "pan_cancer_mesenchymal_signature.txt")))

## ---- Helpers ---------------------------------------------------------
## Raw counts of a snRNA-seq sample (10x matrix directory or H5 file)
read_flex <- function(id) {
    f <- file.path(CELLRANGER_DIR, sample_info[id, "flex_counts"])
    counts <- if (dir.exists(f)) Read10X(f) else Read10X_h5(f)
    CreateSeuratObject(counts, project = id, min.cells = 3, min.features = 200)
}

## Raw counts and image of a Visium sample
read_visium <- function(id) {
    so <- Load10X_Spatial(file.path(SPACERANGER_DIR, sample_info[id, "visium_dir"]),
        filename = "filtered_feature_bc_matrix.h5", filter.matrix = FALSE)
    so$orig.ident <- id
    so
}

## Processed per-sample object written by 01 or 02
read_sample <- function(platform, id) readRDS(file.path(OUTPUT_DIR, platform, paste0(id, ".rds")))

## Signature score: mean log-normalised expression of the signature genes
signature_means <- function(data, sigs) {
    sapply(sigs, function(g) Matrix::colMeans(data[intersect(g, rownames(data)), ]))
}
