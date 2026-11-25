# Hybrid epithelial/mesenchymal states in ovarian carcinosarcoma

This repository contains the code for the single-nucleus RNA-seq and spatial transcriptomics analyses in the following paper:<br>
<br>
**Multi-omics analysis of distinct hybrid epithelial/mesenchymal states in ovarian carcinosarcoma infers functional diversity and therapeutic vulnerabilities**<br>
<br>
Anthony Hadla, Xueyi Dong, Andrew Farrell, Casey J. A. Anttila, Ruvimbo D. Mishi, Ling Ling, Daniela Amann-Zalcenstein, Jiayu Wang, Vineet Vaibhav, Jumana Yousef, Laura F. Gianni, Gayanie Ratnayake, Cassandra J. Vandenberg, Rory Bowden&dagger;, Yunshun Chen&dagger;, Clare L. Scott&dagger; and Holly E. Barker&dagger;<br>

&dagger;These authors jointly supervised this work.

## Analyses

**Detailed setup and analysis instructions are provided in the respective subfolders:**
1) [`EMT_analysis/`](./EMT_analysis): from the Cell Ranger and Space Ranger output to the processed objects and result tables (per-sample QC, OCS-specific epithelial and mesenchymal signatures, merging of the samples, EMT bins, pseudo-bulk analysis of the bin markers).
2) [`Figures/`](./Figures): figure panels regenerated from the output of `EMT_analysis/`, and Fig 4 and Supp Fig 4 from the published KEGG results.

[`annotation/`](./annotation) contains the sample sheet with the per-sample QC thresholds (Supp Table 3) and the pan-cancer epithelial and mesenchymal signatures (Supp Tables 4-5).

## System requirements

The analyses used R 4.4.1 and the package versions listed below. For operating-system support and installation prerequisites, refer to the [R installation documentation](https://cran.r-project.org/doc/manuals/r-release/R-admin.html) and the documentation for the required packages.

Operating system used for the analyses: **Red Hat Enterprise Linux 9.3 (Plow)**.

## Software versions used in the manuscript analyses

| Software/package | Version |
|---|---|
| R | 4.4.1 |
| DropletUtils | 1.26.0 |
| Seurat | 5.3.0 |
| edgeR | 4.5.8 |
| AUCell | 1.30.1 |
| limma | 3.64.3 |
| hdf5r | 1.3.11 |
| ggplot2 | 4.0.0 |
| patchwork | 1.3.2 |
| RColorBrewer | 1.1-3 |
| org.Hs.eg.db | 3.21.0 |
| Matrix | 1.7-3 |
| AnnotationDbi | 1.70.0 |
| BiocManager | 1.30.26 |

Installation time for the required packages and their dependencies varies with the operating system, hardware and network speed. Installation from precompiled binaries is estimated to take approximately 10–30 minutes.

## License

GPL-3; see [`LICENSE`](LICENSE).
