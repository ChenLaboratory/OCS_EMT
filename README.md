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

## License

GPL-3; see [`LICENSE`](LICENSE).
