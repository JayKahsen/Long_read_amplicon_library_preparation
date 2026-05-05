# Long_read_amplicon_library_preparation

Reproducible R code and supporting data for the manuscript:

**A rapid and flexible method for long-read amplicon library preparation and balancing**

This repository contains the current doubleton-focused workflow used to rebuild the rarefied matrices and reproduce Figures 2-8 from the manuscript.

---

## Repository purpose

- Rebuild the ASV-centered matrix layer from the doubleton QIIME export
- Generate rarefied `Phylum`, `Genus`, and `ASV` matrices
- Rebuild reduced metadata and ASV reference tables for the retained samples and ASVs
- Reproduce manuscript figure outputs from the rarefied matrices

---

## Directory structure

```text
long_read_amplicon_balancing/
|- _globalStuff.R
|- helperJ.R
|- 1_QIIME Tables to working files.R
|- Figure_2_taxonomy.R
|- Figure_3_6_alpha_diversity_stats_ANOVA.R
|- Figure_4_7_PCA.R
|- Figure_5_pipeline.R
|- Figure_8_abundance.R
|- LICENSE
|- README.md
|- R_package_citations.txt
|
\- data_tables/
   |- ASV_sequence.csv
   |- ASV_sequence_all.zip
   |- ASV_taxa.csv
   |- ASV_taxa_all.zip
   |- matrix_names.csv
   |- meta.csv
   |- meta.xlsx
   |- meta_ALL.csv
   |- meta_ALL.xlsx
   |- sample_reads.csv
   |
   |- original/
   |  \- Doubleton.zip
   |
   |- raw_matrix/
   |  |- raw_matrix_ASV.csv
   |  |- raw_matrix_Genus.csv
   |  \- raw_matrix_Phylum.csv
   |
   \- rarefied_matrix/
      |- rarefied_matrix_ASV_Set1.csv
      |- rarefied_matrix_Genus_Set1.csv
      \- rarefied_matrix_Phylum_Set1.csv
```

---

## File descriptions

- `_globalStuff.R`: shared libraries, data loads, palettes, factor orders, and plotting wrappers used across scripts.
- `helperJ.R`: shared helper functions used by the run-first and figure scripts.
- `1_QIIME Tables to working files.R`: rebuilds the doubleton-only raw and rarefied matrices, trims metadata to retained samples, trims ASV reference files to retained ASVs, and rewrites matrix registry files. `Phylum` and `Genus` matrices are collapsed from the ASV matrix.
- `Figure_2_taxonomy.R`: creates the Figure 2 taxonomy profile panels from the rarefied `Phylum` matrix.
- `Figure_3_6_alpha_diversity_stats_ANOVA.R`: creates the Figure 3 and Figure 6 alpha-diversity outputs from the rarefied `ASV` matrix.
- `Figure_4_7_PCA.R`: creates the Figure 4 and Figure 7 PCA, PERMDISP, and PERMANOVA outputs from the rarefied `ASV` matrix.
- `Figure_5_pipeline.R`: creates the Figure 5 workflow and read-retention outputs from the rarefied `ASV` matrix.
- `Figure_8_abundance.R`: creates the Figure 8 abundance outputs from the rarefied `Phylum` matrix.
- `data_tables/ASV_taxa.csv`: reduced ASV taxonomy table containing one row per retained rarefied ASV.
- `data_tables/ASV_sequence.csv`: reduced ASV sequence table containing one row per retained rarefied ASV.
- `data_tables/ASV_taxa_all.zip`: zipped full ASV taxonomy source file used to rebuild `ASV_taxa.csv`.
- `data_tables/ASV_sequence_all.zip`: zipped full ASV sequence source file used to rebuild `ASV_sequence.csv`.
- `data_tables/meta.csv`: reduced metadata table containing only the samples retained in the rarefied matrices.
- `data_tables/meta.xlsx`: formatted inspection copy of `meta.csv`.
- `data_tables/meta_ALL.csv`: full metadata source table used by the run-first script.
- `data_tables/meta_ALL.xlsx`: Excel copy of the full metadata source table.
- `data_tables/matrix_names.csv`: matrix registry used by downstream figure scripts.
- `data_tables/sample_reads.csv`: per-sample raw and rarefied read counts.
- `data_tables/original/Doubleton.zip`: zipped QIIME-derived doubleton input tables used to rebuild the matrix layer.
- `data_tables/raw_matrix/*.csv`: raw matrices rebuilt from the doubleton input data.
- `data_tables/rarefied_matrix/*.csv`: rarefied matrices used by the active figure scripts.

---

## Required inputs

- `data_tables/meta_ALL.csv`
- `data_tables/ASV_taxa_all.zip`
- `data_tables/ASV_sequence_all.zip`
- `data_tables/original/Doubleton.zip`

---

## Running from scratch

1. Unzip `data_tables/ASV_taxa_all.zip`.
2. Unzip `data_tables/ASV_sequence_all.zip`.
3. Unzip `data_tables/original/Doubleton.zip` into `data_tables/original/Doubleton/`.
4. Run `1_QIIME Tables to working files.R`.
5. Run the figure scripts as needed.

The active figure scripts use the rarefied matrices listed in `data_tables/matrix_names.csv`.
The current repository already includes the outputs from `1_QIIME Tables to working files.R`, so rerunning it is only needed if you want to rebuild those files.

---

## Figure levels

- `Figure_2_taxonomy.R`: `Phylum`
- `Figure_3_6_alpha_diversity_stats_ANOVA.R`: `ASV`
- `Figure_4_7_PCA.R`: `ASV`
- `Figure_5_pipeline.R`: `ASV`
- `Figure_8_abundance.R`: `Phylum`

---

## Outputs

- Generated figures are written to `output_plot/`.
- Generated tables and cached intermediate results are written to `output_data/`.
- Rebuilt matrices are written to `data_tables/raw_matrix/` and `data_tables/rarefied_matrix/`.

---

## Analysis notes

- The active workflow is doubleton-only.
- `ASV` is the feature-table level.
- Community analyses use CLR-transformed `ASV` data with Euclidean distance.
- `PERMANOVA` uses `vegan::adonis2`.
- `PERMDISP` uses `vegan::betadisper` and `permutest`.

---

## Citation

Wu LYA#, Kahsen J#, Kunstman K, Naqib A, Green SJ.  
**A rapid and flexible method for long-read amplicon library preparation and balancing.**  
(Manuscript in preparation / under review)

---

## Contact

Jeremy Kahsen  
Rush University Medical Center  
Jeremy_Kahsen@Rush.edu
