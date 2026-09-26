# abl_dms: code and processed data for the BCR-ABL imatinib VUDR study

This repository contains the analysis code and processed data for:

> Inam H.\*, Tomaszkiewicz M.\*, Reynolds J.A., Yang Z., Leighow S.M., and Pritchard J.R.
> **Characterizing Variants of Uncertain Drug Response (VUDRs) Using Quantitative Measurements at Clinical Exposures.**
> (\*equal contribution). Manuscript in revision at Cell Press.

The study uses deep mutational scanning (DMS) of the BCR-ABL kinase domain (residues 242–512) in Ba/F3 cells. Mutants were read out by TileSeq and duplex sequencing and screened at 300, 600 and 1200 nM imatinib plus a no-drug control. The resulting net growth rates were converted to concentration–response curves for about 4,900 missense variants. Each variant was then classified as resistant or sensitive at clinical imatinib exposures: 400 mg QD (444 nM), 400 mg BID (760 nM) and 500 mg BID (916 nM).

The repository contains:

- **`Figures/`**: one folder per manuscript figure panel (or group of panels). Each folder holds the script, its input files and the output plots. **Most readers will only need this folder.**
- **`analysis/dms_abl_analysis/`**: the upstream processing pipeline, in three parts: variant calling, net-growth-rate calculation, and the logistic-regression resistance classifier.
- Shared helper code, intermediate data and a rendered analysis website from the project's original [workflowr](https://github.com/workflowr/workflowr) layout (`code/`, `data/`, `output/`, `docs/`).

---

## Data availability

| Resource | Location |
|---|---|
| Raw sequencing reads (DMS screens) | NCBI SRA, BioProject [PRJNA1279649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279649) |
| Processed data (net growth rates, IC50s, classifications) | This repository (see tables below); all supplementary tables are compiled in `Figures/Robustness_Supp_Fig_9/ABL_Supplemental_tables.xlsx` |
| Interactive resistance map | https://amplicomics.shinyapps.io/version2/ (source code in `ShinyApp/`) |
| gnomAD v4.1.0 ABL1 variants | https://gnomad.broadinstitute.org/ (export for ENSG00000097007 included in `Figures/Figure_5/inputs/`) |
| COSMIC / Sanger CML ABL1 mutations | Exports included in `Figures/Figure_4D/`, `Figures/Figure_2B&C/inputs/Cosmic_ABL/` and `Figures/Sanger_Supp_Fig_10/` |

---

## Repository structure

```
abl_dms/
├── README.md                 This file
├── LICENSE                   MIT license
├── .zenodo.json              Metadata used by Zenodo when archiving a GitHub release
├── abl_dms.Rproj             RStudio project file; open this first
├── _workflowr.yml, .Rprofile workflowr configuration (optional; see "Software")
│
├── Figures/                  Code, inputs and outputs for every manuscript figure
│   ├── Figure_1/             Fig. 1E–I, Fig. S3  (region-1 pilot: NGS vs duplex vs TileSeq)
│   ├── Figure_2B&C/          Fig. 2B–C; builds the full-kinase net-growth table (Tables S2–S5)
│   ├── Figure_2D/            Fig. 2D   (1200 nM net-growth heatmap)
│   ├── Figure_3A/            Fig. 3A   (3- vs 10-point IC50s; Table S6)
│   ├── Figure_3C/            Fig. 3C   (DMS vs single-mutant viability, 18 standards)
│   ├── Figure_3D/            Fig. 3D   (DMS-inferred IC50 map)
│   ├── Figure_3E/            Fig. 3E   (DMS-inferred concentration–response curves)
│   ├── Figure_3F/            Fig. 3F   (test set of 17 mutants)
│   ├── Figure_4C/            Fig. 4C   (resistance by highest clinical dose; Table S10)
│   ├── Figure_4D/            Fig. 4D   (COSMIC mutation frequencies; Table S8)
│   ├── Figure_4E/            Fig. 4E   (VUDR classification by dose)
│   ├── Figure_5/             Fig. 5    (gnomAD-enabled variants; Table S13)
│   ├── Robustness_Supp_Fig_9/  Fig. S9 (PK-variability robustness; Table S14) + all supplementary tables
│   ├── Sanger_Supp_Fig_10/     Fig. S10 (cross-validation against Sanger/COSMIC clinical data)
│   └── Asciminib_analysis/     Reviewer-response Figure R1 (asciminib pilot; not in the manuscript)
│
├── analysis/
│   ├── dms_abl_analysis/     Upstream pipeline (see "Upstream processing pipeline")
│   │   ├── step1.generate.variantcalls/
│   │   ├── step2.generate.netgrowthrates/
│   │   └── Clinical_resistance_classification_of_mutations/
│   └── *.Rmd                 workflowr pages for exploratory/earlier analyses (spike-ins, error rates, simulations)
│
├── code/                     Shared R helper functions (used by Figures/Figure_1 and analysis/*.Rmd)
├── data/                     Intermediate inputs (variant-caller outputs, references, IC50 data, COSMIC, Twinstrand)
├── output/                   Intermediate outputs written by the legacy analyses
├── docs/                     Rendered workflowr website (HTML) for analysis/*.Rmd
└── ShinyApp/                 Source of the interactive resistance-map Shiny app
```

Folders named `Archive/` or `archive/` inside a figure folder hold superseded versions of scripts, inputs and plots. They are kept for provenance and are not needed to reproduce the figures.

---

## Figure-to-code map

"Run from" gives the working directory each script expects. The two cases are:

- **Figure folder**: the folder containing the script. This is RStudio's default when you knit an `.Rmd`. For a plain `.R` script, `setwd()` to its folder first.
- **Project root**: the top-level `abl_dms/` folder. Open `abl_dms.Rproj` in RStudio. When these `.Rmd` files are knitted, they set their root directory to the project root automatically.

| Manuscript item | Folder | Script | Run from | Key inputs | Outputs |
|---|---|---|---|---|---|
| Region-1 pilot data parsing (input to Fig. 1, S3) | `Figures/Figure_1/` | `1.ErrorRates_dataparser_region1.Rmd` | Project root | `data/Consensus_Data/novogene_lane18b_rerun/` (copy in `Figures/Figure_1/inputs/`), helpers in `code/` | `output/ABLEnrichmentScreens/ABL_Region1_lane18b/<comparison>/` (copies used downstream are in `Figures/Figure_1/outputs/`) |
| Fig. 1E–F: intended vs. unintended variant counts (NGS, duplex, TileSeq) | `Figures/Figure_1/` | `2.Fig1AEF_TileseqvsNGSvsDuplex.Rmd` | Project root | `Figures/Figure_1/outputs/*/screen_comparison_*.csv` | `fig1A_abl_ngs_untreated.pdf`, `fig1E_abl_tileseq_untreated.pdf`, `fig1F_abl_duplex_untreated.pdf` |
| Fig. S3: null-normalized mean difference, ROC | `Figures/Figure_1/` | `3.SuppFig3_NNMDAnalyses_region1.Rmd` | Project root | as above | `Supp3B.pdf` |
| Fig. 1G–I: mutant standards, pooled vs. single-mutant net growth | `Figures/Figure_1/` | `4.Fig1GHI_mutant_standards.Rmd` | Project root | as above + `inputs/ic50data_all_conc.csv` | `Fig1G_ngs.pdf`, `Fig1H_duplex.pdf`, `Fig1I_tileseq.pdf` |
| Fig. 2B–C; Tables S2–S5 | `Figures/Figure_2B&C/` | `tileseq_dataparser.Rmd` | Figure folder | `inputs/TileSeq_full_kinase_*_renamed.csv` (IL-3/no drug, low, medium, high imatinib), `inputs/ic50data/`, `inputs/Refs/`, helper functions in `inputs/*.R` | `outputs/TileSeq_full_kinase_alldoses_with_lfc_corrected_netgrowths.csv`; `outputs/Fig_2b_*.pdf`, `outputs/Fig_2c_*.pdf`, per-dose correlation and density plots |
| Fig. 2D | `Figures/Figure_2D/` | `High_imatinib_netgrowth_map.Rmd` | Figure folder | `high_imatinib_netgrowth_data.csv` | `High_imatinib_netgrowth_map_simple.{pdf,png}` |
| Fig. 3A; Table S6 | `Figures/Figure_3A/` | `Minimal_number_of_concentrations_profiling.Rmd` | Figure folder | `IC50_3_versus_10_concentrations.csv` | `IC50_correlation_3_versus_10_concentration.pdf` |
| Fig. 3C | `Figures/Figure_3C/` | `mutant_standards_verification.Rmd` | Project root | `inputs/IC50HeatMap_v2.csv` (single-mutant data, Table S0), `inputs/all_corrected_without_poorly_fit_removed_w_drc_11.csv` (**unzip first**) | `clonalvsdms_mutantstandards.pdf` |
| Fig. 3D | `Figures/Figure_3D/` | `Concentration_response_map.Rmd` | Figure folder | `all_mutants_with_IC50_estimates.csv` | `Concentration response map.pdf` |
| Fig. 3E | `Figures/Figure_3E/` | `Imatinib_Concentration_Response_Curves.Rmd` | Figure folder | `all_mutants_with_relative_viability.csv` (**unzip first**) | `Conc_response_curves.{pdf,png}` |
| Fig. 3F | `Figures/Figure_3F/` | `dose_escalation_verification.Rmd` | Project root | `inputs/IC50HeatMap_verification_Elvin.csv`, `inputs/all_corrected_without_poorly_fit_removed_w_drc_11.csv` (**unzip first**) | `IC50s_clonalvsdms.pdf` |
| Fig. 4C; Table S10 | `Figures/Figure_4C/` | `Resistance_map.Rmd` | Figure folder | `all_mutants_with_relative_viability_at_three_clinical_doses.csv` | `Map_by_highest_dose.pdf`, `Highest_resistance_mutation_classification.csv` |
| Fig. 4D; Table S8 | `Figures/Figure_4D/` | `COSMIC_mutations_piechart.Rmd` | Figure folder | `COSMIC_mutations.csv` | `Distribution of COSMIC Mutation Frequencies including VUDRs.pdf` |
| Fig. 4E | `Figures/Figure_4E/` | `VUDR_classification_by_highest_dose.Rmd` | Figure folder | `VUDR_data_classified.csv` | `VUDR classification by highest clinical dose.pdf` |
| Fig. 5A–C, E–G; Table S13 | `Figures/Figure_5/` | `gnomad_analysis.Rmd` | Anywhere in the project (uses `here`) | `inputs/gnomAD_v4.1.0_ENSG00000097007_*.csv`, `inputs/Refs/ABL/`, `inputs/codon_table.csv`, `inputs/TileSeq_full_kinase_alldoses_with_lfc_corrected_netgrowths.csv`, `inputs/Tileseq_clinical_predictions_10.28.25.csv` | `gnomad_col_individuals.pdf`, `gnomad_col_mutants.pdf`, `gnomad_enabled_mutants_penetrance.pdf`, `clin_var_piechart_cleaned.pdf`, `resistance_distributions_cleaned.pdf`, `table_s13.csv` |
| Fig. S9; Table S14 | `Figures/Robustness_Supp_Fig_9/` | `figure4_robustness.Rmd` (`figure4_robustness.R` is an earlier script version) | Figure folder | `ABL_Supplemental_tables.xlsx` (sheets `Table_S10`, `Table_S14`) | `figure_combined_robustness.{pdf,png}` (edited version: `figure_combined_robustness_edited.pdf`), `figure_robustness_by_dose_single.{pdf,png}`; assembled figure `supplemental_fig_9_08.27.26.jpg` |
| Fig. S10 | `Figures/Sanger_Supp_Fig_10/` | `abl_Sanger_enrichment_analysis_updated.Rmd` | Figure folder | `Tileseq_clinical_predictions_10.28.25.csv`, `ABL_Sanger_Gene_mutations_2026.csv`, `hammingdistance_ablmutants.csv`, `known_clinicallyresistant_mutations_betterannotations.csv` | `figure_cosmic_enrichment.{pdf,png}`, `sanger_enrichment_stats.csv`, `sanger_enrichment_tiers.csv`; assembled figure `supplemental_fig_10_08.19.26.jpg` |
| Reviewer Figure R1 (asciminib pilot, 304 variants; not in manuscript) | `Figures/Asciminib_analysis/` | `make_figure_R1.R` | Figure folder | `K5.ABL.v1-01_260813_*_filtered.csv` (unfiltered data in `alldata/`) | `Figure_R1_asciminib.{pdf,png}` |

Schematic panels (e.g., Fig. 1C, Fig. 5D, Fig. S11) were drawn by hand and have no code.

### Supplementary tables

All supplementary tables (S0–S14) are compiled as separate sheets in `Figures/Robustness_Supp_Fig_9/ABL_Supplemental_tables.xlsx`. The code that generates or uses each table:

| Table | Content | Generated by / source file |
|---|---|---|
| S0 | Single-mutant IC50s (21 mutants + WT) | `Figures/Figure_3C/inputs/IC50HeatMap_v2.csv`, `Figures/Figure_1/inputs/ic50data_all_conc.csv` |
| S2–S5 | Net growth rates at 1200, 300 and 600 nM, and no drug | `Figures/Figure_2B&C/tileseq_dataparser.Rmd` → `outputs/TileSeq_full_kinase_alldoses_with_lfc_corrected_netgrowths.csv` |
| S6 | 10- vs. 3-point IC50s | `Figures/Figure_3A/IC50_3_versus_10_concentrations.csv` |
| S7 | Classifier training set | `analysis/dms_abl_analysis/Clinical_resistance_classification_of_mutations/train_logistic_classifier_and_define_resistance_threshold.R` |
| S8 | COSMIC mutations with resistance annotations / VUDRs | `Figures/Figure_4D/COSMIC_mutations.csv`, `Figures/Figure_4E/VUDR_data_classified.csv` |
| S9 | Classifier test set | `analysis/dms_abl_analysis/Clinical_resistance_classification_of_mutations/validate_logistic_classifier_testset.R` |
| S10 | Highest clinically resistant dose for each variant | `Figures/Figure_4C/Resistance_map.Rmd` → `Highest_resistance_mutation_classification.csv` |
| S13 | gnomAD-enabled substitutions | `Figures/Figure_5/gnomad_analysis.Rmd` → `table_s13.csv` |
| S14 | Robustness to simulated PK variability | `Figures/Robustness_Supp_Fig_9/` (read from `ABL_Supplemental_tables.xlsx`; CV sweep recomputed in `figure4_robustness.Rmd`) |

Table S1 (spike-in experiments) relates to the spike-in analyses in `analysis/spikeins_*.Rmd` (rendered in `docs/`). Tables S11 and S12 (TileSeq primers and sequencing summary) are experimental records and are provided only in the workbook.

---

## Upstream processing pipeline (`analysis/dms_abl_analysis/`)

These scripts turn sequencing-derived variant counts into the net growth rates, IC50s and resistance calls that the figure scripts use.

1. **`step1.generate.variantcalls/`**: `variant_caller_2024.Rmd` is the custom R variant caller. It reads consensus-called reads that were aligned with bwa-mem2 and filtered, then produces per-variant counts and per-residue depths. `readme.md` describes the upstream command-line steps (du novo consensus calling, bwa-mem2 alignment, removal of mouse reads).
2. **`step2.generate.netgrowthrates/`**: `dms_abl_analysis.Rmd` and the helpers in its `code/` folder calculate mutant net growth rates. They combine variant allele frequencies with cell densities and dilution factors for each timepoint (Day 0 vs. Day 6). Example consensus data are in `data/Consensus_Data/Novogene_lane18/` and an example output is in `output/`.
3. **`Clinical_resistance_classification_of_mutations/`**: all three scripts read `merged_table_imatinib_full_kinase_ngr_with_resmuts.csv`. Run them from this folder.
   - `predicting_IC50_relviab_from_netgrowth.R`: fits a `dr4pl` concentration–response curve to each variant, giving IC50s and relative viability at 444, 760 and 916 nM.
   - `train_logistic_classifier_and_define_resistance_threshold.R`: trains the logistic-regression classifier on 22 gold-standard mutants. The optimal cutoff comes from the Youden index and is back-calculated to the relative-viability resistance threshold (0.42). This produces Table S7 and the analyses behind Fig. 4B and Figs. S7–S8.
   - `validate_logistic_classifier_testset.R`: checks the classifier on the independent test set of 15 literature-curated resistant mutants (Table S9).

`analysis/dms_abl_analysis/gnomad_analyses.Rmd` is an earlier version of the gnomAD analysis. It has been superseded by `Figures/Figure_5/gnomad_analysis.Rmd`.

Before these steps, the command-line read processing was done in Linux and is described in the STAR Methods. For duplex sequencing, du novo did consensus calling and bwa-mem2 aligned the reads to ABL1 NM_005157.6. For TileSeq, PEAR v0.9.11 merged paired-end reads and bwa-mem2 aligned them. Raw reads are on SRA (PRJNA1279649).

---

## Other folders

- **`code/`**: shared R functions. They cover variant parsing (`variants_parser.R`), merging and comparing samples and screens to compute growth rates (`merge_samples*.R`, `compare_samples.R`, `compare_screens*.R`, `depth_finder.R`), annotation (`resmuts_adder.R`, `res_residues_adder.R`, `is_intended_adder.R`, `cosmic_data_adder.R`, `shortest_codon_finder.R`) and plotting (`plotting/`). The `Figures/Figure_1` scripts `source()` these files from the project root. `Figures/Figure_2B&C/inputs/` holds self-contained copies of the helpers it needs.
- **`data/`**: intermediate inputs for the region-1 pilot and earlier analyses:
  - `Consensus_Data/`: variant-caller outputs
  - `Refs/`: ABL1, EGFR and LTK reference sequences and coordinates
  - `Cosmic_ABL/`: COSMIC ABL1 mutations
  - `ic50data/`: single-mutant IC50 data
  - `codon_table.csv`, `codon_tables/`: codon tables
  - `Twinstrand/`: duplex-sequencing deliverables for the spike-in experiments
- **`output/`**: intermediate results written by the legacy analyses (enrichment screens, spike-in and error-rate figures, enrichment simulations).
- **`analysis/*.Rmd` and `docs/`**: source pages and the rendered HTML site for earlier and exploratory analyses: spike-in mutant experiments, error rates, enrichment simulations, TileSeq parsing and dose–response fitting. Open `docs/index.html` in a browser to explore them.
- **`ShinyApp/`**: `app.R` and its input (`Highest_resistance_mutation_classification.csv`) for the interactive resistance map at https://amplicomics.shinyapps.io/version2/.

---

## Software

The analyses were written in R (≥ 4.0) and run in RStudio. Install the required CRAN packages with:

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "reshape2", "tidyverse", "ggplot2", "ggpubr",
  "ggrepel", "patchwork", "scales", "RColorBrewer", "plotly", "dr4pl", "pROC",
  "binom", "pracma", "readxl", "here", "doParallel", "foreach", "tictoc",
  "rmarkdown", "knitr", "rstudioapi"
))
# Optional: "workflowr" (legacy site in analysis/ and docs/), "shiny" (ShinyApp/)
```

`.Rprofile` loads `workflowr` if it is installed and prints a message if it is not. The figure scripts do not need `workflowr`.

## Reproducing a figure

1. Download this repository (from Zenodo or GitHub) and open `abl_dms.Rproj` in RStudio.
2. Three large inputs (~107 MB each once unzipped) are stored zipped to stay under GitHub's file-size limit. Unzip each in place before running Figures 3C, 3E or 3F:
   - `Figures/Figure_3C/inputs/all_corrected_without_poorly_fit_removed_w_drc_11.csv.zip`
   - `Figures/Figure_3F/inputs/all_corrected_without_poorly_fit_removed_w_drc_11.csv.zip`
   - `Figures/Figure_3E/all_mutants_with_relative_viability.csv.zip`

   (`.gitignore` excludes the unzipped files so they are not committed by accident.)
3. Open the script listed in the figure-to-code map and knit it, or run its chunks from the working directory given in the "Run from" column. Outputs are written back into the same figure folder.

**Note for Linux users:** the scripts were run on macOS, whose file system ignores case. `code/compare_screens_archive.R` (used by `Figures/Figure_1/1.ErrorRates_dataparser_region1.Rmd`) builds paths as `data/Consensus_Data/Novogene_lane18b_rerun/...`, but the folder on disk is `data/Consensus_Data/novogene_lane18b_rerun/`. On a case-sensitive file system, add a symlink before running that script: `ln -s novogene_lane18b_rerun data/Consensus_Data/Novogene_lane18b_rerun`. None of the other figure scripts are affected.

---

## License

The code in this repository is released under the MIT License (see `LICENSE`). Third-party data (gnomAD, COSMIC/Sanger) remain subject to the terms of their original providers.

## Citation

If you use this code or data, please cite the manuscript above and the archived version of this repository on Zenodo (DOI listed on the Zenodo record and in the manuscript's Data and Code Availability section).

## Contact

Lead contact: Justin R. Pritchard (jrp94@psu.edu), Department of Biomedical Engineering, The Pennsylvania State University.
Code questions: Haider Inam and Marta Tomaszkiewicz.
