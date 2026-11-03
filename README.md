# Code for responses to ALK inhibitors in heterogeneous differentiating cultures

## Repository layout
```
ALK_paper_code/
├── README.md
├── treatment/
│   ├── preprocessing/              # QC, cell type annotation, UMAP, sample/cell type summary tables
│   ├── differential_expression/
│   │   └── dreamlet/               # 1v4 DE, mash shrinkage, GSEA (Hallmark + GO:BP), DE/GSEA figures + supp tables
│   └── matrixeqtl/                 # cis-eQTL pipeline (treatment arm) + interaction concordance sensitivity
├── control/                        # Control HDC pipeline (feeds joint mash)
├── fastqtl2mash/                   # Reformat MatrixEQTL output for mash
├── mash/                           # Joint mash fit, response-eQTL definition, cCRE enrichment, supplementary figures/tables
└── chem_perturbation/              # CRISP zero-shot inference, response-eGene fgsea, methods schematic, validation/
```


## Single-cell preprocessing
- Cellranger alignment output is converted to on-disk BPCells matrices in `treatment/preprocessing/01_convert_disk_storage.R`
- Per-batch objects are merged into a single Seurat object in `treatment/preprocessing/02_merge_samples.R`
- QC filtering, TMM normalization, and scaling are applied in `treatment/preprocessing/03_filter_scale_normalize.R`
- The output of this sequence is a filtered, normalized Seurat object used for cell type annotation

### Cell type annotation
- Cell type annotation using CelliD (unguided, cluster-free) with a combined fetal cell atlas + first-trimester brain reference is run in `treatment/preprocessing/04_celltype_annotation.R`
- The same annotation script produces Fig 1C (`celltype_proportion_by_treatment`) and Supp Fig S1 (`celltype_proportion_by_individual`)
- Sample-level, cell-type-level, and per-individual summary tables (Supp Tables S2/S3/S4) are compiled in `treatment/preprocessing/04b_sample_celltype_indiv_summary_tables.R`

### UMAP visualization
- Sketch-based UMAP of 850K cells (Fig 1B) is generated in `treatment/preprocessing/05_umap_sketch.R`

## Differential expression

### 1v4 differential expression (dreamlet)
- The dreamlet pseudobulk SingleCellExperiment `pb_treat` (one assay per cell type) is built from the post-QC annotated Seurat object in `treatment/differential_expression/dreamlet/00_build_pb_treat.R`. This is a one-time preparation step; downstream scripts load the cached RDS
- Canonical 1v4 DE — each ALK inhibitor vs. the mean of the other four, with individual as a random effect — is run in `treatment/differential_expression/dreamlet/01_dreamlet_1v4_DE.R`
- Every gene retained by `dreamlet::processAssays()` is tested (no biotype / MT- / RP[SL] / avg-logCPM pre-filter); cell types with any drug arm below the 4-donor threshold are dropped up-front. Output is `top_all_1v4.csv`, the input to every downstream DE analysis
- A per-arm QC table (`DE_QC_powertier.csv`) is also produced

### Cross-cell-type shrinkage (mash)
- Per-drug `dreamlet::run_mash` fits (30 cell types each) are run via `treatment/differential_expression/dreamlet/02_run_mash.R`
- Post-mash summaries are in `treatment/differential_expression/dreamlet/03_summarize_mash.Rmd`
- FGSEA (Hallmark and GO:BP) on the mash posterior t-statistic is run in `treatment/differential_expression/dreamlet/04_dreamlet_mash_fgsea.R`. Output is `dreamlet_mash_gsea.tsv` (Supp Table S9) and the input to Fig 2C, Fig 2D, and Supp Figs S5/S6
- The pre-mash companion runs the same fgsea on the raw dreamlet moderated t-statistics (no mash shrinkage) in `treatment/differential_expression/dreamlet/04b_dreamlet_pre_mash_fgsea.R` -> `dreamlet_pre_mash_gsea.tsv`. This feeds the pre-vs-post sensitivity analysis, Supp Fig S6 (GO:BP pre-mash), and the Hallmark pre-mash companion to Supp Fig S5

### DE breadth tables
- Per-(drug, gene) DE breadth tables (pre-mash (`adj.P.Val < 0.05`), post-mash (`lfsr < 0.05`), and cell-type-specific (breadth K=1)) are produced in `treatment/differential_expression/dreamlet/05_de_breadth_tables.R` (Supp Tables S6/S7/S8, written as both `.csv` and `.xlsx` with the gene column force-formatted as text to avoid Excel date-mangling)

### DE figure and supplementary table code
- Fig 2A (DEG fraction dotplot), Fig 2B (featured volcanoes), and Supp Fig S2 (all-cell-type volcanoes) are produced together in `treatment/differential_expression/dreamlet/de_volcanoes_and_dotplot.R`
- Fig 2C (Hallmark heatmap): `treatment/differential_expression/dreamlet/gsea_hallmark_heatmap.R`
- Fig 2D (eye-development pathway tile): `treatment/differential_expression/dreamlet/gsea_eye_pathways_tile.R`
- Supp Fig S3 (DEG count vs cell-type abundance) and the DEG count/fraction CSVs consumed by Fig 2A: `treatment/differential_expression/dreamlet/de_counts_vs_abundance.R`
- Supp Fig S4 (cell-type-specific DEG bar chart): `treatment/differential_expression/dreamlet/ctspecific_deg_bar.R`
- Supp Table S5 (top DE calls as Excel with text-formatted gene column): `treatment/differential_expression/dreamlet/de_sig_table.R`
- Supp Fig S5 (all-drug stacked Hallmark heatmap) and Supp Fig S7 (per-drug top-N GO:BP heatmaps, post-mash): `treatment/differential_expression/dreamlet/gsea_per_celltype_heatmaps.R`. The pre-mash companion that produces Supp Fig S6 (GO:BP pre-mash) and the Hallmark pre-mash companion to Supp Fig S5 is `treatment/differential_expression/dreamlet/gsea_per_celltype_heatmaps_premash.R`
- Supp Tables S10 (full pre-mash GSEA) and S11 (pre-vs-post sensitivity summary for the gene sets cited in the Results): `treatment/differential_expression/dreamlet/gsea_pre_vs_post_sensitivity.R` (depends on both `04` and `04b`)

## cis-eQTL calling

### Control HDC pipeline
- The Popp et al. 2024 untreated HDC dataset is split, QC-filtered, cell-type-annotated, and cis-eQTL-scanned in `control/01_split_control_data.py`, `control/02_QC_control.R`, `control/03_celltype_annotation_control.R`, `control/04_matrixeqtl_control.R`, with eigenMT applied in `control/05_eigenmt_control.R`
- The control eQTL output is combined with the treatment eQTL output in the joint mash step below

### Treatment pseudobulk + cis-eQTL scan
- Pseudobulk aggregation, covariate selection, and inverse-normal transformation are prepared in `treatment/matrixeqtl/01_prep_pseudobulk_covariates.Rmd`
- Per-(cell type × treatment) cis-eQTL scans (100 kb window, MatrixEQTL, chromosome-parallelized) are run in `treatment/matrixeqtl/02_matrixeqtl_treatment.R`
- The joint treatment+control expression matrix used by the interaction model is prepared in `treatment/matrixeqtl/03_prep_interaction_input.Rmd`
- The interaction model (`modelLINEAR_CROSS`, genotype × treatment) is run in `treatment/matrixeqtl/04_matrixeqtl_interaction.R`
- eigenMT + Benjamini–Hochberg multiple-testing correction is applied to the interaction scan, and post-eigenMT calls are consolidated in `treatment/matrixeqtl/06_post_eigenmt_interaction.Rmd`

### Reformatting for mash
- MatrixEQTL output is reformatted into mash's expected structure in `fastqtl2mash/01_prep_mash.Rmd` and `fastqtl2mash/02_fastqtl_to_mash.ipynb`. Shared utility functions live in `fastqtl2mash/.sos/utils.py`

### Mash
- The joint mash fit on top variant-gene pairs across all conditions is in `mash/01_run_mash.R`
- Response eQTLs are defined as `lfsr < 0.05` in either treatment or control and a >=3-fold effect-size shift between treatment and matched control
- ENCODE cCRE backgrounds and enrichment tests are built in `mash/02_build_ccre_backgrounds.R` and `mash/03_run_ccre_enrichment.R`
- The main mash post-analysis narrative and Supp Table S15 (all significant eQTLs) live in `mash/04_post_mash_analysis.Rmd`
- Shared helper functions used across mash scripts are in `mash/mash_helpers.R`

### eQTL figure and supplementary table code
- Supp Fig S10 (cis-eQTL count heatmap): `mash/cis_eqtl_count_heatmap.py`
- Supp Figs S11 and S17 (power scatter and response-eQTL vs power) and Supp Table S20 (per-condition power and eQTL counts): `mash/power_vs_eqtl_summary.R`
- Supp Fig S12 (expression-eQTL sharing heatmap): `mash/pairwise_sharing.Rmd`
- Supp Fig S13 (TSS distance density per germ layer): `mash/tss_distance_density.R`
- Supp Fig S14 (response vs non-response cCRE enrichment) and Supp Table S16 (cCRE enrichment TSV): `mash/ccre_enrichment_plots.Rmd`
- Supp Fig S16 (effect-size distribution by drug): `mash/effect_size_distribution.Rmd`
- Supp Tables S18/S19/S21/S22 (response eQTL subset tables including within-cell-type drug-specific loci, condition-specific loci, and the Table 1-anchored candidate set): `mash/response_eqtl_subset_tables.Rmd`

### Interaction-model concordance sensitivity check
- Supp Fig S15 (Pearson r and sign concordance between the mash-derived delta_beta and the independent interaction-model beta_GI, plus an lfc5 sensitivity check) is generated in `treatment/matrixeqtl/interaction_vs_mash_concordance.R`

## Chemical perturbation prediction (CRISP)
- Per-cell-type h5ad slices of the untreated HDC data used as CRISP input are exported in `chem_perturbation/01_export_control_cells.R`
- CRISP model inputs (h5ad + drug SMILES embeddings) are prepared in `chem_perturbation/02_prep_input.py`
- CRISP inference (GPU) is run in `chem_perturbation/03_run_crisp.py`
- Per-gene z-scores (drug vs. mean of the other four drugs) are computed in `chem_perturbation/04_zscore.py`
- FGSEA of CRISP-predicted drug-specific signals against the response eGene set (Supp Table S23) is run in `chem_perturbation/05_crisp_response_eqtl_fgsea.R`
- Supp Fig S18 (CRISP methods schematic) is generated in `chem_perturbation/crisp_methods_schematic.py`
- Shared pipeline config, SMILES validation, and RDKit 2D normalization helpers live in `chem_perturbation/config.py`, `chem_perturbation/validate_smile_df.py`, and `chem_perturbation/make_rdkit2d_normstats.py`
