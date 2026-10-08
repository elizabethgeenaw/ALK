# Code for responses to ALK inhibitors in heterogeneous differentiating cultures

## Repository layout
```
ALK_paper_code/
├── README.md
├── treatment/
│   ├── preprocessing/              # 01–05  QC, annotation, UMAP (Fig 1, Supp S1, Supp Tables S2–S4)
│   ├── differential_expression/
│   │   ├── dreamlet/               # 01–05 analysis (+ 04b pre-mash fgsea) + fig2A/B/S2, fig2C, fig2D, supp_S3–S6 (+ supp_S10_S11 pre-mash sensitivity, supp_S5fig_S6_premash companion)
│   │   └── cormotif/               # 01–04 analysis + 06_regen_s7_figure + fig3A, supp_S8
│   └── matrixeqtl/                 # 01–06 pipeline + supp_S14_interaction_concordance
├── control/                        # 01–05 control HDC pipeline (feeds joint mash)
├── fastqtl2mash/                   # 01–03 reformat MatrixEQTL output for mash
├── mash/                           # 01–04 analysis + supp_S9…supp_S22 (figs + tables)
└── chem_perturbation/              # 01–05 CRISP inference + supp_S17 (methods schematic) + validation/
```


## Single-cell preprocessing
- Cellranger alignment output is converted to on-disk BPCells matrices in `treatment/preprocessing/01_convert_disk_storage.R`
- Per-batch objects are merged into a single Seurat object in `treatment/preprocessing/02_merge_samples.R`
- QC filtering, TMM normalization, and scaling are applied in `treatment/preprocessing/03_filter_scale_normalize.R`
- The output of this sequence is a filtered, normalized Seurat object used for cell type annotation

### Cell type annotation
- Cell type annotation using CelliD (unguided, cluster-free) with a combined fetal cell atlas + first-trimester brain reference is run in `treatment/preprocessing/04_celltype_annotation.R`
- The same annotation script produces Fig 1C (`celltype_proportion_by_treatment`) and Supp Fig S1 (`celltype_proportion_by_individual`)
- Sample-level, cell-type-level, and per-individual summary tables (Supp Tables S2/S3/S4) are compiled in `treatment/preprocessing/04b_supp_tables_S2_S3_S4.R`

### UMAP visualization
- Sketch-based UMAP of 850K cells (Fig 1B) is generated in `treatment/preprocessing/05_fig1B_sketch_umap.R`

## Differential expression

### 1v4 differential expression (dreamlet)
- The dreamlet pseudobulk SingleCellExperiment `pb_treat` (one assay per cell type; consumed by both dreamlet and cormotif) is built from the post-QC annotated Seurat object in `treatment/differential_expression/dreamlet/00_build_pb_treat.R`. This is a one-time preparation step; downstream scripts load the cached RDS
- Canonical 1v4 DE — each ALK inhibitor vs. the mean of the other four, with individual as a random effect — is run in `treatment/differential_expression/dreamlet/01_dreamlet_1v4_DE.R`
- Every gene retained by `dreamlet::processAssays()` is tested (no biotype / MT- / RP[SL] / avg-logCPM pre-filter); cell types with any drug arm below the 4-donor threshold are dropped up-front. Output is `top_all_1v4.csv`, the input to every downstream DE analysis
- A per-arm QC table (`DE_QC_powertier.csv`) is also produced
- The script additionally writes per-cell-type gene universes to `results/treatment/dreamlet_1v4_DE/universes/<cell type>.txt`; Cormotif reads these so the two analyses share the same per-cell-type gene set

### Cross-cell-type shrinkage (mash)
- Per-drug `dreamlet::run_mash` fits (30 cell types each) are launched in parallel via `treatment/differential_expression/dreamlet/02_run_mash_launcher.sh`, which sbatches `02_run_mash_single.sh` per drug; both wrap `02_run_mash.R`
- Post-mash summaries are in `treatment/differential_expression/dreamlet/03_summarize_mash.Rmd`
- FGSEA (Hallmark and GO:BP) on the mash posterior t-statistic is run in `treatment/differential_expression/dreamlet/04_dreamlet_mash_fgsea.R`. Output is `dreamlet_mash_gsea.tsv` (Supp Table S9) and the input to Fig 2C, Fig 2D, and Supp Figs S5/S6
- The pre-mash companion runs the same fgsea on the raw dreamlet moderated t-statistics (no mash shrinkage) in `treatment/differential_expression/dreamlet/04b_dreamlet_pre_mash_fgsea.R` -> `dreamlet_pre_mash_gsea.tsv`. This feeds the pre-vs-post sensitivity analysis, Supp Fig S6 (GO:BP pre-mash), and the Hallmark pre-mash companion to Supp Fig S5

### DE breadth tables
- Per-(drug, gene) DE breadth tables (pre-mash (`adj.P.Val < 0.05`), post-mash (`lfsr < 0.05`), and cell-type-specific (breadth K=1)) are produced in `treatment/differential_expression/dreamlet/05_supp_S6_S7_S8_breadth_tables.R` (Supp Tables S6/S7/S8, written as both `.csv` and `.xlsx` with the gene column force-formatted as text to avoid Excel date-mangling)

### DE figure and supplementary table code
- Fig 2A (DEG fraction dotplot), Fig 2B (featured volcanoes), and Supp Fig S2 (all-cell-type volcanoes) are produced together in `treatment/differential_expression/dreamlet/fig2A_2B_S2_volcanoes_dotplot.R`
- Fig 2C (Hallmark heatmap): `treatment/differential_expression/dreamlet/fig2C_hallmark_heatmap.R`
- Fig 2D (eye-development pathway tile): `treatment/differential_expression/dreamlet/fig2D_eye_pathways.R`
- Supp Fig S3 (DEG count vs cell-type abundance) and the DEG count/fraction CSVs consumed by Fig 2A: `treatment/differential_expression/dreamlet/supp_S3_counts_power.R`
- Supp Fig S4 (cell-type-specific DEG bar chart): `treatment/differential_expression/dreamlet/supp_S4_ctspecific_bar.R`
- Supp Table S5 (top DE calls as Excel with text-formatted gene column): `treatment/differential_expression/dreamlet/supp_S5_de_sig_table.R`
- Supp Fig S5 (all-drug stacked Hallmark heatmap) and Supp Fig S7 (per-drug top-N GO:BP heatmaps, post-mash): `treatment/differential_expression/dreamlet/supp_S5fig_S6_gsea_ct_patterns.R`. The pre-mash companion that produces Supp Fig S6 (GO:BP pre-mash) and the Hallmark pre-mash companion to Supp Fig S5 is `treatment/differential_expression/dreamlet/supp_S5fig_S6_premash.R`
- Supp Tables S10 (full pre-mash GSEA) and S11 (pre-vs-post sensitivity summary for the gene sets cited in the Results): `treatment/differential_expression/dreamlet/supp_S10_S11_pre_mash_sensitivity.R` (depends on both `04` and `04b`)

## Cormotif
- Cell types passing filters are identified in `treatment/differential_expression/cormotif/01_setdiff_check.R`
- Raw per-drug count matrices used as Cormotif input are extracted in `treatment/differential_expression/cormotif/02_prep_raw_counts.R`
- The main Cormotif run (K=2–8, BIC-selected per cell type; input restricted to the per-cell-type gene universe written by `01_dreamlet_1v4_DE.R` — no biotype / MT / RP[SL] / avg-logCPM filter) is in `treatment/differential_expression/cormotif/03_cormotif_run.R`; this produces Supp Fig S8 (Cormotif GO:BP barchart pages) and Supp Tables S10–S14 (motif pattern table, gene lists, GO:BP/CC/MF enrichment)
- The motif pattern table can be regenerated without refitting Cormotif via `treatment/differential_expression/cormotif/04_rerun_pattern_table.R`, which sources `03_cormotif_run.R` with the `CORMOTIF_HELPERS_ONLY` guard so classification stays in one place
- Group-specific motifs (DE in the 6 between-group comparisons for a 2-drug grouping and in none of the 4 within-group comparisons) are classified directly inside `compute_motif_meta()` in `03_cormotif_run.R`; the pattern table carries `motif_type == "group-specific"` and `drug_pattern == "<group> vs <others>"` (e.g. `"Cer+Lor vs Ale+Bri+Cri"`) on the first write, and the S7 panel titles render these as `"<group> vs others"` via `motif_type_display_label()`
- Supp Fig S8 (Cormotif GO BP enrichment pages) can be re-rendered without refitting via `treatment/differential_expression/cormotif/06_regen_s7_figure.R`, which also sources `03_cormotif_run.R` with the helpers-only guard and calls `assemble_go_supplementary_pages()`
- Fig 3A (combined-specificity heatmap for cardiac mesoderm + hepatoblasts): `treatment/differential_expression/cormotif/fig3A_motifs.R`
- Fig 3B/3C (drug-specific cell-type matrix + Lor motif upset) are re-generated by the wrapper `treatment/differential_expression/cormotif/fig3B_3C_refresh.R`, which sources shared definitions from `cormotif_figures.R`. The drug x cell-type matrix in `cormotif_figures.R::fig3B_drug_celltype_matrix()` has an `include_cerlor_column` option that appends an extra "Cer+Lor vs others" column marking cell types with a Cer+Lor group-specific motif
- Fig 3B companion (Lor motif x gene heatmap): `treatment/differential_expression/cormotif/fig3B_lor_heatmap.R`
- Supp Fig S9 (per-drug UpSet plots): `treatment/differential_expression/cormotif/supp_S8_upset.R`

## cis-eQTL calling

### Control HDC pipeline
- The Popp et al. 2024 untreated HDC dataset is split, QC-filtered, cell-type-annotated, and cis-eQTL-scanned in `control/01_split_control_data.py`, `control/02_QC_control.R`, `control/03_celltype_annotation_control.R`, `control/04_matrixeqtl_control.R`, `control/05_eigenmt_control.sh`
- The control eQTL output is combined with the treatment eQTL output in the joint mash step below

### Treatment pseudobulk + cis-eQTL scan
- Pseudobulk aggregation, covariate selection, and inverse-normal transformation are prepared in `treatment/matrixeqtl/01_prep_pseudobulk_covariates.Rmd`
- Per-(cell type × treatment) cis-eQTL scans (100 kb window, MatrixEQTL, chromosome-parallelized as a SLURM array) are run in `treatment/matrixeqtl/02_matrixeqtl_treatment.R` via `02_matrixeqtl_treatment.sh`
- The joint treatment+control expression matrix used by the interaction model is prepared in `treatment/matrixeqtl/03_prep_interaction_input.Rmd`
- The interaction model (`modelLINEAR_CROSS`, genotype × treatment) is run in `treatment/matrixeqtl/04_matrixeqtl_interaction.R` via `04_matrixeqtl_interaction.sh`
- eigenMT + Benjamini–Hochberg multiple-testing correction on the interaction scan is applied via `treatment/matrixeqtl/05_run_eigenmt_interaction.sh`, and post-eigenMT calls are consolidated in `06_post_eigenmt_interaction.Rmd`

### Reformatting for mash
- MatrixEQTL output is reformatted into mash's expected structure in `fastqtl2mash/01_prep_mash.Rmd` and `fastqtl2mash/02_fastqtl_to_mash.ipynb`, run via `03_run_fastqtl2mash.sh`. Shared utility functions live in `fastqtl2mash/.sos/utils.py`

### Mash
- The joint mash fit on top variant-gene pairs across all conditions is in `mash/01_run_mash.R`
- Response eQTLs are defined as `lfsr < 0.05` in either treatment or control and a >=3-fold effect-size shift between treatment and matched control
- ENCODE cCRE backgrounds and enrichment tests are built in `mash/02_build_ccre_backgrounds.R` and `mash/03_run_ccre_enrichment.R` (each with its own `.sbatch`)
- The main mash post-analysis narrative and Supp Table S15 (all significant eQTLs) live in `mash/04_post_mash_analysis_S15.Rmd`
- Shared helper functions used across mash scripts are in `mash/mash_helpers.R`

### eQTL figure and supplementary table code
- Supp Fig S10 (cis-eQTL count heatmap): `mash/supp_S9_cis_eqtl_count_heatmap.py`
- Supp Figs S11 and S17 (power scatter and response-eQTL vs power) and Supp Table S20 (per-condition power and eQTL counts): `mash/supp_S10_S16_S20_power_summary.R`
- Supp Fig S12 (expression-eQTL sharing heatmap): `mash/supp_S11_pairwise_sharing.Rmd`
- Supp Fig S13 (TSS distance density per germ layer): `mash/supp_S12_tss_distance.R` (submitted via `supp_S12_tss_distance.sbatch`)
- Supp Fig S14 (response vs non-response cCRE enrichment) and Supp Table S16 (cCRE enrichment TSV): `mash/supp_S13fig_S16table_ccre_plots.Rmd`
- Supp Fig S16 (effect-size distribution by drug): `mash/supp_S15fig_effect_sizes.Rmd`
- Supp Tables S18/S19/S21/S22 (response eQTL subset tables including within-cell-type drug-specific loci, condition-specific loci, and the Table 1-anchored candidate set): `mash/supp_S18_S19_S21_S22_reQTL_subsets.Rmd`

### Interaction-model concordance sensitivity check
- Supp Fig S15 (Pearson r and sign concordance between the mash-derived delta_beta and the independent interaction-model beta_GI, plus an lfc5 sensitivity check) is generated in `treatment/matrixeqtl/supp_S14_interaction_concordance.R`

## Chemical perturbation prediction (CRISP)
- Per-cell-type h5ad slices of the untreated HDC data used as CRISP input are exported in `chem_perturbation/01_export_control_cells.R`
- CRISP model inputs (h5ad + drug SMILES embeddings) are prepared in `chem_perturbation/02_prep_input.py`, submitted via `02_run_prep_input.sbatch` (or `02_run_prep_input_highmem.sbatch` for high-memory cases)
- CRISP inference (GPU) is run in `chem_perturbation/03_run_crisp.py`; per-cell-type submissions and OOM re-runs are handled by the family of `03_run_crisp_*.sbatch` scripts (and the combined `03_04_run_crisp_and_zscore.sbatch`)
- Per-gene z-scores (drug vs. mean of the other four drugs) are computed in `chem_perturbation/04_zscore.py` (submitted via `04_run_zscore.sbatch`)
- FGSEA of CRISP-predicted drug-specific signals against the response eGene set (Supp Table S23) is run in `chem_perturbation/05_fgsea_S23.R`
- Supp Fig S18 (CRISP methods schematic) is generated in `chem_perturbation/supp_S17_methods_schematic.py`
- Shared pipeline config, SMILES validation, and RDKit 2D normalization helpers live in `chem_perturbation/config.py`, `chem_perturbation/validate_smile_df.py` (with `validate_smile_df.sbatch`), and `chem_perturbation/make_rdkit2d_normstats.py`
