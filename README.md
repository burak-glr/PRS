# Code availability — T2D PRS pipeline (Clinical Exome Sequencing cohort)

Type-2-Diabetes polygenic risk score (PRS) calculation for CES
samples (80 individuals: 50 cases + 30 controls) 
Software and versions: `env/environments.md`.

## Inputs

FASTQ + `sarek_samplesheet.csv`, hg38 reference, Gencode v44 GTF, 
`clinical_data.csv` (per-patient table containing clinical paramter information and miRNA expression), 
cohort Excel, Suzuki 2024 EUR sumstats + hg19→hg38 chain,
PRSice-2, somalier + 1000G reference.

- **core10x** = genes with ≥10× mean depth in **every** one of the 80 samples (per-sample
  intersection), not the cohort average.

## Stages

| Directory | Does | Key output |
|---|---|---|
| `variant_calling/` | Sarek align + germline call → VEP → keep canonical coding/splice/UTR variants → plink2 | `exonic_assoc` (pgen) |
| `coverage_coregenes/` | mosdepth per-gene depth → **core10x** (≥10× mean depth in all 80 samples = 4,947 genes) → restrict variant positions | `core10x_keep_pos.bed` |
| `prs/` | base GWAS QC + liftover hg19→hg38 (Suzuki 2024 EUR) → target QC + core10x restriction + KING → PRSice-2 C+T | `T2D_PRS_core10x.best` |
| `ancestry/` | somalier and 1000G AIMs → interactive 3D PCA | ancestry `.tsv`, `ancestry_pca_3d.html` |
| `correlations_prs/` | PRS/miRNA × clinical Spearman, PRS distribution by group, clinical summary stats (one analysis per file) | figures + CSVs in `results/` |


