# Code availability — T2D PRS pipeline 

Type-2-Diabetes polygenic risk score (PRS) calculation for CES
samples (80 individuals: 50 cases + 30 controls) 
Software and versions: `env/environments.md`.

## Inputs

FASTQ + `sarek_samplesheet.csv`, hg38 reference, Gencode v44 GTF, 
`clinical_data.csv` (per-patient table containing clinical paramter information and miRNA expression), 
Suzuki 2024 EUR sumstats, hg19→hg38 chain,
PRSice-2, somalier and 1000G reference.


## Stages

| Directory | Does | Key output |
|---|---|---|
| `variant_calling/` | sarek align + germline call → VEP → keep canonical coding/splice/UTR variants → plink2 | `exonic_assoc` (pgen) |
| `coverage_coregenes/` | mosdepth per-gene depth → **core10x** (≥10× mean depth in all 80 samples) → restrict variant positions | `core10x_keep_pos.bed` |
| `prs/` | base GWAS QC + liftover hg19→hg38 (Suzuki 2024 EUR) → target QC + core10x restriction + KING → PRSice-2 Clumping + Thresholding | `T2D_PRS_core10x.best` |
| `ancestry/` | somalier and 1000G AIMs → interactive 3D PCA | ancestry `.tsv`, `ancestry_pca_3d.html` |
| `correlations_prs/` | PRS/miRNA × clinical Spearman, PRS distribution by group, clinical summary stats (one analysis per file) | figures + CSVs in `results/` |


