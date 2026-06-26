#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=prsice_core10x
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/prsice_core10x_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/prsice_core10x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00


eval "$(conda shell.bash hook)"
conda activate rplots

PLINK2=/cta/users/burak.guler/miniforge3/bin/plink2
PRSICE_BIN=/cta/users/burak.guler/miRNA/CES_DATALARI/association/PRSice/bin/PRSice
PRSICE_R=/cta/users/burak.guler/miRNA/CES_DATALARI/association/PRSice/PRSice.R

OUT_DIR=/cta/users/burak.guler/miRNA/CES_DATALARI/association/prs_core10x
BASE_OVERLAP="$OUT_DIR/base_exonic_overlap.txt"
BASE_WITH_ID="$OUT_DIR/base_with_snpid.txt"
TARGET_QC="$OUT_DIR/target_qc"
TARGET_NAMED="$OUT_DIR/target_named"
PHENO=/cta/users/burak.guler/miRNA/CES_DATALARI/association/phenotype.txt
COVAR=/cta/users/burak.guler/miRNA/CES_DATALARI/association/covariates.txt
PHENO_FIXED="$OUT_DIR/phenotype_nohash.txt"
COVAR_WITH_PCS="$OUT_DIR/covariates_with_pcs.txt"
KEEP="$OUT_DIR/target_king.king.cutoff.in.id"
OUT="$OUT_DIR/T2D_PRS_core10x"


# PART 1: Add SNP ID column (CHR_HG38:POS_HG38)
awk 'NR==1 {print "SNP\t"$0; next}
         { print $1":"$2"\t"$0 }' "$BASE_OVERLAP" > "$BASE_WITH_ID"
echo "Base variants with SNP ID: $(($(wc -l < "$BASE_WITH_ID") - 1))"

# PART 2: Convert QC-clean pgen → plink1 bed
$PLINK2 \
    --pfile "$TARGET_QC" \
    --set-all-var-ids '@:#' \
    --new-id-max-allele-len 50 \
    --rm-dup exclude-mismatch \
    --make-bed \
    --out "$TARGET_NAMED" \
    --threads 8
echo "Target variants (bed): $(wc -l < "${TARGET_NAMED}.bim")"
echo "Target samples  (bed): $(wc -l < "${TARGET_NAMED}.fam")"

# PART 3: LD pruning + PCA (recomputed on the restricted variant set)
$PLINK2 \
    --bfile "$TARGET_NAMED" \
    --indep-pairwise 200 50 0.25 \
    --out "$OUT_DIR/target_prune" \
    --threads 8
echo "Pruned SNPs for PCA: $(wc -l < "$OUT_DIR/target_prune.prune.in")"

$PLINK2 \
    --bfile "$TARGET_NAMED" \
    --extract "$OUT_DIR/target_prune.prune.in" \
    --pca 6 \
    --out "$OUT_DIR/target_pcs" \
    --threads 8

# PART 4: Merge PCs into covariates
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
Rscript "$SCRIPTDIR/merge_covariates_pcs.R" \
    /cta/users/burak.guler/miRNA/CES_DATALARI/association/covariates.txt \
    "$OUT_DIR/target_pcs.eigenvec" \
    "$OUT_DIR/covariates_with_pcs.txt"

sed 's/^#//' "$PHENO" > "$PHENO_FIXED"

# PART 5: PRSice-2
Rscript "$PRSICE_R" \
    --prsice        "$PRSICE_BIN" \
    --base          "$BASE_WITH_ID" \
    --snp           SNP \
    --chr           CHR_HG38 \
    --bp            POS_HG38 \
    --A1            EffectAllele \
    --A2            NonEffectAllele \
    --stat          Beta \
    --beta \
    --pvalue        Pval \
    --base-maf      EAF:0.01 \
    --target        "$TARGET_NAMED" \
    --pheno         "$PHENO_FIXED" \
    --pheno-col     PHENO \
    --cov           "$COVAR_WITH_PCS" \
    --cov-col       SEX,AGE,PC1,PC2,PC3,PC4,PC5,PC6 \
    --keep          "$KEEP" \
    --ignore-fid \
    --binary-target T \
    --clump-kb      250 \
    --clump-r2      0.1 \
    --bar-levels    5e-8,1e-6,1e-4,1e-3,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --print-snp \
    --out           "$OUT" \
    --thread        8

echo ""
echo " PRSice RESULTS (core10x: 4,947 genes) "
cat "${OUT}.summary"
echo ""
echo "Job finished at $(date)"
