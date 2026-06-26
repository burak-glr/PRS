#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=somalier_ancestry
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/somalier_ancestry_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/somalier_ancestry_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Genetic ancestry inference for the 80-sample cohort by using somalier.

ANC=/cta/users/burak.guler/miRNA/CES_DATALARI/ancestry
SOMALIER=$ANC/somalier
SITES=$ANC/sites.hg38.vcf.gz
FASTA=/cta/users/burak.guler/miRNA/results/genome/GRCh38.p14.genome.fa
JOINT_VCF=/cta/users/burak.guler/miRNA/CES_DATALARI/vep_results/joint_germline_vep.vcf.gz
LABELS=/cta/users/burak.guler/miRNA/CES_DATALARI/ancestry-labels-1kg.tsv.txt
REF_DIR=$ANC/1kg-somalier
OUT_PREFIX=$ANC/cohort

cd "$ANC"
mkdir -p extracted

echo "Job started at $(date)"

# PART 1: extract cohort genotypes at ancestry sites for all samples
echo " somalier extract (full joint VCF, all 80 samples) "
"$SOMALIER" extract -d extracted/ --sites "$SITES" -f "$FASTA" "$JOINT_VCF"
echo "extracted cohort .somalier files: $(ls extracted/*.somalier 2>/dev/null | wc -l)"

# PART 2: ancestry prediction 
echo ""
echo " somalier ancestry (cohort vs 2,504 labeled 1000G samples) "
"$SOMALIER" ancestry --labels "$LABELS" -o "$OUT_PREFIX" \
    "$REF_DIR"/*.somalier ++ extracted/*.somalier

echo ""
echo " outputs "
ls -la "${OUT_PREFIX}".somalier-ancestry.* 2>/dev/null
echo ""
echo " per-sample predicted ancestry (cohort samples only) "
head -1 "${OUT_PREFIX}.somalier-ancestry.tsv" 2>/dev/null
grep -E "^SG" "${OUT_PREFIX}.somalier-ancestry.tsv" 2>/dev/null | head -85

echo ""
echo "Job finished at $(date)"
