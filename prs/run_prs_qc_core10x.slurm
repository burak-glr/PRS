#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=prs_qc_core10x
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/prs_qc_core10x_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/prs_qc_core10x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00


# Target sample QC: remove heterozygosity (+/-3 SD) + KING relatedness (cutoff 0.125);


eval "$(conda shell.bash hook)"
conda activate rplots

PLINK2=/cta/users/burak.guler/miniforge3/bin/plink2
TARGET_FULL=/cta/users/burak.guler/miRNA/CES_DATALARI/association/exonic_assoc
KEEP_BED=/cta/users/burak.guler/miRNA/CES_DATALARI/association/core10x_keep_pos.bed
OUT_DIR=/cta/users/burak.guler/miRNA/CES_DATALARI/association/prs_core10x
BASE_HG38=/cta/users/burak.guler/miRNA/CES_DATALARI/association/prs/base_hg38.txt

mkdir -p "$OUT_DIR"
echo "Job started at $(date)"
[ -f "$BASE_HG38" ] || { echo "ERROR: base_hg38.txt missing"; exit 1; }

# Part 1: restrict exonic target to the curated-gene variant positions

TARGET="$OUT_DIR/exonic_core10x"
$PLINK2 --pfile "$TARGET_FULL" --extract bed0 "$KEEP_BED" \
    --make-pgen --out "$TARGET" --threads 8
echo "Variants after curated-gene restriction: $(grep -vc '^#' "${TARGET}.pvar")"

TARGET_QC="$OUT_DIR/target_qc"

# Part 2: remove ambiguous SNPs (A/T, C/G)

awk 'NR>1 && $1!~/^#/ {ref=$4;alt=$5
    if((ref=="A"&&alt=="T")||(ref=="T"&&alt=="A")) next
    if((ref=="C"&&alt=="G")||(ref=="G"&&alt=="C")) next
    print $3}' "${TARGET}.pvar" > "$OUT_DIR/non_ambiguous_snps.txt"
echo "Non-ambiguous SNPs: $(wc -l < "$OUT_DIR/non_ambiguous_snps.txt")"

# Part 3: heterozygosity outliers (+/-3 SD)

$PLINK2 --pfile "$TARGET" --extract "$OUT_DIR/non_ambiguous_snps.txt" \
    --het --out "$OUT_DIR/target_het" --threads 8
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
Rscript "$SCRIPTDIR/het_outliers.R" \
    "$OUT_DIR/target_het.het" \
    "$OUT_DIR/het_outliers.txt"

# Part 4: KING relatedness (cutoff 0.125)

$PLINK2 --pfile "$TARGET" --extract "$OUT_DIR/non_ambiguous_snps.txt" \
    --king-cutoff 0.125 --out "$OUT_DIR/target_king" --threads 8
KING_REMOVE="$OUT_DIR/target_king.king.cutoff.out.id"
[ -f "$KING_REMOVE" ] || touch "$KING_REMOVE"
echo "Relatedness outliers removed: $(grep -vc '^#' "$KING_REMOVE")"

# Part 5: combine samples to remove heterozygous and relatedness outliers

cat "$OUT_DIR/het_outliers.txt" <(grep -v '^#' "$KING_REMOVE") | sort -u > "$OUT_DIR/samples_to_remove.txt"
echo "Total samples to remove: $(wc -l < "$OUT_DIR/samples_to_remove.txt")"

# Part 6: QC filters, removing flagged samples

$PLINK2 --pfile "$TARGET" \
    --extract "$OUT_DIR/non_ambiguous_snps.txt" \
    --remove "$OUT_DIR/samples_to_remove.txt" \
    --chr 1-22 --maf 0.05 --geno 0.01 --mind 0.01 --hwe 1e-6 \
    --make-pgen --out "$TARGET_QC" --threads 8
echo "QC-clean target variants: $(grep -vc '^#' "${TARGET_QC}.pvar")"
echo "QC-clean target samples:  $(grep -vc '^#' "${TARGET_QC}.psam")"

# PART 7: match base and target data by using CHR and POS columns 

awk 'NR>1 && $1!~/^#/ {gsub("chr","",$1); print $1"_"$2}' "${TARGET_QC}.pvar" \
    > "$OUT_DIR/target_positions.txt"
awk 'NR==FNR {pos[$1]=1; next} FNR==1{print;next}
     {gsub("chr","",$1); k=$1"_"$2; if(k in pos) print}' \
    "$OUT_DIR/target_positions.txt" "$BASE_HG38" > "$OUT_DIR/base_exonic_overlap.txt"
    
echo "Overlap base ∩ target: $(($(wc -l < "$OUT_DIR/base_exonic_overlap.txt") - 1))"
echo "Job finished at $(date)"
