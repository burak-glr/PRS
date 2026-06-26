#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=base_qc_liftover
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/base_qc_liftover_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/base_qc_liftover_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate rplots

LIFTOVER=/cta/users/burak.guler/miniforge3/bin/liftOver

SUZUKI_ZIP=/cta/users/burak.guler/miRNA/CES_DATALARI/association/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip
CHAIN=/cta/users/burak.guler/miRNA/CES_DATALARI/association/hg19ToHg38.over.chain.gz
OUT_DIR=/cta/users/burak.guler/miRNA/CES_DATALARI/association/prs

mkdir -p "$OUT_DIR"

# BASE DATA QC  (Suzuki et al.,2024)

echo "=== BASE DATA QC ==="

BASE_RAW="$OUT_DIR/suzuki_raw.txt"
BASE_QC="$OUT_DIR/base_qc.txt"

if [ ! -f "$BASE_RAW" ]; then
    echo "Unzipping Suzuki sumstats..."
    unzip -p "$SUZUKI_ZIP" > "$BASE_RAW"
fi
echo "Input variants: $(($(wc -l < "$BASE_RAW") - 1))"

# PART 1: Remove ambiguous (A/T, C/G) SNPs, apply MAF filter (EAF 0.01–0.99),
#          remove rows with missing Beta or SE, then deduplicate by CHR:POS
awk 'NR==1 {print; next}
{
    ea=$3; nea=$4; maf=$7
    if ((ea=="A" && nea=="T") || (ea=="T" && nea=="A")) next
    if ((ea=="C" && nea=="G") || (ea=="G" && nea=="C")) next
    if (maf+0 < 0.01 || maf+0 > 0.99) next
    if ($5=="NA" || $6=="NA") next
    print
}' "$BASE_RAW" \
| sort -k1,1n -k2,2n \
| awk 'NR==1 {print; prev=""; next}
       { key=$1"_"$2; if (key != prev) print; prev=key }' \
> "$BASE_QC"

echo "  Variants remaining: $(($(wc -l < "$BASE_QC") - 1))"

# PART 2: Liftover hg19 → hg38
BASE_BED="$OUT_DIR/base_hg19.bed"
BASE_HG38_BED="$OUT_DIR/base_hg38_lifted.bed"
BASE_UNMAPPED="$OUT_DIR/base_hg38_unmapped.bed"
BASE_HG38="$OUT_DIR/base_hg38.txt"

if [ ! -f "$BASE_HG38" ]; then
    echo "Lifting over base data hg19 → hg38..."
    awk 'NR>1 {print "chr"$1"\t"$2-1"\t"$2"\t"NR}' "$BASE_QC" > "$BASE_BED"
    "$LIFTOVER" "$BASE_BED" "$CHAIN" "$BASE_HG38_BED" "$BASE_UNMAPPED"
    echo "  Mapped:   $(wc -l < "$BASE_HG38_BED")"
    echo "  Unmapped: $(grep -vc '^#' "$BASE_UNMAPPED" 2>/dev/null || echo 0)"

    # Attach hg38 CHR/POS as new leading columns (strip "chr" prefix)
    awk '{gsub("chr","",$1); print $4"\t"$1"\t"$3}' "$BASE_HG38_BED" \
        > "$OUT_DIR/liftover_map.txt"

    awk 'NR==FNR {chr[$1]=$2; pos[$1]=$3; next}
         FNR==1  {print "CHR_HG38\tPOS_HG38\t"$0; next}
         FNR in chr {print chr[FNR]"\t"pos[FNR]"\t"$0}
    ' "$OUT_DIR/liftover_map.txt" "$BASE_QC" > "$BASE_HG38"

    echo "  Liftover complete: $(($(wc -l < "$BASE_HG38") - 1)) variants"
fi

echo "Base QC + liftover done at $(date)"
