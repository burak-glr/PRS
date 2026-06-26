#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=plink2_exonic
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/plink2_exonic_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/plink2_exonic_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

eval "$(conda shell.bash hook)"
conda activate rplots

PLINK2=/cta/users/burak.guler/miniforge3/bin/plink2
BCFTOOLS=/cta/users/burak.guler/miniforge3/bin/bcftools
VEP_VCF=/cta/users/burak.guler/miRNA/CES_DATALARI/vep_results/joint_germline_vep.vcf.gz
SEX=/cta/users/burak.guler/miRNA/CES_DATALARI/association/sex.txt

OUT=/cta/users/burak.guler/miRNA/CES_DATALARI/association/exonic_assoc
EXONIC_VCF=/cta/users/burak.guler/miRNA/CES_DATALARI/association/exonic_variants.vcf.gz

echo "Job started at $(date)"

# PART 1: Extract coding + UTR variants from VEP VCF

echo "Extracting coding + UTR variants from VEP VCF..."
$BCFTOOLS view "$VEP_VCF" | \
  awk '
  /^#/ { print; next }
  {
    if (match($8, /CSQ=([^;]+)/, csq)) {
      n = split(csq[1], transcripts, ",")
      for (i=1; i<=n; i++) {
        split(transcripts[i], f, "|")
        if (f[25]=="YES" && (f[3]=="HIGH"     ||
                             f[3]=="MODERATE" ||
                             f[3]=="LOW"      ||
                             f[2] ~ /UTR/)) {
          print; break
        }
      }
    }
  }' | \
  $BCFTOOLS view -O z -o "$EXONIC_VCF"
$BCFTOOLS index -t "$EXONIC_VCF"

echo "Variants extracted by category:"
echo "  Total exonic+UTR: $($BCFTOOLS view -H $EXONIC_VCF | wc -l)"

# PART 2: Convert to PLINK2 binary format

echo ""
echo "Converting to PLINK2 format..."
$PLINK2 \
    --vcf "$EXONIC_VCF" \
    --make-pgen \
    --max-alleles 2 \
    --update-sex "$SEX" \
    --split-par hg38 \
    --out "$OUT" \
    --threads 8

echo "Variants in pgen: $(tail -n +2 ${OUT}.pvar | wc -l)"
echo "Samples in pgen:  $(tail -n +2 ${OUT}.psam | wc -l)"
echo ""
echo "Job finished at $(date)"
