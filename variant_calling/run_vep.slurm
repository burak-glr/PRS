#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=vep_annotation
#SBATCH --output=logs/vep_%j.log
#SBATCH --error=logs/vep_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

eval "$(conda shell.bash hook)"

# Create VEP environment if it doesn't exist
if ! conda env list | grep -q "^vep "; then
    echo "Creating VEP conda environment..."
    conda create -y -n vep -c conda-forge -c bioconda ensembl-vep=110
fi

conda activate vep

VCF=/cta/users/burak.guler/miRNA/CES_DATALARI/sarek_results/variant_calling/haplotypecaller/joint_variant_calling/joint_germline.vcf.gz
CACHE_DIR=/cta/users/burak.guler/miRNA/CES_DATALARI/vep_cache
OUT_DIR=/cta/users/burak.guler/miRNA/CES_DATALARI/vep_results

mkdir -p "$CACHE_DIR" "$OUT_DIR"

# Download cache if not already present
if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
    echo "Downloading VEP cache for GRCh38..."
    vep_install \
        --AUTO cf \
        --SPECIES homo_sapiens \
        --ASSEMBLY GRCh38 \
        --CACHEDIR "$CACHE_DIR" \
        --NO_HTSLIB \
        --NO_UPDATE
fi

echo "Running VEP annotation..."

vep \
    --input_file "$VCF" \
    --output_file "${OUT_DIR}/joint_germline_vep.vcf" \
    --format vcf \
    --vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache \
    --dir_cache "$CACHE_DIR" \
    --cache_version 110 \
    --fork 16 \
    --everything \
    --stats_file "${OUT_DIR}/joint_germline_vep_summary.html" \
    --warning_file "${OUT_DIR}/vep_warnings.txt" \
    --force_overwrite

# Compress and index output
bgzip -f "${OUT_DIR}/joint_germline_vep.vcf"
tabix -f -p vcf "${OUT_DIR}/joint_germline_vep.vcf.gz"

echo "Annotation complete. Output: ${OUT_DIR}/joint_germline_vep.vcf.gz"
