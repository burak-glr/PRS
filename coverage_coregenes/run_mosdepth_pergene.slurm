#!/bin/bash
#SBATCH -A adelab
#SBATCH --qos=adelab
#SBATCH -p genomics
#SBATCH --job-name=mosdepth_pergene
#SBATCH --output=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/mosdepth_pergene_%j.log
#SBATCH --error=/cta/users/burak.guler/miRNA/CES_DATALARI/logs/mosdepth_pergene_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

eval "$(conda shell.bash hook)"
conda activate mapping

WD=/cta/users/burak.guler/miRNA/CES_DATALARI
GTF=/cta/users/burak.guler/miRNA/results/genome/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
OUT_DIR="$WD/mosdepth_pergene"
mkdir -p "$OUT_DIR"

echo "Job started at $(date)"

# Part 1: Create gene BED file by using CDS and UTR coordinates from Gencode vv44 GTF
GENE_BED="$OUT_DIR/coding_genes_exons.bed"

if [ ! -f "$GENE_BED" ]; then
    echo "Building coding gene exon BED from GTF..."
    awk '$3=="CDS" || $3=="UTR" {
        if ($1 !~ /^(chr[0-9]+|chrX|chrY|chrM)$/) next
        match($0, /gene_name "([^"]+)"/, g)
        if (g[1]=="") next
        print $1"\t"$4-1"\t"$5"\t"g[1]
    }' "$GTF" | sort -k1,1 -k2,2n | uniq > "$GENE_BED"
    echo "Coding exon intervals: $(wc -l < "$GENE_BED")"
    echo "Unique genes: $(awk '{print $4}' "$GENE_BED" | sort -u | wc -l)"
else
    echo "Reusing existing $GENE_BED"
fi

# Part 2: Run depth calculation for all samples

echo "Running mosdepth for all samples..."

find "$WD/sarek_bam_all/preprocessing" -name "*.bam" | sort | while read BAM; do
    SAMPLE=$(basename "$BAM" | sed 's/\..*$//')

    # Skip if already done
    if [ -f "$OUT_DIR/${SAMPLE}.regions.bed.gz" ]; then
        echo "  Skipping $SAMPLE (already done)"
        continue
    fi

    echo "  Processing $SAMPLE at $(date)"
    mosdepth \
        --by "$GENE_BED" \
        --no-per-base \
        --threads 8 \
        "$OUT_DIR/$SAMPLE" \
        "$BAM"
done

echo ""
echo " mosdepth output files "
ls "$OUT_DIR/"*.regions.bed.gz 2>/dev/null | wc -l
echo "samples processed"

echo "Job finished at $(date)"
