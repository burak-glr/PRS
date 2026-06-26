#!/bin/bash
# Build core10x_keep_pos.bed : the variant positions used to restrict the PRS.
# Input : exonic_variants.vcf.gz  (from variant_calling/run_plink2_exonic.slurm)
#         core10x_4947genes.txt   (from coverage_coregenes/define_core10x_genes.R)
# Output: core10x_keep_pos.bed    (chr  pos-1  pos ; 58,246 positions; chr prefix stripped)
BCF=/cta/users/burak.guler/miniforge3/bin/bcftools
EXONIC_VCF=/cta/users/burak.guler/miRNA/CES_DATALARI/association/exonic_variants.vcf.gz
CORELIST=/cta/users/burak.guler/miRNA/CES_DATALARI/mosdepth_pergene/core10x_4947genes.txt
OUT=/cta/users/burak.guler/miRNA/CES_DATALARI/association/core10x_keep_pos.bed

$BCF view -H "$EXONIC_VCF" | awk -v cl="$CORELIST" '
  BEGIN{ while((getline g < cl) > 0) core[g]=1 }
  match($8, /CSQ=([^;]+)/, csq){
    n=split(csq[1], tr, ","); sym=""
    for(i=1;i<=n;i++){ split(tr[i], f, "|"); if(f[25]=="YES"){ sym=f[4]; break } }
    if(sym in core){ c=$1; sub(/^chr/,"",c); print c"\t"($2-1)"\t"$2 }
  }' | sort -k1,1 -k2,2n | uniq > "$OUT"

echo "core10x_keep_pos.bed: $(wc -l < "$OUT") positions"
