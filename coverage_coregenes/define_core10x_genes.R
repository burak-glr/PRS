#!/usr/bin/env Rscript

# Find high coverage gene (10x) set to get reliable results

# core10x = genes whose mean exon read-depth is >= 10x in all samples

# Input : output of run_mosdepth_pergene.slurm
#         columns: "CHROM" "START" "END" "GENE" "MEAN_DP"

# Output: core10x_4947genes.txt  — the 4,947 core10x gene symbols
#         core10x_genes.bed       — exon regions restricted to those genes

suppressMessages(library(dplyr))

cov_dir  <- "/cta/users/burak.guler/miRNA/CES_DATALARI/mosdepth_pergene"
gene_bed <- file.path(cov_dir, "coding_genes_exons.bed")
out_list <- file.path(cov_dir, "core10x_4947genes.txt")
out_bed  <- file.path(cov_dir, "core10x_genes.bed")

files <- list.files(cov_dir, pattern="\\.regions\\.bed\\.gz$", full.names=TRUE)
cat("Samples found:", length(files), "\n")

# The average depth of each gene is calculated for every person
per <- lapply(files, function(f) {
  d <- read.table(f, header=FALSE, sep="\t",
                  col.names=c("CHROM","START","END","GENE","MEAN_DP"))
  g <- d %>% group_by(GENE) %>% summarise(dp=mean(MEAN_DP, na.rm=TRUE), .groups="drop")
  setNames(g$dp, g$GENE)
})

genes <- sort(Reduce(union, lapply(per, names)))

# matrix (M) : rows = gene, cols = sample

M <- sapply(per, function(v) v[genes]); rownames(M) <- genes
cat("Genes x samples:", nrow(M), "x", ncol(M), "\n")

# filtering genes that have >=10x depth
core <- genes[ apply(M, 1, function(x) all(!is.na(x) & x >= 10)) ]
cat("core10x genes (>=10x in all", length(files), "samples):", length(core), "\n")

writeLines(core, out_list)

# Create BED file containing the core10x genes 
bed <- read.table(gene_bed, header=FALSE, sep="\t",
                  col.names=c("chr","start","end","gene"))
bed <- bed[bed$gene %in% core, ]
write.table(bed, out_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("Wrote:", out_list, "(", length(core), "genes ) and", out_bed,
    "(", nrow(bed), "exon rows )\n")
