#!/usr/bin/env Rscript
# Build the PRSice covariate file
# #IID SEX AGE PC1..PC6.
# Usage: Rscript merge_covariates_pcs.R <covariates.txt> <pcs.eigenvec> <out.txt>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: merge_covariates_pcs.R <covar> <eigenvec> <out>")

covar <- read.table(args[1], header = TRUE, sep = "\t", comment.char = "")
colnames(covar)[1] <- "IID"

pcs <- read.table(args[2], header = TRUE, comment.char = "")
colnames(pcs)[1] <- "FID"
pcs <- pcs[, c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")]

merged <- merge(covar, pcs, by = "IID")
colnames(merged)[1] <- "#IID"
write.table(merged, args[3],
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat("Covariates + PCs written:", nrow(merged), "samples\n")
