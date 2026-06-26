#!/usr/bin/env Rscript

# Heterozygosity sample-QC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: het_outliers.R <het_file> <out_file>")

het <- read.table(args[1], header = TRUE, comment.char = "")
colnames(het)[1] <- "IID"
het$HET_RATE <- (het$OBS_CT - het$O.HOM.) / het$OBS_CT

m <- mean(het$HET_RATE); s <- sd(het$HET_RATE)
out <- het[abs(het$HET_RATE - m) > 3 * s, , drop = FALSE]
cat("Heterozygosity outliers (+/-3 SD):", nrow(out), "\n")

write.table(if (nrow(out) > 0) out[, "IID", drop = FALSE] else data.frame(IID = character(0)),
            args[2], row.names = FALSE, col.names = FALSE, quote = FALSE)
