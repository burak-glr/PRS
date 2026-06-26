# miRNA x clinical-parameter Spearman correlations

setwd("/cta/users/burak.guler/miRNA/CES_DATALARI")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
dir.create("clinical_corr", showWarnings = FALSE)

# Load clinical data
cl <- read.csv("clinical_data.csv", check.names = FALSE, fileEncoding = "UTF-8")
names(cl)[27] <- "HBA1C_pct"        # was HBA1C%
names(cl)[28] <- "VitaminD"         # was VİTAMİND
cl$id <- gsub("-", "_", cl$NGS)

clin_vars <- c(setdiff(names(cl)[5:28], "Yas"), "Kreatinin", "URE") 
mir_vars  <- names(cl)[29:50]       
for (v in c(clin_vars, mir_vars)) cl[[v]] <- suppressWarnings(as.numeric(cl[[v]]))

for (v in mir_vars) cl[[v]] <- log2(cl[[v]])

prs <- read.table("association/prs_core10x/T2D_PRS_core10x.best", header = TRUE)
prs <- prs[prs$In_Regression == "Yes", ]
prs$id <- sub("_SG.*", "", prs$IID)
m <- merge(prs[, c("id","PRS")], cl, by = "id")
cat("Matched individuals:", nrow(m), "\n")

sp <- function(x, y){ ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(r=NA, p=NA, n=sum(ok)))
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(r = unname(t$estimate), p = t$p.value, n = sum(ok)) }

# miRNA x clinical
B <- expand.grid(miRNA = mir_vars, clinical = clin_vars, stringsAsFactors = FALSE)
B[c("r","p","n")] <- t(mapply(function(mi, cv) sp(m[[mi]], m[[cv]]), B$miRNA, B$clinical))
write.csv(B[order(B$p), c("miRNA","clinical","r","p","n")],
          "clinical_corr/miRNA_clinical_spearman.csv", row.names = FALSE)
cat(sprintf("\n miRNA x clinical: %d pairs, %d significant (raw p<0.05) \n",
            nrow(B), sum(B$p < 0.05, na.rm = TRUE)))
cat("Top 12 by p:\n")
print(head(B[order(B$p), c("miRNA","clinical","r","p","n")], 12), row.names = FALSE)

# heatmap (nominal-p stars) - FDR version: mirna_clinical_heatmap_fdr.py
nlab <- sapply(clin_vars, function(v) sum(is.finite(m[[v]])))
B$clinical_n <- factor(sprintf("%s (n=%d)", B$clinical, nlab[B$clinical]),
                       levels = sprintf("%s (n=%d)", clin_vars, nlab[clin_vars]))
B$star <- ifelse(!is.na(B$p) & B$p < 0.05, "*", "")
pB <- ggplot(B, aes(clinical_n, miRNA, fill = r)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = star), size = 4, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.6, 0.6), na.value = "grey95") +
  labs(title = "miRNA vs clinical markers - Spearman r  (* = p<0.05)",
       subtitle = "n per marker shown on x-axis; pairwise-complete",
       x = NULL, y = NULL, fill = "Spearman r") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9))
ggsave("clinical_corr/miRNA_clinical_heatmap.png", pB, width = 11.5, height = 8, dpi = 150)

cat("\nSaved: miRNA_clinical_spearman.csv, miRNA_clinical_heatmap.png\n")
