# PRS (core10x, EUR weights) vs clinical parameters - Spearman correlation.
setwd("/cta/users/burak.guler/miRNA/CES_DATALARI")
suppressMessages({library(ggplot2)})

prs <- read.table("association/prs_core10x/T2D_PRS_core10x.best", header=TRUE)
prs <- prs[prs$In_Regression=="Yes",]; prs$id <- sub("_SG.*","",prs$IID)
cl  <- read.csv("clinical_data.csv", check.names=FALSE); cl$id <- gsub("-","_",cl$NGS)
names(cl)[27] <- "HBA1C_pct"
clin <- c(setdiff(names(cl)[5:28], "Yas"), "Kreatinin", "URE")   # Yas excl + Kreatinin + URE
for(v in clin) cl[[v]] <- suppressWarnings(as.numeric(cl[[v]]))
d <- merge(prs[,c("id","PRS")], cl[,c("id","Grup",clin)], by="id")
cat("N merged =", nrow(d), "| cases =", sum(d$Grup=="hasta"), "| controls =", sum(d$Grup=="kontrol"), "\n\n")

sp <- function(sub, tag){
  do.call(rbind, lapply(clin, function(v){
    x<-sub$PRS; yv<-sub[[v]]; ok<-!is.na(x)&!is.na(yv); nn<-sum(ok)
    if(nn<8) return(NULL)
    h<-suppressWarnings(cor.test(x[ok],yv[ok],method="spearman"))
    data.frame(set=tag, clinical=v, n=nn, rho=round(unname(h$estimate),3),
               p=signif(h$p.value,3)) }))
}
all <- sp(d, "all(79)")
cas <- sp(d[d$Grup=="hasta",], "cases")
res <- rbind(all, cas)
resA <- all[order(all$p),]
cat("=== PRS vs clinical - WHOLE COHORT (sorted by p) ===\n")
print(resA, row.names=FALSE)
cat(sprintf("\nSignificant (p<0.05), whole cohort: %d / %d\n", sum(resA$p<0.05), nrow(resA)))
write.csv(res, "clinical_corr/prs_vs_clinical.csv", row.names=FALSE)

# barplot of rho, whole cohort, flag significant
resA$sig <- ifelse(resA$p<0.05, "p<0.05", "n.s.")
resA$clinical <- factor(resA$clinical, levels=resA$clinical[order(resA$rho)])
resA$lab <- sprintf("n=%d", resA$n)
g <- ggplot(resA, aes(rho, clinical, fill=sig)) + geom_col(width=0.7) +
  geom_vline(xintercept=0, color="grey40") +
  geom_text(aes(label=lab, hjust=ifelse(rho>=0, -0.15, 1.15)), size=2.9, color="grey20") +
  scale_fill_manual(values=c("p<0.05"="#E74C3C","n.s."="grey75")) +
  scale_x_continuous(expand=expansion(mult=c(0.18,0.18))) +
  labs(title="PRS vs Clinical Parameters",
       x="Spearman rho", y=NULL, fill=NULL) +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold",size=11,hjust=0.5))
ggsave("clinical_corr/prs_vs_clinical.png", g, width=8.5, height=6.5, dpi=150)
cat("\nSaved: prs_vs_clinical.csv, prs_vs_clinical.png\n")
