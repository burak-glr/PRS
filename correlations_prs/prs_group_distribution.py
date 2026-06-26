#!/usr/bin/env python
# PRS distribution by group (T2D cases vs controls)
# Wilcoxon and t-test p-values
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
np.random.seed(0)                     
BD = "/cta/users/burak.guler/miRNA/CES_DATALARI/association/prs_core10x"

best = pd.read_csv(f"{BD}/T2D_PRS_core10x.best", sep=r"\s+"); best = best[best.In_Regression == "Yes"]
ph = pd.read_csv(f"{BD}/phenotype_nohash.txt", sep=r"\s+")
d = best.merge(ph, on="IID")
d["z"] = (d.PRS - d.PRS.mean()) / d.PRS.std()

def build(grp_case, grp_ctrl, title, fname):
    d["Group"] = np.where(d.PHENO == 2, grp_case, grp_ctrl)
    order = [grp_ctrl, grp_case]
    pal = {grp_ctrl: "#AED6F1", grp_case: "#F1948A"}; dotpal = {grp_ctrl: "#2980B9", grp_case: "#C0392B"}
    t = d.z[d.Group == grp_case]; c = d.z[d.Group == grp_ctrl]
    pw = stats.mannwhitneyu(t, c).pvalue; pt = stats.ttest_ind(t, c, equal_var=False).pvalue
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.violinplot(data=d, x="Group", y="z", order=order, hue="Group", palette=pal,
                   legend=False, inner=None, cut=1, linewidth=1.2, ax=ax)
    for v in ax.collections: v.set_alpha(0.55)
    sns.stripplot(data=d, x="Group", y="z", order=order, hue="Group", palette=dotpal,
                  legend=False, size=5, alpha=0.8, jitter=0.12, ax=ax)
    for i, g in enumerate(order):
        x = d.z[d.Group == g]; m = x.mean(); ci = 1.96 * x.std() / np.sqrt(len(x))
        ax.errorbar(i, m, yerr=ci, fmt="o", color="black", ms=8, capsize=5, zorder=10)
    y = d.z.max() + 0.35; ax.plot([0, 0, 1, 1], [y, y + 0.12, y + 0.12, y], color="black", lw=1.2)
    ax.text(0.5, y + 0.18, f"Wilcoxon p = {pw:.3f}   |   t-test p = {pt:.3f}", ha="center", fontsize=11)
    ax.set_ylabel("Standardised PRS (Z-score)", fontsize=12); ax.set_xlabel("")
    ax.set_title(title, fontsize=13, fontweight="bold"); ax.set_ylim(d.z.min() - 0.4, y + 0.6)
    sns.despine(); plt.tight_layout(); plt.savefig(fname, dpi=160, bbox_inches="tight")
    print(f"Saved {fname}  (Wilcoxon p={pw:.3f}, t-test p={pt:.3f})")

# Shapiro-Wilk normality per group (PRS expected ~normal by CLT)
for g, lab in [(2, "T2D"), (1, "Control")]:
    print(f"Shapiro-Wilk {lab}: p={stats.shapiro(d.z[d.PHENO == g]).pvalue:.3f}")

build("T2D Cases", "Controls", "PRS Distribution by Group  (10x depth genes)", f"{BD}/prs_group_violin_seaborn.png")
