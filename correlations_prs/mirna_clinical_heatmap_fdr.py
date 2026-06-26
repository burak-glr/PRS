#!/usr/bin/env python
# miRNA x clinical-parameter Spearman heatmap with FDR (Benjamini-Hochberg)

import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from statsmodels.stats.multitest import multipletests
BD = "/cta/users/burak.guler/miRNA/CES_DATALARI/clinical_corr"

d = pd.read_csv(f"{BD}/miRNA_clinical_spearman.csv").dropna(subset=["p"])
d["q"] = multipletests(d["p"].values, method="fdr_bh")[1]
R = d.pivot(index="miRNA", columns="clinical", values="r")
Q = d.pivot(index="miRNA", columns="clinical", values="q")
# glycemic + kidney markers first, then the rest; miRNAs by mean |rho|
front = ["Glukoz", "HBA1C_pct", "Kreatinin", "URE"]
order_cl = [c for c in front if c in R.columns] + [c for c in R.columns if c not in front]
R = R[order_cl]; Q = Q[order_cl]
R = R.loc[R.abs().mean(axis=1).sort_values(ascending=False).index]; Q = Q.loc[R.index]

def star(q):
    return "" if pd.isna(q) else ("***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else "")))
ann = np.vectorize(star)(Q.values)

def draw(title, fname):
    fig, ax = plt.subplots(figsize=(15, 9))
    sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6, annot=ann, fmt="",
                annot_kws={"size": 11, "weight": "bold"}, linewidths=0.5, linecolor="white",
                cbar_kws={"label": "Spearman rho", "shrink": 0.6}, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title(title + "\n* q<0.05  ** q<0.01  *** q<0.001", fontsize=13, fontweight="bold", pad=12)
    plt.tight_layout(); plt.savefig(fname, dpi=160, bbox_inches="tight"); print("Saved", fname)

draw("miRNA vs Clinical Parameters", f"{BD}/miRNA_clinical_heatmap_fdr.png")
