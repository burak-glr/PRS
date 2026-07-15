#!/usr/bin/env python
# PRS (core10x) vs each miRNA — Spearman correlation, horizontal bar.
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.patches import Patch
BD = "/cta/users/burak.guler/miRNA/CES_DATALARI"

prs = pd.read_csv(f"{BD}/association/prs_core10x/T2D_PRS_core10x.best", sep=r"\s+")
prs = prs[prs.In_Regression == "Yes"].copy()
prs["id"] = prs.IID.str.replace(r"_SG.*", "", regex=True)
cl = pd.read_csv(f"{BD}/clinical_data.csv"); cl["id"] = cl["NGS"].str.replace("-", "_", regex=False)
mir = list(cl.columns[28:50])

rows = []
for m in mir:
    y = np.log2(pd.to_numeric(cl[m], errors="coerce"))
    d = prs[["id", "PRS"]].merge(pd.DataFrame({"id": cl["id"], m: y}), on="id")
    ok = d["PRS"].notna() & d[m].notna()
    r, p = stats.spearmanr(d["PRS"][ok], d[m][ok])
    rows.append((m, round(r, 2), p))
res = pd.DataFrame(rows, columns=["miRNA", "r", "p"]).sort_values("r")
res["sig"] = res["p"].round(3) < 0.05          

def plot(title, fname):
    fig, ax = plt.subplots(figsize=(8.5, 8))
    col = ["#E74C3C" if s else "#B0B0B0" for s in res.sig]
    ax.barh(range(len(res)), res.r, color=col, edgecolor="black", linewidth=0.4)
    for i, (r, p) in enumerate(zip(res.r, res.p)):
        ax.text(r + (0.012 if r >= 0 else -0.012), i, f"{r:.2f}", va="center",
                ha="left" if r >= 0 else "right", fontsize=8)
        ax.text(-0.015 if r >= 0 else 0.015, i, f"p={p:.3f}", va="center",
                ha="right" if r >= 0 else "left", fontsize=7, color="#444")
    ax.set_yticks(range(len(res))); ax.set_yticklabels(res.miRNA, fontsize=8)
    ax.axvline(0, color="black", lw=0.8); ax.set_xlabel("Spearman r", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(handles=[Patch(color="#B0B0B0", label="n.s."), Patch(color="#E74C3C", label="p<0.05")],
              loc="lower right", ncol=2)
    ax.set_xlim(res.r.min() - 0.16, res.r.max() + 0.10)
    plt.tight_layout(); plt.savefig(fname, dpi=160, bbox_inches="tight"); print("Saved", fname)

plot("PRS (X10 genes) vs miRNA - Spearman Correlation", f"{BD}/clinical_corr/PRS_miRNA_spearman.png")
