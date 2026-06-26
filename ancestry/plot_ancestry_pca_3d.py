#!/usr/bin/env python
# 3D ancestry PCA (PC1/PC2/PC3): cohort projected onto 1000G AIM genotypes.
# Interactive rotatable HTML
import pandas as pd
import plotly.graph_objects as go

D="/cta/users/burak.guler/miRNA/CES_DATALARI/ancestry"
a=pd.read_csv(f"{D}/cohort.somalier-ancestry.tsv",sep="\t")
a=a.rename(columns={a.columns[0]:"sample_id"})
pops=["AFR","AMR","EAS","EUR","SAS"]
ref=a[a["given_ancestry"].isin(pops)].copy()
coh=a[a["sample_id"].str.startswith("SG")].copy()
print(f"Reference={len(ref)} Cohort={len(coh)}")

v=ref[["PC1","PC2","PC3","PC4","PC5"]].var(); vp=(100*v/v.sum()).round(1)
ax=lambda p:f"{p} ({vp[p]:.1f}%)"
cols={"AFR":"#E74C3C","AMR":"#F39C12","EAS":"#27AE60","EUR":"#2980B9","SAS":"#8E44AD"}

# interactive 3D HTML 
fig=go.Figure()
for p in pops:
    s=ref[ref["given_ancestry"]==p]
    fig.add_trace(go.Scatter3d(x=s.PC1,y=s.PC2,z=s.PC3,mode="markers",
        name=f"1000G {p}",marker=dict(size=2.5,color=cols[p],opacity=0.45)))
fig.add_trace(go.Scatter3d(x=coh.PC1,y=coh.PC2,z=coh.PC3,mode="markers",
    name="Cohort (n=80)",marker=dict(size=4.5,color="black",
    symbol="diamond",line=dict(color="white",width=1))))
fig.update_layout(title="3D ancestry PCA - cohort projected onto 1000G AIM genotypes",
    scene=dict(xaxis_title=ax("PC1"),yaxis_title=ax("PC2"),zaxis_title=ax("PC3")),
    width=1000,height=800,legend=dict(itemsizing="constant"))
fig.write_html(f"{D}/ancestry_pca_3d.html")
print("Saved: ancestry_pca_3d.html (rotatable)")
