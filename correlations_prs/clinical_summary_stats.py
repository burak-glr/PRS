#!/usr/bin/env python
# Descriptive summary statistics for the clinical parameters, by group

import pandas as pd, numpy as np
from scipy import stats
BD = "/cta/users/burak.guler/miRNA/CES_DATALARI"

cl = pd.read_csv(f"{BD}/clinical_data.csv")
clin = [c for c in cl.columns[4:28] if c != "Yas"] + ["Kreatinin", "URE"]   
grp = cl["Grup"].astype(str).str.lower()

def fps(p):  return "0" if p < 0.0005 else f"{p:.3g}"      
def fpd(p):  return "< 0.001" if p < 0.001 else f"{p:.4f}" 

rows = []
for c in clin:
    x = pd.to_numeric(cl[c], errors="coerce")
    h = x[grp == "hasta"].dropna(); k = x[grp == "kontrol"].dropna()
    if len(h) < 3 or len(k) < 3:
        continue
    shh = stats.shapiro(h).pvalue; shk = stats.shapiro(k).pvalue
    normal = (shh > 0.05) and (shk > 0.05)
    if normal:
        test = "t-test"; p = stats.ttest_ind(h, k, equal_var=False).pvalue
    else:
        test = "Mann-Whitney"; p = stats.mannwhitneyu(h, k).pvalue
    norm = f"G1:{fps(shh)} G2:{fps(shk)}"; pd_ = fpd(p)        # G1=hasta, G2=kontrol
    for g, v in [("hasta", h), ("kontrol", k)]:
        rows.append({"Degisken": c, "Grup": g,
                     "Ortalama ± Standart Sapma (Medyan)": f"{v.mean():.2f}±{v.std():.2f} (Med:{v.median():g})",
                     "Normallik_Testi_p": norm, "p_degeri": pd_, "Test_Tipi": test})

t = pd.DataFrame(rows)
t.to_csv(f"{BD}/clinical_corr/clinical_summary_stats.csv", index=False)
print(f"{len(t)} rows ({len(t)//2} variables x 2 groups) -> clinical_summary_stats.csv")
