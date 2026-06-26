#!/usr/bin/env python

# Input : exonic_assoc.psam
#         clinical_data.csv

# Output: sex.txt

import pandas as pd, re

BASE  = "/cta/users/burak.guler/miRNA"
PSAM  = f"{BASE}/CES_DATALARI/association/exonic_assoc.psam"
EXCEL = f"{BASE}/klinik veriler ve miRNA dataları/DR.Tez grup listesi-ProfDrCananCacına.xlsx"
CLIN  = f"{BASE}/CES_DATALARI/clinical_data.csv"
OUT   = f"{BASE}/CES_DATALARI/association/sex.txt"

# sample IIDs
psam = pd.read_csv(PSAM, sep="\t")
iidcol = [c for c in psam.columns if "IID" in c][0]
iids = psam[iidcol].astype(str).tolist()

# Map miRNA kod to Cinsiyet information
xl = pd.read_excel(EXCEL, sheet_name="dr tez gruplar (80 kişi)", header=1)
xl = xl[xl["miRNA kod"].astype(str).str.match(r"S\d+", na=False)]
sexmap = dict(zip(xl["miRNA kod"].astype(str),
                  xl["Cinsiyet"].astype(str).str.upper().str.strip().map({"ERKEK": 1, "KADIN": 2})))

# NGS key -> kod (clinical_data)
cl = pd.read_csv(CLIN)
key2kod = dict(zip(cl["NGS"].str.replace("-", "_", regex=False), cl["kod"]))

with open(OUT, "w") as f:
    n = 0
    for iid in iids:
        key  = re.sub(r"_SG\d.*", "", iid)
        kod  = key2kod.get(key)
        code = sexmap.get(kod)
        if code is None:
            raise ValueError(f"no sex for {iid} (key={key}, kod={kod})")
        f.write(f"{iid}\t{int(code)}\n"); n += 1
print(f"sex.txt: {n} samples written")
