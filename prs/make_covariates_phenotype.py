#!/usr/bin/env python
#Build covariates and phenotype data
# Input: SEX  <- Excel "Cinsiyet" (ERKEK=1, KADIN=2)
#        AGE  <- clinical_data "Yas"
#        PHENO<- clinical_data "Grup" (hasta=2, kontrol=1)

# Build the PRSice covariate and phenotype files from cohort metadata.
#   covariates.txt : #IID  SEX  AGE     
#   phenotype.txt  : #IID  PHENO        

import pandas as pd, re

BASE  = "/cta/users/burak.guler/miRNA"
PSAM  = f"{BASE}/CES_DATALARI/association/exonic_assoc.psam"
EXCEL = f"{BASE}/klinik veriler ve miRNA dataları/DR.Tez grup listesi-ProfDrCananCacına.xlsx"
CLIN  = f"{BASE}/CES_DATALARI/clinical_data.csv"
COV   = f"{BASE}/CES_DATALARI/association/covariates.txt"
PHENO = f"{BASE}/CES_DATALARI/association/phenotype.txt"

psam = pd.read_csv(PSAM, sep="\t")
iidcol = [c for c in psam.columns if "IID" in c][0]
iids = psam[iidcol].astype(str).tolist()

xl = pd.read_excel(EXCEL, sheet_name="dr tez gruplar (80 kişi)", header=1)
xl = xl[xl["miRNA kod"].astype(str).str.match(r"S\d+", na=False)]
sexmap = dict(zip(xl["miRNA kod"].astype(str),
                  xl["Cinsiyet"].astype(str).str.upper().str.strip().map({"ERKEK": 1, "KADIN": 2})))

cl = pd.read_csv(CLIN)
cl["key"] = cl["NGS"].str.replace("-", "_", regex=False)
key2kod = dict(zip(cl["key"], cl["kod"]))
kod2age = dict(zip(cl["kod"], cl["Yas"]))
kod2grp = dict(zip(cl["kod"], cl["Grup"]))
phmap = {"hasta": 2, "kontrol": 1}

with open(COV, "w") as fc, open(PHENO, "w") as fp:
    fc.write("#IID\tSEX\tAGE\n"); fp.write("#IID\tPHENO\n")
    for iid in iids:
        kod = key2kod[re.sub(r"_SG\d.*", "", iid)]
        fc.write(f"{iid}\t{int(sexmap[kod])}\t{int(kod2age[kod])}\n")
        fp.write(f"{iid}\t{phmap[str(kod2grp[kod]).strip().lower()]}\n")
print("Wrote covariates.txt and phenotype.txt (80 samples)")
