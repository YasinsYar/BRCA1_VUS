"""Combine ClinVar, VEP, pLDDT, and UniProt domains into a master table."""
import re
import numpy as np
import pandas as pd
from io import StringIO
from pathlib import Path

CLINVAR = Path("data_int/brca1_clinvar_missense.tsv")
VEP_TSV = Path("results/tables/vep_brca1.tsv")
PLDDT_NPY = Path("data_int/brca1_plddt_by_residue.npy")
DOMAINS = Path("data_int/uniprot_brca1_domains.tsv")
OUT = Path("results/tables/master.tsv")


def plddt_bin(x):
    if np.isnan(x):
        return "NA"
    if x >= 90:
        return ">=90"
    if x >= 70:
        return "70-89.999"
    if x >= 50:
        return "50-69.999"
    return "<50"


def add_domain(pos, domdf):
    hits = domdf[(domdf.start <= pos) & (pos <= domdf.end)]
    if len(hits) == 0:
        return "None"
    return "|".join(hits.domain_group.unique())


def parse_location(loc):
    if pd.isna(loc):
        return None, None
    s = str(loc)
    parts = s.split(":")
    if len(parts) != 2:
        return None, None
    chrom = parts[0]
    pos_part = parts[1]
    pos = pos_part.split("-")[0]
    try:
        return chrom, int(pos)
    except:
        return None, None


def normalize_chrom(chrom):
    if chrom is None:
        return None
    s = str(chrom)
    if s.lower().startswith("chr"):
        return s[3:]
    return s


def main():
    clin = pd.read_csv(CLINVAR, sep="\t")
    with open(VEP_TSV, "r", encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith("##")]
    vep = pd.read_csv(StringIO("".join(lines)), sep="\t")
    if vep.columns[0].startswith("#"):
        vep.columns = [c.lstrip("#") for c in vep.columns]
    if "CANONICAL" in vep.columns:
        vep = vep[vep["CANONICAL"] == "YES"].copy()
    loc_parsed = vep["Location"].apply(parse_location)
    vep["chrom"] = loc_parsed.apply(lambda x: normalize_chrom(x[0]))
    vep["pos"] = loc_parsed.apply(lambda x: x[1])
    vep = vep[vep["chrom"].notnull() & vep["pos"].notnull()].copy()
    vep["alt"] = vep["Allele"].astype(str)
    clin["chrom"] = clin["chrom"].apply(normalize_chrom)
    m = clin.merge(vep, on=["chrom", "pos", "alt"], how="inner", suffixes=("", "_vep"))
    if "Protein_position" in m.columns:
        m["prot_pos"] = pd.to_numeric(m["Protein_position"].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    else:
        rx = re.compile(r"p\..*?(\d+)")
        m["prot_pos"] = m["HGVSp"].astype(str).str.extract(rx)[0].astype(float)
    plddt = np.load(PLDDT_NPY)
    def pick_plddt(p):
        if pd.isna(p):
            return np.nan
        i = int(p)
        if i <= 0 or i >= len(plddt):
            return np.nan
        return float(plddt[i])
    m["plddt"] = m["prot_pos"].apply(pick_plddt)
    m["plddt_bin"] = m["plddt"].apply(plddt_bin)
    dom = pd.read_csv(DOMAINS, sep="\t")
    if len(dom) == 0:
        m["uniprot_domain_group"] = "None"
    else:
        m["uniprot_domain_group"] = m["prot_pos"].apply(lambda p: add_domain(int(p), dom) if pd.notnull(p) else "NA")
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    m.to_csv(OUT, sep="\t", index=False)
    print("Saved:", OUT, "N=", len(m))


if __name__ == "__main__":
    main()
