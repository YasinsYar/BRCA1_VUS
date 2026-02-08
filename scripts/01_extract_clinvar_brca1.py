"""Extract BRCA1 missense SNVs from ClinVar VCF into a clean TSV."""
import re
import pandas as pd
import vcfpy
from pathlib import Path

IN_VCF = Path("data_raw/clinvar/clinvar.vcf.gz")
OUT_TSV = Path("data_int/brca1_clinvar_missense.tsv")
SO_MISSENSE = "SO:0001583|missense_variant"


def normalize_clnsig(raw):
    if not raw:
        return []
    s = str(raw).replace(" ", "_")
    s = re.sub(r"[\\|,;/]", "|", s)
    parts = [p for p in s.split("|") if p]
    return parts


def pick_class(tokens):
    if not tokens:
        return ""
    if "Conflicting_interpretations_of_pathogenicity" in tokens:
        return ""
    if "Pathogenic" in tokens:
        return "Pathogenic"
    if "Likely_pathogenic" in tokens:
        return "Likely_pathogenic"
    if "Benign" in tokens:
        return "Benign"
    if "Likely_benign" in tokens:
        return "Likely_benign"
    if "Uncertain_significance" in tokens:
        return "Uncertain_significance"
    return ""


def to_str_list(v):
    if v is None:
        return ""
    if isinstance(v, list):
        return "|".join(str(x) for x in v)
    return str(v)


def main():
    Path("data_int").mkdir(parents=True, exist_ok=True)
    reader = vcfpy.Reader.from_path(str(IN_VCF))
    rows = []
    for rec in reader:
        geneinfo = to_str_list(rec.INFO.get("GENEINFO"))
        if "BRCA1:" not in geneinfo:
            continue
        mc = to_str_list(rec.INFO.get("MC"))
        if SO_MISSENSE not in mc:
            continue
        if len(rec.REF) != 1:
            continue
        alts = [a.value for a in rec.ALT]
        if not alts:
            continue
        clnsig_raw = to_str_list(rec.INFO.get("CLNSIG"))
        tokens = normalize_clnsig(clnsig_raw)
        cls = pick_class(tokens)
        if not cls:
            continue
        revstat = to_str_list(rec.INFO.get("CLNREVSTAT")).replace(" ", "_")
        for alt in alts:
            if len(alt) != 1:
                continue
            rows.append({
                "chrom": str(rec.CHROM),
                "pos": int(rec.POS),
                "ref": rec.REF,
                "alt": alt,
                "variation_id": rec.ID,
                "clnsig_raw": clnsig_raw,
                "class": cls,
                "revstat": revstat,
                "mc": mc
            })
    reader.close()
    df = pd.DataFrame(rows).drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print("Saved:", OUT_TSV, "N=", len(df))


if __name__ == "__main__":
    main()
