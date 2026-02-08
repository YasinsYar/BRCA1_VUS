"""Fetch BRCA1 domain ranges from UniProt and save a domain table."""
import requests
import re
import pandas as pd
from pathlib import Path

ACC = "P38398"
URL = f"https://rest.uniprot.org/uniprotkb/{ACC}.json"
OUT = Path("data_int/uniprot_brca1_domains.tsv")

KEEP_PATTERNS = [
    ("RING", re.compile(r"RING|zinc finger", re.I)),
    ("BRCT", re.compile(r"BRCT", re.I))
]


def main():
    Path("data_int").mkdir(parents=True, exist_ok=True)
    r = requests.get(URL, timeout=60)
    r.raise_for_status()
    j = r.json()
    rows = []
    for ft in j.get("features", []):
        if ft.get("type") != "Domain":
            continue
        desc = ft.get("description", "") or ""
        loc = ft.get("location", {}) or {}
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        label = None
        for name, rx in KEEP_PATTERNS:
            if rx.search(desc):
                label = name
                break
        if label is None:
            continue
        rows.append({
            "domain_group": label,
            "description": desc,
            "start": int(start),
            "end": int(end)
        })
    df = pd.DataFrame(rows)
    if len(df) == 0:
        df = pd.DataFrame(columns=["domain_group", "description", "start", "end"])
    df.sort_values(["domain_group", "start"], inplace=True)
    df.to_csv(OUT, sep="\t", index=False)
    print("Saved:", OUT, "N=", len(df))


if __name__ == "__main__":
    main()
