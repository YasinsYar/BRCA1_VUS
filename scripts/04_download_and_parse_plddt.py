"""Download AlphaFold BRCA1 model and extract per-residue pLDDT values."""
import requests
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser

UNIPROT = "P38398"
API_URL = f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}"
OUT_PDB = Path("data_raw/alphafold") / f"AF-{UNIPROT}-F1.pdb"
OUT_NPY = Path("data_int/brca1_plddt_by_residue.npy")


def fetch_urls():
    r = requests.get(API_URL, timeout=60)
    r.raise_for_status()
    data = r.json()
    if not data:
        raise RuntimeError("AlphaFold API returned empty response")
    entry = data[0]
    urls = []
    for k in ["pdbUrl", "pdb_url", "cifUrl", "cif_url"]:
        if k in entry and entry[k]:
            urls.append(entry[k])
    if not urls:
        raise RuntimeError("No pdb/cif URL found in AlphaFold API response")
    return urls


def download(urls, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        print("File exists:", path)
        return
    last_err = None
    for url in urls:
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            path.write_bytes(r.content)
            print("Downloaded:", path)
            return
        except Exception as e:
            last_err = e
    raise RuntimeError(f"All download attempts failed: {last_err}")


def choose_chain(structure):
    chains = list(structure.get_chains())
    if not chains:
        raise RuntimeError("No chains found in structure")
    def count_res(chain):
        return sum(1 for r in chain if r.id[0] == " ")
    chains.sort(key=count_res, reverse=True)
    return chains[0]


def main():
    Path("data_int").mkdir(parents=True, exist_ok=True)
    urls = fetch_urls()
    download(urls, OUT_PDB)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("brca1", str(OUT_PDB))
    chain = choose_chain(structure)
    residues = [r for r in chain if r.id[0] == " "]
    if not residues:
        raise RuntimeError("No standard residues found")
    max_res = max(r.id[1] for r in residues)
    plddt = np.full(max_res + 1, np.nan, dtype=float)
    for res in residues:
        i = res.id[1]
        if "CA" in res:
            plddt[i] = float(res["CA"].bfactor)
        else:
            vals = [a.bfactor for a in res.get_atoms()]
            plddt[i] = float(np.mean(vals)) if vals else np.nan
    np.save(OUT_NPY, plddt)
    print("Saved:", OUT_NPY, "L=", max_res)


if __name__ == "__main__":
    main()
