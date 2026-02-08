"""Download ClinVar GRCh38 VCF and index into data_raw/clinvar."""
import requests
from pathlib import Path

VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
TBI_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
OUT_DIR = Path("data_raw/clinvar")


def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        print("File exists:", path)
        return
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    print("Downloaded:", path)


def main():
    download(VCF_URL, OUT_DIR / "clinvar.vcf.gz")
    download(TBI_URL, OUT_DIR / "clinvar.vcf.gz.tbi")


if __name__ == "__main__":
    main()
