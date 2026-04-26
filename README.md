# AlphaFold-derived local structural confidence in BRCA1 missense variant interpretation: a ClinVar-based computational analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19789707.svg)](https://doi.org/10.5281/zenodo.19789707)

This repository contains a reproducible pipeline for analysing ClinVar BRCA1 missense variants together with Ensembl VEP annotations, AlphaFold per-residue pLDDT values, and UniProt domain annotations.

Version `v1.0.0` is archived on Zenodo with DOI [`10.5281/zenodo.19789707`](https://doi.org/10.5281/zenodo.19789707).

## Repository Contents

- `scripts/00_download_clinvar.py`: downloads the ClinVar GRCh38 VCF and index.
- `scripts/01_extract_clinvar_brca1.py`: filters BRCA1 missense SNVs and selected ClinVar classes.
- `scripts/02_make_vep_input_vcf.py`: builds a minimal VCF for Ensembl VEP.
- `scripts/03_run_vep.ps1`: runs VEP through Docker.
- `scripts/04_download_and_parse_plddt.py`: downloads the AlphaFold BRCA1 model and extracts pLDDT by residue.
- `scripts/05_uniprot_domains.py`: downloads UniProt BRCA1 domain intervals.
- `scripts/06_build_master_table.py`: merges ClinVar, VEP, pLDDT, and domain annotations.
- `scripts/07_stats_and_figures.py`: generates summary tables, statistical tests, and figures.
- `scripts/run_pipeline.py`: orchestrates the full workflow.

## Data Sources

- ClinVar GRCh38 VCF from the official NCBI FTP endpoint.
- Ensembl VEP Docker image: `ensemblorg/ensembl-vep:release_115.2`.
- AlphaFold BRCA1 model: `AF-P38398-F1`.
- UniProt BRCA1 accession: `P38398`.

For the result files currently present in the local workspace, the methods note reports ClinVar `fileDate=2026-02-01`, with AlphaFold and UniProt accessed on `2026-02-04`.

## Requirements

- Python 3.13.5 was used in the local environment.
- Docker is required for VEP.
- PowerShell is used by the setup scripts.
- Python dependencies are listed in `requirements.txt`.

## Setup

From the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/00_setup_python_env.ps1
powershell -ExecutionPolicy Bypass -File scripts/00_setup_vep_cache.ps1
```

Alternatively, install Python dependencies into an existing environment:

```powershell
python -m pip install -r requirements.txt
```

## Running the Pipeline

List available steps:

```powershell
python scripts/run_pipeline.py --list
```

Run the full workflow:

```powershell
python scripts/run_pipeline.py
```

Run a single step:

```powershell
python scripts/run_pipeline.py --only stats
```

Run a range of steps:

```powershell
python scripts/run_pipeline.py --start build_master --end stats
```

## Outputs

The pipeline writes intermediate files to `data_int/`, raw inputs to `data_raw/`, and final tables and figures to `results/`.

Large raw inputs, VEP cache files, and the large derived VEP table are intentionally ignored by git because they are reproducible from public sources and local setup scripts.

The DOI release scope is source code plus `results/tables/master.tsv`. Additional local result tables and figures can be regenerated from the pipeline and are not required for the first Zenodo release.

## License

Code is released under the MIT License. Generated result files, including `results/tables/master.tsv`, are released under CC0 1.0 Universal; see `LICENSE-DATA.md`.

## Citation and DOI

Citation metadata are provided in `CITATION.cff` and Zenodo release metadata are provided in `.zenodo.json`.

Please cite the archived release:

Yasinskyi, Y. (2026). AlphaFold-derived local structural confidence in BRCA1 missense variant interpretation: a ClinVar-based computational analysis (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19789707
