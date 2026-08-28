# Data provenance

| Resource | Identifier/release | Access date | Role in study |
| --- | --- | --- | --- |
| NCBI ClinVar | GRCh38 VCF `clinvar_20260201.vcf.gz` (`fileDate=2026-02-01`); https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260201.vcf.gz | 2026-02-04 | Source of BRCA1 missense variants and ClinVar Variation IDs |
| UniProtKB | Accession `P38398` (BRCA1_HUMAN); https://www.uniprot.org/uniprotkb/P38398 | 2026-02-04 | BRCA1 sequence and domain interval annotations |
| AlphaFold Protein Structure Database | Model `AF-P38398-F1`, file `AF-P38398-F1-model_v6.pdb` (HEADER date 01-AUG-25); https://alphafold.ebi.ac.uk/entry/P38398 | 2026-02-04 | Per-residue pLDDT from the PDB B-factor field |
| Ensembl VEP | Docker image `ensemblorg/ensembl-vep:release_115.2` | 2026-02-05 | Missense consequence, protein position, SIFT, and PolyPhen annotation |
| Zenodo | DOI `10.5281/zenodo.19789707` | n/a | Deposited pipeline and analysis-ready `master.tsv` |

Individual ClinVar Variation IDs are retained in the deposited analysis table. No novel genetic variation, protein sequence, or experimentally determined macromolecular structure data were generated in this study.
