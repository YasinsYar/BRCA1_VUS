# Zenodo Release Checklist

This repository is connected to GitHub at:

`https://github.com/YasinsYar/BRCA1_VUS`

Zenodo archives GitHub releases, not untracked local files. The release archive will contain only files committed to the release tag.

## Current Published Scope

The current `origin/main` commit contains only the pipeline scripts under `scripts/`.

Local files that were not published in the initial GitHub commit include:

- `.gitignore`
- manuscript `.docx` drafts
- `data_raw/`
- `data_int/`
- `results/`
- `outsource/`
- local Python/VEP environments

## v1.0.0 Release Scope

For the first DOI release, publish the repository as software plus the generated master table:

- source scripts under `scripts/`
- `README.md`
- `requirements.txt`
- `.gitignore`
- `CITATION.cff`
- `.zenodo.json`
- `LICENSE`
- `LICENSE-DATA.md`
- `CHANGELOG.md`
- `results/README.md`
- `results/tables/master.tsv`

Exclude:

- `.venv/`
- `env/`
- `data_raw/`
- `data_int/`
- manuscript drafts
- `outsource/`
- `results/tables/vep_brca1.tsv`
- `results/tables/vep_brca1.tsv_summary.html`
- generated summary tables and figures not selected for the first release

## Release Metadata

- Title: `AlphaFold-derived local structural confidence in BRCA1 missense variant interpretation: a ClinVar-based computational analysis`
- Creator: `Yasinskyi, Yaroslav`
- ORCID: `0000-0002-4801-8678`
- Version: `1.0.0`
- Git tag: `v1.0.0`
- Upload type: `software`
- Repository: `https://github.com/YasinsYar/BRCA1_VUS`
- Code license: MIT
- Generated result license: CC0 1.0 Universal
- Related manuscript: planned, not yet submitted; no related identifier is included.

## Files To Commit

The finalized metadata files are `CITATION.cff` and `.zenodo.json`. Zenodo will use `.zenodo.json` for the GitHub release metadata.

## Release Steps

1. Enable the GitHub repository in Zenodo.
2. Commit the finalized release files and `results/tables/master.tsv`.
3. Push `main` to GitHub.
4. Create an annotated tag:

```powershell
git tag -a v1.0.0 -m "Release v1.0.0"
git push origin v1.0.0
```

5. Create a GitHub release from the tag.
6. Wait for Zenodo to archive the release and assign a DOI.
7. In the Zenodo record, verify that the rights information reflects both the MIT code license and the CC0 license for generated result files where Zenodo's mixed-license interface is available.
8. Add the Zenodo DOI badge to `README.md` after the DOI is available.

## Post-release

After Zenodo assigns the DOI, add the DOI badge and DOI citation to `README.md`. If the manuscript is later published, add its DOI as a related identifier in the next release.
