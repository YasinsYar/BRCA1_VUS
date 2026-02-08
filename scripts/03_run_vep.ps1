Write-Host "Run VEP annotation in Docker for BRCA1 missense variants."
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$VepData = Join-Path $Root "env\vep_data"
$Input = "/work/data_int/brca1_missense_for_vep.vcf"
$Output = "/work/results/tables/vep_brca1.tsv"
$Image = "ensemblorg/ensembl-vep:release_115.2"

docker run --rm `
  -v "${VepData}:/data" `
  -v "${Root}:/work" `
  $Image `
  vep `
    --cache --offline --dir /data `
    --species homo_sapiens --assembly GRCh38 `
    --input_file $Input `
    --output_file $Output `
    --tab --force_overwrite `
    --symbol --canonical --protein --hgvs --numbers `
    --sift b --polyphen b --domains
