Write-Host "Download VEP cache and FASTA into env/vep_data using Docker."
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$VepData = Join-Path $Root "env\vep_data"
New-Item -ItemType Directory -Force -Path $VepData | Out-Null
$Image = "ensemblorg/ensembl-vep:release_115.2"
$DockerRoot = docker info --format '{{.DockerRootDir}}'
docker rm -f vep_cache_install 2>$null | Out-Null
docker pull $Image
docker run --name vep_cache_install --rm -v "${VepData}:/data" $Image INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data
