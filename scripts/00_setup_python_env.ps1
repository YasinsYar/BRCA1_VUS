Write-Host "Create a local Python virtual environment and install dependencies."
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$Venv = Join-Path $Root ".venv"
$PipCache = Join-Path $Root "env\pip_cache"
New-Item -ItemType Directory -Force -Path $PipCache | Out-Null
$env:PIP_CACHE_DIR = $PipCache
python -m venv $Venv
& "$Venv\Scripts\python.exe" -m pip install -U pip
& "$Venv\Scripts\pip.exe" install pandas numpy vcfpy requests biopython scipy statsmodels scikit-learn matplotlib tqdm
