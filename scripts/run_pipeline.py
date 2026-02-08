"""Run the BRCA1 VUS pipeline"""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_cmd(cmd, cwd=ROOT):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def run_py(script):
    run_cmd([sys.executable, str(ROOT / script)])


def run_vep():
    vep_data = (ROOT / "env" / "vep_data").resolve()
    if not vep_data.exists():
        raise RuntimeError("Missing env/vep_data. Run the VEP cache installer first.")
    image = "ensemblorg/ensembl-vep:release_115.2"
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{vep_data}:/data",
        "-v", f"{ROOT.resolve()}:/work",
        image,
        "vep",
        "--cache", "--offline", "--dir", "/data",
        "--species", "homo_sapiens", "--assembly", "GRCh38",
        "--input_file", "/work/data_int/brca1_missense_for_vep.vcf",
        "--output_file", "/work/results/tables/vep_brca1.tsv",
        "--tab", "--force_overwrite",
        "--symbol", "--canonical", "--protein", "--hgvs", "--numbers",
        "--sift", "b", "--polyphen", "b", "--domains"
    ]
    run_cmd(cmd)


def pause_if_needed(pause):
    if pause:
        input("Press Enter to continue...")


def main():
    steps = [
        ("download_clinvar", lambda: run_py("scripts/00_download_clinvar.py")),
        ("extract_brca1", lambda: run_py("scripts/01_extract_clinvar_brca1.py")),
        ("make_vep_vcf", lambda: run_py("scripts/02_make_vep_input_vcf.py")),
        ("run_vep", run_vep),
        ("plddt", lambda: run_py("scripts/04_download_and_parse_plddt.py")),
        ("uniprot_domains", lambda: run_py("scripts/05_uniprot_domains.py")),
        ("build_master", lambda: run_py("scripts/06_build_master_table.py")),
        ("stats", lambda: run_py("scripts/07_stats_and_figures.py"))
    ]
    step_names = [s[0] for s in steps]

    parser = argparse.ArgumentParser()
    parser.add_argument("--list", action="store_true", help="List steps")
    parser.add_argument("--only", choices=step_names, help="Run only one step")
    parser.add_argument("--start", choices=step_names, help="Start from this step")
    parser.add_argument("--end", choices=step_names, help="End at this step")
    parser.add_argument("--pause", action="store_true", help="Pause after each step")
    args = parser.parse_args()

    if args.list:
        for name in step_names:
            print(name)
        return

    if args.only:
        idx = step_names.index(args.only)
        steps[idx][1]()
        return

    start_idx = step_names.index(args.start) if args.start else 0
    end_idx = step_names.index(args.end) if args.end else len(steps) - 1
    if start_idx > end_idx:
        raise RuntimeError("start must be before end")

    for i in range(start_idx, end_idx + 1):
        name, fn = steps[i]
        print("Step:", name)
        fn()
        pause_if_needed(args.pause)


if __name__ == "__main__":
    main()
