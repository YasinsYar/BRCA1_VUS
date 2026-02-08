"""Create a VCF with BRCA1 missense SNVs for VEP annotation."""
import pandas as pd
import vcfpy
from pathlib import Path

IN_TSV = Path("data_int/brca1_clinvar_missense.tsv")
CLINVAR_VCF = Path("data_raw/clinvar/clinvar.vcf.gz")
OUT_VCF = Path("data_int/brca1_missense_for_vep.vcf")


def main():
    df = pd.read_csv(IN_TSV, sep="\t")
    wanted = set(zip(df.chrom.astype(str), df.pos.astype(int), df.ref, df.alt))
    reader = vcfpy.Reader.from_path(str(CLINVAR_VCF))
    Path("data_int").mkdir(parents=True, exist_ok=True)
    writer = vcfpy.Writer.from_path(str(OUT_VCF), reader.header)
    count = 0
    for rec in reader:
        for alt in rec.ALT:
            key = (str(rec.CHROM), int(rec.POS), rec.REF, alt.value)
            if key in wanted:
                new = vcfpy.Record(
                    CHROM=rec.CHROM,
                    POS=rec.POS,
                    ID=rec.ID,
                    REF=rec.REF,
                    ALT=[alt],
                    QUAL=rec.QUAL,
                    FILTER=rec.FILTER,
                    INFO=dict(rec.INFO),
                    FORMAT=rec.FORMAT,
                    calls=rec.calls
                )
                writer.write_record(new)
                count += 1
    writer.close()
    reader.close()
    print("VCF for VEP:", OUT_VCF, "records=", count)


if __name__ == "__main__":
    main()
