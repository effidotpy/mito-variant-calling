#!/usr/bin/env python3

import argparse
from pathlib import Path

from reporting.tsv import TsvVariants


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_in")
    parser.add_argument("--tsv_out")
    args = parser.parse_args()

    parsed_records = TsvVariants.parse_vcf(vcf_path=Path(args.vcf_in))
    report = TsvVariants(records=parsed_records)
    report.to_tsv(output_tsv=Path(args.tsv_out))


