#!/usr/bin/env python3

import argparse
from reporting.report import Report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--var", nargs="+", required=True)
    parser.add_argument("--hc", nargs="+", required=True)
    parser.add_argument("--hg", nargs="+", required=True)
    parser.add_argument("--dp", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    report = Report(tsv_variants=args.var, tsv_haplocheck=args.hc, tsv_haplogrep=args.hg,
                    tsv_depth=args.dp, html_output=args.out)

    report.generate_report()





