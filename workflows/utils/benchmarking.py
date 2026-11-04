#!/usr/bin/env python
"""
Combine Snakemake benchmark files into one run-level table.

Each rule writes <bench_dir>/<rule>_<wildcards>.tsv; this collects them into a
single TSV with one row per (rule, job), wall time as hh:mm:ss.

Adapted from PhanstielLab/bagPipes workflows/utils/benchmarking.py.

Usage:
    python workflows/utils/benchmarking.py <bench_dir> <out.tsv>
"""
import glob
import os
import sys

import pandas as pd


def hms(seconds):
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    return "%d:%02d:%02d" % (hours, minutes, sec)


def main(bench_dir, out_path):
    tables = []
    for f in sorted(glob.glob(os.path.join(bench_dir, "*.tsv"))):
        table = pd.read_csv(f, sep="\t")
        job = os.path.splitext(os.path.basename(f))[0]
        table["rule"] = job.split("_")[0]
        table["job"] = "_".join(job.split("_")[1:])
        tables.append(table)
    if not tables:
        sys.exit(f"No benchmark files found in {bench_dir}")
    df = pd.concat(tables, ignore_index=True)
    agg = df.groupby(["rule", "job"]).agg("max", numeric_only=True)
    agg["s"] = agg["s"].apply(hms)
    agg.to_csv(out_path, sep="\t")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
