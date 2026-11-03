#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path
from collections import defaultdict

def parse_region(region: str):
    # "chr4:1000-6000"
    chrom, rest = region.split(":")
    start_s, end_s = rest.replace(",", "").split("-")
    start = int(start_s)
    end = int(end_s)
    if end < start:
        raise ValueError("Region end must be >= start")
    return chrom, start, end

def open_cov(path: str):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

def main():
    ap = argparse.ArgumentParser(
        description="Merge multiple Bismark .cov(.gz) files within a region by summing m/u per site, then recompute %."
    )
    ap.add_argument("--region", default="chr4:1000-6000",
                    help='Region like "chr4:1000-6000" (default: chr4:1000-6000)')
    ap.add_argument("-o", "--out", default="merged.chr4_1000_6000.cov.gz",
                    help="Output .cov.gz path (default: merged.chr4_1000_6000.cov.gz)")
    ap.add_argument("--mode", choices=["within", "overlap"], default="within",
                    help=(
                        'Filtering mode. "within": keep rows fully inside region (start>=Rstart and end<=Rend). '
                        '"overlap": keep rows overlapping region (not (end<Rstart or start>Rend)). '
                        "Default: within."
                    ))
    ap.add_argument("--keep-zero", action="store_true",
                    help="Keep sites with total coverage 0 (normally none).")
    ap.add_argument("inputs", nargs="+",
                    help="Input Bismark coverage files (*.cov or *.cov.gz), 6 columns: chrom start end meth% m u")
    args = ap.parse_args()

    chrom, rstart, rend = parse_region(args.region)

    # key: (chrom, start, end) -> [m_sum, u_sum]
    sums = defaultdict(lambda: [0, 0])

    def in_region(s, e):
        if args.mode == "within":
            return (s >= rstart) and (e <= rend)
        else:  # overlap
            return not (e < rstart or s > rend)

    for cov_path in args.inputs:
        with open_cov(cov_path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                c = parts[0]
                if c != chrom:
                    continue
                try:
                    s = int(parts[1])
                    e = int(parts[2])
                    m = int(parts[4])
                    u = int(parts[5])
                except ValueError:
                    continue

                if not in_region(s, e):
                    continue

                key = (c, s, e)
                sums[key][0] += m
                sums[key][1] += u

    # write output (sorted by position)
    out_path = args.out
    with gzip.open(out_path, "wt") as out:
        for (c, s, e) in sorted(sums.keys(), key=lambda x: (x[0], x[1], x[2])):
            m_sum, u_sum = sums[(c, s, e)]
            n = m_sum + u_sum
            if n == 0 and not args.keep_zero:
                continue
            meth_pct = 0.0 if n == 0 else (100.0 * m_sum / n)
            # Format matches Bismark cov: chrom start end meth% m u
            out.write(f"{c}\t{s}\t{e}\t{meth_pct:.2f}\t{m_sum}\t{u_sum}\n")

if __name__ == "__main__":
    main()
