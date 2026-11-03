#!/usr/bin/env python3
import argparse
import gzip
import math
from pathlib import Path

import matplotlib.pyplot as plt

def parse_region(region: str):
    chrom, rest = region.split(":")
    start_s, end_s = rest.replace(",", "").split("-")
    start = int(start_s)
    end = int(end_s)
    if end < start:
        raise ValueError("Region end must be >= start")
    return chrom, start, end

def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if Path(path).suffix == ".gz" else open(path, "rt")

def wilson_ci(k: int, n: int, z: float = 1.96):
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + (z*z)/n
    center = (p + (z*z)/(2*n)) / denom
    half = (z * math.sqrt((p*(1-p) + (z*z)/(4*n))/n)) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return lo, hi

def size_map(n: int, min_size: float, max_size: float, scale: float):
    s = math.log1p(max(0, n)) * scale
    return max(min_size, min(max_size, s))

def read_cov_region(cov_path: str, chrom: str, start: int, end: int, min_cov: int):
    # Bismark cov: chrom start end meth% m u
    xs, p_hat, cov, lo, hi = [], [], [], [], []
    with open_maybe_gzip(cov_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            if parts[0] != chrom:
                continue
            try:
                s = int(parts[1])
                e = int(parts[2])
                m = int(parts[4])
                u = int(parts[5])
            except ValueError:
                continue

            # overlap with region
            if e < start or s > end:
                continue

            n = m + u
            if n < min_cov:
                continue

            p = m / n
            ci_lo, ci_hi = wilson_ci(m, n)

            xs.append(s)
            p_hat.append(p)
            cov.append(n)
            lo.append(ci_lo)
            hi.append(ci_hi)

    order = sorted(range(len(xs)), key=lambda i: xs[i])
    xs   = [xs[i] for i in order]
    p_hat= [p_hat[i] for i in order]
    cov  = [cov[i] for i in order]
    lo   = [lo[i] for i in order]
    hi   = [hi[i] for i in order]
    return xs, p_hat, cov, lo, hi

def short_name(path: str):
    name = Path(path).name
    for suf in [".bismark.cov.gz", ".bismark.cov", ".cov.gz", ".cov"]:
        if name.endswith(suf):
            return name[:-len(suf)]
    return name

def main():
    ap = argparse.ArgumentParser(
        description="Overlay CpG methylation across samples with fixed legend panel on the right."
    )
    ap.add_argument("--region", required=True, help='Region like "chr4:1000-6000"')
    ap.add_argument("-o", "--out", default="methylation_overlay.png", help="Output png/pdf path.")
    ap.add_argument("--min-cov", type=int, default=1, help="Minimum coverage to plot (default: 1).")

    ap.add_argument("--alpha", type=float, default=0.45, help="Point alpha (default: 0.45).")
    ap.add_argument("--size-scale", type=float, default=18.0, help="Dot size scale (log1p) (default: 18).")
    ap.add_argument("--min-size", type=float, default=10.0, help="Minimum dot size (default: 10).")
    ap.add_argument("--max-size", type=float, default=70.0, help="Maximum dot size cap (default: 70).")

    ap.add_argument("--show-ci", action="store_true",
                    help="Show 95% CI vertical bars (can get busy).")
    ap.add_argument("--legend-cov-values", default="5,10,20,50,100",
                    help="Comma-separated coverage values for dot-size legend.")
    ap.add_argument("--legend-cols", type=int, default=0,
                    help="Columns for sample legend (0=auto).")
    ap.add_argument("--legend-width", type=float, default=1.6,
                    help="Relative width of legend panel (bigger if sample names are long).")

    ap.add_argument("inputs", nargs="+", help="Input Bismark coverage files (*.cov or *.cov.gz).")
    args = ap.parse_args()

    chrom, rstart, rend = parse_region(args.region)

    # ---- Fixed layout: main plot + legend panel ----
    fig_w = 12.0 + max(0.0, args.legend_width - 1.6) * 2.0
    fig = plt.figure(figsize=(fig_w, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.6, args.legend_width], wspace=0.05)

    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.axis("off")

    # plot
    for path in args.inputs:
        name = short_name(path)
        xs, p_hat, cov, lo, hi = read_cov_region(path, chrom, rstart, rend, args.min_cov)

        color = ax._get_lines.get_next_color()

        if args.show_ci:
            for x, l, h in zip(xs, lo, hi):
                ax.plot([x, x], [l, h], linewidth=1, alpha=0.15, color=color)

        sizes = [size_map(n, args.min_size, args.max_size, args.size_scale) for n in cov]
        ax.scatter(xs, p_hat, s=sizes, alpha=args.alpha, color=color, label=name)

    ax.set_xlim(rstart, rend)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(f"Genomic position in {chrom}:{rstart}-{rend} (bp)")
    ax.set_ylabel("Methylation proportion (m / (m + u))")
    ax.set_title("CpG-level methylation overlay (dot size = log(coverage), capped)")
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    # ---- Legends in right panel ----
    sample_handles, sample_labels = ax.get_legend_handles_labels()
    if args.legend_cols and args.legend_cols > 0:
        ncol = args.legend_cols
    else:
        ncol = 2 if len(sample_labels) >= 10 else 1

    leg1 = ax_leg.legend(
        sample_handles, sample_labels,
        title="Samples (color)",
        loc="upper left",
        frameon=True,
        ncol=ncol,
        fontsize=8,
        title_fontsize=9,
        borderaxespad=0.0
    )
    ax_leg.add_artist(leg1)

    # size legend
    try:
        cov_values = [int(x.strip()) for x in args.legend_cov_values.split(",") if x.strip()]
    except ValueError:
        cov_values = [5, 10, 20, 50, 100]

    size_handles = []
    size_labels = []
    for n in cov_values:
        s = size_map(n, args.min_size, args.max_size, args.size_scale)
        h = ax_leg.scatter([], [], s=s, alpha=args.alpha, color="black")
        size_handles.append(h)
        size_labels.append(f"cov={n}")

    leg2 = ax_leg.legend(
        size_handles, size_labels,
        title="Dot size (coverage)",
        loc="lower left",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        borderaxespad=0.0
    )
    ax_leg.add_artist(leg2)

    fig.savefig(args.out, dpi=200)

if __name__ == "__main__":
    main()
