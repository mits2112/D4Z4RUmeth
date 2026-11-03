#!/usr/bin/env python3
import sys
import re
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Read FASTA (1 chr) and SAM from stdin, output BED with optional filters."
    )
    p.add_argument("fasta", help="reference fasta (single chromosome)")
    # --- filter options (default = previous hard-coded behavior) ---
    p.add_argument("--min-as", type=int, default=1000,
                   help="minimum AS score to keep (default: 1000)")
    p.add_argument("--max-as", type=int, default=None,
                   help="maximum AS score to keep (default: unset)")
    p.add_argument("--min-nm", type=int, default=None,
                   help="minimum NM to keep (default: unset)")
    p.add_argument("--max-nm", type=int, default=10000,
                   help="maximum NM to keep (default: 10000)")
    # optional: filter by strand
    p.add_argument("--strand", choices=["+", "-", "both"], default="both",
                   help="keep only alignments on this strand (default: both)")
    return p.parse_args()

args = parse_args()
fasta_path = args.fasta

# --- 1) Read chromosome name and length from FASTA ---
chrom_name = None
chrom_len = 0
with open(fasta_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if chrom_name is None:
                chrom_name = line[1:].split()[0]
        else:
            chrom_len += len(line)

if chrom_name is None or chrom_len == 0:
    sys.stderr.write("Error: Could not read chromosome name or length from FASTA.\n")
    sys.exit(1)

# --- CIGAR regex ---
cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

def get_hangs(cigar):
    """Return (left_hang, right_hang) from leading/trailing H/S in CIGAR."""
    ops = cigar_re.findall(cigar)
    if not ops:
        return 0, 0
    left_hang = 0
    right_hang = 0
    first_len, first_op = ops[0]
    last_len, last_op = ops[-1]
    if first_op in ('H', 'S'):
        left_hang = int(first_len)
    if last_op in ('H', 'S'):
        right_hang = int(last_len)
    return left_hang, right_hang

def pass_filters(ascore, nm, strand):
    # AS
    if args.min_as is not None and ascore < args.min_as:
        return False
    if args.max_as is not None and ascore > args.max_as:
        return False
    # NM
    if args.min_nm is not None and (nm < args.min_nm):
        return False
    if args.max_nm is not None and (nm > args.max_nm):
        return False
    # strand
    if args.strand != "both" and strand != args.strand:
        return False
    return True

# --- 2) Process SAM ---
for line in sys.stdin:
    line = line.rstrip("\n")
    if not line or line.startswith("@"):
        continue

    fields = line.split("\t")
    if len(fields) < 11:
        continue

    qname = fields[0]
    flag  = int(fields[1])
    cigar = fields[5]
    tags  = fields[11:]

    # --- Extract NM and AS tags ---
    nm = None
    ascore = None
    for t in tags:
        if t.startswith("NM:i:"):
            try:
                nm = int(t.split(":", 2)[2])
            except ValueError:
                nm = None
        elif t.startswith("AS:i:"):
            try:
                ascore = int(t.split(":", 2)[2])
            except ValueError:
                ascore = None

    if nm is None or ascore is None:
        continue

    # --- Compute strand ---
    strand = "-" if (flag & 16) else "+"

    # --- Apply filters from CLI ---
    if not pass_filters(ascore, nm, strand):
        continue

    # --- Compute hangs from CIGAR ---
    left_hang, right_hang = get_hangs(cigar)

    # --- Compute genomic coordinates depending on strand ---
    if strand == "+":
        # 1-based inclusive
        start1 = left_hang + 1
        end1   = chrom_len - right_hang
    else:
        # minus strand: left/right hangs are swapped on genome coordinates
        start1 = right_hang + 1
        end1   = chrom_len - left_hang

    # sanity check
    if end1 < start1:
        continue

    # BED (0-based start, half-open end)
    bed_start = start1 - 1
    bed_end   = end1

    # --- Output BED: chrom, start, end, name, score, strand, NM ---
    print("\t".join(map(str, [
        chrom_name,
        bed_start,
        bed_end,
        qname,
        ascore,
        strand,
        nm
    ])))
