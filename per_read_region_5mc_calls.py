#!/usr/bin/env python3
"""
Per-read, per-BED-region CpG methylation summary from an aligned modBAM,
using `modkit extract calls`.

Rules:
  * mapped primary alignments only
  * a read is analyzed for a BED region only when its primary alignment
    completely spans the BED interval:
        alignment_start <= BED_start
        alignment_end   >= BED_end
  * CpG calls are obtained with `modkit extract calls --cpg`
  * calls with fail=true are excluded from methylation fractions
  * call_code:
        m  = 5mC
        h  = 5hmC
        -  = canonical/unmodified C
  * per-read/per-region 5mC fraction:
        n_5mC / n_informative
    where n_informative is the number of passing CpG calls
  * a second metric is also reported:
        (n_5mC + n_5hmC) / n_informative
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pysam


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str
    order: int


@dataclass(frozen=True)
class CpGCall:
    read_id: str
    chrom: str
    cpg_pos: int
    call_code: str
    call_prob: Optional[float]
    fail: bool
    ref_strand: str
    ref_mod_strand: str


def log(msg: str, handle) -> None:
    print(msg, flush=True)
    handle.write(msg + "\n")
    handle.flush()


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def parse_bool(value: object) -> bool:
    return clean(value).lower() in {"true", "1", "yes", "y"}


def parse_probability(value: object) -> Optional[float]:
    text = clean(value)
    if not text:
        return None
    try:
        x = float(text)
    except ValueError:
        return None
    if x < 0:
        return None
    if x <= 1.0:
        return x
    if x <= 255.0:
        return x / 255.0
    return None


def read_bed(path: Path) -> List[Region]:
    regions: List[Region] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {lineno} has fewer than 3 columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                raise ValueError(f"BED line {lineno}: end <= start")
            name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}:{start}-{end}"
            regions.append(Region(chrom, start, end, name, len(regions)))
    if not regions:
        raise ValueError("No BED regions found")
    return regions


def load_reference(fasta: Path, regions: Sequence[Region]) -> Dict[str, str]:
    with pysam.FastaFile(str(fasta)) as ref:
        needed = sorted({r.chrom for r in regions})
        missing = [chrom for chrom in needed if chrom not in ref.references]
        if missing:
            raise ValueError("BED contig(s) not present in FASTA: " + ", ".join(missing))
        return {chrom: ref.fetch(chrom).upper() for chrom in needed}


def cpg_positions_in_region(seq: str, region: Region) -> List[int]:
    positions: List[int] = []
    stop = min(region.end, len(seq) - 1)
    for pos in range(region.start, stop):
        if seq[pos:pos + 2] == "CG":
            positions.append(pos)
    return positions


def normalize_cpg_site(seq: str, ref_pos: int) -> Optional[int]:
    if ref_pos < 0 or ref_pos >= len(seq):
        return None
    if ref_pos + 1 < len(seq) and seq[ref_pos:ref_pos + 2] == "CG":
        return ref_pos
    if ref_pos > 0 and seq[ref_pos - 1:ref_pos + 1] == "CG":
        return ref_pos - 1
    return None


def ensure_bam_index(bam: Path, threads: int, log_handle) -> None:
    candidates = [Path(str(bam) + ".bai"), bam.with_suffix(".bai"), Path(str(bam) + ".csi")]
    if any(path.exists() for path in candidates):
        return
    log(f"INFO: BAM index not found; creating: {bam}", log_handle)
    pysam.index("-@", str(threads), str(bam))


def collect_full_span_reads(
    bam: Path,
    regions: Sequence[Region],
    threads: int,
) -> Dict[Tuple[str, str], Dict[str, object]]:
    eligible: Dict[Tuple[str, str], Dict[str, object]] = {}
    with pysam.AlignmentFile(str(bam), "rb", threads=threads) as alignment_file:
        bam_contigs = set(alignment_file.references)
        for region in regions:
            if region.chrom not in bam_contigs:
                continue
            for read in alignment_file.fetch(region.chrom, region.start, region.end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.reference_start is None or read.reference_end is None:
                    continue
                if read.reference_start <= region.start and read.reference_end >= region.end:
                    eligible[(read.query_name, region.name)] = {
                        "read_id": read.query_name,
                        "region": region.name,
                        "chrom": region.chrom,
                        "region_start": region.start,
                        "region_end": region.end,
                        "alignment_start": read.reference_start,
                        "alignment_end": read.reference_end,
                        "alignment_strand": "-" if read.is_reverse else "+",
                        "mapq": read.mapping_quality,
                    }
    return eligible


def run_modkit_extract_calls(
    bam: Path,
    bed: Path,
    fasta: Path,
    out_tsv: Path,
    threads: int,
    modkit: str,
    modkit_log: Path,
    log_handle,
) -> None:
    cmd = [
        modkit, "extract", "calls",
        str(bam), str(out_tsv),
        "--reference", str(fasta),
        "--include-bed", str(bed),
        "--mapped-only",
        "--cpg",
        "--threads", str(threads),
        "--log-filepath", str(modkit_log),
        "--force",
    ]
    log("INFO: running: " + " ".join(cmd), log_handle)
    subprocess.run(cmd, check=True)


def find_regions_for_cpg(
    chrom: str,
    cpg_pos: int,
    regions_by_chrom: Mapping[str, Sequence[Region]],
) -> List[Region]:
    return [
        region
        for region in regions_by_chrom.get(chrom, [])
        if region.start <= cpg_pos < region.end
    ]


def parse_modkit_calls(
    tsv: Path,
    seqs: Mapping[str, str],
    regions_by_chrom: Mapping[str, Sequence[Region]],
    eligible: Mapping[Tuple[str, str], Mapping[str, object]],
) -> Dict[Tuple[str, str], Dict[int, List[CpGCall]]]:
    calls: Dict[Tuple[str, str], Dict[int, List[CpGCall]]] = defaultdict(
        lambda: defaultdict(list)
    )

    with tsv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return calls

        required = {"read_id", "ref_position", "chrom", "call_prob", "call_code", "fail"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                "modkit extract calls output lacks columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            read_id = clean(row.get("read_id"))
            chrom = clean(row.get("chrom"))
            if chrom not in seqs:
                continue

            try:
                ref_pos = int(row["ref_position"])
            except (TypeError, ValueError):
                continue

            cpg_pos = normalize_cpg_site(seqs[chrom], ref_pos)
            if cpg_pos is None:
                continue

            call = CpGCall(
                read_id=read_id,
                chrom=chrom,
                cpg_pos=cpg_pos,
                call_code=clean(row.get("call_code")),
                call_prob=parse_probability(row.get("call_prob")),
                fail=parse_bool(row.get("fail")),
                ref_strand=clean(row.get("ref_strand")),
                ref_mod_strand=clean(row.get("ref_mod_strand")),
            )

            for region in find_regions_for_cpg(chrom, cpg_pos, regions_by_chrom):
                key = (read_id, region.name)
                if key in eligible:
                    calls[key][cpg_pos].append(call)

    return calls


def select_one_call(calls: Sequence[CpGCall]) -> Tuple[Optional[CpGCall], int]:
    if not calls:
        return None, 0

    def rank(call: CpGCall):
        prob = call.call_prob if call.call_prob is not None else -1.0
        return (0 if call.fail else 1, prob)

    selected = max(calls, key=rank)
    return selected, max(0, len(calls) - 1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--modkit", default="modkit")
    parser.add_argument("--keep-modkit-tsv", action="store_true")
    args = parser.parse_args()

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_csv = Path(str(args.out_prefix) + ".per_read_region.csv")
    per_cpg_csv = Path(str(args.out_prefix) + ".per_cpg.csv")
    main_log = Path(str(args.out_prefix) + ".log")
    modkit_log = Path(str(args.out_prefix) + ".modkit.log")

    with main_log.open("w", encoding="utf-8") as log_handle:
        for path, label in [(args.bam, "BAM"), (args.bed, "BED"), (args.fasta, "FASTA")]:
            if not path.exists():
                raise FileNotFoundError(f"{label} does not exist: {path}")

        if shutil.which(args.modkit) is None:
            raise RuntimeError(f"modkit command not found: {args.modkit}")

        regions = read_bed(args.bed)
        seqs = load_reference(args.fasta, regions)

        regions_by_chrom: Dict[str, List[Region]] = defaultdict(list)
        for region in regions:
            regions_by_chrom[region.chrom].append(region)

        cpg_by_region = {
            region.name: cpg_positions_in_region(seqs[region.chrom], region)
            for region in regions
        }

        log(f"INFO: regions: {len(regions)}", log_handle)
        log(
            "INFO: reference CpGs across BED regions: "
            + str(sum(len(x) for x in cpg_by_region.values())),
            log_handle,
        )

        ensure_bam_index(args.bam, args.threads, log_handle)

        eligible = collect_full_span_reads(args.bam, regions, args.threads)
        log(
            f"INFO: fully spanning primary read-region pairs: {len(eligible)}",
            log_handle,
        )

        if args.keep_modkit_tsv:
            modkit_tsv = Path(str(args.out_prefix) + ".modkit_extract_calls.tsv")
            cleanup = False
        else:
            tmp = tempfile.NamedTemporaryFile(
                prefix="modkit_extract_calls_",
                suffix=".tsv",
                delete=False,
            )
            tmp.close()
            modkit_tsv = Path(tmp.name)
            cleanup = True

        run_modkit_extract_calls(
            bam=args.bam,
            bed=args.bed,
            fasta=args.fasta,
            out_tsv=modkit_tsv,
            threads=args.threads,
            modkit=args.modkit,
            modkit_log=modkit_log,
            log_handle=log_handle,
        )

        try:
            calls = parse_modkit_calls(
                tsv=modkit_tsv,
                seqs=seqs,
                regions_by_chrom=regions_by_chrom,
                eligible=eligible,
            )

            region_order = {region.name: region.order for region in regions}
            sorted_keys = sorted(
                eligible,
                key=lambda key: (region_order.get(key[1], 10**9), key[0]),
            )

            summary_fields = [
                "read_id",
                "region",
                "chrom",
                "region_start",
                "region_end",
                "region_length",
                "alignment_start",
                "alignment_end",
                "alignment_strand",
                "mapq",
                "fully_spans_region",
                "n_cpg_reference",
                "n_cpg_with_any_call",
                "n_cpg_informative",
                "n_cpg_5mC",
                "n_cpg_5hmC",
                "n_cpg_canonical",
                "n_cpg_other_modified",
                "n_cpg_failed",
                "n_duplicate_call_rows",
                "informative_fraction_of_reference_CpGs",
                "methylation_fraction_5mC",
                "methylation_percent_5mC",
                "modified_cytosine_fraction_5mC_plus_5hmC",
                "modified_cytosine_percent_5mC_plus_5hmC",
                "cpg_calls",
            ]

            cpg_fields = [
                "read_id",
                "region",
                "chrom",
                "region_start",
                "region_end",
                "cpg_position_0based",
                "cpg_position_1based",
                "call_code",
                "call_label",
                "call_prob",
                "fail",
                "included_in_fraction",
                "ref_strand",
                "ref_mod_strand",
            ]

            with summary_csv.open("w", newline="", encoding="utf-8") as summary_handle, \
                 per_cpg_csv.open("w", newline="", encoding="utf-8") as cpg_handle:

                summary_writer = csv.DictWriter(summary_handle, fieldnames=summary_fields)
                cpg_writer = csv.DictWriter(cpg_handle, fieldnames=cpg_fields)

                summary_writer.writeheader()
                cpg_writer.writeheader()

                for key in sorted_keys:
                    read_id, region_name = key
                    meta = dict(eligible[key])
                    site_rows = calls.get(key, {})

                    n_any = 0
                    n_informative = 0
                    n_5mc = 0
                    n_5hmc = 0
                    n_canonical = 0
                    n_other = 0
                    n_failed = 0
                    n_duplicates = 0
                    compact_calls: List[str] = []

                    for cpg_pos in sorted(site_rows):
                        selected, duplicates = select_one_call(site_rows[cpg_pos])
                        n_duplicates += duplicates

                        if selected is None:
                            continue

                        n_any += 1

                        if selected.fail:
                            label = "fail"
                            included = False
                            n_failed += 1
                        elif selected.call_code == "m":
                            label = "5mC"
                            included = True
                            n_informative += 1
                            n_5mc += 1
                        elif selected.call_code == "h":
                            label = "5hmC"
                            included = True
                            n_informative += 1
                            n_5hmc += 1
                        elif selected.call_code == "-":
                            label = "canonical"
                            included = True
                            n_informative += 1
                            n_canonical += 1
                        else:
                            label = "modified:" + selected.call_code
                            included = True
                            n_informative += 1
                            n_other += 1

                        prob_text = (
                            ""
                            if selected.call_prob is None
                            else f"{selected.call_prob:.6f}"
                        )

                        compact_calls.append(
                            f"{cpg_pos}:{label}:{prob_text}"
                        )

                        cpg_writer.writerow({
                            "read_id": read_id,
                            "region": region_name,
                            "chrom": meta["chrom"],
                            "region_start": meta["region_start"],
                            "region_end": meta["region_end"],
                            "cpg_position_0based": cpg_pos,
                            "cpg_position_1based": cpg_pos + 1,
                            "call_code": selected.call_code,
                            "call_label": label,
                            "call_prob": prob_text,
                            "fail": str(selected.fail).upper(),
                            "included_in_fraction": str(included).upper(),
                            "ref_strand": selected.ref_strand,
                            "ref_mod_strand": selected.ref_mod_strand,
                        })

                    n_ref = len(cpg_by_region[region_name])

                    informative_fraction = (
                        n_informative / n_ref if n_ref > 0 else None
                    )

                    fraction_5mc = (
                        n_5mc / n_informative if n_informative > 0 else None
                    )

                    fraction_mod_c = (
                        (n_5mc + n_5hmc) / n_informative
                        if n_informative > 0
                        else None
                    )

                    summary_writer.writerow({
                        **meta,
                        "region_length":
                            int(meta["region_end"]) - int(meta["region_start"]),
                        "fully_spans_region": "TRUE",
                        "n_cpg_reference": n_ref,
                        "n_cpg_with_any_call": n_any,
                        "n_cpg_informative": n_informative,
                        "n_cpg_5mC": n_5mc,
                        "n_cpg_5hmC": n_5hmc,
                        "n_cpg_canonical": n_canonical,
                        "n_cpg_other_modified": n_other,
                        "n_cpg_failed": n_failed,
                        "n_duplicate_call_rows": n_duplicates,
                        "informative_fraction_of_reference_CpGs":
                            "" if informative_fraction is None else f"{informative_fraction:.6f}",
                        "methylation_fraction_5mC":
                            "" if fraction_5mc is None else f"{fraction_5mc:.6f}",
                        "methylation_percent_5mC":
                            "" if fraction_5mc is None else f"{100 * fraction_5mc:.3f}",
                        "modified_cytosine_fraction_5mC_plus_5hmC":
                            "" if fraction_mod_c is None else f"{fraction_mod_c:.6f}",
                        "modified_cytosine_percent_5mC_plus_5hmC":
                            "" if fraction_mod_c is None else f"{100 * fraction_mod_c:.3f}",
                        "cpg_calls": ";".join(compact_calls),
                    })

            log(f"RESULT: {summary_csv}", log_handle)
            log(f"RESULT: {per_cpg_csv}", log_handle)
            log(f"RESULT: {main_log}", log_handle)

        finally:
            if cleanup and modkit_tsv.exists():
                modkit_tsv.unlink()

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
