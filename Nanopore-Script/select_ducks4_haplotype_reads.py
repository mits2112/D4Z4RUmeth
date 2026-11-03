#!/usr/bin/env python3
"""
Select DUCKS4 reads into haplotype-specific groups.

Outputs:
  <OUTDIR>/<haplotype>_selected_rows.csv
  <OUTDIR>/<haplotype>_read_ids_aroundD4Z4.txt
  <OUTDIR>/10qA_PAS_undefined_selected_rows.csv
  <OUTDIR>/10qA_PAS_undefined_read_ids_aroundD4Z4.txt
  <OUTDIR>/selection_summary.tsv
  <OUTDIR>/selection.log

The selected_rows CSV files retain all columns available in the source CSV,
plus metadata columns beginning with "_selection_".
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple


GAP_RE = re.compile(r"gap_\d+bp", re.IGNORECASE)
DEFAULT_REQUIREMENTS = {
    "4qA": ">0",
    "4qB": ">0",
    "10q_long": ">12",
    "10q_shrt": "=12",
}
HAPLOTYPES = ("4qA", "4qB", "10q_long", "10q_shrt")
UNDEFINED_GROUP = "10qA_PAS_undefined"


@dataclass(frozen=True)
class AlleleRule:
    haplotype: str
    ru: int
    fasta: str
    ru_requirement: str


@dataclass
class SelectedRow:
    row: Dict[str, str]
    sources: List[str]
    reasons: List[str]


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().strip('"').strip("“”")


def parse_int(value: object) -> Optional[int]:
    text = clean(value)
    if text in {"", "NA", "NaN", "nan"}:
        return None
    try:
        return int(float(text.replace(",", ".")))
    except ValueError:
        return None


def normalize_header(name: str) -> str:
    return clean(name).lstrip("\ufeff")


def read_delimited(path: Path, delimiter: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header found: {path}")

        original_fields = [normalize_header(x) for x in reader.fieldnames]
        rows: List[Dict[str, str]] = []

        for raw in reader:
            row: Dict[str, str] = {}
            for old_name, new_name in zip(reader.fieldnames, original_fields):
                row[new_name] = clean(raw.get(old_name))
            rows.append(row)

    return original_fields, rows


def load_allele_rules(path: Path, log) -> "OrderedDict[str, AlleleRule]":
    fields, rows = read_delimited(path, "\t")
    required = {"haplotype", "RU", "fasta"}
    missing = required.difference(fields)
    if missing:
        raise ValueError(
            f"Allele file is missing columns {sorted(missing)}: {path}"
        )

    has_requirement = "RU_requirement" in fields
    if not has_requirement:
        log(
            "WARNING: RU_requirement column is absent from the allele file. "
            "Using built-in defaults: "
            + ", ".join(f"{k}:{v}" for k, v in DEFAULT_REQUIREMENTS.items())
        )

    rules: "OrderedDict[str, AlleleRule]" = OrderedDict()
    for row in rows:
        hap = clean(row.get("haplotype"))
        if not hap:
            continue
        if hap not in HAPLOTYPES:
            log(f"WARNING: unsupported haplotype in allele file; skipped: {hap}")
            continue

        ru = parse_int(row.get("RU"))
        if ru is None:
            raise ValueError(f"Invalid RU for {hap}: {row.get('RU')}")

        requirement = clean(row.get("RU_requirement"))
        if not requirement:
            requirement = DEFAULT_REQUIREMENTS[hap]

        # Normalize typographic quotation marks and == to =
        requirement = requirement.replace(" ", "").replace("==", "=")
        if not re.fullmatch(r"(?:>|>=|=|<=|<)\d+", requirement):
            raise ValueError(
                f"Unsupported RU_requirement for {hap}: {requirement!r}"
            )

        rules[hap] = AlleleRule(
            haplotype=hap,
            ru=ru,
            fasta=clean(row.get("fasta")),
            ru_requirement=requirement,
        )

    missing_haps = [h for h in HAPLOTYPES if h not in rules]
    if missing_haps:
        raise ValueError(f"Allele file lacks haplotypes: {', '.join(missing_haps)}")
    return rules


def ru_matches(value: Optional[int], requirement: str) -> bool:
    if value is None:
        return False
    match = re.fullmatch(r"(>=|<=|>|<|=)(\d+)", requirement)
    if not match:
        raise ValueError(f"Invalid RU requirement: {requirement}")
    op, rhs_text = match.groups()
    rhs = int(rhs_text)
    return {
        ">": value > rhs,
        ">=": value >= rhs,
        "<": value < rhs,
        "<=": value <= rhs,
        "=": value == rhs,
    }[op]


def has_gap(row: Mapping[str, str]) -> bool:
    return GAP_RE.search(clean(row.get("repeat_sequence"))) is not None


def exact_tokens(sequence: str) -> set[str]:
    # DUCKS4 repeat_sequence uses " - " as a token separator.
    return {part.strip() for part in clean(sequence).split(" - ") if part.strip()}


def add_selected(
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
    group: str,
    row: Mapping[str, str],
    source: str,
    reason: str,
) -> None:
    read_id = clean(row.get("read.id"))
    if not read_id:
        return

    copied = {normalize_header(k): clean(v) for k, v in row.items()}
    existing = groups[group].get(read_id)
    if existing is None:
        groups[group][read_id] = SelectedRow(
            row=copied,
            sources=[source],
            reasons=[reason],
        )
    else:
        if source not in existing.sources:
            existing.sources.append(source)
        if reason not in existing.reasons:
            existing.reasons.append(reason)


def load_profile(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        profile = json.load(handle)
    return profile


def select_standard_haplotype(
    hap: str,
    rule: AlleleRule,
    rows: Sequence[Mapping[str, str]],
    source: str,
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
) -> None:
    for row in rows:
        status = clean(row.get("status")).lower()
        ru = parse_int(row.get("RU_count"))

        if status == "complete":
            # Complete reads bypass the gap filter, but RU is still needed
            # to separate the two chr10 alleles.
            if hap.startswith("10q"):
                if ru_matches(ru, rule.ru_requirement):
                    add_selected(
                        groups, hap, row, source,
                        f"complete; RU_count {rule.ru_requirement}"
                    )
            else:
                add_selected(groups, hap, row, source, "complete")
            continue

        if status == "partial":
            if has_gap(row):
                continue
            if ru_matches(ru, rule.ru_requirement):
                add_selected(
                    groups, hap, row, source,
                    f"partial; no gap; RU_count {rule.ru_requirement}"
                )


def select_chr10(
    rules: Mapping[str, AlleleRule],
    rows: Sequence[Mapping[str, str]],
    source: str,
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
) -> None:
    long_rule = rules["10q_long"]
    short_rule = rules["10q_shrt"]
    short_ru = short_rule.ru

    for row in rows:
        status = clean(row.get("status")).lower()
        ru = parse_int(row.get("RU_count"))
        pas_type = clean(row.get("PAS.type"))

        if status not in {"complete", "partial"}:
            continue
        if status == "partial" and has_gap(row):
            continue

        status_reason = "complete" if status == "complete" else "partial; no gap"

        if ru_matches(ru, long_rule.ru_requirement):
            add_selected(
                groups, "10q_long", row, source,
                f"{status_reason}; RU_count {long_rule.ru_requirement}"
            )
            continue

        if ru_matches(ru, short_rule.ru_requirement):
            add_selected(
                groups, "10q_shrt", row, source,
                f"{status_reason}; RU_count {short_rule.ru_requirement}"
            )
            continue

        if (
            ru is not None
            and 1 < ru < short_ru
            and pas_type == "10qA_PAS"
        ):
            add_selected(
                groups, UNDEFINED_GROUP, row, source,
                f"{status_reason}; 1 < RU_count < {short_ru}; PAS.type=10qA_PAS"
            )


def select_chr4_undefined(
    rows: Sequence[Mapping[str, str]],
    source: str,
    profile: Mapping[str, object],
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
) -> None:
    cfg = profile.get("chr4_undefined", {})
    if not isinstance(cfg, Mapping):
        raise ValueError("profile.chr4_undefined must be an object")

    min_ru = int(cfg.get("min_ru_exclusive", 1))
    exclude_gap = bool(cfg.get("exclude_gap", True))

    for row in rows:
        ru = parse_int(row.get("RU_count"))
        if ru is None or ru <= min_ru:
            continue
        if exclude_gap and has_gap(row):
            continue

        tokens = exact_tokens(clean(row.get("repeat_sequence")))

        for hap in ("4qA", "4qB"):
            hap_cfg = cfg.get(hap, {})
            if not isinstance(hap_cfg, Mapping):
                continue
            required_tokens = {
                clean(x) for x in hap_cfg.get("required_exact_tokens", [])
            }
            if required_tokens and required_tokens.issubset(tokens):
                add_selected(
                    groups,
                    hap,
                    row,
                    source,
                    "chr4_undefined; "
                    f"RU_count > {min_ru}; no gap; exact token(s): "
                    + ",".join(sorted(required_tokens)),
                )


def add_chimeric_rows(
    rows: Sequence[Mapping[str, str]],
    source: str,
    targets: Sequence[str],
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
) -> None:
    for target in targets:
        if target not in groups:
            raise ValueError(f"Unsupported chimeric target: {target}")
        for row in rows:
            add_selected(
                groups, target, row, source,
                f"chimeric file explicitly assigned to {target}"
            )


def load_exclusions(path: Optional[Path], log) -> set[str]:
    if path is None:
        return set()
    if not path.exists():
        log(f"INFO: exclusion file does not exist; none applied: {path}")
        return set()

    excluded: set[str] = set()
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            read_id = clean(line)
            if read_id and not read_id.startswith("#"):
                excluded.add(read_id)
    return excluded


def apply_exclusions(
    groups: MutableMapping[str, "OrderedDict[str, SelectedRow]"],
    excluded: set[str],
) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for group, records in groups.items():
        hits = [read_id for read_id in records if read_id in excluded]
        for read_id in hits:
            del records[read_id]
        counts[group] = len(hits)
    return counts


def write_group_csv(
    path: Path,
    group: str,
    records: "OrderedDict[str, SelectedRow]",
) -> None:
    metadata_fields = [
        "_selection_haplotype",
        "_selection_source",
        "_selection_reason",
    ]

    field_order: List[str] = []
    seen: set[str] = set()
    for selected in records.values():
        for field in selected.row:
            if field not in seen:
                seen.add(field)
                field_order.append(field)

    # Keep read.id first when available.
    if "read.id" in field_order:
        field_order.remove("read.id")
        field_order.insert(0, "read.id")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=field_order + metadata_fields,
            extrasaction="ignore",
        )
        writer.writeheader()
        for selected in records.values():
            out = {field: selected.row.get(field, "") for field in field_order}
            out["_selection_haplotype"] = group
            out["_selection_source"] = " | ".join(selected.sources)
            out["_selection_reason"] = " | ".join(selected.reasons)
            writer.writerow(out)


def write_ids(path: Path, records: "OrderedDict[str, SelectedRow]") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for read_id in records:
            handle.write(read_id + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select DUCKS4 reads into haplotype-specific groups."
    )
    parser.add_argument("--alleles", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--profile", required=True, type=Path)

    parser.add_argument("--4qa-csv", dest="csv_4qa", type=Path)
    parser.add_argument("--4qb-csv", dest="csv_4qb", type=Path)
    parser.add_argument("--chr10-csv", type=Path)
    parser.add_argument("--chr4-undefined-csv", type=Path)
    parser.add_argument("--chimeric-csv", type=Path)
    parser.add_argument(
        "--chimeric-target",
        action="append",
        default=[],
        choices=list(HAPLOTYPES),
        help=(
            "Assign every row of chimeric_reads.csv to this haplotype. "
            "Repeat the option only when intentional."
        ),
    )
    parser.add_argument("--exclude-reads", type=Path)

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "selection.log"

    def log(message: str) -> None:
        print(message)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")

    log_path.write_text("", encoding="utf-8")

    rules = load_allele_rules(args.alleles, log)
    profile = load_profile(args.profile)

    groups: Dict[str, "OrderedDict[str, SelectedRow]"] = {
        group: OrderedDict()
        for group in (*HAPLOTYPES, UNDEFINED_GROUP)
    }

    def optional_rows(path: Optional[Path], delimiter: str, label: str):
        if path is None:
            log(f"INFO: {label} was not specified; skipped.")
            return []
        if not path.exists():
            log(f"INFO: {label} does not exist; skipped: {path}")
            return []
        if path.stat().st_size == 0:
            log(f"INFO: {label} is empty; skipped: {path}")
            return []
        fields, rows = read_delimited(path, delimiter)
        log(f"INFO: loaded {len(rows)} rows from {label}: {path}")
        return rows

    rows_4qa = optional_rows(args.csv_4qa, ";", "A_4qA_all-reads.csv")
    rows_4qb = optional_rows(args.csv_4qb, ";", "B_4qB_all-reads.csv")
    rows_chr10 = optional_rows(args.chr10_csv, ";", "chr10_all-reads.csv")
    rows_chr4u = optional_rows(
        args.chr4_undefined_csv, ";", "chr4_undefined_all-reads.csv"
    )
    rows_chim = optional_rows(args.chimeric_csv, ";", "chimeric_reads.csv")

    if rows_4qa:
        select_standard_haplotype(
            "4qA", rules["4qA"], rows_4qa,
            str(args.csv_4qa), groups
        )
    if rows_4qb:
        select_standard_haplotype(
            "4qB", rules["4qB"], rows_4qb,
            str(args.csv_4qb), groups
        )
    if rows_chr10:
        select_chr10(rules, rows_chr10, str(args.chr10_csv), groups)
    if rows_chr4u:
        select_chr4_undefined(
            rows_chr4u, str(args.chr4_undefined_csv), profile, groups
        )
    if rows_chim:
        if args.chimeric_target:
            add_chimeric_rows(
                rows_chim,
                str(args.chimeric_csv),
                args.chimeric_target,
                groups,
            )
        else:
            log(
                "INFO: chimeric_reads.csv exists, but no --chimeric-target "
                "was supplied; no chimeric reads were added."
            )

    excluded = load_exclusions(args.exclude_reads, log)
    removed = apply_exclusions(groups, excluded)

    summary_path = args.outdir / "selection_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "selected_reads", "excluded_reads"])
        for group, records in groups.items():
            writer.writerow([group, len(records), removed[group]])

    for group, records in groups.items():
        if group == UNDEFINED_GROUP:
            csv_name = "10qA_PAS_undefined_selected_rows.csv"
            txt_name = "10qA_PAS_undefined_read_ids_aroundD4Z4.txt"
        else:
            csv_name = f"{group}_selected_rows.csv"
            txt_name = f"{group}_read_ids_aroundD4Z4.txt"

        write_group_csv(args.outdir / csv_name, group, records)
        write_ids(args.outdir / txt_name, records)
        log(
            f"RESULT: {group}: {len(records)} selected; "
            f"{removed[group]} removed by exclusion list."
        )

    log(f"RESULT: summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
