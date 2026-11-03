#!/usr/bin/env python3
import csv
import os
import re
import sys
from collections import Counter, defaultdict


def normalize_value(value: str) -> str:
    """Join feature values safely for source_name output."""
    if value is None:
        value = ""
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    return value


def get_samplefeatures() -> list[str]:
    if len(sys.argv) > 1:
        return [x for x in sys.argv[1:] if x]

    env = os.environ.get("samplefeatures", "")
    if env:
        return [x for x in env.split() if x]

    sys.exit(
        "[ERROR] samplefeatures is empty. "
        "Set bash variable/environment variable or pass columns as arguments.\n"
        "Example:\n"
        "  samplefeatures=\"source_name cell_state protocol\" ./extract_cols_dict_v3.py < SraRunTable.csv\n"
        "  export samplefeatures\n"
        "or\n"
        "  ./extract_cols_dict_v3.py source_name cell_state protocol < SraRunTable.csv"
    )


def get_assay_type_override() -> str:
    """Return user-specified Assay_Type override from environment."""
    return os.environ.get("Assay_Type", "")


reader = csv.DictReader(sys.stdin)
header = reader.fieldnames if reader.fieldnames else []

if not header:
    sys.exit("[ERROR] Input CSV header is empty")

samplefeatures = get_samplefeatures()
assay_type_override = get_assay_type_override()

missing = [col for col in samplefeatures if col not in header]
if missing:
    sys.exit(f"[ERROR] Columns not found: {', '.join(missing)}")

required_cols = ["Run", "LibraryLayout"]
if not assay_type_override:
    required_cols.append("Assay Type")

missing_required = [col for col in required_cols if col not in header]
if missing_required:
    sys.exit(f"[ERROR] Required columns not found: {', '.join(missing_required)}")

rows = []
src_list = []

for row in reader:
    run = row.get("Run", "")
    assay = assay_type_override if assay_type_override else row.get("Assay Type", "")
    layout = row.get("LibraryLayout", "")

    src_parts = [normalize_value(row.get(col, "")) for col in samplefeatures]
    src = "_".join(src_parts)

    rows.append({
        "run": run,
        "assay": assay,
        "layout": layout,
        "src": src,
    })
    src_list.append(src)

counts = Counter(src_list)

rows_sorted = sorted(rows, key=lambda r: r["run"])
rep_counters = defaultdict(int)
renamed_by_run = {}

for r in rows_sorted:
    src = r["src"]
    if counts[src] > 1:
        rep_counters[src] += 1
        new_src = f"{src}_rep{rep_counters[src]}"
    else:
        new_src = src
    renamed_by_run[r["run"]] = new_src

print("Run\tsource_name\tAssay_Type\tLibraryLayout")
for r in rows:
    print(f"{r['run']}\t{renamed_by_run[r['run']]}\t{r['assay']}\t{r['layout']}")
