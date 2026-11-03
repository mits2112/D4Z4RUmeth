#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  cat >&2 <<'EOF'
Usage:
  bash make_condition_list_from_conditions.sh <output.tsv> "${CONDITIONS}"

Arguments:
  output.tsv   Output condition list TSV path
  CONDITIONS   Space-separated full names.
               Each value will be written to both full_name and short_name.

Output columns:
  full_name    Used by R_methyl_dot_v7.r to find <INDIR>/<full_name>.local_roi.cov.gz
  short_name   Used as the plot label

Example:
  CONDITIONS="5iLAF_SSEA4+ 5iLAF_SSEA4- Primed_SSEA4+"
  bash make_condition_list_from_conditions.sh condition_list.tsv "${CONDITIONS}"
EOF
  exit 1
fi

OUT=$1
CONDITIONS_STR=$2

read -r -a conditions <<< "${CONDITIONS_STR}"

if [[ ${#conditions[@]} -eq 0 ]]; then
  echo "[ERROR] CONDITIONS is empty." >&2
  exit 1
fi

out_dir=$(dirname "${OUT}")
if [[ "${out_dir}" != "." ]]; then
  mkdir -p "${out_dir}"
fi

printf "full_name\tshort_name\n" > "${OUT}"

for condition in "${conditions[@]}"; do
  if [[ -z "${condition}" ]]; then
    echo "[ERROR] Empty condition detected." >&2
    exit 1
  fi
  printf "%s\t%s\n" "${condition}" "${condition}" >> "${OUT}"
done

echo "[OK] Wrote ${OUT} with ${#conditions[@]} rows."
echo "[INFO] R script will look for files like: <INDIR>/<full_name>.local_roi.cov.gz"
