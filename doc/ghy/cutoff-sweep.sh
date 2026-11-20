#!/usr/bin/env bash
# cutoff-sweep.sh -- locate test zeros with method B (zzz --ghy) for a
# log-spaced range of prime counts k, to study how the located zero
# gamma(X) oscillates with the cutoff X = p_k.
#
# Output TSV (one row per run): tag <TAB> k <TAB> gamma
# Companion analysis: cutoff-osc-analysis.wls (reads cutoff-sweep.tsv).
#
# Usage: doc/ghy/cutoff-sweep.sh [outfile]    (default: doc/ghy/cutoff-sweep.tsv)
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export ZZZ="$ROOT/build/zzz"
OUT="${1:-$ROOT/doc/ghy/cutoff-sweep.tsv}"

KS=$(awk 'BEGIN{k=200; while(k<=10000){print k; k=int(k*1.09)+1}}')
: > "$OUT"
for spec in "n1000:1000:0" "n10000:10000:0" "n99950:99950:0" "n100000:100000:0" "n1e12:1e12:50"; do
  tag=${spec%%:*}; rest=${spec#*:}; base=${rest%%:*}; off=${rest#*:}
  for k in $KS; do echo "$tag $k $base $off"; done
done | xargs -P8 -n4 sh -c 'g=$("$ZZZ" --ghy -k $1 -d 8 $2 $3 2>/dev/null); printf "%s\t%s\t%s\n" "$0" "$1" "$g"' >> "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT") rows)" >&2
