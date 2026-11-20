#!/usr/bin/env bash
# boot-validate.sh -- locate batches of zeros with method B (--ghy) and the
# self-consistent hybrid bootstrap (--boot W) at several heights, for
# comparison against the Odlyzko tables (see boot-validate.wls).
#
# Output TSV: tag <TAB> offset <TAB> method <TAB> gamma
# Usage: doc/ghy/boot-validate.sh [outfile]   (default doc/ghy/boot-validate.tsv)
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export ZZZ="$ROOT/build/zzz"
OUT="${1:-$ROOT/doc/ghy/boot-validate.tsv}"

# tag : ordinal base : first offset : count : k : W
SPECS="
low:1000:0:20:1000:32
mid:99900:0:20:1000:32
e12:1e12:30:20:10000:32
e21:1e21:30:20:10000:32
e22:1e22:30:10:10000:32
"
: > "$OUT"
for spec in $SPECS; do
  tag=$(echo "$spec" | cut -d: -f1); base=$(echo "$spec" | cut -d: -f2)
  off0=$(echo "$spec" | cut -d: -f3); cnt=$(echo "$spec" | cut -d: -f4)
  k=$(echo "$spec" | cut -d: -f5);   W=$(echo "$spec" | cut -d: -f6)
  for i in $(seq 0 $((cnt-1))); do
    echo "$tag $base $((off0+i)) $k $W"
  done
done | xargs -P8 -n5 sh -c '
  gb=$("$ZZZ" --ghy -k $3 -d 8 $1 $2 2>/dev/null)
  gc=$("$ZZZ" --boot $4 -k $3 -d 8 $1 $2 2>/dev/null)
  printf "%s\t%s\tB\t%s\n%s\t%s\tboot\t%s\n" "$0" "$2" "$gb" "$0" "$2" "$gc"
' >> "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT") rows)" >&2
