#!/bin/sh
# Gap-scale regression guard, run by ctest:  smoke.sh <case> <bindir> <zeros_file>
# Expected ordinates are Odlyzko's; tolerances sit ~4-30x above today's errors,
# so a real accuracy regression trips them while FLINT/compiler jitter does not.
set -u
case_=$1; bin=$2; zeros=$3
G1=14.134725141734693
G1000=1419.4224809459956

# check <label> <got> <expected> <tol>
check() {
    awk -v l="$1" -v g="$2" -v e="$3" -v t="$4" 'BEGIN {
        if (g == "") { printf("FAIL %s: no output\n", l); exit 1 }
        d = g - e; if (d < 0) d = -d
        if (!(d <= t)) { printf("FAIL %s: got %s, expected %s +- %s\n", l, g, e, t); exit 1 }
        printf("ok   %s: %s  (|err| %.1e <= %s)\n", l, g, d, t) }'
}

# rows <label> <expected_rows> <tsv_file>: right row count, every field numeric
rows() {
    awk -v l="$1" -v n="$2" '
        NF == 0 || /^#/ { next }
        { r++; for (i = 1; i <= NF; i++)
            if ($i !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
                printf("FAIL %s: non-numeric field \"%s\" on row %d\n", l, $i, r); exit 1 } }
        END { if (r != n) { printf("FAIL %s: %d rows, expected %d\n", l, r, n); exit 1 }
              printf("ok   %s: %d numeric rows\n", l, r) }' "$3"
}

case $case_ in
    zero1) check "zzz 1 (method A, k=100)"     "$("$bin/zzz" 1)"                     $G1    0.01 ;;
    ghy)   check "zzz --ghy -k 300 1000"        "$("$bin/zzz" --ghy -k 300 1000)"     $G1000 0.01 ;;
    boot)  check "zzz --boot 2 -k 300 1000"     "$("$bin/zzz" --boot 2 -k 300 1000)"  $G1000 0.02 ;;
    weil)  check "zzz --weil -k 300 1000 0 3"   "$("$bin/zzz" --weil -k 300 1000 0 3 | head -n 1)" $G1000 0.001 ;;
    grids)
        tmp=${TMPDIR:-/tmp}/zzz-smoke.$$; trap 'rm -f "$tmp"' EXIT
        "$bin/zghy"    10 60 100 1000               > "$tmp" && rows "zghy grid"    101 "$tmp" || exit 1
        "$bin/zproxy"  10 60 100 50                 > "$tmp" && rows "zproxy grid"  101 "$tmp" || exit 1
        "$bin/zhad"    "$zeros" 10 60 100 1000 20   > "$tmp" && rows "zhad grid"    101 "$tmp" || exit 1
        "$bin/zhybrid" "$zeros" 10 60 100 1000 20   > "$tmp" && rows "zhybrid grid" 101 "$tmp" || exit 1 ;;
    *) echo "unknown case: $case_"; exit 2 ;;
esac
