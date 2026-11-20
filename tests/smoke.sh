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

# primes <label> <min_largest> <file>: the --loop stream is exactly the primes up to its largest
primes() {
    awk -v l="$1" -v need="$2" '
        /^#/ || NF == 0 { next }
        { n++; got[$1] = 1; if ($1 + 0 > mx) mx = $1 + 0 }
        END { for (i = 2; i <= mx; i++) if (!c[i]) { np++; if (!got[i]) { printf("FAIL %s: prime %d missing\n", l, i); exit 1 }
                                                     for (j = 2 * i; j <= mx; j += i) c[j] = 1 }
                                        else if (got[i]) { printf("FAIL %s: composite %d emitted\n", l, i); exit 1 }
              if (n != np) { printf("FAIL %s: %d lines for %d primes\n", l, n, np); exit 1 }
              if (mx < need) { printf("FAIL %s: reached only %d, expected >= %d\n", l, mx, need); exit 1 }
              printf("ok   %s: exactly the %d primes up to %d\n", l, np, mx) }' "$3"
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
    loop)   # an earlier small-x branch emitted 18, 22, 24, ... from the 6-zero seed
        tmp=${TMPDIR:-/tmp}/zzz-smoke.$$; trap 'rm -f "$tmp" "$tmp.st"' EXIT
        L="--loop --loop-fresh --loop-iters 6 --loop-state $tmp.st"
        "$bin/zzz" $L                 > "$tmp" 2>/dev/null && primes "zzz --loop"                 1000 "$tmp" || exit 1
        "$bin/zzz" $L --loop-boxc 0   > "$tmp" 2>/dev/null && primes "zzz --loop --loop-boxc 0"   1000 "$tmp" || exit 1
        "$bin/zzz" $L --loop-seed 6   > "$tmp" 2>/dev/null && primes "zzz --loop --loop-seed 6"    500 "$tmp" || exit 1 ;;
    *) echo "unknown case: $case_"; exit 2 ;;
esac
