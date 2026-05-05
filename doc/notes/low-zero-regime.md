# `--ghy` low-zero regime

`zzz --ghy` is a hard truncation of the prime sum at `X = p_k`. GHY Theorem 1
gives a useful error bound only while `X ≲ T log T`; outside that regime
`P_X(1/2+iT)` no longer tracks `ζ` closely and B drifts.

The implementation `ghy_log_px` in `ghy.c` matches GHY eq. 6 — verified by
hand against the reference. The behaviour below is the math, not a bug.

## Empirical drift at zero #1

`T ≈ 14.13`, true value `14.134725...`. Run as `./build/zzz [-G] -k <k> 1 0`.

| k    | X = p_k | A finds | B finds (`--ghy`) |
|------|---------|---------|-------------------|
| 5    | 11      | 14.203  | 14.192            |
| 10   | 29      | 14.128  | 14.154            |
| 50   | 229     | 14.141  | **14.133**        |
| 100  | 541     | 14.132  | 14.103            |
| 300  | 1987    | 14.130  | 14.082            |
| 1000 | 7919    | 14.136  | 14.076            |
| 3000 | 27449   | 14.131  | 14.057            |

A converges to `|err| ~ 10⁻³` and stays there. B has a sweet spot near
`k ≈ 50` (`X ≈ T log T` ≈ 50) and drifts roughly linearly with k thereafter,
ultimately bounded by Goldston's `(log T)/(π log X)` projected through
`1/(dN/dT) = 2π/log(T/2π)`.

The pattern persists at zeros #2, #5, #10 — the cliff is sharpest at #1
and softens as T grows past `~50`.

## Mechanism

A's damping `1 − exp(−√(T/p))` suppresses primes with `p ≳ T`, so A behaves
as if it had a soft self-limiting cutoff at `p ≈ T`. B has no such
self-limiting — every prime in the truncation contributes with weight 1, and
once `X ≳ T log T` the GHY identity error term `X^{K+2}/(T log X)^K` blows up.

## Practical guidance

- **Production regime** (high-T zeros, `T ≥ 10⁶`, `k ≤ 10⁴`): `X / T ≪ 1`,
  B is firmly inside its rigor regime, and A and B agree to within ~1.5×
  (see `doc/ghy/ab-scan.wls` and `doc/notes/rigor-bound-b.md`).
- **Very low zeros** (`T ≲ 100`): prefer the heuristic A. With `--ghy`,
  raising k past `~π(T log T)` actively makes B worse.
- A clamp `k_max ≈ π(T log T)` would make `--ghy` safe across regimes; not
  yet implemented.

The earlier `doc/ghy/ab-scan.wls` scan started at ordinal #10 and missed
this regime; it is visible only at the very lowest zeros.
