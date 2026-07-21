# High-height verification of the damped locator (t = 1e18, t ≈ 1e26)

2026-07-21. Empirical check of the plain damped counter (`zzz -k 10000`)
far above the Odlyzko tables, against an independent referee: full
Riemann–Siegel Z(t) sweeps from a GPU-accelerated port of Bober's
zetacalc, on a δ = 0.01 grid, with actual zeros located by linear
interpolation of the sign changes.

## Protocol

* pick the center index n_c so the asymptotic location lands at the
  target height, via N(T) ≈ (T/2π)(log(T/2π) − 1) + 7/8;
* `zzz -k 10000 <n_c> <-w> <2w+1>` names 2w+1 consecutive zeros;
* the referee sweeps one window bracketing all predictions; each
  prediction is matched to the nearest sign change; offsets are reported
  in mean-gap units (2π/log(T/2π)).

## Results

| height | zeros | matched | mean abs offset | median | max |
|---|---|---|---|---|---|
| 1e18 | 11 | 11/11, distinct + ordered | 0.201 | 0.164 | 0.435 |
| ≈ 1e26 | 91 | 91/91, distinct + ordered | 0.198 | 0.155 | 0.640 |

Offsets in mean-gap units; k = 10000, default damping, tolerance 1e-6.

The 1e26 run: indices 907640000000000000000000000 ± 45
(`zzz -k 10000 907640000000000000000000000 -45 91`, ~3 s per zero),
heights near t = 9.99989e25, mean gap 0.10828. Referee window
t₀ = 99998910500341654095708115.978, 1100 points, span 11.0
(~1.2 h wall on a consumer GPU box). The window held 102 sign changes
vs ~101.5 expected from the zero density — the count closes, so no
prediction was matched against a spurious or missing crossing.

## Reading

* **Accuracy is flat in height**: mean |offset| ≈ 0.2 gaps at both
  1e18 and 1e26 — eight orders apart. In this regime the damped
  counter's error is governed by k, not by t (consistent with the
  1/log X law discussed in `doc/ghy.md`).
* **Indexing is exact**: every prediction matched a *distinct* actual
  zero, in order — zzz is naming the n-th zero correctly, not merely
  landing near some zero.
* The previous validation ceiling was γ ~ 4e5 (Odlyzko tables); this
  extends the empirically checked range by ~21 orders of magnitude in t.
