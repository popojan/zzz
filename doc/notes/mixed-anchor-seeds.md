# Mixed exact/self-seeded anchors: how many exact zeros help, and where

*Session 2026-06-15. Follow-up to [`band-saturation.md`](band-saturation.md),
testing the doc's "an oracle buys the E1 gains" bullet quantitatively: if you
inject a **few exact zeros** among self-seeded ones, does the hybrid improve,
and how many do you need? All numbers measured against `doc/refs/odlyzko/zeros3`
(#10¹²+1…+10⁴, accurate to 1e-8). Driver: `/tmp/mixed_seed_exp.py`,
`/tmp/single_anchor.py` (ephemeral).*

## Setup

Vehicle: `zzz --seeds FILE -R 0` — seed a window of consecutive ordinates,
do **one** leave-one-out hybrid relocation of the centre (`-R 0` = relocate
target only, neighbours held exactly as supplied). This is the E1 oracle
experiment; the only variable is how many of the centre's neighbours are
**exact Odlyzko** vs **self-seeded method-B**. Metric: `|relocated − truth|`,
ensemble-averaged.

Below-threshold regime (κ = gap·log X < π) reached by lowering `k` at fixed
height 10¹² — the same proxy trick `band-saturation.md` uses (10²²/k=600 for
10³⁶). At 10¹², gap ≈ 0.2567, so κ(k=30)=1.21 stands in for "10³⁶ at any
feasible k".

## Result 1 — a handful of exact anchors recovers most of the oracle gain

Gain over plain B (16 centres), nearest `a` neighbours each side exact:

| k | κ | all-B (self-seeded) | ±1 | ±2 | ±4 | ±8 | ALL (oracle) |
|---|---|---|---|---|---|---|---|
| 30   | 1.21 | 0.86× | 0.76× | 1.49× | 2.14× | 2.17× | **2.74×** |
| 100  | 1.62 | 1.10× | 1.29× | 1.57× | 1.94× | 2.65× | **3.36×** |
| 300  | 1.95 | 1.00× | 2.20× | 3.26× | 2.63× | 2.86× | **3.44×** |
| 1000 | 2.30 | 1.21× | 1.89× | 2.27× | 2.10× | 3.06× | **3.04×** |
| 3000 | 2.62 | 1.23× | 2.47× | 4.10× | 6.28× | 4.54× | **4.87×** |
| 10⁵  | 3.61 | 1.23× | 4.61× | 4.43× | 3.66× |  —   | 3.74× (n=4) |

- **Self-seeding alone is inert below threshold** (all-B ≈ 1×, harmful at the
  bottom) — the band-saturation dead zone, reconfirmed.
- **Exact anchors unlock a bounded gain that saturates after a few.** The
  leave-one-out Z_X subtraction is dominated by the nearest zeros (kernel
  decays ≈ e^(−κ·d)); cleaning the nearest few captures the dominant error,
  far neighbours barely matter. Anchors needed each side ≈ kernel width ≈ 1/κ:
  ~±4 deep below (κ≈1.2), ~±2 mid-band, **±1 above threshold**.
- **The ceiling is the E1 value (2.7–4.9×), not a wall-crossing.** Even the
  all-exact oracle is a bounded constant-factor discount over B. Matches E1's
  reported 2.7–5.9× — the harness reproduces the known result.

## Result 2 — one *single* exact zero is the fragile regime

`±1` above means *two* zeros (one per side). A literally single one-sided
anchor (10 centres):

| | k=30 (κ=1.21) | k=300 (κ=1.95) | k=3000 (κ=2.62) |
|---|---|---|---|
| one zero (−1) | 0.86× | 1.04× | 1.45× |
| one zero (+1) | 1.41× | 1.20× | 1.60× |
| pair (±1)     | 1.19× | 1.73× | 2.00× |
| quad (±2)     | 1.57× | 1.89× | 3.34× |

A single exact zero is **noise** deep below threshold (sign flips with side at
κ=1.21) and only a reliable ~1.5× near threshold. The **cluster always beats
the singleton.** Practical rule: spend an exact-zero budget on a small
**contiguous cluster** (~±4) near the targets, not on isolated anchors.

## Result 3 — primes and exact zeros are substitutes; the exchange rate is brutal at height

κ = gap·log X, so raising κ costs **k = π(X) ≈ e^(κ/gap)/(κ/gap)** primes per
evaluation. At the showcase zero #10³⁶+… (T ≈ 8.10×10³⁴, gap ≈ 0.0800,
RS = √(T/2π) ≈ 1.14×10¹⁷):

| κ | k = π(X)/eval | k / RS | one zero buys |
|---|---|---|---|
| 1.21 | 2.4×10⁵  | 2×10⁻¹² | ≤0.76× (can mislead); cluster needed |
| 2.62 | 5.1×10¹² | 5×10⁻⁵  | a pair starts paying (~2.5×) |
| 3.14=π | 2.8×10¹⁵ | 0.025 | self-seeding works, exact zero redundant |
| 3.60 | 7.8×10¹⁷ | **>RS** | one zero gives full 4.6× — but past direct eval |

To make a lone zero pay by primes alone you must climb to κ≳2.6 (≈2×10⁷× more
primes), and the κ where one zero gives the *full* gain is already past
Riemann–Siegel — where self-seeding needs no oracle. **Raising κ until one
zero suffices converges to raising κ until you don't need it.** The efficient
direction at height is the opposite: keep primes cheap (k~10⁵–10⁶, κ~1.2–1.5)
and acquire a small cluster of exact zeros.

## Code: `--weil --seeds` and `--boot --seeds`

`--seeds FILE` now also feeds the Weil window fit (`main.c::weil_locate`):
odd count of consecutive ordinates, window centre on the middle line, marching
B pass skipped, in-window zeros refined toward the (possibly exact) prior,
ring held fixed. Verified: full ±20-t-unit all-truth seeds → centre error
**2.8×10⁻⁶** (data floor); B-seeds → B-quality — the fit holds whatever you
seed, the degeneracy made literal. **Caveat:** the file must span the full
±20 t-units (≈ ±⌈20/gap⌉ zeros), far more than `--boot`'s window; a short file
starves the fixed kernel and the fit diverges.

## Bottom line for the high-height hybrid

Inject exact zeros as a **small contiguous cluster** near the targets; you do
*not* need all seeds precise, and one isolated zero is the worst spend. The
payoff is a bounded ~3–5× over B, local to each cluster's kernel (~±2–4 zeros),
greatest just below threshold. It is the E1 oracle discount, confirmed — not a
route across the √T wall.
