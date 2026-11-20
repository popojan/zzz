# Snap-margin detectors: can a majorant window widen the zeros→primes tolerance?

*Synthesis of the 2026-06-21 experiment. The `--loop` recovers **exact** primes
from **approximate** zeros (gap-scale, median error ~0.036) because the backward
step is a quantizer onto $\mathbb{Z}$: zero error below the lattice resolution
snaps away (see [`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md), the
integer-snapping crack). This note tests a concrete way to **enlarge** that
tolerance — taper the zero-sum with a majorant window to kill the Dirichlet
sidelobes that drive false positives near the reach edge. Script:
`doc/ghy/maxsnap-detector.wls`; tables `doc/ghy/maxsnap-detector-nz{500,1200,4500}.txt`.*

## The idea under test

The default detector reads prime powers from the truncated
$$\psi'_T(x)\approx 1+\frac1{x-x^3}-\frac4{\sqrt x}\sum_{|\gamma_k|<T}\cos(\gamma_k\log x),$$
a **hard cutoff** $|\gamma_k|<T$ — i.e. a Dirichlet (sinc) kernel with large Gibbs
sidelobes. Near the reach edge those *deterministic* sidelobes, not Gaussian zero
jitter, are what manufacture a spurious peak at a composite (false prime) or a
ringing dip at a prime (miss). Hypothesis: a non-negative / tapered kernel (Fejér,
Hann, Gaussian) suppresses the sidelobes and so tolerates larger zero noise — at
the cost of a wider main lobe (less reach). Test it.

## Setup

Clean baseline = first $n_z$ Odlyzko zeros (`doc/refs/odlyzko/zeros1`, $\sim10^{-9}$).
Inject i.i.d. Gaussian jitter of size $\varepsilon$ (worst case — the real loop
error is partly common-mode, which cancels more). Detector = the project's
$\psi'/\mathrm{env}$ with a window $w(\gamma)$ on the cos-sum, **DC-normalized**
($\sum w=n_z$) so every window has the *same* prime-power peak height — only the
sidelobe structure differs, making the same $\tau=0.30$ threshold a fair test.
Scoring (external, loop-faithful): a failure is a **missed prime** ($x$ prime,
$r\le\tau$) or a **false positive at a genuine composite** ($x$ not a prime power,
$r>\tau$); higher prime powers $p^m$ ($m\ge2$, peak $\sim1/m$) are don't-cares, as
in `loop.c`. *Clean reach* = largest $X$ whose prefix $[2,X]$ is failure-free.

## Result

Median clean reach vs injected noise, frontier config ($n_z=500$, scan $2..220$,
20 reps; "220" = scan ceiling, i.e. detector still over-resolves):

| $\varepsilon$ | Dirichlet | Fejér | Hann | Gauss |
|---|---|---|---|---|
| 0.00 – 0.05 | **220** | **220** | 191 | 179 |
| 0.07 | **220** | 204 | 191 | 179 |
| **0.10** | 109 | **192** | 179 | 166 |
| 0.15 | 100 | 98 | 100 | 105 |
| 0.20 | 62 | 66 | 66 | 66 |

Snap margin $\varepsilon^*(X)$ = largest $\varepsilon$ holding reach $\ge X$:
at $X=120$ and $160$, Dirichlet $0.07$ vs Fejér/Hann/Gauss $0.10$.

Three regimes, consistent across $n_z=500/1200/4500$:

1. **The loop's operating regime — $\varepsilon\lesssim0.07$ (native $0.036$):**
   every window is flawless to the resolution limit, **with no noise sensitivity
   up to $\sim2\times$ the native error.** This quantifies the user's observation:
   approximate zeros ($\le7\%$ of a gap) give exact primes with comfortable margin,
   and *the detector choice is irrelevant here.* The readout sits far below its
   collapse threshold.

2. **Resolution cost of tapering ($\varepsilon=0$):** Dirichlet and Fejér reach the
   ceiling; **Hann and Gauss lose ~15–20% reach** to main-lobe widening *even at
   zero noise*. Aggressive tapers are strictly worse on resolution.

3. **Stress regime — $\varepsilon\approx0.10$–$0.15$ ($3$–$4\times$ native), sparse
   zeros:** the **Fejér majorant clearly widens the margin** — at $n_z=500$,
   $\varepsilon=0.10$ it holds reach 192 where Dirichlet collapses to 109 ($1.8\times$);
   $\varepsilon^*$ improves $0.07\!\to\!0.10$ ($\sim1.4\times$). This confirms the
   sidelobe mechanism: once Dirichlet's Gibbs ringing (relatively larger with few
   zeros) starts forging false positives, the non-negative Fejér kernel removes
   them. With *many* zeros ($n_z=4500$) the sidelobes are relatively small and the
   flat window's lower noise amplification ($\sum w^2$ is minimized by $w\equiv1$ at
   fixed peak) wins back — Dirichlet ties or leads at $\varepsilon\ge0.10$.

## Mechanism

At equal peak height, two effects compete. Flat ($w\equiv1$) **minimizes noise
amplification** $\propto\sum w_k^2$ and has the **narrowest main lobe** (best
reach) — it is optimal when the failure mode is Gaussian jitter or resolution.
Tapering **suppresses deterministic sidelobes** — it wins only when those
sidelobes are the failure mode, i.e. few zeros (big relative sidelobes) *and*
noise large enough to push the threshold into the ringing. The crossover is near
$\varepsilon\approx0.1$, well above the loop's native $0.036$.

## Verdict

The majorant idea is **partially confirmed but not a lever where the loop lives.**
Fejér (the Cesàro / non-negative majorant) is a **Pareto-safe drop-in**: never
worse than the hard cutoff on resolution, and $1.4$–$1.8\times$ better snap margin
under stress — so if the loop is ever run with a *thin* zero supply near collapse
(e.g. pushing $\kappa\to\pi$ on a sparse seed), swapping the hard cutoff for a
Fejér taper buys free noise headroom. But at the actual operating point the
readout is already flawless with $\sim2\times$ margin for **every** window, so the
dominant lever stays the $\kappa$-margin / zero count, exactly as
[`band-saturation.md`](band-saturation.md) predicts. Aggressive tapers (Hann,
Gauss) are counterproductive — they pay reach for a sidelobe benefit that only
matters past collapse.

This is itself a band-saturation instance: the window is an **in-band functional**;
reorganizing it buys a capped, constant-factor robustness, and that gain vanishes
in the regime where the loop operates. Integer rigidity makes the backward step
*robust*, and a majorant makes it *slightly more robust off-regime* — neither
synthesizes the out-of-band information that the $\sqrt T$ wall withholds.

## Artifact index

`doc/ghy/maxsnap-detector.wls` (env-vars `MSNZ`, `MSX` set $n_z$, scan width);
result tables `doc/ghy/maxsnap-detector-nz{500,1200,4500}.txt`; TSV of the last
run `doc/ghy/maxsnap-detector.tsv`. Companions:
[`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md) (the integer-snapping
crack the window tries to widen), [`band-saturation.md`](band-saturation.md)
(why in-band reorganization is capped).
